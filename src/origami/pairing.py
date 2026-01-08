from typing import List, Dict, Any, Tuple
from math import log
from origami.utils import gc_at_content


def linear_edge_distance(d: Dict[str, Any], o: Dict[str, Any]) -> int:
    """Compute the edge-to-edge genomic distance between dnaA and oriC."""
    d_start, d_end = sorted((d["qstart"], d["qend"]))
    o_start, o_end = sorted((o["qstart"], o["qend"]))
    if o_end < d_start:
        return d_start - o_end
    elif d_end < o_start:
        return o_start - d_end
    else:
        return 0  # overlap


def adjust_for_extension(seq: str, chosen: Dict[str, Any], extension_len: int = 4000
                         ) -> Tuple[str, Dict[str, Any]]:
    """Trim genome extension if hits fall in the appended region and adjust coords."""
    genome_len = len(seq) - extension_len
    dnaa_coords = chosen.get("dnaA_coords")
    oric_coords = chosen.get("oriC_coords")

    if not dnaa_coords and not oric_coords:
        return seq, chosen

    start_positions = []
    if dnaa_coords:
        start_positions.append(dnaa_coords[0])
    if oric_coords:
        start_positions.append(oric_coords[0])

    # Default: no trimming
    seq_adj = seq
    offset = 0  # coordinate shift

    # Case 1: Hits fall in the *end extension* (we drop the first extension)
    if any(pos >= genome_len for pos in start_positions):
        seq_adj = seq[extension_len:]
        offset = extension_len * -1  # shift coordinates left by extension_len

    # Case 2: Hits fall in the *start extension* (we drop the tail extension)
    elif any(pos < extension_len for pos in start_positions):
        seq_adj = seq[:-extension_len]
        offset = 0  # coords stay the same

    else:
        seq_adj = seq[:-extension_len]
        offset = 0

    # --- Adjusted coordinates ---
    chosen["dnaA_coords_adj"] = None
    chosen["oriC_coords_adj"] = None

    if dnaa_coords:
        sA, eA = dnaa_coords
        chosen["dnaA_coords_adj"] = (sA + offset, eA + offset)
    if oric_coords:
        sO, eO = oric_coords
        chosen["oriC_coords_adj"] = (sO + offset, eO + offset)

    return seq_adj, chosen



def pair_and_choose(header: str, seq: str,
                    dnaA_hits: List[Dict[str, Any]],
                    oriC_hits: List[Dict[str, Any]],
                    extension_len: int = 4000
                    ) -> Tuple[List[Dict[str, Any]], Dict[str, Any], str]:
    """
    Pair dnaA and oriC hits based on:
        Edge-to-edge distance (not midpoint)
        Composite score favoring longer, higher-identity hits
        Tiered distance thresholds: 1%, then 5%, then whole genome without distance penalty
        Prefer longest oriC, and if equal, one in original genome region
    """
    genome_len = len(seq) - extension_len
    pairs: List[Dict[str, Any]] = []

    # Deduplicate dnaA hits (top 10 by evalue)
    unique_dnaa_hits = []
    seen_coords_dnaa = set()
    for d in dnaA_hits:
        coords = (d["qstart"], d["qend"])
        if coords not in seen_coords_dnaa:
            unique_dnaa_hits.append(d)
            seen_coords_dnaa.add(coords)
    unique_dnaa_hits = sorted(unique_dnaa_hits, key=lambda x: x.get("evalue", 0))[:10]

    # Filter and deduplicate oriC hits (>70bp)
    filtered_oric_hits = [o for o in oriC_hits if abs(o["qend"] - o["qstart"]) + 1 > 70]
    unique_oric_hits = []
    seen_coords_oric = set()
    for o in filtered_oric_hits:
        coords = (o["qstart"], o["qend"])
        if coords not in seen_coords_oric:
            unique_oric_hits.append(o)
            seen_coords_oric.add(coords)

    # Deduplicate by sseqid, prefer longest oriC, then prefer original region
    chosen_orics = {}
    for o in unique_oric_hits:
        sid = o["sseqid"]
        length = abs(o["qend"] - o["qstart"]) + 1
        in_original = o["qstart"] < extension_len

        if sid not in chosen_orics:
            chosen_orics[sid] = o
        else:
            prev = chosen_orics[sid]
            prev_length = abs(prev["qend"] - prev["qstart"]) + 1
            prev_in_original = prev["qstart"] < extension_len

            if length > prev_length or (length == prev_length and in_original and not prev_in_original):
                chosen_orics[sid] = o

    unique_oric_hits = list(chosen_orics.values())

    # Helper score functions
    def gc_at_for_hit(o):
        sO, eO = sorted((o["qstart"], o["qend"]))
        return gc_at_content(seq[sO:eO + 1])

    def composite_score(o: Dict[str, Any]) -> float:
        pid = o.get("pident", 0)
        length = abs(o["qend"] - o["qstart"]) + 1
        evalue = o.get("evalue", 1e-9)
        return (pid * log(length + 1)) / (1 + evalue)

    def overlaps(d: Dict[str, Any], o: Dict[str, Any]) -> bool:
        d_start, d_end = sorted((d["qstart"], d["qend"]))
        o_start, o_end = sorted((o["qstart"], o["qend"]))
        return d_start <= o_end and o_start <= d_end

    # ---------------------------
    # Tiered pairing logic
    # ---------------------------
    tier_limits = [
        (0.01 * genome_len, True),    # Tier 1: <= 1%, penalize distance
        (0.05 * genome_len, True),    # Tier 2: <= 5%, penalize distance
        (float("inf"), False)         # Tier 3: any distance, NO penalty
    ]

    for d in unique_dnaa_hits:
        best_ori = None
        best_score = float("-inf")
        best_dist = None

        for max_dist, penalize in tier_limits:
            candidates = []

            for o in unique_oric_hits:
                if overlaps(d, o):
                    continue

                dist = linear_edge_distance(d, o)
                if dist > max_dist:
                    continue

                score = composite_score(o)
                if penalize:
                    score -= dist * 0.001

                candidates.append((score, o, dist))

            if candidates:
                best_score, best_ori, best_dist = max(candidates, key=lambda x: x[0])
                break  # stop at the first tier that yields results

        if best_ori:
            at, gc = gc_at_for_hit(best_ori)
            pairs.append({
                "dnaA_id": d["sseqid"],
                "dnaA_coords": (d["qstart"], d["qend"]),
                "oriC_id": best_ori["sseqid"],
                "oriC_coords": (best_ori["qstart"], best_ori["qend"]),
                "AT%": at,
                "GC%": gc,
                "distance": best_dist,
                "pair_score": best_score,
            })

    # Choose final best pair
    chosen = None
    best_score = None
    for p in pairs:
        combined_score = p["pair_score"] - p["distance"] * 0.001
        if best_score is None or combined_score > best_score:
            best_score = combined_score
            chosen = p

    # Fallback if no matches
    if not chosen:
        chosen = {
            "dnaA_id": None,
            "dnaA_coords": None,
            "oriC_id": None,
            "oriC_coords": None,
            "AT%": 0.0,
            "GC%": 0.0,
            "distance": 0,
            "pair_score": 0.0,
        }

    # Wrap info if oriC extends past genome length
    if chosen.get("oriC_coords"):
        sO, eO = chosen["oriC_coords"]
        chosen["wrap_info"] = {
            "is_wrapped": True,
            "qstart": sO,
            "remaining": eO - genome_len,
            "genome_len": genome_len
        } if eO > genome_len else None

    if chosen.get("dnaA_coords"):
        sA, eA = chosen["dnaA_coords"]
        chosen["dnaa_wrap_info"] = {
            "is_wrapped": True,
            "qstart": sA ,
            "remaining": eA - genome_len,
            "genome_len": genome_len
        } if eA > genome_len else None

    # Adjust for extension
    seq_adjusted, chosen = adjust_for_extension(seq, chosen, extension_len)
    return pairs, chosen, seq_adjusted
