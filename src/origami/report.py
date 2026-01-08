from typing import List, Dict, Any

def write_section_header(report_path: str, title: str, append: bool = False) -> None:
    mode = "a" if append else "w"
    with open(report_path, mode) as f:
        f.write(f"{title}\n")

def _fmt_hit_line(h: Dict[str, Any]) -> str:
    return f"{h['sseqid']} {h['qstart']}…{h['qend']}"

def write_trim_record_block(report_path: str, header: str, trim_msg: str) -> None:
    with open(report_path, "a") as f:
        f.write(f"{header}\n")
        f.write(f"Trim: {trim_msg}\n\n")

def _pair_display_dict(pd: Dict[str, Any]) -> Dict[str, Any]:
    """
    Build the literal-style dict the spec wants for the 'pairs' lines:
    { <dnaA_id>: (sA,eA), <oriC_id>: (sO,eO), AT%: xx.x, GC%: yy.y }
    (Omit the None entries if a site is missing.)
    """
    disp = {"AT%": pd["AT%"], "GC%": pd["GC%"]}
    if pd.get("dnaA_id") and pd.get("dnaA_coords"):
        disp[pd["dnaA_id"]] = tuple(pd["dnaA_coords"])
    if pd.get("oriC_id") and pd.get("oriC_coords"):
        disp[pd["oriC_id"]] = tuple(pd["oriC_coords"])
    # Reorder to put IDs first, then AT/GC (purely cosmetic)
    id_items = [(k, v) for k, v in disp.items() if k not in ("AT%", "GC%")]
    return {**{k: v for k, v in id_items}, "AT%": disp["AT%"], "GC%": disp["GC%"]}

def write_oric_record_block(report_path: str,
                            header: str,
                            top_dnaa: List[Dict[str, Any]],
                            top_oric: List[Dict[str, Any]],
                            pairs: List[Dict[str, Any]],
                            chosen: Dict[str, Any]) -> None:
    with open(report_path, "a") as f:
        f.write(f"{header}\n")

        # Top hits
        f.write("Top dnaA hits (up to 3):\n")
        for h in top_dnaa:
            f.write(_fmt_hit_line(h) + "\n")
        f.write("Top oriC hits (up to 7):\n")
        for h in top_oric:
            f.write(_fmt_hit_line(h) + "\n")

        # Pairs list (header as you requested)
        f.write("\ndnaA - OriC pairs:\n")
        if pairs:
            for pd in pairs:
                f.write(str(_pair_display_dict(pd)) + "\n")
        else:
            f.write("(none)\n")

        # Chosen
        if chosen.get("dnaA_id") and chosen.get("dnaA_coords"):
            sA, eA = chosen["dnaA_coords"]
            f.write(f"\nChosen dnaA: {chosen['dnaA_id']} {sA}…{eA}\n")
        if chosen.get("oriC_id") and chosen.get("oriC_coords"):
            sO, eO = chosen["oriC_coords"]
            f.write(f"Chosen oriC: {chosen['oriC_id']} {sO}…{eO}\n")

        f.write(f"AT%: {chosen['AT%']}\n")
        f.write(f"GC%: {chosen['GC%']}\n")
        f.write("Reason: highest AT/GC ratio\n\n")

def write_extended_record_block(report_path: str, header: str,
                                chosen_dnaa: dict = None, chosen_oric: dict = None) -> None:
    """
    Write a block for the extended report. Handles wrap-around coordinates for oriC if needed.
    chosen_dnaa: dict with keys sseqid, qstart, qend
    chosen_oric: dict with keys sseqid and optionally wrap_info (is_wrapped, qstart, remaining, genome_len)
    """
    with open(report_path, "a") as f:
        f.write(f">{header}\n")

        # dnaA
        if chosen_dnaa:
            f.write("Chosen dnaA:\n")
            f.write(f"  ID: {chosen_dnaa['sseqid']}\n")
            f.write(f"  Coords: {chosen_dnaa['qstart']} - {chosen_dnaa['qend']}\n")
            f.write(f"  Score: {chosen_dnaa.get('bitscore', 'NA')}\n")
            f.write(f"  E-value: {chosen_dnaa.get('evalue', 'NA')}\n")
        else:
            f.write("No dnaA found\n")

        # oriC
        if chosen_oric:
            f.write("Chosen oriC:\n")
            f.write(f"  ID: {chosen_oric['sseqid']}\n")

            # Check for wrap-around info
            wrap_info = chosen_oric.get("wrap_info")
            if wrap_info and wrap_info.get("is_wrapped"):
                f.write(f"  Coords: {wrap_info['qstart']} - {wrap_info['genome_len']} and 1 - {wrap_info['remaining']}\n")
            else:
                qstart = chosen_oric.get("qstart")
                qend = chosen_oric.get("qend")
                if qstart is not None and qend is not None:
                    f.write(f"  Coords: {qstart} - {qend}\n")

            f.write(f"  Score: {chosen_oric.get('bitscore', 'NA')}\n")
            f.write(f"  E-value: {chosen_oric.get('evalue', 'NA')}\n")
        else:
            f.write("No oriC found\n")

        f.write("\n")


def write_original_record_block(report_path: str, header: str,
                                chosen_dnaa: dict, chosen_oric: dict,
                                genome_len: int):
    with open(report_path, "a") as f:
        f.write(f">{header}\n")

        # ---------------- dnaA ----------------
        if chosen_dnaa:
            f.write("Chosen dnaA:\n")
            f.write(f"  ID: {chosen_dnaa['sseqid']}\n")

            qstart, qend = chosen_dnaa["qstart"], chosen_dnaa["qend"]
            wrap_info = chosen_dnaa.get("wrap_info")

            # Case 1: dnaA completely in extension
            if qstart > genome_len and qend > genome_len:
                adj_start = qstart - genome_len
                adj_end = qend - genome_len
                f.write(f"  Coords: {adj_start} - {adj_end}\n")

            # Case 2: dnaA wraps across genome end
            elif wrap_info and wrap_info.get("is_wrapped"):
                qstart = wrap_info["qstart"]
                remaining = wrap_info["remaining"]
                genome_len = wrap_info["genome_len"]
                f.write(f"  Coords: {qstart} - {genome_len} and 1 - {remaining}\n")

            # Case 3: dnaA within normal genome
            else:
                f.write(f"  Coords: {qstart} - {qend}\n")

            f.write(f"  Score: {chosen_dnaa.get('bitscore', 'NA')}\n")
            f.write(f"  E-value: {chosen_dnaa.get('evalue', 'NA')}\n")
        else:
            f.write("No dnaA found\n")

        # ---------------- oriC ----------------
        if chosen_oric:
            f.write("Chosen oriC:\n")
            f.write(f"  ID: {chosen_oric['sseqid']}\n")

            qstart, qend = chosen_oric["qstart"], chosen_oric["qend"]
            wrap_info = chosen_oric.get("wrap_info")

            # Case 1: oriC completely in extension
            if qstart > genome_len and qend > genome_len:
                adj_start = qstart - genome_len
                adj_end = qend - genome_len
                f.write(f"  Coords: {adj_start} - {adj_end}\n")

            # Case 2: oriC wraps across genome end
            elif wrap_info and wrap_info.get("is_wrapped"):
                qstart = wrap_info["qstart"]
                remaining = wrap_info["remaining"]
                genome_len = wrap_info["genome_len"]
                f.write(f"  Coords: {qstart} - {genome_len} and 1 - {remaining}\n")

            # Case 3: oriC within normal genome
            else:
                f.write(f"  Coords: {qstart} - {qend}\n")

            f.write(f"  Score: {chosen_oric.get('bitscore', 'NA')}\n")
            f.write(f"  E-value: {chosen_oric.get('evalue', 'NA')}\n")
        else:
            f.write("No oriC found\n")

        f.write("\n")
