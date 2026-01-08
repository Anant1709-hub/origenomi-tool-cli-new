import os
import subprocess
from subprocess import DEVNULL
from origami.blast import run_blast_trim, parse_blast_trim

def detect_terminal_overlap(header: str, seq: str, paths, logger=None):
    """
    Split into halves; BLAST half2 (query) vs DB(half1).
    If a credible terminal overlap is found, trim exactly the aligned query span from the front.
    """
    if not seq:
        return seq, "no terminal overlap"

    # --- split into halves ---
    half = len(seq) // 2
    half1, half2 = seq[:half], seq[half:]

    # --- temp files under trim_dir ---
    h1_path   = os.path.join(paths.trim_dir, f"{paths.prefix}_first_half.fna")
    h2_path   = os.path.join(paths.trim_dir, f"{paths.prefix}_second_half.fna")
    db_prefix = os.path.join(paths.trim_dir, f"{paths.prefix}_half1_db")
    trim_out  = os.path.join(paths.trim_dir, f"{paths.prefix}_trim.blast")

    with open(h1_path, "w") as f:
        f.write(">half1\n" + half1 + "\n")
    with open(h2_path, "w") as f:
        f.write(">half2\n" + half2 + "\n")

    # --- build tiny DB from half1 (quiet) ---
    subprocess.run(
        ["makeblastdb", "-in", h1_path, "-dbtype", "nucl", "-out", db_prefix],
        check=True, stdout=DEVNULL, stderr=DEVNULL
    )
    if logger: logger.debug(f"makeblastdb → {db_prefix}")

    # --- BLAST half2 → DB(half1) with trim-specific flags/outfmt (quiet) ---
    run_blast_trim(h2_path, db_prefix, trim_out, evalue="1e-10", quiet=True)
    if logger: logger.debug(f"blastn (trim) → {trim_out}")

    # --- parse and pick the best terminal overlap ---
    all_hits = [h for h in parse_blast_trim(trim_out) if h["qseqid"] == "half2"]

    SLACK = 10  # allow slight slack at the ends
    candidates = []
    for h in all_hits:
        if h["pident"] < 95.0 or h["length"] < 100:
            continue

        # Accept either terminal configuration:
        #  (A) query near START & subject near END
        #  (B) query near END   & subject near START
        q_near_start = h["qstart"] <= SLACK
        q_near_end   = (h["qlen"] - h["qend"]) <= SLACK
        s_near_start = h["sstart"] <= SLACK
        s_near_end   = (h["slen"] - h["send"]) <= SLACK

        if (q_near_start and s_near_end) or (q_near_end and s_near_start):
            candidates.append(h)

    if not candidates:
        return seq, "no terminal overlap"

    # Longest overlap wins; tie-break on identity, then evalue
    candidates.sort(key=lambda h: (h["length"], h["pident"], -h["evalue"]), reverse=True)
    best = candidates[0]

    # Trim exactly the aligned query span (handle any orientation safely)
    overlap_len = abs(best["qend"] - best["qstart"]) + 1
    trimmed_seq = seq[overlap_len:]

    if logger:
        logger.debug(
            f"trim chosen: qstart={best['qstart']} qend={best['qend']} "
            f"len={best['length']} pident={best['pident']:.2f} → trim {overlap_len} bp"
        )

    return trimmed_seq, f"trimmed {overlap_len} bp from front"
