def select_top_hits(hits, k):
    """
    Rank hits by:
    1. Percent identity (descending)
    2. E-value (ascending)
    3. Coverage (descending)
    4. Bit score (descending)

    Do NOT deduplicate here — duplicates will be handled later in pairing.
    Return top-k hits.
    """
    def _coverage(h):
        return (h["qend"] - h["qstart"] + 1) / max(1, h.get("qlen", (h["qend"] - h["qstart"] + 1)))

    # Sort by pident desc, evalue asc, coverage desc, bitscore desc
    ordered = sorted(
        hits,
        key=lambda h: (
            -h["pident"],
            h["evalue"],
            -_coverage(h),
            -h["bitscore"],
        ),
    )

    return ordered[:k]
