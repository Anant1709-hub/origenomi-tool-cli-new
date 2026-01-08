
def prepare_chromosome_sequence(header: str, sequence: str, duplicate_bp: int = 4000) -> str | None:
    """
    Prepares a bacterial chromosome sequence for BLAST:
    - Skips plasmids
    - Adds the first `duplicate_bp` bases at the end to handle circular overlap
    """
    # checking for plasmid in the header
    if "plasmid" in header.lower():
        return sequence  # skip plasmids

    # Make sure we don't exceed sequence length
    bp_to_duplicate = min(len(sequence), duplicate_bp)
    extended_seq = sequence + sequence[:bp_to_duplicate]
    return extended_seq


def trim_extended_if_needed(seq: str, hit: dict, extension_len: int = 4000):
    """
    Adjust sequence and coordinates:
      - If hit is in the extension (last extension_len bp), 
        drop first extension_len bp.
      - Otherwise, drop artificial tail.
    """
    seq_len = len(seq)
    orig_len = seq_len - extension_len

    if hit["qend"] > orig_len:
        # Hit was in extension → drop front
        return seq[extension_len:], "front"
    else:
        # No hit in extension → drop tail
        return seq[:-extension_len], "tail"
