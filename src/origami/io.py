from typing import Iterator, Tuple

def stream_fasta(path: str) -> Iterator[Tuple[str, str]]:
    header = None
    seq_chunks = []
    with open(path) as f:
        for line in f:
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    yield header.rstrip(), "".join(seq_chunks)
                header = line.strip()
                seq_chunks = []
            else:
                seq_chunks.append(line.strip())
    if header is not None:
        yield header.rstrip(), "".join(seq_chunks)

def append_fasta_record(path: str, header: str, seq_with_info, wrap: int = 60) -> None:
    # If a tuple is passed, unpack sequence and info
    if isinstance(seq_with_info, tuple):
        seq, info = seq_with_info
        # Append info to the header
        header = f"{header} | trimmed={info}"
    else:
        seq = seq_with_info

    # Write to file
    with open(path, "a") as out:
        out.write(f"{header}\n")
        for i in range(0, len(seq), wrap):
            out.write(seq[i:i+wrap] + "\n")
