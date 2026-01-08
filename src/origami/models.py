from dataclasses import dataclass
from typing import Optional, Tuple

@dataclass
class Hit:
    qseqid: str
    sseqid: str
    qstart: int
    qend: int
    sstart: int
    send: int
    qlen: int
    slen: int
    pident: float
    evalue: float
    bitscore: float
    mismatch: int
    gaps: int
    length: int

@dataclass
class PairingDecision:
    dnaA_id: Optional[str]
    dnaA_coords: Optional[Tuple[int, int]]
    oriC_id: Optional[str]
    oriC_coords: Optional[Tuple[int, int]]
    at_pct: float
    gc_pct: float
    start_coord: int