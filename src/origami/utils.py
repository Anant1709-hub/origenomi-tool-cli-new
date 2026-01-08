import os
import sys
import time

def midpoint(start: int, end: int) -> int:
    return (start + end) // 2

def gc_at_content(seq: str):
    at = sum(1 for b in seq if b in "AaTt")
    gc = sum(1 for b in seq if b in "GgCc")
    total = at + gc
    if total == 0:
        return 0.0, 0.0
    return round(at * 100.0 / total, 1), round(gc * 100.0 / total, 1)

def rotate_sequence_to_start(seq: str, start: int) -> str:
    start = max(0, min(len(seq), start))
    return seq[start:] + seq[:start]

def ensure_blast_db(prefix: str) -> str:
    """
    Ensure a BLAST nucl DB exists for the given prefix (expects .nin/.nsq/.nhr).
    Returns the prefix if valid; raises with a helpful message otherwise.
    """
    required = [prefix + ext for ext in (".nin", ".nsq", ".nhr")]
    missing = [p for p in required if not os.path.exists(p)]
    if missing:
        raise FileNotFoundError(
            f"BLAST DB missing for prefix '{prefix}'. "
            f"Missing: {', '.join(os.path.basename(m) for m in missing)}"
        )
    return prefix

class GlobalProgress:
    """
    One progress bar for the whole run.
    Weighted phases (defaults): Trim 40%, dnaA 20%, oriC 20%, Rotate 20%.
    Trimming and rotate scale with #records; BLAST phases flip when done.
    """
    def __init__(self, n_records: int, weights=(0.40, 0.20, 0.20, 0.20), label="origami"):
        self.n = max(1, n_records)
        self.w_trim, self.w_dnaa, self.w_oric, self.w_rotate = weights
        self.trim_i = 0
        self.rotate_i = 0
        self.dnaa_done = False
        self.oric_done = False
        self._t0 = None
        self.label = label

    def start(self):
        self._t0 = time.perf_counter()
        self._render(0.0)

    def update_trim(self, i: int):
        self.trim_i = min(i, self.n)
        self._render()

    def mark_dnaa_done(self):
        self.dnaa_done = True
        self._render()

    def mark_oric_done(self):
        self.oric_done = True
        self._render()

    def update_rotate(self, i: int):
        self.rotate_i = min(i, self.n)
        self._render()

    def finish(self):
        self._render(1.0, done=True)

    # ----- internals -----
    def _fraction(self):
        f_trim   = (self.trim_i / self.n) * self.w_trim
        f_dnaa   = (1.0 if self.dnaa_done else 0.0) * self.w_dnaa
        f_oric   = (1.0 if self.oric_done else 0.0) * self.w_oric
        f_rotate = (self.rotate_i / self.n) * self.w_rotate
        return min(1.0, f_trim + f_dnaa + f_oric + f_rotate)

    def _render(self, force_frac=None, done=False):
        frac = self._fraction() if force_frac is None else force_frac
        pct = int(frac * 100)
        width = 28
        filled = int(frac * width)
        bar = "[" + "=" * filled + " " * (width - filled) + "]"
        elapsed = 0.0 if self._t0 is None else (time.perf_counter() - self._t0)
        msg = f"\r{self.label} {bar} {pct:3d}%  elapsed: {elapsed:6.1f}s"
        if done:
            msg += "  ✓\n"
        sys.stdout.write(msg)
        sys.stdout.flush()
