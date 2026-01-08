# Origami

Command-line tool that:
- Trims duplicated terminal overlaps in circular assemblies (half2 → DB(half1) via BLAST).
- Runs **one BLAST per DB** (dnaA and oriC) on the concatenated FASTA.
- Keeps **top 3 dnaA** and **top 7 oriC** per record (dedup by `sseqid`, ordered by evalue→identity→coverage→bitscore).
- Pairs dnaA–oriC via midpoint proximity (1% then 5% of record length, else closest), computes AT/GC on the oriC segment, and chooses the pair with **highest AT/GC ratio** (tie → higher AT%).
- Rotates so the **earlier** of dnaA/oriC starts the sequence (or the single site if only one found).
- Writes exactly two outputs: `origami_<prefix>.fna` and `origami_<prefix>_report.txt`.

## Install
```bash
pip install -e .
```

## Usage

### Full pipeline
```bash
origami run -i genome.fna --dnaA-db /path/new_dnaA_DB --oric-db /path/DoriC_DB
```

### Trim only
```bash
origami trim -i genome.fna
```

### OriC only (no trim)
```bash
origami oric -i genome.fna --dnaA-db /path/new_dnaA_DB --oric-db /path/DoriC_DB
```

## Outputs
1. origami_<prefix>.fna
2. origami_<prefix>_report.txt
3. Temps in ./temp/<prefix>/, cleaned unless --keep-temp.