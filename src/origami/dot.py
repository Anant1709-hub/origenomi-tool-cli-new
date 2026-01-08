#!/usr/bin/env python3
import os
import re
import subprocess
import matplotlib.pyplot as plt
from Bio import SeqIO


# ------------------------------------------------------------
# EXTRACT dnaA + oriC FROM REPORT
# ------------------------------------------------------------
def extract_coords(report_path):
    dnaA = None
    oriC = None

    with open(report_path) as f:
        lines = f.readlines()

    for i, line in enumerate(lines):
        if "Chosen dnaA:" in line:
            for j in range(i, min(i + 8, len(lines))):
                if "Coords:" in lines[j]:
                    nums = re.findall(r"\d+", lines[j])
                    if nums:
                        dnaA = int(nums[0])
                        break

        if "Chosen oriC:" in line:
            for j in range(i, min(i + 8, len(lines))):
                if "Coords:" in lines[j]:
                    nums = re.findall(r"\d+", lines[j])
                    if nums:
                        oriC = int(nums[0])
                        break

    if dnaA is None or oriC is None:
        raise ValueError("dnaA or oriC coordinates missing in report.")

    return dnaA, oriC


# ------------------------------------------------------------
# ROTATION CALCULATION (MATCH CIRCOS LOGIC)
# ------------------------------------------------------------
def compute_rotated_positions(dnaA, oriC, genome_len):
    # shortest-arc rule (same as circos)
    a, b = sorted([dnaA, oriC])
    d = b - a
    shift_bp = b if d > genome_len / 2 else a

    dnaA_rot = (dnaA - shift_bp) % genome_len
    oriC_rot = (oriC - shift_bp) % genome_len

    return shift_bp, dnaA_rot, oriC_rot


# ------------------------------------------------------------
# PARSE NUCMER COORDS
# ------------------------------------------------------------
def parse_blocks(coords_path):
    blocks = []
    with open(coords_path) as f:
        for line in f:
            parts = line.split()
            if len(parts) >= 4 and parts[0].isdigit():
                s1, e1, s2, e2 = map(int, parts[:4])
                blocks.append((s1, e1, s2, e2))

    return blocks


# ------------------------------------------------------------
# RUN NUCMER + SHOW-COORDS
# ------------------------------------------------------------
def run_nucmer(original_fna, rotated_fna, work_dir):
    prefix = os.path.join(work_dir, "dotplot_output")

    subprocess.run(
        ["nucmer", "--mum", "-g", "1000", "-p", prefix, original_fna, rotated_fna],
        check=True
    )

    coords_path = os.path.join(work_dir, "coords.txt")

    with open(coords_path, "w") as out:
        subprocess.run(
            ["show-coords", "-THrd", f"{prefix}.delta"],
            stdout=out,
            check=True
        )

    return coords_path


# ------------------------------------------------------------
# GENERATE DOT PLOT
# ------------------------------------------------------------
def generate_dotplot(blocks, dnaA, oriC, dnaA_rot, oriC_rot, out_path):
    plt.figure(figsize=(8, 8))

    # alignment blocks
    for s1, e1, s2, e2 in blocks:
        same_direction = (e2 - s2) * (e1 - s1) > 0
        color = "blue" if same_direction else "orange"
        plt.plot([s1, e1], [s2, e2], color=color, linewidth=1)

    # vertical markers (original genome)
    plt.axvline(dnaA, color="red", linewidth=1.8)
    plt.axvline(oriC, color="green", linestyle=":", linewidth=1.8)

    # horizontal markers (rotated genome)
    plt.axhline(dnaA_rot, color="red", linewidth=1.8)
    plt.axhline(oriC_rot, color="green", linestyle=":", linewidth=1.8)

    plt.title("Dot Plot (dnaA = red, oriC = green)", fontsize=11, fontweight="bold")
    plt.xlabel("Original Genome Position (bp)")
    plt.ylabel("Rotated Genome Position (bp)")
    plt.tight_layout()

    # Save
    png = out_path + ".png"
    svg = out_path + ".svg"

    plt.savefig(png, dpi=300)
    plt.savefig(svg)
    plt.close()

    return png


# ------------------------------------------------------------
# MAIN EXPORTED FUNCTION FOR FastAPI
# ------------------------------------------------------------
def plot_dotplot(original_fna, rotated_fna, report_file, output_dir):
    """
    FastAPI expects:
        original_fna  = basename string
        rotated_fna   = basename string
        report_file   = basename string
        output_dir    = absolute path

    Returns:
        PNG file path
    """
    original_path = os.path.join(output_dir, original_fna)
    rotated_path = os.path.join(output_dir, rotated_fna)
    report_path = os.path.join(output_dir, report_file)

    # extract dnaA/oriC
    dnaA, oriC = extract_coords(report_path)

    # genome length
    genome_len = len(next(SeqIO.parse(rotated_path, "fasta")).seq)

    # compute rotated alignment positions
    _, dnaA_rot, oriC_rot = compute_rotated_positions(dnaA, oriC, genome_len)

    # run nucmer + parse coords
    coords_path = run_nucmer(original_path, rotated_path, output_dir)
    blocks = parse_blocks(coords_path)

    # output prefix
    out_path = os.path.join(output_dir, "dotplot_highlighted_circos_logic")

    # create plot
    return generate_dotplot(blocks, dnaA, oriC, dnaA_rot, oriC_rot, out_path)
