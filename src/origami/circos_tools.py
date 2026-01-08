import os
import re
import numpy as np
import matplotlib.pyplot as plt
from Bio import SeqIO

#  SHARED FUNCTION
def fast_gc_skew(seq, win=1000):
    """Compute GC skew using cumulative NumPy operations."""
    seq = np.frombuffer(seq.upper().encode(), dtype='S1')
    g = (seq == b'G').astype(int)
    c = (seq == b'C').astype(int)

    cum_g = np.cumsum(g)
    cum_c = np.cumsum(c)

    skew = np.empty(len(seq) // win)

    for i in range(len(skew)):
        start = i * win
        end = start + win

        g_win = cum_g[end - 1] - (cum_g[start - 1] if start > 0 else 0)
        c_win = cum_c[end - 1] - (cum_c[start - 1] if start > 0 else 0)

        denom = g_win + c_win
        skew[i] = (g_win - c_win) / denom if denom > 0 else 0

    return skew


def read_positions(report):
    """Extract dnaA and oriC coordinates from Origami report."""
    dnaA = None
    oriC_list = []

    with open(report) as f:
        for line in f:
            # dnaA coordinates
            if "Chosen dnaA:" in line:
                next(f)
                coords_line = next(f).strip()
                m = re.search(r"Coords:\s*(\d+)\s*-\s*(\d+)", coords_line)
                if m:
                    dnaA = int(m.group(1))

            # oriC coordinates (can be multiple)
            if "Chosen oriC:" in line:
                next(f)
                coords_line = next(f).strip()
                m = re.findall(r"(\d+)\s*-\s*(\d+)", coords_line)
                for a, b in m:
                    oriC_list.append((int(a), int(b)))

    return dnaA, oriC_list



def plot_circos(fna_path, report_path, out_dir):
    """
    This is the NON-rotated circos plot (Version B).
    Used for original genome orientation.
    """
    record = next(SeqIO.parse(fna_path, "fasta"))
    seq = str(record.seq)
    genome_len = len(seq)
    # print(f"Genome Length in non rotated circos:{genome_len}")
    dnaA, oriCs = read_positions(report_path)
    oriC = oriCs[0][0] if oriCs else None

    # Compute GC skew
    skew = fast_gc_skew(seq, win=1000)
    angles = np.linspace(0, 2*np.pi, len(skew), endpoint=False)

    # Polar plot setup
    fig, ax = plt.subplots(figsize=(10, 10), subplot_kw={'projection': 'polar'})
    ax.set_theta_direction(-1)
    ax.set_theta_offset(np.pi / 2)
    ax.set_axis_off()

    # Outer genome ring
    outer_r = 1.0
    ax.plot(np.linspace(0, 2*np.pi, 2000),
            [outer_r] * 2000,
            color='navy', lw=4)

    # GC skew bars
    bar_r = 0.8
    bar_scale = 0.25

    ax.bar(
        angles[skew > 0],
        skew[skew > 0] * bar_scale,
        width=2*np.pi/len(skew),
        bottom=bar_r,
        color='lightblue',
        linewidth=0
    )
    ax.bar(
        angles[skew < 0],
        -skew[skew < 0] * bar_scale,
        width=2*np.pi/len(skew),
        bottom=bar_r - bar_scale * np.abs(skew[skew < 0]),
        color='orange',
        linewidth=0
    )

    # Ticks (always 28)
    num_ticks = 28
    tick_spacing = genome_len / num_ticks
    label_offset = 0.10

    for i in range(num_ticks):
        pos = tick_spacing * i
        theta = 2 * np.pi * (pos / genome_len)

        ax.plot([theta, theta], [outer_r, outer_r + 0.03],
                color='navy', lw=1.5)
        ax.text(theta, outer_r + label_offset,
                f"{pos/1e6:.2f} Mb",
                ha='center', va='center', fontsize=9)

    # dnaA / oriC markers
    if dnaA:
        theta_dnaA = 2 * np.pi * (dnaA / genome_len)
        ax.plot([theta_dnaA], [outer_r], marker='s',
                color='blue', markersize=10)

    if oriC:
        theta_oriC = 2 * np.pi * (oriC / genome_len)
        ax.plot([theta_oriC], [outer_r], marker='o',
                color='red', markersize=10)

    # Legend
    fig.text(0.05, 0.93, u'\u25A0 dnaA', color='blue', fontsize=13, fontweight='bold')
    fig.text(0.05, 0.89, u'\u25CF oriC', color='red', fontsize=13, fontweight='bold')

    # Title text
    desc = record.description.split(",")[0]
    parts = desc.split(" ", 1)
    line1 = parts[0]
    line2 = parts[1] if len(parts) > 1 else ""

    ax.text(0, 0,
            f"{line1}\n{line2}\n({genome_len:,} bp)",
            ha='center', va='center',
            fontsize=10, fontweight='bold')

    # Save
    base = os.path.splitext(os.path.basename(fna_path))[0]
    png = os.path.join(out_dir, f"{base}_final_circos_v9.png")
    svg = os.path.join(out_dir, f"{base}_final_circos_v9.svg")

    plt.savefig(png, dpi=400, bbox_inches='tight')
    plt.savefig(svg, bbox_inches='tight')
    plt.close()

    return png  # main.py expects a PNG path



def plot_rotated_circos(fna_path, report_path, out_dir):
    """
    ROTATED circos plot (Version A).
    Adjusts genome so the origin (dnaA or oriC) starts at 0°.
    """
    record = next(SeqIO.parse(fna_path, "fasta"))
    seq = str(record.seq)
    genome_len = len(seq)
    # print(f"Genome Length in rotated circos:{genome_len}")
    dnaA, oriC_list = read_positions(report_path)
    oriC = oriC_list[0][0] if oriC_list else None

    # Determine shift point: whichever (dnaA or oriC) is further.
    a, b = sorted([dnaA, oriC])
    d = b - a

    if d > genome_len / 2:
        shift_bp = b
        rotated_from = "oriC" if b == oriC else "dnaA"
    else:
        shift_bp = a
        rotated_from = "oriC" if a == oriC else "dnaA"

    # Compute GC skew and rotate
    skew = fast_gc_skew(seq, win=1000)
    bins = len(skew)
    bin_size = genome_len / bins
    # shift_bins = int(shift_bp / bin_size)
    # skew = np.roll(skew, -shift_bins)

    # Polar plot
    fig, ax = plt.subplots(figsize=(10, 10), subplot_kw={'projection': 'polar'})
    ax.set_theta_direction(-1)
    ax.set_theta_offset(np.pi / 2)
    ax.set_axis_off()

    # Outer ring
    outer_r = 1.0
    ax.plot(np.linspace(0, 2*np.pi, 2000),
            [outer_r] * 2000,
            color='navy', lw=4)

    # GC skew bars
    r_inner = 0.8
    bar_scale = 0.25
    angles = np.linspace(0, 2*np.pi, bins, endpoint=False)

    ax.bar(angles[skew > 0], skew[skew > 0] * bar_scale,
           width=2*np.pi/bins, bottom=r_inner,
           color='lightblue', linewidth=0)

    ax.bar(angles[skew < 0], -skew[skew < 0] * bar_scale,
           width=2*np.pi/bins,
           bottom=r_inner - bar_scale*np.abs(skew[skew < 0]),
           color='orange', linewidth=0)

    # Ticks (28)
    num_ticks = 28
    tick_spacing = genome_len / num_ticks
    label_offset = 0.08

    for i in range(num_ticks):
        pos = tick_spacing * i
        theta = 2 * np.pi * i / num_ticks
        ax.plot([theta, theta], [outer_r, outer_r + 0.03],
                color='navy', lw=1.5)
        ax.text(theta, outer_r + label_offset,
                f"{pos/1e6:.2f} Mb",
                ha='center', va='center', fontsize=9)

    # Rotated marker positions
    dnaA_rot = (dnaA - shift_bp) % genome_len
    oriC_rot = (oriC - shift_bp) % genome_len

    ax.plot([2*np.pi*(dnaA_rot/genome_len)], [outer_r],
            marker='s', color='blue', markersize=10)
    ax.plot([2*np.pi*(oriC_rot/genome_len)], [outer_r],
            marker='o', color='red', markersize=10)

    # Center label
    desc = record.description.split(",")[0]
    parts = desc.split(" ", 1)
    line1 = parts[0]
    line2 = parts[1] if len(parts) > 1 else ""

    ax.text(0, 0,
            f"{line1}\n{line2}\n({genome_len:,} bp)",
            ha='center', va='center',
            fontsize=10, fontweight='bold')

    # Legend
    fig.text(0.05, 0.93, u'\u25A0 dnaA', color='blue', fontsize=13, fontweight='bold')
    fig.text(0.05, 0.89, u'\u25CF oriC', color='red', fontsize=13, fontweight='bold')

    # Save
    base = os.path.splitext(os.path.basename(fna_path))[0]
    png = os.path.join(out_dir, f"{base}_rotated_circular_gcskew.png")
    svg = os.path.join(out_dir, f"{base}_rotated_circular_gcskew.svg")

    plt.savefig(png, dpi=300, bbox_inches='tight')
    plt.savefig(svg, bbox_inches='tight')
    plt.close()

    # Return the PNG, the chosen marker, and the coordinate
    return png, rotated_from, shift_bp
