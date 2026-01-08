#!/usr/bin/env python3
import os
import glob
import shutil
import subprocess


def run_linear_synteny_plot(work_dir):

    # FIND THE TWO INPUT FILES
    # print(f"Work Directory in the linear.py:{work_dir}")
    temp_work_dir = os.path.dirname(work_dir)
    fna_files = sorted(glob.glob(os.path.join(temp_work_dir, "*.fna")))

    for f in fna_files:
        if f.endswith(".fna"):

            if os.path.basename(f).startswith("GCF"):
                original_fna = f
            elif "origenomi" in f and "original" not in f:
                rotated_fna = f

    if not original_fna or not rotated_fna:
        # print("[ERROR] Could not find required origami FASTA files.")
        # print("Found:", fna_files)
        return None

    # print("[INFO] Using ORIGINAL :", original_fna)
    # print("[INFO] Using ROTATED  :", rotated_fna)

    # ------------------------------------------------------------
    # 3. COPY TO TEMPORARY NAMES
    # ------------------------------------------------------------
    tmp_original = os.path.join(work_dir, "Original.fna")
    tmp_rotated  = os.path.join(work_dir, "Origenomi.fna")

    shutil.copy(original_fna, tmp_original)
    shutil.copy(rotated_fna, tmp_rotated)

    # ------------------------------------------------------------
    # 4. RUN PMAUVE
    # ------------------------------------------------------------
    outdir = os.path.join(work_dir, "pgv_out")
    os.makedirs(outdir, exist_ok=True)

    cmd = [
        "pgv-pmauve",
        tmp_original, tmp_rotated,
        "-o", outdir,
        "--block_plotstyle", "bigbox",
        "--fig_width", "25",
        "--fig_track_height", "1.2",
        "--track_align_type", "left",
        "--block_cmap", "viridis",
        "--curve",
        "--show_scale_xticks",
        "--reuse",
        "--quiet"
    ]
    
#     subprocess.run(
#     cmd,
#     stdout=subprocess.DEVNULL,
#     stderr=subprocess.DEVNULL,
#     check=True
#    )

    # print("\n[INFO] Running PMAUVE:", " ".join(cmd))
    try:
        subprocess.run(cmd, check=True)
    except Exception as e:
        # print("[ERROR] pgv-pmauve failed:", e)
        return None

    # ------------------------------------------------------------
    # 5. SAVE PNG
    # ------------------------------------------------------------
    src = os.path.join(outdir, "result.png")
    dst = os.path.join(work_dir, "pgv.png")

    if os.path.exists(src):
        shutil.copy(src, dst)
        # print("Synteny plot saved as:", dst)
        return dst

    # print("result.png missing — pgv-pmauve failed")
    return None
