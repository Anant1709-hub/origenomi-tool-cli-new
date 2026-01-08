import os
import shutil
import subprocess
from Bio import SeqIO

from .circos_tools import plot_circos, plot_rotated_circos
from .dot import plot_dotplot
from .linear import run_linear_synteny_plot
import re


# HELPERS
# def extract_fna_from_gbff(gbff_path: Path, output_fna: Path):
#     version_line = None
#     definition_line = None
#     sequence_lines = []
#     in_origin = False

#     with open(gbff_path, "r") as f:
#         for line in f:
#             if line.startswith("VERSION"):
#                 version_line = ">" + line.split("VERSION", 1)[1].strip()
#             elif line.startswith("DEFINITION"):
#                 definition_line = line.split("DEFINITION", 1)[1].strip()
#             elif line.startswith("ORIGIN"):
#                 in_origin = True
#             elif in_origin:
#                 if line.startswith("//"):
#                     break
#                 seq = "".join(filter(str.isalpha, line)).upper()
#                 sequence_lines.append(seq)

#     if not version_line or not definition_line or not sequence_lines:
#         raise ValueError(f"Invalid or incomplete GBFF file: {gbff_path}")

#     with open(output_fna, "w") as out:
#         out.write(version_line + " " + definition_line + "\n")
#         out.write("".join(sequence_lines) + "\n")

#     print(f"[INFO] Extracted .fna saved to {output_fna}")


def parse_report(report_path):
    """
    Parses Origami original_report.txt and extracts dnaA + oriC coordinates
    exactly the same way the ORIGINAL working code does inside read_positions():

        dnaA = first match of:  Coords: A - B
        oriC = all matches of:  A - B   (even split cases)

    Supports wrap-around e.g.:
        "2309025 - 2309262 and 1 - 308"
    → oriC_list = [(2309025,2309262), (1,308)]
    """

    results = {}
    current_seq = None

    with open(report_path, "r") as f:
        for line in f:
            line = line.strip()

            if line.startswith(">>"):
                seq_id = line[2:].split()[0]
                results[seq_id] = {
                    "dnaA_found": False,
                    "oriC_found": False,
                    "dnaA_coords": None,       # raw string
                    "oriC_coords": None,       # raw string
                    "dnaA_start": None,        # numeric for plotting
                    "oriC_list": [],           # list of (start,end)
                }
                current_seq = seq_id
                continue

            # ────────────────────────────────────────────────
            # dnaA region
            # ────────────────────────────────────────────────
            if "Chosen dnaA:" in line:
                next_line = next(f).strip()      # ID line
                coords_line = next(f).strip()    # Coords: ...
                # Save raw coords string
                results[current_seq]["dnaA_coords"] = coords_line.replace("Coords:", "").strip()

                # Extract numeric A-B using regex (same as original code)
                m = re.search(r"(\d+)\s*-\s*(\d+)", coords_line)
                if m:
                    results[current_seq]["dnaA_found"] = True
                    results[current_seq]["dnaA_start"] = int(m.group(1))
                continue

            # ────────────────────────────────────────────────
            # oriC region (may include wrap-around)
            # ────────────────────────────────────────────────
            if "Chosen oriC:" in line:
                next_line = next(f).strip()      
                coords_line = next(f).strip()    

                raw_coords = coords_line.replace("Coords:", "").strip()
                results[current_seq]["oriC_coords"] = raw_coords

                pairs = re.findall(r"(\d+)\s*-\s*(\d+)", raw_coords)
                oriC_intervals = [(int(a), int(b)) for a, b in pairs]

                if oriC_intervals:
                    results[current_seq]["oriC_found"] = True
                    results[current_seq]["oriC_list"] = oriC_intervals

                continue

    return results




def generate_all_plots(fasta_output_path, report_path, rotated_fasta_path, work_dir, gcf_fasta_path):

    plot_dir = os.path.join(work_dir, "plots")
    os.makedirs(plot_dir, exist_ok=True)

    # 1. Parse report to get per-sequence results
    report_info = parse_report(report_path)

    # print("[LOG] Parsed report results:")
    # for seq_id, data in report_info.items():
    #     print("  ", seq_id, data)

    # 2. Loop over sequences
    for seq_id, data in report_info.items():

        if not (data["dnaA_found"] and data["oriC_found"]):
            # print(f"[LOG] {seq_id}: dnaA or oriC missing — skipping plots")
            continue

        # print(f"[LOG] {seq_id}: Generating plots")

        orig_plot = plot_circos(gcf_fasta_path, report_path,out_dir=work_dir)
        shutil.copy(orig_plot, os.path.join(plot_dir, f"{seq_id}_original_circos.png"))

        rot_result = plot_rotated_circos(rotated_fasta_path, report_path,out_dir=work_dir)
        if not rot_result:
            continue

        rot_png = rot_result[0] if isinstance(rot_result, (tuple, list)) else rot_result
        shutil.copy(rot_png, os.path.join(plot_dir, f"{seq_id}_rotated_circos.png"))

        dot_png = plot_dotplot(
            os.path.basename(fasta_output_path),
            os.path.basename(rotated_fasta_path),
            os.path.basename(report_path),
            output_dir=work_dir,
        )
        if dot_png:
            shutil.copy(dot_png, os.path.join(plot_dir, f"{seq_id}_dotplot.png"))

        # ===== LINEAR SYNTENY =====
        synteny_dir = os.path.join(work_dir, "linear_synteny")
        os.makedirs(synteny_dir, exist_ok=True)

        synteny_png = run_linear_synteny_plot(synteny_dir)
        if synteny_png:
            shutil.copy(synteny_png, os.path.join(plot_dir, f"{seq_id}_linear_synteny.png"))

    # print("[LOG] All multi-sequence plot generation complete.")
