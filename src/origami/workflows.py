import os
import shutil
from .main import generate_all_plots
from origami.circularize import prepare_chromosome_sequence, trim_extended_if_needed
from origami.paths import resolve_paths
from origami.io import stream_fasta, append_fasta_record
from origami.trim import detect_terminal_overlap
from origami.blast import run_blast_annotation, parse_blast_annotation
from origami.hits import select_top_hits
from origami.pairing import pair_and_choose
from origami.rotate import rotate_sequence
from origami.report import (
    write_section_header,
    write_trim_record_block,
    write_oric_record_block,
    write_extended_record_block,
    write_original_record_block,   
)
from origami.utils import GlobalProgress
import glob


# Helper functions
def _finalize(paths, keep_temp: bool, logger):
    if not keep_temp and os.path.isdir(paths.temp_dir):
        shutil.rmtree(paths.temp_dir, ignore_errors=True)
        if logger:
            logger.debug(f"Cleaned temp dir {paths.temp_dir}")


def _record_id(header: str) -> str:
    """Normalize FASTA header to qseqid used by BLAST."""
    return header[1:].split()[0] if header.startswith(">") else header.split()[0]


def _count_records(fasta_path: str) -> int:
    n = 0
    with open(fasta_path) as f:
        for line in f:
            if line.startswith(">"):
                n += 1
    return n


# run_trim_only()
def run_trim_only(args, logger=None):
    paths = resolve_paths(args.input, args.out_dir)
    write_section_header(paths.report, "Trim Analysis")

    n = _count_records(args.input)
    gp = GlobalProgress(n, weights=(1.0, 0.0, 0.0, 0.0), label="origami")
    gp.start()

    i = 0
    for header, seq in stream_fasta(args.input):
        trimmed_seq, trim_msg = detect_terminal_overlap(header, seq, paths, logger=logger)
        append_fasta_record(paths.final_fasta, header, trimmed_seq)
        write_trim_record_block(paths.report, header, trim_msg)
        i += 1
        gp.update_trim(i)

    gp.finish()
    _finalize(paths, args.keep_temp, logger)


# run_oric_only()
def run_oric_only(args, logger=None):
    paths = resolve_paths(args.input, args.out_dir)
    shutil.copyfile(args.input, paths.final_fasta)
    if logger:
        logger.debug(f"Copied input → {paths.final_fasta}")

    n = _count_records(paths.final_fasta)
    gp = GlobalProgress(n, weights=(0.0, 0.40, 0.40, 0.20), label="origami")
    gp.start()

    run_blast_annotation(paths.final_fasta, args.dnaA_db, paths.blast_dnaa, evalue="1e-5", quiet=True)
    gp.mark_dnaa_done()

    run_blast_annotation(paths.final_fasta, args.oric_db, paths.blast_oric, evalue="1e-5", quiet=True)
    gp.mark_oric_done()

    dnaA_hits_all = parse_blast_annotation(paths.blast_dnaa)
    oriC_hits_all = parse_blast_annotation(paths.blast_oric)

    write_section_header(paths.report, "OriC Analysis")

    rotated_records = []
    j = 0
    for header, seq in stream_fasta(paths.final_fasta):
        rid = _record_id(header)
        top_dnaa = select_top_hits(dnaA_hits_all.get(rid, []), k=3)
        top_oric = select_top_hits(oriC_hits_all.get(rid, []), k=10)
        pairs, chosen = pair_and_choose(rid, seq, top_dnaa, top_oric)

        rotated = rotate_sequence(seq, chosen["start_coord"])
        rotated_records.append((header, rotated))
        write_oric_record_block(paths.report, header, top_dnaa, top_oric, pairs, chosen)
        j += 1
        gp.update_rotate(j)

    # Rewrite FASTA
    if os.path.exists(paths.final_fasta):
        os.remove(paths.final_fasta)
    for header, seq in rotated_records:
        append_fasta_record(paths.final_fasta, header, seq)

    gp.finish()
    _finalize(paths, args.keep_temp, logger)

def ensure_input_in_outdir(input_fasta, out_dir, logger=None):
    """
    Ensure input FASTA exists inside out_dir.
    If not, copy it there and return the new path.
    """
    input_fasta = os.path.abspath(input_fasta)
    out_dir = os.path.abspath(out_dir)

    os.makedirs(out_dir, exist_ok=True)

    dest_path = os.path.join(out_dir, os.path.basename(input_fasta))

    if input_fasta != dest_path:
        if not os.path.exists(dest_path):
            msg = f"[LOG] Copying input FASTA to output directory: {dest_path}"
            print(msg)
            if logger:
                logger.info(msg)

            shutil.copy2(input_fasta, dest_path)
        else:
            print(f"[LOG] Input FASTA already exists in output dir: {dest_path}")

        return dest_path

    return input_fasta

def cleanup_unwanted_outputs(out_dir, logger=None):
    file_patterns = [
        "*_genomic_final_circos_v9.svg",
        "*_genomic_final_circos_v9.png",
        "dotplot_highlighted_circos_logic.*",
        "origami_*_genomic_rotated_circular_gcskew.*",
        "coords.txt",
        "dotplot_output.delta",
    ]

    dir_patterns = [
        "linear_synteny"
    ]

    # Remove matching files
    for pattern in file_patterns:
        for fpath in glob.glob(os.path.join(out_dir, pattern)):
            if os.path.isfile(fpath):
                os.remove(fpath)
                msg = f"[LOG] Removed file: {fpath}"
                print(msg)
                if logger:
                    logger.info(msg)

    # Remove matching directories
    for pattern in dir_patterns:
        for dpath in glob.glob(os.path.join(out_dir, pattern)):
            if os.path.isdir(dpath):
                shutil.rmtree(dpath)
                msg = f"[LOG] Removed directory: {dpath}"
                print(msg)
                if logger:
                    logger.info(msg)


def run_full(args, logger=None):

    args.input = ensure_input_in_outdir(args.input, args.out_dir, logger)

    paths = resolve_paths(args.input, args.out_dir)
    write_section_header(paths.report, "Trim Analysis")

    n = _count_records(args.input)
    gp = GlobalProgress(n, weights=(0.40, 0.20, 0.20, 0.20), label="origami")
    gp.start()

    # Trim + Extend
    adjusted_records = []
    for header, seq in stream_fasta(args.input):
        trimmed_seq, trim_msg = detect_terminal_overlap(header, seq, paths, logger=logger)
        append_fasta_record(paths.original_fasta, header, trimmed_seq)

        extended_seq = prepare_chromosome_sequence(header, trimmed_seq)
        append_fasta_record(paths.final_fasta, header, extended_seq)

        write_trim_record_block(paths.report, header, trim_msg)
        gp.update_trim(1)

    # BLAST on extended sequence only
    run_blast_annotation(paths.final_fasta, args.dnaA_db, paths.blast_dnaa_ext, evalue="1e-5", quiet=True)
    run_blast_annotation(paths.final_fasta, args.oric_db, paths.blast_oric_ext, evalue="1e-5", quiet=True)

    dnaA_hits_ext = parse_blast_annotation(paths.blast_dnaa_ext)
    oriC_hits_ext = parse_blast_annotation(paths.blast_oric_ext)

    # Extended + Original reports (from extended hits)
    for header, seq in stream_fasta(paths.final_fasta):
        rid = _record_id(header)
        top_dnaa_ext = select_top_hits(dnaA_hits_ext.get(rid, []), k=3)
        top_oric_ext = select_top_hits(oriC_hits_ext.get(rid, []), k=7)

        # Pairing using extended sequence
        pairs, chosen, seq_adjusted = pair_and_choose(rid, seq, top_dnaa_ext, top_oric_ext)

        adjusted_records.append((header, seq_adjusted, pairs, chosen))

        # Build chosen oriC dict
        chosen_oric_dict = None
        if chosen.get("oriC_id") and chosen.get("oriC_coords"):
            qstart, qend = chosen["oriC_coords"]
            wrap_info = chosen.get("wrap_info")
            chosen_oric_dict = {
                "sseqid": chosen["oriC_id"],
                "qstart": qstart,
                "qend": qend,
                "wrap_info": wrap_info
            }

        # Build chosen dnaA dict
        chosen_dnaa_dict = None
        if chosen.get("dnaA_id") and chosen.get("dnaA_coords"):
            sA, eA = chosen["dnaA_coords"]
            wrap_info = chosen.get("dnaa_wrap_info")
            chosen_dnaa_dict = {
                "sseqid": chosen["dnaA_id"],
                "qstart": sA,
                "qend": eA,
                "bitscore": chosen.get("dnaA_bitscore", "NA"),
                "evalue": chosen.get("dnaA_evalue", "NA"),
                "wrap_info": wrap_info
            }

        # Write extended report
        write_extended_record_block(paths.extended_report, header,
                                    chosen_dnaa=chosen_dnaa_dict,
                                    chosen_oric=chosen_oric_dict)

        # Write original report (coordinates mapped back using wrap_info)
        genome_len = len(seq) - 4000  # original genome length
        write_original_record_block(paths.original_report, header,
                                    chosen_dnaa=chosen_dnaa_dict,
                                    chosen_oric=chosen_oric_dict,
                                    genome_len=genome_len)

    #  Overwrite final FASTA with adjusted sequences
    if os.path.exists(paths.final_fasta):
        os.remove(paths.final_fasta)
    for header, seq_adjusted, _, _ in adjusted_records:
        append_fasta_record(paths.final_fasta, header, seq_adjusted)

    # Rotation + final OriC report
    write_section_header(paths.report, "OriC Analysis", append=True)
    rotated_records = []

    for header, seq_adjusted, pairs, chosen in adjusted_records:
        rid = _record_id(header)
        top_dnaa = select_top_hits(dnaA_hits_ext.get(rid, []), k=3)
        top_oric = select_top_hits(oriC_hits_ext.get(rid, []), k=7)

        # Always use adjusted coords (they're identical to originals in normal cases)
        dnaa_coords = chosen.get("dnaA_coords_adj")
        oric_coords = chosen.get("oriC_coords_adj")

        if not dnaa_coords and not oric_coords:
            print(f"[WARN] No dnaA or oriC found for {header}, skipping rotation.")
            continue

        dnaa_start, dnaa_end = dnaa_coords if dnaa_coords else (None, None)
        oric_start, oric_end = oric_coords if oric_coords else (None, None)

        genome_len = len(seq_adjusted)

        # ---- CASE: both markers present ----
        if dnaa_start is not None and oric_start is not None:

            dnaa_start = int(dnaa_start)
            oric_start = int(oric_start)

            # Compute circular distances
            dist_dnaa_to_oric = (oric_start - dnaa_start) % genome_len
            dist_oric_to_dnaa = (dnaa_start - oric_start) % genome_len

            # Choose anchor giving minimal separation
            if dist_dnaa_to_oric < dist_oric_to_dnaa:
                # dnaA → oriC forward distance is smaller
                start_coord = dnaa_start
                chosen_feature = "dnaA"
            else:
                # oriC → dnaA forward distance is smaller
                start_coord = oric_start
                chosen_feature = "oriC"

        # ---- CASE: only dnaA exists ----
        elif dnaa_start is not None:
            start_coord = int(dnaa_start)
            chosen_feature = "dnaA"

        # ---- CASE: only oriC exists ----
        elif oric_start is not None:
            start_coord = int(oric_start)
            chosen_feature = "oriC"

        # ---- Perform rotation ----
        rotated = rotate_sequence(seq_adjusted, start_coord - 1)
        rotated_records.append((header, rotated))
        print(f"Rotated {header} using minimal-distance anchor: {chosen_feature} at {start_coord}")

        write_oric_record_block(paths.report, header, top_dnaa, top_oric, pairs, chosen)


    # Rewrite final.fasta with rotated sequences
    if os.path.exists(paths.final_fasta):
        os.remove(paths.final_fasta)
    for header, seq in rotated_records:
        append_fasta_record(paths.final_fasta, header, seq)
    
    if args.keep_plot:
        print("[LOG] --keep-plot enabled → generating multi-sequence plots...")
        generate_all_plots(
            fasta_output_path=paths.original_fasta,
            report_path=paths.original_report,
            rotated_fasta_path=paths.final_fasta,
            work_dir=args.out_dir,
            gcf_fasta_path=args.input
        )
    else:
        print("[LOG] --keep-plot NOT enabled → skipping plot generation.")

    gp.finish()

    cleanup_unwanted_outputs(args.out_dir, logger)

    _finalize(paths, args.keep_temp, logger)

