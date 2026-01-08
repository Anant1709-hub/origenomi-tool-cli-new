import os
from pathlib import Path

class Paths:
    def __init__(self, input_file: str, out_dir: str):
        self.prefix = Path(input_file).stem
        self.out_dir = out_dir

        self.temp_dir = os.path.join(out_dir, "temp_files")
        self.annotation_dir = os.path.join(self.temp_dir, "annotation")
        self.trim_dir = os.path.join(self.temp_dir, "trim")  
        os.makedirs(self.annotation_dir, exist_ok=True)
        os.makedirs(self.trim_dir, exist_ok=True)

        # original input fasta (copied for clarity)
        self.original_fasta = os.path.join(out_dir, f"origenomi_{self.prefix}_original.fna")

        # extended fasta after trim + extension
        self.final_fasta = os.path.join(out_dir, f"origenomi_{self.prefix}.fna")

        self.report = os.path.join(out_dir, f"origenomi_{self.prefix}_report.txt")              
        self.original_report = os.path.join(out_dir, f"origenomi_{self.prefix}_original_report.txt")
        self.extended_report = os.path.join(out_dir, f"origenomi_{self.prefix}_extended_report.txt")

        self.blast_dnaa_ext = os.path.join(self.annotation_dir, f"{self.prefix}_dnaa_extended.blast")
        self.blast_oric_ext = os.path.join(self.annotation_dir, f"{self.prefix}_oric_extended.blast")
        self.blast_dnaa_orig = os.path.join(self.annotation_dir, f"{self.prefix}_dnaa_original.blast")
        self.blast_oric_orig = os.path.join(self.annotation_dir, f"{self.prefix}_oric_original.blast")


def resolve_paths(input_file: str, out_dir: str) -> Paths:
    """
    Creates a Paths object with organized file structure for the pipeline.
    """
    return Paths(input_file, out_dir)
