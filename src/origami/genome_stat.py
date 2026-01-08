from Bio import SeqIO

def compute_genome_stats(fasta_path: str):
    genome_size = 0
    contig_count = 0
    plasmid_count = 0
    gc_count = 0
    total_bases = 0

    for record in SeqIO.parse(fasta_path, "fasta"):
        contig_count += 1
        seq = record.seq.upper()
        length = len(seq)
        
        genome_size += length
        total_bases += length
        gc_count += seq.count("G") + seq.count("C")

        # detect plasmid by name
        if "plasmid" in record.description.lower():
            plasmid_count += 1

    gc_percent = round((gc_count / total_bases * 100), 2) if total_bases > 0 else 0

    return {
        "genome_size_bp": genome_size,
        "contig_count": contig_count,
        "plasmid_count": plasmid_count,
        "gc_percent": gc_percent
    }
