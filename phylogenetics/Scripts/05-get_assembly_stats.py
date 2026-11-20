from pathlib import Path
from Bio import SeqIO
import pandas as pd
import re

BASE = Path("getorganelle/all_runs")

rows = []

for sample_dir in BASE.glob("*_getorganelle"):

    sample_id = sample_dir.name.replace("_getorganelle", "")

    fasta = sample_dir / "extended_spades" / "K115" / "final_contigs.fasta"
    logfile = sample_dir / "get_org.log.txt"

    row = {
        "sample_id": sample_id,
        "fasta_present": fasta.exists(),
        "n_contigs": None,
        "total_length": None,
        "largest_contig": None,
        "N50": None,
        "GC_percent": None,
        "circular_path": False
    }

    # Check if fasta exists
    if fasta.exists():

        records = list(SeqIO.parse(fasta, "fasta"))

        lengths = [len(r.seq) for r in records]

        if lengths:

            total_length = sum(lengths)
            largest = max(lengths)

            # N50 calculation
            lengths_sorted = sorted(lengths, reverse=True)

            cumulative = 0
            half = total_length / 2

            for l in lengths_sorted:
                cumulative += l
                if cumulative >= half:
                    n50 = l
                    break

            # GC calculation
            seq_concat = "".join(str(r.seq).upper() for r in records)

            gc = (
                (seq_concat.count("G") + seq_concat.count("C"))
                / len(seq_concat)
            ) * 100

            row.update({
                "n_contigs": len(records),
                "total_length": total_length,
                "largest_contig": largest,
                "N50": n50,
                "GC_percent": round(gc, 2)
            })

    # Look for evidence of circularisation
    if logfile.exists():

        with open(logfile) as f:
            logtext = f.read().lower()

            if "circular" in logtext:
                row["circular_path"] = True

    rows.append(row)

df = pd.DataFrame(rows)

df = df.sort_values("sample_id")

df.to_csv("getorganelle_assembly_stats.tsv", sep="\t", index=False)

print(df.head())
print(f"\nProcessed {len(df)} samples")