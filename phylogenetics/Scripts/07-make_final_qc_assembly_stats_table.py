import pandas as pd
import re

# ==========================================
# Load existing master table
# ==========================================

master = pd.read_csv(
    "getorganelle_assembly_stats_with_metadata.tsv",
    sep="\t"
)

# ==========================================
# Load MultiQC general stats
# ==========================================

qc = pd.read_csv(
    "qc/multiqc/multiqc_data/multiqc_general_stats.txt",
    sep="\t"
)

# ==========================================
# Parse sample names
# ==========================================

pattern = re.compile(
    r'^illumina_'
    r'(?P<batch>LIMS\d+(?:_x)?)_'
    r'Sample_'
    r'(?P<sample_no>\d+)-'
    r'(?P<lifestyle>FL|s)-'
    r'(?P<strain>\d+)_.*_'
    r'(?P<read>R[12])_001$'
)

parsed = qc["Sample"].str.extract(pattern)

# Construct sample_id matching master table
qc["sample_id"] = (
    parsed["batch"] + "_"
    + parsed["sample_no"] + "-"
    + parsed["lifestyle"] + "-"
    + parsed["strain"]
)

qc["read"] = parsed["read"]

# ==========================================
# Keep only useful QC columns
# ==========================================

qc = qc[[
    "sample_id",
    "read",
    "FastQC_mqc-generalstats-fastqc-percent_duplicates",
    "FastQC_mqc-generalstats-fastqc-percent_gc",
    "FastQC_mqc-generalstats-fastqc-avg_sequence_length",
    "FastQC_mqc-generalstats-fastqc-percent_fails",
    "FastQC_mqc-generalstats-fastqc-total_sequences"
]]

# ==========================================
# Rename columns to simpler names
# ==========================================

qc = qc.rename(columns={
    "FastQC_mqc-generalstats-fastqc-percent_duplicates": "percent_duplicates",
    "FastQC_mqc-generalstats-fastqc-percent_gc": "percent_gc",
    "FastQC_mqc-generalstats-fastqc-avg_sequence_length": "avg_sequence_length",
    "FastQC_mqc-generalstats-fastqc-percent_fails": "percent_fails",
    "FastQC_mqc-generalstats-fastqc-total_sequences": "total_sequences"
})

# ==========================================
# Pivot R1/R2 into separate columns
# ==========================================

dupes = qc.groupby(["sample_id", "read"]).size()

print(dupes[dupes > 1])

qc_wide = qc.pivot_table(
    index="sample_id",
    columns="read",
    aggfunc="first"
)

# Flatten column names
qc_wide.columns = [
    f"{metric}_{read}"
    for metric, read in qc_wide.columns
]

qc_wide = qc_wide.reset_index()

# ==========================================
# Merge with master table
# ==========================================

final = master.merge(
    qc_wide,
    on="sample_id",
    how="left"
)

# ==========================================
# Save final table
# ==========================================

final.to_csv(
    "master_assembly_qc_table.tsv",
    sep="\t",
    index=False
)

print(final.head())
print(f"\nFinal merged table contains {len(final)} rows")
