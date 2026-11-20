import pandas as pd
import re

# =========================
# Load assembly stats table
# =========================

df = pd.read_csv(
    "getorganelle_assembly_stats.tsv",
    sep="\t"
)

# ==========================================
# Extract metadata from sample_id
# ==========================================

# Example sample_ids:
# LIMS317701_38-FL-318
# LIMS317701_40-s-13
# LIMS32034_x_11-FL-62

pattern = re.compile(
    r'^(?P<batch>LIMS\d+(?:_x)?)_'
    r'(?P<sample_no>\d+)-'
    r'(?P<lifestyle>FL|s)-'
    r'(?P<strain>\d+)$'
)

# Extract into temporary dataframe
meta = df["sample_id"].str.extract(pattern)

# Create strain_no column
meta["strain_no"] = meta.apply(
    lambda row:
        row["strain"]
        if row["lifestyle"] == "FL"
        else f's{row["strain"]}',
    axis=1
)

# Add extracted columns back into dataframe
df["batch"] = meta["batch"]
df["sample_no"] = meta["sample_no"]
df["lifestyle"] = meta["lifestyle"]
df["strain_no"] = meta["strain_no"]

# ==========================================
# Merge pond information
# ==========================================

pond_df = pd.read_csv(
    "pond_info.tsv",
    sep="\t"
)

# Merge on strain_no
df = df.merge(
    pond_df,
    on="strain_no",
    how="left"
)

# ==========================================
# Save updated table
# ==========================================

df.to_csv(
    "getorganelle_assembly_stats_with_metadata.tsv",
    sep="\t",
    index=False
)

print(df.head())
print(f"\nFinal table contains {len(df)} rows")