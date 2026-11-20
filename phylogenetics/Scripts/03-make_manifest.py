import re
from pathlib import Path
from collections import defaultdict

files = list(Path(".").glob("*.fastq.gz"))

samples = defaultdict(dict)

pattern = re.compile(
    r"illumina_(LIMS\d+(?:_x)?)_Sample_([0-9]+-(?:FL|s)-\d+)_.*_(R[12])_001\.fastq\.gz"
)

for f in files:
    m = pattern.match(f.name)
    if not m:
        continue

    lims = m.group(1)
    strain = m.group(2)
    read = m.group(3)

    sample_id = f"{lims}_{strain}"
    samples[sample_id][read] = str(f.resolve())

with open("samples.tsv", "w") as out:
    out.write("sample_id\tR1\tR2\n")
    for sample_id, reads in sorted(samples.items()):
        r1 = reads.get("R1", "")
        r2 = reads.get("R2", "")
        if r1 and r2:
            out.write(f"{sample_id}\t{r1}\t{r2}\n")