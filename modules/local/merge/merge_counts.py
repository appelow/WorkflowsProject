#!/usr/bin/env python3
import pandas as pd
import os
import sys

# Input arguments: all except the last are input files, last is output file
if len(sys.argv) < 3:
    sys.exit("Usage: merge_counts.py <file1> <file2> ... <output_file>")

files = sys.argv[1:-1]
output_file = sys.argv[-1]

dfs = []
for f in files:
    if not os.path.exists(f):
        sys.exit(f"❌ Datei nicht gefunden: {f}")
    sample = os.path.basename(f).replace(".featureCounts.tsv", "")
    df = pd.read_csv(f, sep="\t", comment="#")
    if "Geneid" not in df.columns:
        print(f"⚠️ Datei {f} enthält keine Geneid-Spalte, wird übersprungen.")
        continue
    last_col = df.columns[-1]
    df = df[["Geneid", last_col]]
    df.columns = ["Geneid", sample]
    dfs.append(df)

if not dfs:
    sys.exit("❌ Keine gültigen featureCounts-Dateien gefunden.")

merged = dfs[0]
for df in dfs[1:]:
    merged = pd.merge(merged, df, on="Geneid", how="outer")

merged = merged.fillna(0)
merged = merged.sort_values("Geneid")

merged.to_csv(output_file, sep="\t", index=False)
print(f"✅ Merge abgeschlossen: {output_file} erstellt.")
