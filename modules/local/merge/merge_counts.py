#!/usr/bin/env python3
import pandas as pd
import glob
import os
import sys

# Argumente
input_folder = sys.argv[1]
output_file = sys.argv[2]

# Alle TSV-Dateien finden
files = glob.glob(os.path.join(input_folder, "*.tsv"))

if not files:
    sys.exit(f"Keine TSV-Dateien in '{input_folder}' gefunden.")

dfs = []
for f in files:
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
