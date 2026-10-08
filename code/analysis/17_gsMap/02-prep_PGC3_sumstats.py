"""Prepare PGC3 SCZ (European, Trubetskoy et al. 2022) summary statistics for
`gsmap format_sumstats`.

The raw PGC3 file has `##` header lines and reports case/control counts per
SNP. Here we:
  * drop the `##` header lines,
  * compute the per-SNP effective sample size, N_eff = 4 / (1/NCAS + 1/NCON),
    explicitly rather than relying on the PGC3 `NEFF` column, whose
    definition varies across PGC releases (some report half of N_eff),
  * rename to the column names passed to `gsmap format_sumstats`.

BETA in PGC3 is log(OR) for the A1 allele; FCON (A1 frequency in controls) is
used as the allele frequency for MAF filtering.

Usage: python 02-prep_PGC3_sumstats.py <raw_pgc3.tsv.gz> <out.tsv.gz>
"""

import sys

import pandas as pd

raw_file, out_file = sys.argv[1], sys.argv[2]

required_col = ["ID", "A1", "A2", "BETA", "SE", "PVAL", "NCAS", "NCON", "IMPINFO", "FCON"]

# `##` meta lines are skipped by comment="#"; the column header line
# ("CHROM ID POS ...") has no '#', and data lines contain no '#'.
gwas = pd.read_csv(raw_file, sep="\t", comment="#", low_memory=False)

missing_col = [c for c in required_col if c not in gwas.columns]
if missing_col:
    sys.exit(f"Missing expected PGC3 columns: {missing_col}; found {list(gwas.columns)}")

gwas = gwas.assign(N=4 / (1 / gwas["NCAS"] + 1 / gwas["NCON"]))

out = gwas.rename(
    columns={
        "ID": "SNP",
        "PVAL": "P",
        "IMPINFO": "INFO",
        "FCON": "FRQ",
    }
)[["SNP", "A1", "A2", "BETA", "SE", "P", "N", "INFO", "FRQ"]]

# error prevention
assert out["SNP"].str.startswith("rs").mean() > 0.9, "SNP column is not mostly rsIDs"
assert out[["BETA", "SE", "P", "N"]].notna().all().all(), "NA in BETA/SE/P/N"

print(f"Read {len(gwas):,} SNPs; median N_eff = {out['N'].median():,.0f}")
out.to_csv(out_file, sep="\t", index=False, compression="gzip")
print(f"Wrote {out_file}")
