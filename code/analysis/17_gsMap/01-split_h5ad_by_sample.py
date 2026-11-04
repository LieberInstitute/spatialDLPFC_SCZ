"""Split the existing full-dataset h5ad into one gsMap-ready h5ad per sample.

Source: processed-data/samui/spg_kept_only.h5ad, written from
fnl_spe_kept_spots_only.rds by code/samui/spe_to_anndata.R (zellkonverter;
X = raw `counts`, layers['logcounts']; obs_names = spe$key). Checked to match
the current spe: 279,806 spots x 36,601 genes, 63 samples.

gsMap treats each h5ad as one tissue section and builds spatial neighbourhoods
from its coordinates, so each Visium capture area (sample_id) must be its own
file. For each sample, writes {sample_id}.h5ad with:
  X, layers['count']  raw counts (spots x genes, CSR)
  var_names           gene symbols (gsMap matches symbols to its GTF); the
                      first occurrence is kept for duplicated symbols
  obs                 key, sample_id, brnum, dx, age, sex, spd_label
                      (obs_names = key, so results join back to spe)
  obsm['spatial']     spot coordinates
Also writes the sample lists used by the SLURM array jobs:
  sample_list.txt (line N = SLURM_ARRAY_TASK_ID N), sample_list_ntc.txt,
  sample_list_scz.txt, h5ad_list.txt, sample_dx.csv

Usage: python 01-split_h5ad_by_sample.py <source_h5ad> <h5ad_dir>
"""

import sys
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import scipy.sparse as sp

if len(sys.argv) != 3 or not all(sys.argv[1:]):
    sys.exit(f"Usage: python {sys.argv[0]} <source_h5ad> <h5ad_dir>; got {sys.argv[1:]}")

source_h5ad, h5ad_dir = Path(sys.argv[1]), Path(sys.argv[2])
if not source_h5ad.is_file():
    sys.exit(f"Source h5ad is not a file: '{source_h5ad}'")
h5ad_dir.mkdir(parents=True, exist_ok=True)

keep_col = ["key", "sample_id", "brnum", "dx", "age", "sex", "spd_label"]

adata_all = ad.read_h5ad(source_h5ad, backed="r")
print(adata_all)

# error prevention: source has what gsMap needs
assert adata_all.uns.get("X_name", "counts") == "counts", adata_all.uns.get("X_name")
assert all(c in adata_all.obs for c in keep_col), keep_col
assert "spatial" in adata_all.obsm and "gene_name" in adata_all.var
assert adata_all.obs["spd_label"].notna().all()
assert (adata_all.obs_names == adata_all.obs["key"].astype(str)).all()


# Gene symbols ----
gene_name = adata_all.var["gene_name"].astype(str)
keep_gene = (~gene_name.duplicated()).to_numpy()
print(f"Dropping {(~keep_gene).sum()} genes with duplicated symbols; keeping {keep_gene.sum()}")
gene_symbols = pd.Index(gene_name[keep_gene].to_numpy(), dtype=str)


# Sample lists ----
sample_dx = (
    adata_all.obs[["sample_id", "brnum", "dx"]]
    .astype(str)
    .drop_duplicates()
    .sort_values("sample_id")
)

# error prevention: 63 samples, one diagnosis per sample
assert sample_dx["sample_id"].is_unique and len(sample_dx) == 63, len(sample_dx)

sample_ids = sample_dx["sample_id"].tolist()
h5ad_files = [str(h5ad_dir / f"{s}.h5ad") for s in sample_ids]

(h5ad_dir / "sample_list.txt").write_text("\n".join(sample_ids) + "\n")
(h5ad_dir / "h5ad_list.txt").write_text("\n".join(h5ad_files) + "\n")
for dx in ["ntc", "scz"]:
    ids = sample_dx.loc[sample_dx["dx"] == dx, "sample_id"]
    (h5ad_dir / f"sample_list_{dx}.txt").write_text("\n".join(ids) + "\n")
sample_dx.to_csv(h5ad_dir / "sample_dx.csv", index=False)


# Write one h5ad per sample ----
# Row (spot) slices of a backed CSR matrix read only those rows from disk. If X
# is stored as CSC, each row slice would scan the whole matrix, so load it into
# memory once as CSR instead (~ a few GB for raw counts).
X_all = adata_all.X
if getattr(X_all, "format", "csr") != "csr":
    print(f"X is stored as {X_all.format}; loading into memory as CSR")
    X_all = sp.csr_matrix(X_all[:, :])

# zellkonverter stores obsm['spatial'] as a DataFrame (spatialCoords has column
# names); convert to a plain (n_spots x 2) array so it can be row-indexed.
spatial_all = adata_all.obsm["spatial"]
if isinstance(spatial_all, pd.DataFrame):
    print(f"obsm['spatial'] columns: {list(spatial_all.columns)}")
    spatial_all = spatial_all.to_numpy()
spatial_all = np.asarray(spatial_all, dtype=np.float64)

# error prevention
assert spatial_all.shape == (adata_all.n_obs, 2), spatial_all.shape
assert np.isfinite(spatial_all).all()

sample_of_spot = adata_all.obs["sample_id"].astype(str).to_numpy()
n_spots_written = 0

for sample_id, out_file in zip(sample_ids, h5ad_files):
    idx = np.flatnonzero(sample_of_spot == sample_id)

    # Read only this sample's rows of X (raw counts) from the backed file
    counts = sp.csr_matrix(X_all[idx, :], dtype=np.float32)[:, keep_gene]

    # error prevention: X holds raw integer counts, not log-normalized values
    assert np.array_equal(counts.data, np.round(counts.data)), f"{sample_id}: X is not integer counts"

    obs = adata_all.obs.iloc[idx][keep_col].copy()
    for col in ["key", "sample_id", "brnum", "dx", "sex"]:
        obs[col] = obs[col].astype(str)
    obs["spd_label"] = obs["spd_label"].astype(str).astype("category")
    obs.index = pd.Index(obs["key"].to_numpy(), dtype=str)

    adata = ad.AnnData(
        X=counts,
        obs=obs,
        var=pd.DataFrame(index=gene_symbols),
        layers={"count": counts.copy()},
        obsm={"spatial": spatial_all[idx]},
    )
    adata.write_h5ad(out_file, compression="gzip")
    n_spots_written += adata.n_obs
    print(f"{sample_id}: {adata.n_obs} spots x {adata.n_vars} genes -> {out_file}")

# error prevention: every spot written exactly once
assert n_spots_written == adata_all.n_obs, (n_spots_written, adata_all.n_obs)

print("Finished writing per-sample h5ad files for gsMap")
