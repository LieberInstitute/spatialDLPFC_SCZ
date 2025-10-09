import os
os.chdir('/dcs04/lieber/marmaypag/spatialDLPFC_SCZ_LIBD4100/')

from pyhere import here
from pathlib import Path
import session_info

import numpy as np
import pandas as pd
import json
import sys
from loopy.sample import Sample
import tifffile
from PIL import Image
import re
import matplotlib.pyplot as plt

import scanpy as sc
from rasterio import Affine
from loopy.utils.utils import remove_dupes, Url
import re
from scipy.spatial import KDTree
import glob


this_donor = sys.argv[1:]
this_donor = ''.join(this_donor)

spot_diameter_m = 55e-6 # 5-micrometer diameter for Visium spot
m_per_px = 4.20390369911423e-07

################################################################################
#   Gather gene-expression data into a DataFrame to later as a feature
################################################################################
spg_path = here("processed-data", "samui", "spg_kept_only.h5ad")
spg = sc.read(spg_path)

#path_groups = spg.obs['path_groups'].cat.categories
spgP = spg[spg.obs['BrNumbr'] == this_donor, :]
unique_sample_id = spgP.obs['sample_id'].unique()[0]
samui_dir = Path(here('processed-data', 'samui', f"{this_donor}_{unique_sample_id}"))
samui_dir.mkdir(parents = True, exist_ok = True)

spgP.obs.index = spgP.obs.index.str.replace('_'+unique_sample_id , '')
#   Convert the sparse gene-expression matrix to pandas DataFrame, with the
#   gene symbols as column names
gene_df = pd.DataFrame(
    spgP.X.toarray(),
    index = spgP.obs.index,
    columns = spgP.var['gene_name']
)
#   Some gene symbols are actually duplicated. Just take the first column in
#   any duplicated cases
gene_df = gene_df.loc[: , ~gene_df.columns.duplicated()].copy()
gene_df.index.name = None

precast_columns = spgP.obs.filter(like="PRECAST")
precast_columns = spgP.obs[["spd_label"]].join(precast_columns)
precast_df = pd.DataFrame(precast_columns)

sample_df = spgP.obs[['sample_id']].copy()
spotCalling_df = spgP.obs[['pnn_pos', 'neuropil_pos', 'neun_pos', 'vasc_pos']]
#spotCalling_metrics = spgP.obs[['spg_NDAPI', 'spg_PDAPI', 'spg_IDAPI', 'spg_CNDAPI',
#       'spg_NNeuN', 'spg_PNeuN', 'spg_INeuN', 'spg_CNNeuN', 'spg_NWFA', 'spg_PWFA', 'spg_IWFA', 'spg_CNWFA',
#       'spg_NClaudin5', 'spg_PClaudin5', 'spg_IClaudin5']]

################################################################################
#   Use the Samui API to create the importable directory for this combined "sample"
################################################################################
img_channels = ['DAPI', 'NeuN', 'WFA', 'Claudin5', 'segDAPI', 'segNeuN', 'segWFA', 'segClaudin5']
default_channels = {'blue': 'DAPI', 'green': 'NeuN', 'yellow': 'Claudin5', 'red': 'WFA', 'white':'segDAPI', 'white':'segNeuN', 'white':'segWFA', 'white':'segClaudin5'}
img_path = Path(here('processed-data', 'samui',f'{unique_sample_id}.tif'))
tissue_positions_path = Path(here("processed-data", "spaceranger", unique_sample_id, "outs", "spatial", "tissue_positions.csv"))
tissue_positions = pd.read_csv(tissue_positions_path ,index_col = 0).rename({'pxl_row_in_fullres': 'y', 'pxl_col_in_fullres': 'x'},axis = 1)
tissue_positions.index.name = None
tissue_positions = tissue_positions[['x', 'y']].astype(int)

default_gene = 'SNAP25'
assert default_gene in gene_df.columns, "Default gene not in AnnData"

#notes_md_url = Url('/dcs04/lieber/lcolladotor/spatialHPC_LIBD4035/spatial_hpc/code/VSPG_image_stitching/feature_notes.md')
#this_sample = Sample(name = samui_dir.name, path = samui_dir, notesMd = notes_md_url)
this_sample = Sample(name = samui_dir.name, path = samui_dir)
this_sample.add_coords(tissue_positions, name = "coords", mPerPx = m_per_px, size = spot_diameter_m)
this_sample.add_image(tiff = img_path,channels = img_channels,defaultChannels = default_channels, scale = m_per_px)

this_sample.add_csv_feature(precast_df, name = "Domains", coordName = "coords", dataType = "categorical")
this_sample.add_csv_feature(spotCalling_df, name = "Spot_Calling", coordName = "coords", dataType = "categorical")
#this_sample.add_csv_feature(spotCalling_metrics, name = "Spot_Calling_metrics", coordName = "coords", dataType = "quantitative")
this_sample.add_chunked_feature(gene_df, name = "Genes", coordName = "coords", dataType = "quantitative")
this_sample.set_default_feature(group = "Genes", feature = default_gene)

this_sample.write()


session_info.show()
