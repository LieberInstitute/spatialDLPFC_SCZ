# first module load loopy/1.0.0-next.24

from pathlib import Path
from pyhere import here
import json
import os
import scanpy as sc

import numpy as np
import pandas as pd
from rasterio import Affine

from loopy.sample import Sample
from loopy.utils.utils import remove_dupes, Url

spot_diameter_m = 55e-6 # 55-micrometer diameter for Visium spot
img_channels = [
    'DAPI', 'Claudin-5', 'NeuN', 'WFA', 'Lipofuscin', 'segmented_DAPI',
    'segmented_Claudin', 'segmented_NeuN', 'segmented_WFA'
]
default_channels = {'blue': 'DAPI', 'red': 'WFA'}

# default_gene = 'SNAP25' # we don't have one yet

#   Names of continuous features expected to be columns in the observation
#   data (colData) of the AnnData
# spe_cont_features = ['PpTau', 'PAbeta']

#   Diagnosis by brain number (not included in the sample_info sheet)
# sample_dx = {
#     'Br3854': 'AD', 'Br3873': 'AD', 'Br3880': 'AD', 'Br3874': 'control'
# }

