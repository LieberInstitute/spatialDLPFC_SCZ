# first -- module load loopy/1.0.0-next.24

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

# file paths
img_channels = ['Lipofuscin', 'Claudin5', 'DAPI', 'NeuN', 'WFA']
                # 'segmented_DAPI',
    # 'segmented_Claudin', 'segmented_NeuN', 'segmented_WFA', 'segmented_Lipofuscin']
default_channels = {'blue': 'DAPI', 'red': 'WFA'}

master_excel_path = here('raw-data', 'experiment_info', 'VisiumSPG_PNN_Master.xlsx')
img_path = here('processed-data', 'Samui', 'section_053_A1V12F14-053_A1.tif')
json_path = here('processed-data', 'spaceranger', 'V12F14-053_A1', 'outs', 'spatial', 'scalefactors_json.json')
tissue_positions_path = here('processed-data', 'spaceranger', 'V12F14-053_A1', 'outs', 'spatial', 'tissue_positions.csv')
spot_counts_path = here('processed-data', 'spaceranger', 'V12F14-053_A1', 'outs', 'spatial', 'tissue_spot_counts_correct_counts.csv')


# out_dir = here('processed-data', 'Samui', 'section_053_A1')
out_dir = '/dcs04/lieber/marmaypag/spatialDLPFC_SCZ_LIBD4100/processed-data/Samui/section_053_A1_/'
