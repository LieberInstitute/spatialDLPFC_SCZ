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
img_channels = ['DAPI', 'Claudin-5', 'NeuN', 'WFA', 'Lipofuscin', 'segmented_DAPI',
    'segmented_Claudin', 'segmented_NeuN', 'segmented_WFA', 'segmented_Lipofuscin']
default_channels = {'blue': 'DAPI', 'red': 'WFA'}

master_excel_path = here('raw-data', 'experiment_info', 'VisiumSPG_PNN_Master.xlsx')

img_path = here('processed-data', 'Samui', 'V12F14-053_A1_.tif')
json_path = here('processed-data', 'spaceranger', 'V12F14-053_A1', 'outs', 'spatial', 'scalefactors_json.json')

out_dir = here('processed-data', 'Samui')

#   Read in sample info, subset to relevant columns, and clean
sample_info = (pd.read_excel(master_excel_path)
    .query('`Will be Sequenced? ` == "Yes"')
    .filter(["BrNumbr", "Slide #", "Array #", "Sample #"])
    #   Clean up column names
    .rename(
        columns = {
            "Br####": "br_num",
            "Slide SN #": "sample_id",
            "Array #": "array_num",
            "Sample #": "sample_num"
        }
    )
)



