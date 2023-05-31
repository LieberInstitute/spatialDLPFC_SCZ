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

sample_info_path = here(
    'raw-data', 'experiment_info', 'VisiumSPG_PNN_Master.xlsx'
) # /dcs04/lieber/marmaypag/spatialDLPFC_SCZ_LIBD4100/raw-data/experiment_info/VisiumSPG_PNN_Master.xlsx

# spe_path = here('processed-data', '16_samui', 'spe.h5ad')
# notes_path = str(Path(here('code', '16_samui', 'feature_notes.md')).resolve())



img_path = here('processed-data', 'RealPNN', 'all_channels_segemented', 'Test2', 'V12F14-053_A1.tif') #/dcs04/lieber/marmaypag/spatialDLPFC_SCZ_LIBD4100/processed-data/RealPNN/all_channels_segemented/Test2/V12F14-053_A1.tif
json_path = here(
    'processed-data', 'spaceranger', 'V12F14-053_A1', 'outs', 'spatial',
    'scalefactors_json.json'
) # /dcs04/lieber/marmaypag/spatialDLPFC_SCZ_LIBD4100/processed-data/spaceranger/V12F14-053_A1/outs/spatial/scalefactors_json.json

out_dir = here('processed-data', 'Samui', 'V12F14-053_A1') # /dcs04/lieber/marmaypag/spatialDLPFC_SCZ_LIBD4100/processed-data/Samui/V12F14-053_A1/


################################################################################
#   Read in sample info and clean
################################################################################

#   Read in sample info, subset to relevant columns, and clean
sample_info = (pd.read_excel(sample_info_path)
    .query('`Sequenced? ` == "Yes"')
    .filter(["Br####", "Slide SN #", "Array #", "Sample #"])
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

#   Prepend "Br" tp brain number and make a string
sample_info['br_num'] = "Br" + sample_info['br_num'].astype(str)

#   Add diagnosis using brain number
sample_info['diagnosis'] = sample_info['br_num'].replace(sample_dx)

#   Fix the experiment number column (use strings of integers)
sample_info['experiment_num'] = ((sample_info['sample_num'] - 1) // 4  + 1).astype(str)

#   Different forms of sample IDs appear to be used for spaceranger outputs
#   and raw images
sample_info = (sample_info
    .assign(
        spaceranger_id = sample_info['sample_id'].transform(lambda x: x.replace('-', '')) +
            '_' + sample_info['array_num'] + '_' + sample_info['br_num'],
        image_id = 'VIFAD' + sample_info['experiment_num'] + '_' + sample_info['sample_id'] + '_' + sample_info['array_num']
    )
)

