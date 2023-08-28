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
# img_path = here('processed-data', 'Samui', 'section_053_A1V12F14-053_A1.tif')
img_path = '/dcs04/lieber/marmaypag/spatialDLPFC_SCZ_LIBD4100/processed-data/VistoSeg/captureAreas/V13M06-279_A1.tif'
# json_path = here('processed-data', 'spaceranger', 'V12F14-053_A1', 'outs', 'spatial', 'scalefactors_json.json')
# tissue_positions_path = here('processed-data', 'spaceranger', 'V12F14-053_A1', 'outs', 'spatial', 'tissue_positions.csv')
# spot_counts_path = here('processed-data', 'spaceranger', 'V12F14-053_A1', 'outs', 'spatial', 'tissue_spot_counts_correct_counts.csv')


# out_dir = here('processed-data', 'Samui', 'section_053_A1')
out_dir = '/dcs04/lieber/marmaypag/spatialDLPFC_SCZ_LIBD4100/processed-data/Samui/section_V13M06-279_A1/'

# Read in sample info, subset to relevant columns, and clean
sample_info = (pd.read_excel(master_excel_path)
    .query('`Will be Sequenced? ` == "Yes"')
    .filter(["BrNumbr", "Slide #", "Array #", "Sample #"])
    #   Clean up column names
    .rename(
        columns = {
            "BrNumbr": "br_num",
            "Slide #": "sample_id",
            "Array #": "array_num",
            "Sample #": "sample_num"
        }
    )
)

# extract the experiment number based on the sample number
sample_info['experiment_num'] = (((sample_info['sample_num'].str.extract(r'(\d+)').astype(int)) -1) // 8 + 1).astype(str)

# create the data structure needed for Samui using the image data
test_sample = Sample(name = 'V12F14-053_A1', path = out_dir)

# add the IF image for test sample
test_sample.add_image(tiff = img_path, channels = img_channels, scale = m_per_px,
    defaultChannels = default_channels)

# write to the output directory
test_sample.write()
