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
img_channels = ['Lipofuscin', 'Claudin5', 'DAPI', 'NeuN', 'WFA'
                'segmented_DAPI',
    'segmented_Claudin', 'segmented_NeuN', 'segmented_WFA', 'segmented_Lipofuscin']
default_channels = {'blue': 'DAPI', 'red': 'WFA'}

master_excel_path = here('raw-data', 'experiment_info', 'VisiumSPG_PNN_Master.xlsx')
img_path = here('processed-data', 'Samui', 'section_053_A1V12F14-053_A1.tif')
json_path = here('processed-data', 'spaceranger', 'V12F14-053_A1', 'outs', 'spatial', 'scalefactors_json.json')
tissue_positions_path = here('processed-data', 'spaceranger', 'V12F14-053_A1', 'outs', 'spatial', 'tissue_positions.csv')
spot_counts_path = here('processed-data', 'spaceranger', 'V12F14-053_A1', 'outs', 'spatial', 'tissue_spot_counts_correct_counts.csv')


# out_dir = here('processed-data', 'Samui', 'section_053_A1')
out_dir = '/dcs04/lieber/marmaypag/spatialDLPFC_SCZ_LIBD4100/processed-data/Samui/section_053_A1_/'

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

# Different forms of sample IDs appear to be used for spaceranger outputs and raw images
sample_info = (sample_info
    .assign(
        spaceranger_id = sample_info['sample_id'].apply(lambda x: x.replace('-', '')) +
            '_' + sample_info['array_num'].astype(str) + '_' + sample_info['br_num'].astype(str),
        image_id = 'VIFAD' + sample_info['sample_num'].str.extract(r'(\d+)').astype(str).squeeze() + '_' + sample_info['sample_id'] + '_' + sample_info['array_num'].astype(str)
    )
)

# read and pre-process the spaceranger json file
with open(json_path, 'r') as f:
    spaceranger_json = json.load(f)

m_per_px = spot_diameter_m / spaceranger_json['spot_diameter_fullres']

# create the data structure needed for Samui using the image data
test_sample = Sample(name = 'V12F14-053_A1', path = out_dir)

# get the tissue positions file from spaceranger and extract the coordinates
tissue_positions = pd.read_csv(tissue_positions_path,
    header = 1,
    names = ["in_tissue", "row", "col", "y", "x"],
    index_col = 0)
tissue_positions.index.name = "barcode"
tissue_positions.index = tissue_positions.index.astype(str)
tissue_positions.index.is_unique

# read the spot counts/positions file #"imagerow", "imagecol",
columns_to_read = ["barcode",  "NAF", "PAF", "CNAF", "NClaudin5", "PClaudin5", "CNClaudin5", "NDAPI", "PDAPI", "CNDAPI", "NNeuN", "PNeuN", "CNNeuN", "NWFA", "PWFA", "CNWFA"]
spot_positions = pd.read_csv(spot_counts_path, #"y", "x",
                             # header = 1,
                             # names = ["barcode",  "NAF", "PAF", "NClaudin5", "PClaudin5", "NDAPI", "PDAPI", "NNeuN", "PNeuN", "NWFA", "PWFA"],
                             usecols=columns_to_read)

spot_positions.set_index("barcode", inplace=True) # Set 'barcode' as the index column
spot_positions.index = spot_positions.index.astype(str)
spot_positions.index.is_unique  # Verify if the index is unique


# add the tissue positions for test sample
test_sample.add_coords(
    tissue_positions, name="coords", mPerPx=m_per_px, size=spot_diameter_m)

# add csv feature for coordinates
test_sample.add_csv_feature(
    tissue_positions, name = "Tissue position",
    coordName = "coords", dataType = "quantitative")

# add csv feature for spot coverage data
test_sample.add_csv_feature(spot_positions, name = "Spot coverage",
                            coordName = "coords", dataType = "quantitative")

# add the IF image for test sample
test_sample.add_image(tiff = img_path, channels = img_channels, scale = m_per_px,
    defaultChannels = default_channels)

# write to the output directory
test_sample.write()


