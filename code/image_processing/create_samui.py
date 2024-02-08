from pathlib import Path
import os
os.chdir('/dcs04/lieber/lcolladotor/spatialHPC_LIBD4035/spatial_hpc/')
from pyhere import here
import json

import scanpy as sc

import numpy as np
import pandas as pd
from rasterio import Affine

from loopy.sample import Sample
from loopy.utils.utils import remove_dupes, Url


spot_diameter_m = 55e-6 # 5-micrometer diameter for Visium spot
img_channels = ['DAPI', 'Alexa_488', 'Alexa_555', 'Alexa_594', 'Alexa_647', 'Autofluorescence', 'binaryWFA']
default_channels = {'blue': 'DAPI', 'green': 'Alexa_488', 'yellow': 'Alexa_555', 'red': 'Alexa_594', 'magenta': 'Alexa647', 'cyan': 'Autofluorescence', 'white': 'binaryWFA'}

sample_id = ''
IMG_path = here('processed-data', 'image_processing', 'samui_test', sample_id+'.tif')

this_sample = Sample(name = sample_id, path = out_dir)
this_sample.add_image( tiff = img_path, channels = img_channels, scale = m_per_px, defaultChannels = default_channels)

this_sample.write()