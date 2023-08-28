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
