from pathlib import Path
import os
os.chdir('/dcs04/lieber/marmaypag/spatialDLPFC_SCZ_LIBD4100/')
from pyhere import here
import json

from loopy.sample import Sample
from loopy.utils.utils import remove_dupes, Url

import glob
files = glob.glob(os.path.join('processed-data', 'image_processing', 'samui', '*.tif'))

#spot_diameter_m = 55e-6 # 5-micrometer diameter for Visium spot
#m_per_px = 4.971263040387764e-07
    
for file in files[5:]:    
    sample_id = os.path.basename(file)[:-4]
    samui_dir = Path(here('processed-data', 'image_processing', 'samui', sample_id))
    samui_dir.mkdir(parents = True, exist_ok = True)
    img_channels = ['DAPI', 'NeuN', 'Claudin5', 'WFA', 'Autofluorescence', 'binaryWFA']
    default_channels = {'blue': 'DAPI', 'yellow': 'NeuN', 'green': 'Claudin5', 'red': 'WFA', 'cyan': 'Autofluorescence', 'white': 'binaryWFA'}
    img_path = here('processed-data', 'image_processing', 'samui', sample_id+'.tif')
    this_sample = Sample(name = sample_id, path = samui_dir)
    this_sample.add_image(tiff = img_path, channels = img_channels, scale = m_per_px, defaultChannels = default_channels)
    this_sample.write()
    print(sample_id)