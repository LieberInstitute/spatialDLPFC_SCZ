from PIL import Image
import numpy as np
import pathlib
import pyhere
import os
import pyhere

raw_img_A1 = pyhere.here('processed-data', 'VistoSeg', 'captureAreas', 'V12F14-057_A1.tif')
seg_img_A1 = pyhere.here('processed-data', 'RealPNN', 'all_channels_segemented', 'Test2', 'V12F14-057_A1.tif')

raw_A1 = 

new_image = np.zeros((10, dapi_057_A1.shape[0], dapi_057_A1.shape[1]), dtype=np.uint8)

# Assign existing images to the first four channels of the new numpy array
new_image[0] = dapi_057_D1
new_image[1] = cla_057_D1
new_image[2] = neun_057_D1
new_image[3] = wfa_057_D1
new_image[4] = af_057_D1


# Set metadata for the channels
channel_names = ['DAPI', 'Claudin-5', 'NeuN', 'WFA', 'AF']
metadata = {'axes': 'YXZ', 'channel_names': channel_names} #ZYX right ones


# Save the new image with metadata using tifffile module
tifffile.imwrite('/dcs04/lieber/marmaypag/spatialDLPFC_SCZ_LIBD4100/processed-data/RealPNN/all_channels_segemented/Test2/V12F14-057_D1.tif', new_image, metadata=metadata)

