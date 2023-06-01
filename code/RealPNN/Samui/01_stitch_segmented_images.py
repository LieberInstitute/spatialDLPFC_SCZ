from PIL import Image
import numpy as np
import pathlib
import pyhere
import os
import pyhere

raw_img_A1 = pyhere.here('processed-data', 'VistoSeg', 'captureAreas', 'V12F14-057_A1.tif')
seg_img_A1 = pyhere.here('processed-data', 'RealPNN', 'all_channels_segemented', 'Test2', 'V12F14-057_A1.tif')

Image.MAX_IMAGE_PIXELS = None
raw_A1 = Image.open(raw_img_A1)
seg_A1 = Image.open(seg_img_A1)

# seek the channels in a loop for raw and segmented image

# find the right shape parameter
new_image = np.zeros((10, dapi_057_A1.shape[0], dapi_057_A1.shape[1]), dtype=np.uint8)

# write a loop for filling up all the channels of the new image
# Assign existing images to the first four channels of the new numpy array
new_image[0] = dapi_057_D1
new_image[1] = cla_057_D1
new_image[2] = neun_057_D1
new_image[3] = wfa_057_D1
new_image[4] = af_057_D1


# Set metadata for the channels -- change the name of the channels
channel_names = ['DAPI', 'Claudin-5', 'NeuN', 'WFA', 'AF']
metadata = {'axes': 'YXZ', 'channel_names': channel_names} #ZYX right ones

# create a new folder to store all the new stitched images
# Save the new image with metadata using tifffile module
tifffile.imwrite('/dcs04/lieber/marmaypag/spatialDLPFC_SCZ_LIBD4100/processed-data/RealPNN/all_channels_segemented/Test2/V12F14-057_D1.tif', new_image, metadata=metadata)

