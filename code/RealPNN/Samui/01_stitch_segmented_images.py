'''
For Stitched Visium-IF tissue sections from VistoSeg SplitSlide output
Channel0 = AF
Channel1 = Claudin - 5 (Alex 488),
Channel2 = DAPI,
Channel3 = NeuN,
Channel4 = WFA
'''

from PIL import Image
import numpy as np
import pathlib
import pyhere
import os
import pyhere
import cv2
import matplotlib.pyplot as plt

raw_img_A1 = pyhere.here('processed-data', 'VistoSeg', 'captureAreas', 'V12F14-057_A1.tif')
seg_img_A1 = pyhere.here('processed-data', 'RealPNN', 'all_channels_segemented', 'Test2', 'V12F14-057_A1.tif')

# find the right shape parameter
Image.MAX_IMAGE_PIXELS = None
raw_A1_arr = np.array(Image.open(raw_img_A1))
seg_A1_arr = np.array(Image.open(seg_img_A1))
new_image = np.zeros((10, raw_A1_arr.shape[0], raw_A1_arr.shape[1]), dtype=np.uint8)

# seek the channels in a loop for raw and segmented image
raw_A1 = Image.open(raw_img_A1)
seg_A1 = Image.open(seg_img_A1)
for ch_num_raw in range(5):
    for ch_num_new in range(10):
        print("first",ch_num_raw)
        print("second",ch_num_new)
    # raw_A1.seek(ch_num)
    # raw_A1_ch0 = np.array(raw_A1, dtype = 'uint8')
    # print("second",ch_num)
    # new_image[ch_num] = raw_A1_ch0

# each image has 5 channels, there are 2 images, the new image has 10 channels
# seek one channel, convert to array, save it in the new image
new_im_raw = np.zeros((5, (np.array(Image.open(raw_img_A1))).shape[0], (np.array(Image.open(raw_img_A1))).shape[1]),  dtype=np.uint8)
for ch_num in range(5):
    print("first", ch_num)
    raw_A1.seek(ch_num)
    print("second", ch_num)
    new_im_raw[ch_num] = raw_A1
    print("last", ch_num)

tifffile.imwrite('/dcs04/lieber/marmaypag/spatialDLPFC_SCZ_LIBD4100/processed-data/Samui/V12F14-053_A1/new_im1.tif', new_im_raw)

new_im_seg = np.zeros((5, (np.array(Image.open(raw_img_A1))).shape[0], (np.array(Image.open(raw_img_A1))).shape[1]),  dtype=np.uint8)
for ch_num in range(5):
    print("first", ch_num)
    raw_A1.seek(ch_num)
    print("second", ch_num)
    new_im_raw[ch_num] = raw_A1
    print("last", ch_num)

tifffile.imwrite('/dcs04/lieber/marmaypag/spatialDLPFC_SCZ_LIBD4100/processed-data/Samui/V12F14-053_A1/new_im2.tif', new_im_seg)

out_path = '/dcs04/lieber/marmaypag/spatialDLPFC_SCZ_LIBD4100/processed-data/Samui/V12F14-053_A1/'

full_img = np.concatenate((new_im_raw, new_im_seg), axis = 0)
with tifffile.TiffWriter(out_path) as tiff:
    tiff.write(full_img)

fig,ax = plt.subplots(figsize = (20,20))
ax.imshow(new_im, cmap='gray')
fig.show()


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

