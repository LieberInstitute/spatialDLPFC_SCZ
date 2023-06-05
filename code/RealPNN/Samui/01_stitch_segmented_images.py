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
import tifffile
import matplotlib.pyplot as plt

raw_img_A1 = pyhere.here('processed-data', 'VistoSeg', 'captureAreas', 'V12F14-053_A1.tif')
seg_img_A1 = pyhere.here('processed-data', 'RealPNN', 'all_channels_segemented', 'Test2', 'V12F14-053_A1.tif')
dst_dir = pyhere.here('processed-data', 'Samui', 'V12F14-053_A1')

# find the right shape parameter
Image.MAX_IMAGE_PIXELS = None
raw_A1 = Image.open(raw_img_A1)
seg_A1 = Image.open(seg_img_A1)

# raw image added to new image
new_im_raw = np.zeros((5, (np.array(Image.open(raw_img_A1))).shape[0], (np.array(Image.open(raw_img_A1))).shape[1]),  dtype=np.uint8)
for ch_num in range(5):
    raw_A1.seek(ch_num)
    new_im_raw[ch_num] = raw_A1

# segmented image added to another new image
new_im_seg = np.zeros((5, (np.array(Image.open(seg_img_A1))).shape[0], (np.array(Image.open(seg_img_A1))).shape[1]),  dtype=np.uint8)
for ch_num in range(5):
    seg_A1.seek(ch_num)
    new_im_seg[ch_num] = seg_A1

samui_img = np.concatenate((new_im_raw, new_im_seg), axis = 0)

tifffile.imwrite(dst_dir + os.path.basename(seg_img_A1), samui_img)
