'''
For Stitched Visium-IF tissue sections from VistoSeg SplitSlide output
Channel0 = AF
Channel1 = Claudin - 5 (Alex 488),
Channel2 = DAPI,
Channel3 = NeuN,
Channel4 = WFA
'''

from __future__ import print_function
from skimage.feature import peak_local_max
from skimage.segmentation import find_boundaries, watershed
from scipy import ndimage
import argparse
from argparse import ArgumentParser
import imutils
import numpy as np
import pyhere
from pyhere import here
from pylab import xticks
from pathlib import Path
import pandas as pd
import PIL
from PIL import Image, ImageFont, ImageDraw, ImageEnhance, ImageSequence
import os
import matplotlib
import matplotlib.pyplot as plt
import sys
import cv2
import math
import scipy
from scipy.spatial.distance import *
import skimage
from skimage import *
from skimage import feature, segmentation, draw, measure, morphology
from skimage.filters import threshold_otsu
from skimage.morphology import (erosion,dilation,opening,closing,white_tophat,black_tophat,skeletonize,convex_hull_image)
from skimage.draw import polygon_perimeter
from itertools import product
from collections import defaultdict
from shapely.geometry import Point
from shapely.geometry.polygon import Polygon


import numpy as np
import tifffile

# Load the images
Image.MAX_IMAGE_PIXELS = None
image_DAPI = np.array(Image.open(pyhere.here('processed-data', 'RealPNN', 'capture_area_segmentations', 'DAPI', 'DAPI_binarized', 'V12F14-053_A1_dapi_binarized.tif'))) #'/dcs04/lieber/marmaypag/spatialDLPFC_SCZ_LIBD4100/processed-data/RealPNN/capture_area_segmentations/DAPI/DAPI_binarized/V12F14-053_A1_dapi_binarized.tif'))
image_claudin = np.array(Image.open('/dcs04/lieber/marmaypag/spatialDLPFC_SCZ_LIBD4100/processed-data/RealPNN/capture_area_segmentations/Claudin/claudin_binarized/V12F14-053_A1_claudin_binarized.tif'))
image_NeuN = np.array(Image.open('/dcs04/lieber/marmaypag/spatialDLPFC_SCZ_LIBD4100/processed-data/RealPNN/capture_area_segmentations/NeuN/NeuN_binarized/V12F14-053_A1_neun_binarized.tif'))
image_WFA = np.array(Image.open('/dcs04/lieber/marmaypag/spatialDLPFC_SCZ_LIBD4100/processed-data/RealPNN/capture_area_segmentations/WFA/WFA_binarized/V12F14-053_A1_wfa_binarized.tif'))
image_AF = np.array(Image.open('/dcs04/lieber/marmaypag/spatialDLPFC_SCZ_LIBD4100/processed-data/RealPNN/capture_area_segmentations/AF/AF_segmented_binary/V12F14-053_A1_af_contours_segmented.tif'))

new_image = np.zeros((5, image_DAPI.shape[0], image_DAPI.shape[1]), dtype=np.uint8)

# Assign existing images to the first four channels of the new numpy array
new_image[0] = image_DAPI
new_image[1] = image_claudin
new_image[2] = image_NeuN
new_image[3] = image_WFA
new_image[4] = image_AF


# Set metadata for the channels
channel_names = ['DAPI', 'Claudin-5', 'NeuN', 'WFA', 'AF']
metadata = {'axes': 'YXZ', 'channel_names': channel_names} #ZYX right ones

# Save the new image with metadata using tifffile module
tifffile.imwrite('/dcs04/lieber/marmaypag/spatialDLPFC_SCZ_LIBD4100/processed-data/RealPNN/capture_area_segmentations/all_channels_segemented/Test/A1_new_image_yxz.tif', new_image, imagej=True, metadata=metadata)

# automate this for all the images in the directory
Image.MAX_IMAGE_PIXELS = None
dapi_seg_dir = pyhere.here('processed-data', 'RealPNN', 'capture_area_segmentations', 'DAPI', 'DAPI_binarized') # '/dcs04/lieber/marmaypag/spatialDLPFC_SCZ_LIBD4100/processed-data/RealPNN/capture_area_segmentations/DAPI/DAPI_binarized/'
neun_seg_dir = pyhere.here('processed-data', 'RealPNN', 'capture_area_segmentations', 'NeuN', 'NeuN_binarized') # '/dcs04/lieber/marmaypag/spatialDLPFC_SCZ_LIBD4100/processed-data/RealPNN/capture_area_segmentations/NeuN/NeuN_binarized/'
claudin_seg_dir = pyhere.here('processed-data', 'RealPNN', 'capture_area_segmentations','Claudin','claudin_binarized') # '/dcs04/lieber/marmaypag/spatialDLPFC_SCZ_LIBD4100/processed-data/RealPNN/capture_area_segmentations/Claudin/claudin_binarized/'
wfa_seg_dir = pyhere.here('processed-data', 'RealPNN', 'capture_area_segmentations','WFA', 'WFA_binarized') # '/dcs04/lieber/marmaypag/spatialDLPFC_SCZ_LIBD4100/processed-data/RealPNN/capture_area_segmentations/WFA/WFA_binarized/'
af_seg_dir = pyhere.here('processed-data', 'RealPNN', 'capture_area_segmentations','AF','AF_segmented_binary') # '/dcs04/lieber/marmaypag/spatialDLPFC_SCZ_LIBD4100/processed-data/RealPNN/capture_area_segmentations/AF/AF_segmented_binary/'
stitched_dir = pyhere.here('processed-data', 'RealPNN', 'capture_area_segmentations', 'all_channels_segemented')


dir_list = [dapi_seg_dir, claudin_seg_dir, neun_seg_dir, wfa_seg_dir, af_seg_dir]
for dir in dir_list:
    print("Reading image from: ", dir.split(',')[0].split('_')[0])
    for img_path in os.listdir(dapi_seg_dir):
        img = np.array(Image.open(os.path.join(dapi_seg_dir, img_path)))
        print("Stitching image:", img_path.split('.')[0])
        new_image = np.zeros((5, img.shape[0], img.shape[1]), dtype=np.uint8)
        new_image[0] = image_DAPI
        new_image[1] = image_claudin
        new_image[2] = image_NeuN
        new_image[3] = image_WFA
        new_image[4] = image_AF
    tifffile.imwrite(stitched_dir + '')


import pathlib

dapi_seg_dir = pathlib.Path(pyhere.here('processed-data', 'RealPNN', 'capture_area_segmentations', 'DAPI', 'DAPI_binarized'))
neun_seg_dir = pathlib.Path(pyhere.here('processed-data', 'RealPNN', 'capture_area_segmentations', 'NeuN', 'NeuN_binarized'))
claudin_seg_dir = pathlib.Path(pyhere.here('processed-data', 'RealPNN', 'capture_area_segmentations','Claudin','claudin_binarized'))
wfa_seg_dir = pathlib.Path(pyhere.here('processed-data', 'RealPNN', 'capture_area_segmentations','WFA', 'WFA_binarized'))
af_seg_dir = pathlib.Path(pyhere.here('processed-data', 'RealPNN', 'capture_area_segmentations','AF','AF_segmented_binary'))
stitched_dir = pathlib.Path(pyhere.here('processed-data', 'RealPNN', 'capture_area_segmentations', 'all_channels_segemented'))

dir_list = [dapi_seg_dir, claudin_seg_dir, neun_seg_dir, wfa_seg_dir, af_seg_dir]
for dir in dir_list:
    print("Reading image from:", dir.name)
    for file in dir.glob('*'):  # iterate over all files in the directory
        if file.suffix.lower() in ['.tif']:
            print("  -", file.name)  # display the filename
            image = Image.open(file)
            # img = np.array(Image.open(os.path.join(dapi_seg_dir, img_path)))
            print("Stitching image:", file)
            new_image = np.zeros((5, img.shape[0], img.shape[1]), dtype=np.uint8)
            new_image[0] = image_DAPI
            new_image[1] = image_claudin
            new_image[2] = image_NeuN
            new_image[3] = image_WFA
            new_image[4] = image_AF


# will open 1 image from each directory
for files in zip(*[dir.glob('*') for dir in dir_list]):
    print("Reading images from:", os.path.basename(dir))
    for file in files:
        if file.suffix.lower() in ['.tif']:
            print("  -", file.name)
            # print("Reading image:", file.name)
            # image = Image.open(file)

for dir, files in zip(dir_list, (dir.glob('*') for dir in dir_list)):
    print("Reading images from:", os.path.basename(dir))
    for file in files:
        if file.suffix.lower() in ['.tif']:
            print("  -", file.name)




# testing for channels last to match countnuclei
# Create a new numpy array with shape (x, y, 5)
new_image = np.zeros((image_DAPI.shape[0], image_DAPI.shape[1], 5), dtype=np.uint8)

# Assign existing images to the first four channels of the new numpy array
new_image[..., 0] = image_DAPI
new_image[..., 1] = image_claudin
new_image[..., 2] = image_NeuN
new_image[..., 3] = image_WFA
new_image[..., 4] = image_AF


# Set metadata for the channels
channel_names = ['DAPI', 'Claudin-5', 'NeuN', 'WFA', 'AF']
metadata = {'axes': 'YXC', 'channel_names': channel_names}

# Move the axis to have channels as the last dimension
new_image = np.moveaxis(new_image, 2, -1)

# Save the new image with metadata using tifffile module
tifffile.imwrite('/dcs04/lieber/marmaypag/spatialDLPFC_SCZ_LIBD4100/processed-data/RealPNN/capture_area_segmentations/all_channels_segemented/Test/A1_new_image_yxc_.tif', new_image, metadata=metadata)


###trying to find the right format
import cv2

# Load the four tiff images
images = [cv2.imread("image1.tif"), cv2.imread("image2.tif"), cv2.imread("image3.tif"), cv2.imread("image4.tif")]

# Create a new multichannel image
multichannel_image = cv2.merge([images[0], images[1], images[2], images[3]])

# Save the multichannel image
cv2.imwrite("multichannel_image.tif", multichannel_image)







OMJGSRJH


from PIL import Image
import numpy as np
import pathlib
import pyhere
import os


# Define the directories
dapi_seg_dir = pathlib.Path(pyhere.here('processed-data', 'RealPNN', 'capture_area_segmentations', 'DAPI', 'DAPI_binarized'))
neun_seg_dir = pathlib.Path(pyhere.here('processed-data', 'RealPNN', 'capture_area_segmentations', 'NeuN', 'NeuN_binarized'))
claudin_seg_dir = pathlib.Path(pyhere.here('processed-data', 'RealPNN', 'capture_area_segmentations','Claudin','claudin_binarized'))
wfa_seg_dir = pathlib.Path(pyhere.here('processed-data', 'RealPNN', 'capture_area_segmentations','WFA', 'WFA_binarized'))
af_seg_dir = pathlib.Path(pyhere.here('processed-data', 'RealPNN', 'capture_area_segmentations','AF','AF_segmented_binary'))
stitched_dir = pathlib.Path(pyhere.here('processed-data', 'RealPNN', 'capture_area_segmentations', 'all_channels_segemented'))

# Create a list of directories -- ChatGPT code
dir_list = [dapi_seg_dir, neun_seg_dir, claudin_seg_dir, wfa_seg_dir, af_seg_dir]

# Create a list to store the first images from each directory
first_images = []

# Open the first image from each directory and store them in the list
images = []
for directory in dir_list:
    # get the list of image files in the directory that end with ".tif"
    image_files = [f for f in directory.glob('*') if f.is_file() and f.name.endswith('.tif')]
    print(os.path.basename(image_files))
    # open the first image in the directory
    if image_files:
        image = Image.open(image_files[0])
        images.append(image)

########################################


dapi_057_A1 = np.array(Image.open('/dcs04/lieber/marmaypag/spatialDLPFC_SCZ_LIBD4100/processed-data/RealPNN/capture_area_segmentations/DAPI/DAPI_binarized/V12F14-057_A1_dapi_binarized.tif'))
dapi_057_B1 = np.array(Image.open('/dcs04/lieber/marmaypag/spatialDLPFC_SCZ_LIBD4100/processed-data/RealPNN/capture_area_segmentations/DAPI/DAPI_binarized/V12F14-057_B1_dapi_binarized.tif'))
dapi_057_C1 = np.array(Image.open('/dcs04/lieber/marmaypag/spatialDLPFC_SCZ_LIBD4100/processed-data/RealPNN/capture_area_segmentations/DAPI/DAPI_binarized/V12F14-057_C1_dapi_binarized.tif'))
dapi_057_D1 =  np.array(Image.open('/dcs04/lieber/marmaypag/spatialDLPFC_SCZ_LIBD4100/processed-data/RealPNN/capture_area_segmentations/DAPI/DAPI_binarized/V12F14-057_D1_dapi_binarized.tif'))

cla_057_A1 = np.array(Image.open('/dcs04/lieber/marmaypag/spatialDLPFC_SCZ_LIBD4100/processed-data/RealPNN/capture_area_segmentations/Claudin/claudin_binarized/V12F14-057_A1_claudin_binarized.tif'))
cla_057_B1 = np.array(Image.open('/dcs04/lieber/marmaypag/spatialDLPFC_SCZ_LIBD4100/processed-data/RealPNN/capture_area_segmentations/Claudin/claudin_binarized/V12F14-057_B1_claudin_binarized.tif'))
cla_057_C1 = np.array(Image.open('/dcs04/lieber/marmaypag/spatialDLPFC_SCZ_LIBD4100/processed-data/RealPNN/capture_area_segmentations/Claudin/claudin_binarized/V12F14-057_C1_claudin_binarized.tif'))
cla_057_D1 = np.array(Image.open('/dcs04/lieber/marmaypag/spatialDLPFC_SCZ_LIBD4100/processed-data/RealPNN/capture_area_segmentations/Claudin/claudin_binarized/V12F14-057_D1_claudin_binarized.tif'))

neun_057_A1 = np.array(Image.open('/dcs04/lieber/marmaypag/spatialDLPFC_SCZ_LIBD4100/processed-data/RealPNN/capture_area_segmentations/NeuN/NeuN_binarized/V12F14-057_A1_neun_binarized.tif'))
neun_057_B1 = np.array(Image.open('/dcs04/lieber/marmaypag/spatialDLPFC_SCZ_LIBD4100/processed-data/RealPNN/capture_area_segmentations/NeuN/NeuN_binarized/V12F14-057_B1_neun_binarized.tif'))
neun_057_C1 = np.array(Image.open('/dcs04/lieber/marmaypag/spatialDLPFC_SCZ_LIBD4100/processed-data/RealPNN/capture_area_segmentations/NeuN/NeuN_binarized/V12F14-057_C1_neun_binarized.tif'))
neun_057_D1 = np.array(Image.open('/dcs04/lieber/marmaypag/spatialDLPFC_SCZ_LIBD4100/processed-data/RealPNN/capture_area_segmentations/NeuN/NeuN_binarized/V12F14-057_D1_neun_binarized.tif'))

wfa_057_A1 = np.array(Image.open('/dcs04/lieber/marmaypag/spatialDLPFC_SCZ_LIBD4100/processed-data/RealPNN/capture_area_segmentations/WFA/WFA_binarized/V12F14-057_A1_wfa_binarized.tif'))
wfa_057_B1 = np.array(Image.open('/dcs04/lieber/marmaypag/spatialDLPFC_SCZ_LIBD4100/processed-data/RealPNN/capture_area_segmentations/WFA/WFA_binarized/V12F14-057_B1_wfa_binarized.tif'))
wfa_057_C1 = np.array(Image.open('/dcs04/lieber/marmaypag/spatialDLPFC_SCZ_LIBD4100/processed-data/RealPNN/capture_area_segmentations/WFA/WFA_binarized/V12F14-057_C1_wfa_binarized.tif'))
wfa_057_D1 = np.array(Image.open('/dcs04/lieber/marmaypag/spatialDLPFC_SCZ_LIBD4100/processed-data/RealPNN/capture_area_segmentations/WFA/WFA_binarized/V12F14-057_D1_wfa_binarized.tif'))

af_057_A1 = np.array(Image.open('/dcs04/lieber/marmaypag/spatialDLPFC_SCZ_LIBD4100/processed-data/RealPNN/capture_area_segmentations/AF/AF_segmented_binary/V12F14-057_A1_af_contours_segmented.tif'))
af_057_B1 = np.array(Image.open('/dcs04/lieber/marmaypag/spatialDLPFC_SCZ_LIBD4100/processed-data/RealPNN/capture_area_segmentations/AF/AF_segmented_binary/V12F14-057_B1_af_contours_segmented.tif'))
af_057_C1 = np.array(Image.open('/dcs04/lieber/marmaypag/spatialDLPFC_SCZ_LIBD4100/processed-data/RealPNN/capture_area_segmentations/AF/AF_segmented_binary/V12F14-057_C1_af_contours_segmented.tif'))
af_057_D1 = np.array(Image.open('/dcs04/lieber/marmaypag/spatialDLPFC_SCZ_LIBD4100/processed-data/RealPNN/capture_area_segmentations/AF/AF_segmented_binary/V12F14-057_D1_af_contours_segmented.tif'))

###########################################

dapi_053_A1 = np.array(Image.open('/dcs04/lieber/marmaypag/spatialDLPFC_SCZ_LIBD4100/processed-data/RealPNN/capture_area_segmentations/DAPI/DAPI_binarized/V12F14-053_A1_dapi_binarized.tif'))
dapi_053_B1 = np.array(Image.open('/dcs04/lieber/marmaypag/spatialDLPFC_SCZ_LIBD4100/processed-data/RealPNN/capture_area_segmentations/DAPI/DAPI_binarized/V12F14-053_B1_dapi_binarized.tif'))
dapi_053_C1 = np.array(Image.open('/dcs04/lieber/marmaypag/spatialDLPFC_SCZ_LIBD4100/processed-data/RealPNN/capture_area_segmentations/DAPI/DAPI_binarized/V12F14-053_C1_dapi_binarized.tif'))
dapi_053_D1 = np.array(Image.open('/dcs04/lieber/marmaypag/spatialDLPFC_SCZ_LIBD4100/processed-data/RealPNN/capture_area_segmentations/DAPI/DAPI_binarized/V12F14-053_D1_dapi_binarized.tif'))

cla_053_A1 = np.array(Image.open('/dcs04/lieber/marmaypag/spatialDLPFC_SCZ_LIBD4100/processed-data/RealPNN/capture_area_segmentations/Claudin/claudin_binarized/V12F14-053_A1_claudin_binarized.tif'))
cla_053_B1 = np.array(Image.open('/dcs04/lieber/marmaypag/spatialDLPFC_SCZ_LIBD4100/processed-data/RealPNN/capture_area_segmentations/Claudin/claudin_binarized/V12F14-053_B1_claudin_binarized.tif'))
cla_053_C1 = np.array(Image.open('/dcs04/lieber/marmaypag/spatialDLPFC_SCZ_LIBD4100/processed-data/RealPNN/capture_area_segmentations/Claudin/claudin_binarized/V12F14-053_C1_claudin_binarized.tif'))
cla_053_D1 = np.array(Image.open('/dcs04/lieber/marmaypag/spatialDLPFC_SCZ_LIBD4100/processed-data/RealPNN/capture_area_segmentations/Claudin/claudin_binarized/V12F14-053_D1_claudin_binarized.tif'))

neun_053_A1 = np.array(Image.open('/dcs04/lieber/marmaypag/spatialDLPFC_SCZ_LIBD4100/processed-data/RealPNN/capture_area_segmentations/NeuN/NeuN_binarized/V12F14-053_A1_neun_binarized.tif'))
neun_053_B1 = np.array(Image.open('/dcs04/lieber/marmaypag/spatialDLPFC_SCZ_LIBD4100/processed-data/RealPNN/capture_area_segmentations/NeuN/NeuN_binarized/V12F14-053_B1_neun_binarized.tif'))
neun_053_C1 = np.array(Image.open('/dcs04/lieber/marmaypag/spatialDLPFC_SCZ_LIBD4100/processed-data/RealPNN/capture_area_segmentations/NeuN/NeuN_binarized/V12F14-053_C1_neun_binarized.tif'))
neun_053_D1 = np.array(Image.open('/dcs04/lieber/marmaypag/spatialDLPFC_SCZ_LIBD4100/processed-data/RealPNN/capture_area_segmentations/NeuN/NeuN_binarized/V12F14-053_D1_neun_binarized.tif'))

wfa_053_A1 = np.array(Image.open('/dcs04/lieber/marmaypag/spatialDLPFC_SCZ_LIBD4100/processed-data/RealPNN/capture_area_segmentations/WFA/WFA_binarized/V12F14-053_A1_wfa_binarized.tif'))
wfa_053_B1 = np.array(Image.open('/dcs04/lieber/marmaypag/spatialDLPFC_SCZ_LIBD4100/processed-data/RealPNN/capture_area_segmentations/WFA/WFA_binarized/V12F14-053_B1_wfa_binarized.tif'))
wfa_053_C1 = np.array(Image.open('/dcs04/lieber/marmaypag/spatialDLPFC_SCZ_LIBD4100/processed-data/RealPNN/capture_area_segmentations/WFA/WFA_binarized/V12F14-053_C1_wfa_binarized.tif'))
wfa_053_D1 = np.array(Image.open('/dcs04/lieber/marmaypag/spatialDLPFC_SCZ_LIBD4100/processed-data/RealPNN/capture_area_segmentations/WFA/WFA_binarized/V12F14-053_D1_wfa_binarized.tif'))

af_053_A1 = np.array(Image.open('/dcs04/lieber/marmaypag/spatialDLPFC_SCZ_LIBD4100/processed-data/RealPNN/capture_area_segmentations/AF/AF_segmented_binary/V12F14-053_A1_af_contours_segmented.tif'))
af_053_B1 = np.array(Image.open('/dcs04/lieber/marmaypag/spatialDLPFC_SCZ_LIBD4100/processed-data/RealPNN/capture_area_segmentations/AF/AF_segmented_binary/V12F14-053_B1_af_contours_segmented.tif'))
af_053_C1 = np.array(Image.open('/dcs04/lieber/marmaypag/spatialDLPFC_SCZ_LIBD4100/processed-data/RealPNN/capture_area_segmentations/AF/AF_segmented_binary/V12F14-053_C1_af_contours_segmented.tif'))
af_053_D1 = np.array(Image.open('/dcs04/lieber/marmaypag/spatialDLPFC_SCZ_LIBD4100/processed-data/RealPNN/capture_area_segmentations/AF/AF_segmented_binary/V12F14-053_D1_af_contours_segmented.tif'))


##### Works!!! ----re-run this again to change the wfa channel images

new_image = np.zeros((5, dapi_057_A1.shape[0], dapi_057_A1.shape[1]), dtype=np.uint8)

# Assign existing images to the first four channels of the new numpy array
new_image[0] = dapi_053_D1
new_image[1] = cla_053_D1
new_image[2] = neun_053_D1
new_image[3] = wfa_053_D1
new_image[4] = af_053_D1


# Set metadata for the channels
channel_names = ['DAPI', 'Claudin-5', 'NeuN', 'WFA', 'AF']
metadata = {'axes': 'YXZ', 'channel_names': channel_names} #ZYX right ones


# Save the new image with metadata using tifffile module
tifffile.imwrite('/dcs04/lieber/marmaypag/spatialDLPFC_SCZ_LIBD4100/processed-data/RealPNN/all_channels_segemented/Test2/V12F14-053_D1.tif', new_image, metadata=metadata)


### code to do this for multiple images in onr slide
