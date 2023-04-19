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


dir_list = [dapi_seg_dir, claudin_seg_dir, neun_seg_dir, wfa_seg_dir, af_seg_dir]
for dir in dir_list:
    print("Reading image from: ", dir_list.split[','][0].split['_'][0])
    for img_path in os.listdir(dapi_seg_dir):
        img = np.array(Image.open(os.path.join(dapi_seg_dir, img_path)))
        print("Stitching image:", img_path.split('.')[0])
        new_image = np.zeros((5, img.shape[0], img.shape[1]), dtype=np.uint8)
        new_image[0] = image_DAPI
        new_image[1] = image_claudin
        new_image[2] = image_NeuN
        new_image[3] = image_WFA
        new_image[4] = image_AF















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
