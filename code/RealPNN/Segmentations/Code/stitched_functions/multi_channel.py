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
image_DAPI = np.array(Image.open('/dcs04/lieber/marmaypag/spatialDLPFC_SCZ_LIBD4100/processed-data/RealPNN/capture_area_segmentations/DAPI/DAPI_binarized/V12F14-053_A1_dapi_binarized.tif'))
image_claudin = np.array(Image.open('/dcs04/lieber/marmaypag/spatialDLPFC_SCZ_LIBD4100/processed-data/RealPNN/capture_area_segmentations/Claudin/claudin_binarized/V12F14-053_A1_claudin_binarized.tif'))
image_NeuN = np.array(Image.open('/dcs04/lieber/marmaypag/spatialDLPFC_SCZ_LIBD4100/processed-data/RealPNN/capture_area_segmentations/NeuN/NeuN_binarized/V12F14-053_A1_neun_binarized.tif'))
image_WFA = np.array(Image.open('/dcs04/lieber/marmaypag/spatialDLPFC_SCZ_LIBD4100/processed-data/RealPNN/capture_area_segmentations/WFA/WFA_binarized/V12F14-053_A1_wfa_binarized.tif'))
image_AF = np.array(Image.open('/dcs04/lieber/marmaypag/spatialDLPFC_SCZ_LIBD4100/processed-data/RealPNN/capture_area_segmentations/AF/AF_segmented_binary/V12F14-053_A1_af_contours_segmented.tif'))

# Combine the images to create a multi-channel image
multi_channel_image = np.stack((image_DAPI, image_claudin, image_NeuN, image_WFA, image_AF), axis = 0) # image_NeuN, image_WFA], axis=0
# multi_channel_image = np.concatenate((image_DAPI[..., np.newaxis],
#                                        image_claudin[..., np.newaxis],
#                                        image_NeuN[..., np.newaxis],
#                                        image_WFA[..., np.newaxis],
#                                       image_AF[..., np.newaxis]), axis=-1)


# Create a new numpy array with shape (x, y, 5)
new_image = np.zeros((image_DAPI.shape[0], image_DAPI.shape[1], 5), dtype=np.uint8)

# Assign existing images to the first four channels of the new numpy array
new_image[:, :, 0] = image_DAPI
new_image[:, :, 1] = image_claudin
new_image[:, :, 2] = image_NeuN
new_image[:, :, 3] = image_WFA


# Set metadata for the channels
channel_names = ['DAPI', 'Claudin-5', 'NeuN', 'WFA', 'Channel 5']
metadata = {'axes': 'ZYX', 'channel_names': channel_names}

# Save the new image with metadata using tifffile module
tifffile.imwrite('new_image.tif', new_image, imagej=True, metadata=metadata)
