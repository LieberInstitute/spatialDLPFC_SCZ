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
from skimage.morphology import (erosion,dilation,opening,closing,white_tophat,black_tophat,skeletonize,convex_hull_image)
from skimage.draw import polygon_perimeter
from itertools import product
from collections import defaultdict
from shapely.geometry import Point
from shapely.geometry.polygon import Polygon
from stitched_functions import read_img
from stitched_functions import watershed_segmentation
from stitched_functions import *

# directory path
source_dir = '/dcs04/lieber/marmaypag/spatialDLPFC_SCZ_LIBD4100/processed-data/VistoSeg/captureAreas/'
dst_dir = '/dcs04/lieber/marmaypag/spatialDLPFC_SCZ_LIBD4100/processed-data/RealPNN/capture_area_segmentations/DAPI/'



print("packages imported")
# file paths
img_C1 = pyhere.here('processed-data', 'VistoSeg', 'captureAreas','V12F14-057_C1.tif') # /dcs04/lieber/marmaypag/spatialDLPFC_SCZ_LIBD4100/processed-data/VistoSeg/captureAreas/V12F14-057_C1.tif
img_D1 = pyhere.here('processed-data', 'VistoSeg', 'captureAreas','V12F14-057_D1.tif') # /dcs04/lieber/marmaypag/spatialDLPFC_SCZ_LIBD4100/processed-data/VistoSeg/captureAreas/V12F14-057_D1.tif
print("filepaths imported")

# parser = argparse.ArgumentParser()
# parser.add_argument('--filepath', type=dir_path, required = True)
# args = parser.parse_args()

img_dapi, dapi_shifted, dapi_gray, dapi_thresh = read_img.read_and_preprocess(img_D1, 2)
# plot_imgs(img_dapi, dapi_thresh)
# fig,ax = plt.subplots(figsize = (20,20))
# ax.imshow(dapi_thresh, cmap = 'gray')
# fig.show()
print("dapi preprocessed")
dapi_labels, dapi_localmax = find_labels(dapi_thresh) #watershed_segmentation.
print("dapi segmented")
dpx, dpy, dpw, dph, dp_area, dapi_segmented = draw_rect_dapi(dapi_labels, dapi_gray, img_dapi)
print("dapi segments drawn and saved")
cv2.imwrite('/users/ukaipa/PNN/One_img/dapi_stitched_segmented_D1_1148569.tif', dapi_segmented)
print("segmented image saved")
# dapi_df = save_coordinates.create_df(dpx, dpy, dpw, dph, dp_area, img_dapi, 'dapi')
# print("dapi coordinates saved")

# detect contours on test image
Image.MAX_IMAGE_PIXELS = None
dapi_img = Image.open(img_A1)
dapi_img.seek(2)
dapi = np.array(dapi_img, dtype = 'uint8') # (17799, 16740)
dapi_c = cv2.cvtColor(dapi,cv2.COLOR_BGR2RGB)
gray = cv2.cvtColor(dapi_c,cv2.COLOR_RGB2GRAY)
contours,_ = cv2.findContours(thresh, cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_SIMPLE)
print(len(contours))
dp_cnt = cv2.drawContours(dapi_c, contours, -1, (0, 255, 0), 2)
fig,ax = plt.subplots(figsize = (20,20))
ax.imshow(dp_cnt) #
fig.show()

# find contours for all images in the dir
Image.MAX_IMAGE_PIXELS = None
for img_path in os.listdir(source_dir):
    if img_path.endswith(".tif"):
        dapi_img = Image.open(os.path.join(source_dir, img_path))
        dapi_img.seek(2)
        dapi = np.array(dapi_img, dtype = 'uint8')
        dapi_c = cv2.cvtColor(dapi,cv2.COLOR_BGR2RGB)
        gray = cv2.cvtColor(dapi_c,cv2.COLOR_RGB2GRAY)
        _,thresh = cv2.threshold(gray, np.mean(gray), 255, cv2.THRESH_BINARY | cv2.THRESH_OTSU) #_INV
        contours,_ = cv2.findContours(thresh, cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_SIMPLE)
        print("found", len(contours), "in", img_path)
        dp_cnt = cv2.drawContours(dapi_c, contours, -1, (0, 255, 0), 2)
        cv2.imwrite(dst_dir + img_path + '_dapi_contours_segmented.tif', dp_cnt)






# dapi segmentations by detecting contours for all images in the directory
# for img_path in os.listdir(img_dir):
#     if img_path.endswith(".tif"):
#         im_dapi, dapi_clr = read_img.read_and_preprocess(img_path, 2)
#         print("read", os.path.basename(img_path))
#         # plot_im(im_claudin)
#         dapi_contours = detect_contours.return_contours(im_dapi)
#         dpx, dpy, dpw, dph, dp_area, dapi_segmented = draw_contours.draw_detected_contours(im_dapi, 2, dapi_contours , (255,0,0), 2)
#         img_info_dapi = save_coordinates.create_df(dpx, dpy, dpw, dph, dp_area, dapi_segmented, im_dapi, 'DAPI')

# dapi edge detection for all images in the directory
edges_all_images(source_dir, dst_dir)
