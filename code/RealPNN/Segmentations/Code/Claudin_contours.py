'''
For Visium-IF
Channel0 = DAPI, DAPI
Channel1 = Claudin5 (Alex 488),
Channel2 = NeuN (Alexa 555),
Channel3 = WFA (Alexa 647),
Channel4 = AF (Autofluorescence), sample AF
Channel5 = Thumbnail
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
from stitched_functions import read_img, watershed_segmentation, save_coordinates
from stitched_functions import *


# image paths
img_C1 = pyhere.here('processed-data', 'VistoSeg', 'captureAreas','V12F14-057_C1.tif') # /dcs04/lieber/marmaypag/spatialDLPFC_SCZ_LIBD4100/processed-data/VistoSeg/captureAreas/V12F14-057_C1.tif
img_D1 = pyhere.here('processed-data', 'VistoSeg', 'captureAreas','V12F14-057_D1.tif') # /dcs04/lieber/marmaypag/spatialDLPFC_SCZ_LIBD4100/processed-data/VistoSeg/captureAreas/V12F14-057_D1.tif
img_dir = pyhere.here('processed-data', 'VistoSeg', 'captureAreas')

# claudin segmentations by detecting contours for 1 image
# im_claudin = read_img.read_and_preprocess(img_C1, 1)
# plot_im(im_claudin)
# cla_contours = detect_contours.return_contours(im_claudin)
# clx,cly,clw,clh, cl_area, seg_cla = draw_contours.draw_detected_contours(im_claudin, 1, cla_contours , (255,0,0), 2)
# claudin_df = save_coordinates.create_df(clx,cly,clw,clh, cl_area, im_claudin, 'claudin')

# # claudin segmentations by detecting contours for all images in the directory
# for img_path in os.listdir(img_dir):
#     if img_path.endswith(".tif"):
#         im_claudin = read_img.read_and_preprocess(img_path, 1)
#         print("read", os.path.basename(img_path))
#         # plot_im(im_claudin)
#         cla_contours = detect_contours.return_contours(im_claudin)
#         clx,cly,clw,clh, cl_area, seg_cla = draw_contours.draw_detected_contours(im_claudin, 1, cla_contours , (255,0,0), 2)
#         img_info_claudin = save_coordinates.create_df(clx,cly,clw,clh, cl_area, im_claudin, 'claudin')

# watershed segmentations claudin
img_claudin, claudin_shifted, claudin_gray, claudin_thresh = read_img.read_and_preprocess(img_D1, 3)
# plot_im(img_claudin)
# fig,ax = plt.subplots(figsize = (20,20))
# ax.imshow(neun_thresh, cmap = 'gray')
# fig.show()
claudin_labels, claudin_localmax = watershed_segmentation.find_labels(claudin_thresh)
clx, cly, clw, clh, cl_area, claudin_segmented = watershed_segmentation.draw_rect_dapi(claudin_labels, claudin_gray, claudin_neun)
cv2.imwrite('/users/ukaipa/PNN/One_img/claudin_stitched_segmented_D1_run1.tif', claudin_segmented)
print("segmented image saved")
# claudin_df = save_coordinates.create_df(clx, cly, clw, clh, cl_area, img_claudin, 'Claudin')

