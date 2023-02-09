OMJGJSRJH
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

# file paths
img_C1 = pyhere.here('processed-data', 'VistoSeg', 'captureAreas','V12F14-057_C1.tif') # /dcs04/lieber/marmaypag/spatialDLPFC_SCZ_LIBD4100/processed-data/VistoSeg/captureAreas/V12F14-057_C1.tif
img_D1 = pyhere.here('processed-data', 'VistoSeg', 'captureAreas','V12F14-057_D1.tif') # /dcs04/lieber/marmaypag/spatialDLPFC_SCZ_LIBD4100/processed-data/VistoSeg/captureAreas/V12F14-057_D1.tif

# directory path
img_dir = pyhere.here('processed-data', 'VistoSeg', 'captureAreas')

img_dapi, dapi_shifted, dapi_gray, dapi_thresh = read_img.read_and_preprocess(img_D1, 2)
plot_im(dapi_thresh)
fig,ax = plt.subplots(figsize = (20,20))
ax.imshow(dapi_thresh, cmap = 'gray')
fig.show()
dapi_labels, dapi_localmax = watershed_segmentation.find_labels(dapi_thresh)
dpx, dpy, dpw, dph, dp_area, dapi_segmented = watershed_segmentation.draw_rect_dapi(dapi_labels, dapi_gray, img_dapi)
dapi_df = save_coordinates.create_df(dpx, dpy, dpw, dph, dp_area, im_claudin, 'claudin')

# claudin segmentations by detecting contours for all images in the directory
for img_path in os.listdir(img_dir):
    if img_path.endswith(".tif"):
        im_dapi = read_img.read_and_preprocess(img_path, 2)
        print("read", os.path.basename(img_path))
        # plot_im(im_claudin)
        dapi_contours = detect_contours.return_contours(im_dapi)
        dpx, dpy, dpw, dph, dp_area, dapi_segmented = draw_contours.draw_detected_contours(im_dapi, 2, dapi_contours , (255,0,0), 2)
        img_info_claudin = save_coordinates.create_df(dpx, dpy, dpw, dph, dp_area, dapi_segmented, im_dapi, 'DAPI')


