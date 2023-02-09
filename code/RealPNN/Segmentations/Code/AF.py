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

img_af, af_shifted, af_gray, af_thresh = read_img.read_and_preprocess(img_D1, 0)
plot_im(af_thresh)
fig,ax = plt.subplots(figsize = (20,20))
ax.imshow(dapi_thresh, cmap = 'gray')
fig.show()
af_labels, af_localmax = watershed_segmentation.find_labels(af_thresh)
afx, afy, afw, afh, af_area, af_segmented = watershed_segmentation.draw_rect_dapi(af_labels, af_gray, img_af)
af_df = save_coordinates.create_df(afx, afy, afw, afh, af_area, img_af, 'autofluorescence')

# claudin segmentations by detecting contours for all images in the directory
for img_path in os.listdir(img_dir):
    if img_path.endswith(".tif"):
        im_af = read_img.read_and_preprocess(img_path, 0)
        print("read", os.path.basename(img_path))
        # plot_im(im_claudin)
        af_labels, af_localmax = watershed_segmentation.find_labels(af_thresh)
        afx, afy, afw, afh, af_area, af_segmented = watershed_segmentation.draw_rect_dapi(af_labels, af_gray, img_af)
        af_df = save_coordinates.create_df(afx, afy, afw, afh, af_area, img_af, 'autofluorescence')
