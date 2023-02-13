
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
from stitched_functions import read_img, watershed_segmentation, save_coordinates
from stitched_functions import *


img_neun, neun_shifted, neun_gray, neun_thresh = read_img.read_and_preprocess(img_D1, 3)
plot_im(img_neun)
fig,ax = plt.subplots(figsize = (20,20))
ax.imshow(neun_thresh, cmap = 'gray')
fig.show()
neun_labels, neun_localmax = watershed_segmentation.find_labels(neun_thresh)
nnx, nny, nnw, nnh, nn_area, neun_segmented = watershed_segmentation.draw_rect_dapi(neun_labels, neun_gray, img_neun)
# dapi_df = save_coordinates.create_df(nnx, nny, nnw, nnh, nn_area, img_neun, 'NeuN')

# neun segmentations by detecting contours for all images in the directory
for img_path in os.listdir(img_dir):
    if img_path.endswith(".tif"):
        im_neun = read_img.read_and_preprocess(img_path, 3)
        print("read", os.path.basename(img_path))
        # plot_im(im_claudin)
        neun_contours = detect_contours.return_contours(im_neun)
        nnx, nny, nnw, nnh, nn_area, neun_segmented = draw_contours.draw_detected_contours(im_neun, 3, neun_contours , (255,0,0), 2)
        img_info_claudin = save_coordinates.create_df(nnx, nny, nnw, nnh, nn_area, im_neun, 'NeuN')
