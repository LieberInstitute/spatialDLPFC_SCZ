'''
For Stitched Visium-IF tissue sections from VistoSeg SplitSlide output
Channel0 = AF
Channel1 = Claudin - 5 (Alex 488),
Channel2 = DAPI,
Channel3 = NeuN,
Channel4 = WFA
'''



import numpy as np
import helper_functions
import pyhere
import pylab
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
import itertools
from itertools import product
from collections import defaultdict
from shapely.geometry import Point
from shapely.geometry.polygon import Polygon


img_C1 = pyhere.here('processed-data', 'VistoSeg', 'captureAreas','V12F14-057_C1.tif')
img_D1 = pyhere.here('processed-data', 'VistoSeg', 'captureAreas','V12F14-057_D1.tif')

Image.MAX_IMAGE_PIXELS = None
dapi = Image.open(img_C1)
dapi.seek(2)
im_dapi = np.array(dapi, dtype = 'uint8')

dapi, dapi_clr = read_norm(img_C1, 2)
fig,ax = plt.subplots(figsize = (20,20))
ax.imshow(dapi, cmap = 'gray')
fig.show()
dpx, dpy, dpw, dph, area, ws_img_bb = draw_contours(dapi, 2, contours = None,  color = None, thickness = None, dapi_clr = dapi_clr)
fig,ax = plt.subplots(figsize = (20,20))
ax.imshow(ws_img_bb)
fig.show()
