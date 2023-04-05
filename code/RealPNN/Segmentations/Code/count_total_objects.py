
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
import re
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

img_dir_NTC = '/dcs04/lieber/marmaypag/spatialDLPFC_SCZ_LIBD4100/raw-data/images/RealPNN/round1/20220814_VIF_PNN_S1_NTC/'
img_dir_SCZ = '/dcs04/lieber/marmaypag/spatialDLPFC_SCZ_LIBD4100/raw-data/images/RealPNN/round1/20220814_VIF_PNN_S2_SCZ/'


# gets all the x,y coordinates of the slide
cnt_x, cnt_y = [], []
for img in os.listdir(img_dir_NTC):
    if img.endswith('.tif'):
        print(int(re.split(',', img.split('_')[6])[1].split(']')[0]))
        cnt_x.append(int(re.split(',', img.split('_')[6])[0].split('[')[1]))
        cnt_y.append(int(re.split(',', img.split('_')[6])[1].split(']')[0]))

np.unique(cnt_x)
np.unique(cnt_y)




##### two ways to do this ######
#1) blow up the neun channel
######## look for the fiducial frames by detecting perfect circles
### 2) in VistoSeg tutorial:
# find the size of the image (x,y)
# split x: x/4 = the size of each individual tissue section
# if there are any offsets, adjust them manually
# maybe try this on the NTC section because the tissues are not so wonky/tricky as compared to the SCZ slide
### first step would be to read in the qptiff and then convert this into a numpy array
# or maybe check if you can read the mat file in the NTC folder and convert that to np array
