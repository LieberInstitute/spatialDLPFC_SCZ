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


# morphological transformations
def morph_transform(original_img):
    image_clr = skimage.color.gray2rgb(original_img)
    shifted = cv2.pyrMeanShiftFiltering(image_clr, 21, 51) #dapi_clr
    print("shifted")
    gray = cv2.cvtColor(shifted, cv2.COLOR_BGR2GRAY)
    print("grayed")
    thresh = cv2.threshold(gray, 0, 255, cv2.THRESH_BINARY | cv2.THRESH_OTSU)[1]
    print("thresholded")
    # fig,ax = plt.subplots(figsize = (20,20))
    # ax.imshow(image_clr)
    # fig.show()
    return shifted, thresh, gray


# read and normalise the image
Image.MAX_IMAGE_PIXELS = None
def preprocess(filepath, ch_num):
    img = Image.open(filepath)
    img.seek(ch_num)
    if ch_num == 1: # claudin
        # img_claudin = cv2.normalize(np.array(img, dtype = 'uint8'), np.zeros(np.array(img, dtype = 'uint8').shape, np.double), 1.0, 0.0, cv2.NORM_MINMAX)
        img_claudin = np.array(img, dtype = 'uint8')
        img_claudin[img_claudin <= img_claudin.mean()] = 0
        img_claudin[img_claudin >= img_claudin.mean()] = 255
        return img_claudin
    if ch_num == 2: # DAPI
        # img_dapi = cv2.normalize(np.array(img, dtype = 'uint8'), np.zeros(np.array(img, dtype = 'uint8').shape, np.double), 1.0, 0.0, cv2.NORM_MINMAX)
        img_dapi = np.array(img, dtype = 'uint8')
        dapi_shifted, dapi_gray, dapi_thresh = morph_transform(img_dapi)
        return img_dapi, dapi_shifted, dapi_gray, dapi_thresh
    if ch_num == 3: #NeuN
        img_neun = np.array(img, dtype = 'uint8')
        neun_shifted, neun_gray, neun_thresh = morph_transform(img_neun)
        # img_neun = cv2.normalize(img_neun, np.zeros(img_neun.shape, np.double), 1.0, 0.0, cv2.NORM_MINMAX)
        return img_neun, neun_shifted, neun_gray, neun_thresh
    else: # wfa
        img_wfa = np.array(img, dtype = 'uint8')
        # img_wfa = cv2.normalize(img_arr, np.zeros(img_arr.shape, np.double), 1.0, 0.0, cv2.NORM_MINMAX)
        return img_wfa
