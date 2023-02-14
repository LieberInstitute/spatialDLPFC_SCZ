
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


def find_labels(threshold):
    D = ndimage.distance_transform_edt(threshold) # Euclidean distance from binary to nearest 0-pixel
    print("distance measured")
    localMax = peak_local_max(D, indices=False, min_distance=1, labels=threshold) # find the local maxima for all the individual objects
    print("local max found")
    markers = ndimage.label(localMax, structure=np.ones((3, 3)))[0] # 8-connectivity connected component analysis
    print("local max markers found")
    labels = watershed(-D, markers, mask=threshold)
    print("{} unique segments found".format(len(np.unique(labels)) - 1))
    return labels, localMax


# extract the watershed algorithm labels
def draw_rect_dapi(labels, gray, img):
    x_, y_, w_, h_, area_ = [], [], [], [], []
    # print("1) entering the label loop")
    for counter, label in enumerate(np.unique(labels)):
        if label == 0: # label marked 0 are background
            continue
        mask = np.zeros(gray.shape, dtype="uint8") # create masks that only have the detected labels as foreground and 0 as background
        mask[labels == label] = 255
        # print("2) found the masks", counter)
        # detect contours in the mask and grab the largest one
        cnts = cv2.findContours(mask.copy(), cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_SIMPLE) # detect the watershed contours
        # print("3) found contours", counter)
        cnts = imutils.grab_contours(cnts) # extract only the contours
        # print("4) grabbed contours",counter)
        c = max(cnts, key=cv2.contourArea) # get the area
        x,y,w,h = cv2.boundingRect(c) # BB coordinates
        area1 = cv2.contourArea(c)
        # print("5) found BB coordinates and area", counter)
        if area1 >= 2: # adding shape filter to filter out the smaller contours and the biggest cluter contours 2 >= area1 <= 1000
            # print("6) appending BB coordinates", counter)
            x_.append(x)
            y_.append(y)
            w_.append(w)
            h_.append(h)
            area_.append(cv2.contourArea(c))
            # print("7) appended")
            ws_img_bb = cv2.rectangle(skimage.color.gray2rgb(img), (x,y), (x+w, y+h), (0,255,0), 1) # if a colored BB is not required then, change color to (0,0,0) and thickness to 1
            # print("8) drawing rectangles")
            print("drawing rectangle number",counter, "with area", area1)
    # cv2.imwrite('/users/ukaipa/PNN/One_img/dapi_stitched_segmented_D1_run1.tif', ws_img_bb)
    return x_, y_, w_, h_, area_, ws_img_bb


# # testing the segmentation on a single tile to check if the draw tile function is working
# filepath = '/dcs04/lieber/marmaypag/spatialDLPFC_SCZ_LIBD4100/raw-data/images/RealPNN/round1/20220814_VIF_PNN_S1_NTC/20220814_VIF_PNN_S1_NTC_Scan1_[12056,34215]_component_data.tif'
# dp, dp_clr = read_norm(filepath, 0)
# fig,ax = plt.subplots(figsize = (20,20))
# ax.imshow(dp, cmap = 'gray')
# fig.show()
# dapi_shifted, dapi_gray, dapi_thresh = morph_transform(dp_clr)
# fig,ax = plt.subplots(figsize = (20,20))
# ax.imshow(dapi_thresh, cmap = 'gray')
# fig.show()
# dp_lbls, dp_max = find_labels(dapi_thresh)
# dpx, dpy, dpw, dph, area, ws_img_bb=draw_rect_dapi(dp_lbls, dapi_gray, dp)
#
# fig,ax = plt.subplots(figsize = (20,20))
# ax.imshow(ws_img_bb)
# fig.show()
