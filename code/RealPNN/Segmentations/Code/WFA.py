
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

# directory path
Image.MAX_IMAGE_PIXELS = None
source_dir = '/dcs04/lieber/marmaypag/spatialDLPFC_SCZ_LIBD4100/processed-data/VistoSeg/captureAreas/'
dst_dir = '/dcs04/lieber/marmaypag/spatialDLPFC_SCZ_LIBD4100/processed-data/RealPNN/capture_area_segmentations/WFA/'

# file paths for test
Image.MAX_IMAGE_PIXELS = None
img_A1 = '/dcs04/lieber/marmaypag/spatialDLPFC_SCZ_LIBD4100/processed-data/VistoSeg/captureAreas/V12F14-053_A1.tif'


# functions
def detect_contours(normalised_img, lower_thresh): ### create a separate function for shape detection and run it through this loop for contour detection
    print("1) starting")
    hierachy, img_threshold = cv2.threshold(normalised_img, lower_thresh, 255, cv2.THRESH_BINARY)
    print("thresholded")
    contours,_ = cv2.findContours(img_threshold, cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_SIMPLE)
    print("contours detected", len(contours))
    return contours, img_threshold


def draw_contours(contours, color_img, rect_clr, rect_thick, cont_clr, cont_thick): #NeuN
    # color_img = skimage.color.gray2rgb(normalised_img)
    for count, cnt in enumerate(contours):
        x_, y_, w_, h_ = cv2.boundingRect(cnt)
        area_ = cv2.contourArea(cnt)
        if area_ >= 100 and area_ <= 1000:
            bb_img = cv2.rectangle(color_img, (x_,y_), (x_+w_+10, y_+h_+10), rect_clr, rect_thick,) #(255,0,0), 2-- to draw colored boxes
            box = np.int0(cv2.boxPoints(cv2.minAreaRect(cnt)))
            contour_img = cv2.drawContours(bb_img,[box],0,cont_clr, cont_thick) # change the color and thickness here if contours need to be visible
            print("Segmenting", count, "with area", area_)
    return contour_img

# claudin read and preprocess
claudin_img = Image.open(img_A1)
claudin_img.seek(1)
claudin = np.array(claudin_img, dtype = 'uint8')
clac = cv2.cvtColor(claudin,cv2.COLOR_BGR2RGB)
# read and preprocess wfa
wfa_img = Image.open(img_A1)
wfa_img.seek(4)
wfa = np.array(wfa_img, dtype = 'uint8')
fig,ax = plt.subplots(figsize = (20,20))
ax.imshow(wfa, cmap = 'gray')
fig.show()
wfac = cv2.cvtColor(wfa, cv2.COLOR_BGR2RGB)

# detect claudin contours
claudin_contours, cla_thresh = detect_contours(claudin, 50)
cont_claudin = draw_contours(claudin_contours, wfac, (0,0,0), -1, (0,0,0), -1)
fig,ax = plt.subplots(figsize = (20,20))
ax.imshow(wfa_cont)
fig.show()

# detect contours on wfa after claudin masking
wfa_contours, wfa_thresh = detect_contours(wfa, 100)
wfa_img = Image.open(img_A1)
wfa_img.seek(4)
wfa = np.array(wfa_img, dtype = 'uint8')
# hierachy, img_threshold = cv2.threshold(wfa, 200, 255, cv2.THRESH_BINARY | cv2.THRESH_OTSU)
img_threshold = cv2.adaptiveThreshold(wfa, 255, cv2.ADAPTIVE_THRESH_MEAN_C, cv2.THRESH_BINARY_INV, 21, 10)
# wfa_gray = cv2.cvtColor(wfac, cv2.COLOR_RGB2GRAY)
# thresh_wfa = threshold_otsu(wfa_gray)
# img_otsu  = wfa_gray < thresh_wfa

fig,ax = plt.subplots(figsize = (20,20))
ax.imshow(img_threshold, cmap = 'gray')
fig.show()


contours,_ = cv2.findContours(img_threshold, cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_SIMPLE)

wfa_cont = draw_contours(wfa_contours, wfac, (0,255,0), 1, (0,0,255), 2)
fig,ax = plt.subplots(figsize = (20,20))
ax.imshow(wfa_thresh, cmap = 'gray')
fig.show()


cv2.imwrite('/dcs04/lieber/marmaypag/spatialDLPFC_SCZ_LIBD4100/processed-data/RealPNN/capture_area_segmentations/Claudin/A1_Claudin_bounding_box_test.tif', cont_claudin)

# WFA threshold
wfa_img = Image.open(img_B1)
wfa_img.seek(4)
wfa = np.array(wfa_img, dtype = 'uint8')
wfa_c = cv2.cvtColor(wfa,cv2.COLOR_BGR2RGB)
hierachy, img_threshold = cv2.threshold(wfa,  100, 150, cv2.THRESH_BINARY)
img_th_c = cv2.cvtColor(img_threshold,cv2.COLOR_BGR2RGB)
fig,ax = plt.subplots(figsize = (20,20))
ax.imshow(img_threshold, cmap = 'gray')
fig.show()

# detecting claudin contours
claudin_img = Image.open(img_A1)
claudin_img.seek(1)
claudin = np.array(claudin_img, dtype = 'uint8')
claudin_c = cv2.cvtColor(claudin,cv2.COLOR_BGR2RGB)
gray = cv2.cvtColor(claudin_c,cv2.COLOR_RGB2GRAY)
_,thresh = cv2.threshold(gray, np.mean(gray), 255, cv2.THRESH_BINARY | cv2.THRESH_OTSU) #_INV
claudin_contours,_ = cv2.findContours(thresh, cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_SIMPLE)
cla_cnt = cv2.drawContours(img_th_c, claudin_contours, -1, (255, 153, 255), 2) #pink
# print("found", len(claudin_contours))

# drawing claudin contours on the color thresholded image of WFA
# cla_cnt = cv2.drawContours(img_th_c, claudin_contours, -1, (0, 255, 0), 2)
fig,ax = plt.subplots(figsize = (20,20))
ax.imshow(cla_cnt)
fig.show()

# drawing DAPI on the color thresholded image of WFA
dapi_img = Image.open(img_A1)
dapi_img.seek(2)
dapi = np.array(dapi_img, dtype = 'uint8')
dapi_c = cv2.cvtColor(dapi,cv2.COLOR_BGR2RGB)
gray = cv2.cvtColor(dapi_c,cv2.COLOR_RGB2GRAY)
_,thresh = cv2.threshold(gray, np.mean(gray), 255, cv2.THRESH_BINARY | cv2.THRESH_OTSU) #_INV
dapi_contours,_ = cv2.findContours(thresh, cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_SIMPLE)
print("found", len(dapi_contours), "DAPI")
dapi_cnt = cv2.drawContours(cla_cnt, dapi_contours, -1, (0, 255, 0), 2) #green
fig,ax = plt.subplots(figsize = (20,20))
ax.imshow(dapi_cnt)
fig.show()


# drawing neun contours
neun_img = Image.open(img_A1)
neun_img.seek(3)
neun = np.array(neun_img, dtype = 'uint8')
neun_c = cv2.cvtColor(neun,cv2.COLOR_BGR2RGB)
gray = cv2.cvtColor(neun_c,cv2.COLOR_RGB2GRAY)
_,thresh = cv2.threshold(gray, np.mean(gray), 255, cv2.THRESH_BINARY | cv2.THRESH_OTSU) #_INV
neun_contours,_ = cv2.findContours(thresh, cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_SIMPLE)
print("found", len(neun_contours), "in", img_path)
final = cv2.drawContours(dapi_cnt, neun_contours, -1, (255, 255, 0), 2) #yellow
fig,ax = plt.subplots(figsize = (20,20))
ax.imshow(final)
fig.show()
