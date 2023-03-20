
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
dst_dir_wfa = '/dcs04/lieber/marmaypag/spatialDLPFC_SCZ_LIBD4100/processed-data/RealPNN/capture_area_segmentations/WFA/Segmented_images_binary/'

# file paths for test
Image.MAX_IMAGE_PIXELS = None
img_A1 = '/dcs04/lieber/marmaypag/spatialDLPFC_SCZ_LIBD4100/processed-data/VistoSeg/captureAreas/V12F14-053_A1.tif'
img_B1 = '/dcs04/lieber/marmaypag/spatialDLPFC_SCZ_LIBD4100/processed-data/VistoSeg/captureAreas/V12F14-053_B1.tif'


# WFA threshold
Image.MAX_IMAGE_PIXELS = None
for img_path in os.listdir(source_dir):
    if img_path.endswith(".tif"):
        wfa_img = Image.open(os.path.join(source_dir, img_path))
        wfa_img.seek(4)
        wfa = np.array(wfa_img, dtype = 'uint8')
        wfa_c = cv2.cvtColor(wfa,cv2.COLOR_BGR2RGB)
        hierachy, img_threshold = cv2.threshold(wfa,  100, 150, cv2.THRESH_BINARY)
        img_th_c = cv2.cvtColor(img_threshold,cv2.COLOR_BGR2RGB)
        fig,ax = plt.subplots(figsize = (20,20))
        ax.imshow(img_th_c)
        fig.show()
        wfa_contours,_ = cv2.findContours(img_threshold, cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_SIMPLE)
        print("found", len(wfa_contours), "in", img_path)
        wfa_cnt = cv2.drawContours(img_th_c, wfa_contours, -1, (0, 0, 0), 1) #(255, 153, 255)pink
        gray_segmented_wfa = cv2.cvtColor(wfa_cnt,cv2.COLOR_RGB2GRAY)
        thresh_segmented_wfa = cv2.threshold(gray_segmented_wfa, np.mean(gray_segmented_wfa), 255, cv2.THRESH_BINARY | cv2.THRESH_OTSU)[1] #_INV
        binary_segmented_wfa = cv2.normalize(np.array(thresh_segmented_wfa, dtype = 'uint8'), np.zeros(np.array(thresh_segmented_wfa, dtype = 'uint8').shape, np.double), 1.0, 0.0, cv2.NORM_MINMAX)
        # fig,ax = plt.subplots(figsize = (20,20))
        # ax.imshow(binary_segmented_wfa, cmap = 'gray')
        # fig.show()
        cv2.imwrite(dst_dir_wfa + img_path.split('.')[0] + '_wfa_binarized.tif', binary_segmented_wfa)
        # approx, contours, shape, contour_img = detect_shape_pnns(img_th_c, wfa_contours)



# wfa contours for 1 image
wfa_img = Image.open(img_A1)
wfa_img.seek(4)
wfa = np.array(wfa_img, dtype = 'uint8')
wfa_c = cv2.cvtColor(wfa,cv2.COLOR_BGR2RGB)
hierachy, img_threshold = cv2.threshold(wfa,  100, 150, cv2.THRESH_BINARY)
img_th_c = cv2.cvtColor(img_threshold,cv2.COLOR_BGR2RGB)
# fig,ax = plt.subplots(figsize = (20,20))
# ax.imshow(img_th_c)
# fig.show()
wfa_contours,_ = cv2.findContours(img_threshold, cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_SIMPLE)
print("found", len(wfa_contours), "in", os.path.basename(img_A1))
wfa_cnt = cv2.drawContours(img_th_c, wfa_contours, -1, (0, 255, 0), 2) #(255, 153, 255)pink
fig,ax = plt.subplots(figsize = (20,20))
ax.imshow(img_th_c)
fig.show()
gray_segmented_wfa = cv2.cvtColor(wfa_cnt,cv2.COLOR_RGB2GRAY)
thresh_segmented_wfa = cv2.threshold(gray_segmented_wfa, np.mean(gray_segmented_wfa), 255, cv2.THRESH_BINARY | cv2.THRESH_OTSU)[1] #_INV
# binary_segmented_wfa = cv2.normalize(np.array(thresh_segmented_wfa, dtype = 'uint8'), np.zeros(np.array(thresh_segmented_wfa, dtype = 'uint8').shape, np.double), 1.0, 0.0, cv2.NORM_MINMAX)
# fig,ax = plt.subplots(figsize = (20,20))
# ax.imshow(thresh_segmented_wfa, cmap = 'gray')
# fig.show()



# detecting claudin contours
claudin_img = Image.open(img_A1)
claudin_img.seek(1)
claudin = np.array(claudin_img, dtype = 'uint8')
claudin_c = cv2.cvtColor(claudin,cv2.COLOR_BGR2RGB)
gray = cv2.cvtColor(claudin_c,cv2.COLOR_RGB2GRAY)
_,thresh = cv2.threshold(gray, np.mean(gray), 255, cv2.THRESH_BINARY | cv2.THRESH_OTSU) #_INV
claudin_contours,_ = cv2.findContours(thresh, cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_SIMPLE)
# cla_cnt = cv2.drawContours(img_th_c, claudin_contours, -1, (255, 153, 255), 2) #pink
# print("found", len(claudin_contours))
area_ = []
for cnt in claudin_contours:
    x,y,w,h = cv2.boundingRect(cnt)
    area = cv2.contourArea(cnt)
    area_.append(area)
    if area<1000:
        cla_rect = cv2.rectangle(img_th_c, (x,y), (x+w, y+h), (0,0,255), 2) #green
    elif area>=10000:
        cla_rect = cv2.rectangle(img_th_c, (x,y), (x+50+w+100, y+50+h+100), (255,255,0), 2) #yellow


fig,ax = plt.subplots(figsize = (20,20))
ax.imshow(cla_rect)
fig.show()


cla_area = np.array(area_)

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
dapi_cnt = cv2.drawContours(img_th_c, dapi_contours, -1, (255, 153, 255), 2) #pink
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


cv2.imwrite('/dcs04/lieber/marmaypag/spatialDLPFC_SCZ_LIBD4100/processed-data/RealPNN/capture_area_segmentations/Claudin/A1_Claudin_bounding_box_test.tif', cont_claudin)
