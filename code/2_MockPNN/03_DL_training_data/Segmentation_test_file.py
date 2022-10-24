'''
For Visium-IF
Channel0 = DAPI, DAPI
Channel1 = Claudin5 (Alex 488),
Channel2 = NeuN (Alexa 555),
Channel3 = WFA (Alexa 647),
Channel4 = AF (Autofluorescence), sample AF
Channel5 = Thumbnail
'''


import numpy as np
import numpy.ma as ma
import pyhere
from pathlib import Path
import pandas as pd
import PIL
from PIL import Image, ImageFont, ImageDraw, ImageEnhance, ImageSequence
import os
import matplotlib
import matplotlib.pyplot as plt
import mpldatacursor
import glob
import sys
import cv2
import math
import scipy
from scipy.spatial.distance import *
# from scipy import ndimage, distance
import imageio
import skimage
from skimage import *
from skimage import feature, segmentation, draw, measure, morphology
from skimage.morphology import (erosion,dilation,opening,closing,white_tophat,black_tophat,skeletonize,convex_hull_image)
from skimage.draw import polygon_perimeter
import tifffile as tif
import imagecodecs
import itertools
from itertools import product
from collections import defaultdict
from shapely.geometry import Point
from shapely.geometry.polygon import Polygon

# cd dcs04/lieber/marmaypag/spatialDLPFC_SCZ_LIBD4100/code
img_dir = pyhere.here('raw-data', 'images', '2_MockPNN', 'Training_tiles')
img_dir_NTC = pyhere.here('raw-data', 'images', '2_MockPNN', '20220712_VIF_MockPNN_Strong_NTC_C1_Br5182_MLtraining')
img_NTC = pyhere.here('raw-data', 'images', '2_MockPNN', '20220712_VIF_MockPNN_Strong_NTC_C1_Br5182_MLtraining', '20220712_VIF_MockPNN_Strong_NTC_Scan1_[11013,50974]_component_data.tif')
img_SCZ = pyhere.here('raw-data', 'images', '2_MockPNN', '20220712_VIF_MockPNN_Strong_SCZ_C1_Br2039_MLtraining', '20220712_VIF_MockPNN_Strong_SCZ_Scan1_[10629,49106]_component_data.tif')
img_test = pyhere.here('raw-data', 'images', '2_MockPNN', 'Training_tiles', '20220712_VIF_MockPNN_Strong_Scan1_[6384,53057]_component_data_11.tif')
                       # '20220712_VIF_MockPNN_Strong_SCZ_Scan1_[6925,49106]_component_data_24.tif')

img_test = '/dcs04/lieber/marmaypag/spatialDLPFC_SCZ_LIBD4100/raw-data/images/2_MockPNN/Training_tiles/20220712_VIF_MockPNN_Strong_Scan1_[6384,53057]_component_data_11.tif'
# WFA channel
img_wfa = Image.open(img_test)
img_wfa.seek(3) # channel 1 = Claudin 5
wfa = cv2.normalize(np.array(img_wfa, dtype = 'float32'), np.zeros(np.array(img_wfa, dtype = 'float32').shape, np.double), 1.0, 0.0, cv2.NORM_MINMAX)

# Claudin channel
img_claudin = Image.open(img_test)
img_claudin.seek(1) # channel 1 = Claudin 5
claudin = cv2.normalize(np.array(img_claudin, dtype = 'float32'), np.zeros(np.array(img_claudin, dtype = 'float32').shape, np.double), 1.0, 0.0, cv2.NORM_MINMAX)

# Increasing the contrast
claudin[claudin <= claudin.mean()] = 0.0
claudin[claudin >= claudin.mean()] = 1.0

# Plot the normalized/pre-processed image
fig,ax = plt.subplots(nrows = 1, ncols = 2,figsize = (20,20))
ax[0].imshow(claudin)
ax[1].imshow(wfa)
fig.show()

# Detecting contours
wfa255 = np.array(wfa * 255, dtype = np.uint8)
wfac = skimage.color.gray2rgb(wfa)
claudin_clr = skimage.color.gray2rgb((np.array((claudin * 255), dtype = np.uint8))) # convert to color to draw colored bb
hierachy, img_threshold = cv2.threshold(np.array(claudin * 255, dtype = np.uint8), 10, 255, cv2.THRESH_BINARY)
contours,_ = cv2.findContours(img_threshold, cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_SIMPLE)
for cnt in contours:
    x,y,w,h = cv2.boundingRect(cnt)
    print("CLAUDIN CONTOURS", x,y,w,h)
    out_img = cv2.rectangle(wfac, (x-10,y-10), (x+w+10, y+h+10), (0,0,0), -1)
# out_img now has black colored on most of the blood vessels; out has changed contrast pixels where PNNs are highlighted
# so now we find contours for the PNNs using the out_img

fig,ax = plt.subplots(nrows = 1, ncols = 2,figsize = (20,20))
ax[0].imshow(claudin)
ax[1].imshow(out_img)
fig.show()


out_img_gry = skimage.color.rgb2gray(out_img) # convert to gray to find contours and increase contrast
# out_img_gry[out_img_gry <= 0.2] = 0.0 # decrease contrast of background
# out_img_gry[out_img_gry >= 0.3] = 1.0 # increase the contrast of PNNs
fig,ax = plt.subplots(nrows = 1, ncols = 2, figsize = (20,20))
ax[0].imshow(out_img)
ax[1].imshow(out_img_gry)
fig.show()

out_img255 = np.array(out_img_gry * 255, dtype = np.uint8) # change scale to 0-255 for find contours
out_img_clr = skimage.color.gray2rgb(np.array(out_img_gry * 255, dtype = np.uint8)) # convert to color to draw colored bb
hierachy1, img_threshold1 = cv2.threshold(np.array(out_img_gry * 255, dtype = np.uint8), 100, 255, cv2.THRESH_BINARY) # ADAPTIVE_THRESH_MEAN_C
contours1,_ = cv2.findContours(img_threshold1, cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_SIMPLE)
# cv2.drawContours(out_imgc, contours1, -1, (0, 255, 0), 2, cv2.LINE_AA) # color scheme: BGR
for cnt in contours1:
    coords = cv2.boundingRect(cnt) # x,y,w,h
    x1,y1,w1,h1 = cv2.boundingRect(cnt)
    if (w1*h1) >= 300:
        print((x1-10,y1-10), (x1+w1+10, y1+h1+10))
        out_img1 = cv2.rectangle(out_img_clr, (x1-10,y1-10), (x1+w1+10, y1+h1+10), (0,255,0), 1) # change the color to black (0,0,0) if bb is not needed
        rect = cv2.minAreaRect(cnt)
        box = cv2.boxPoints(rect)
        box = np.int0(box)
        out_img1 = cv2.drawContours(out_img1,[box],0,(0,0,255),3) # comment out if contour box is not needed


fig,ax = plt.subplots(figsize = (20,20))
ax.imshow(out_img1)
fig.show()

# print(os.path.basename(file_name).split('.')[0].split('_')[6])
