import numpy as np
import numpy.ma as ma
import pandas as pd
import csv
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
import re
from itertools import product
from collections import defaultdict
from shapely.geometry import Point
from shapely.geometry.polygon import Polygon

csv_dir = '/dcs04/lieber/marmaypag/spatialDLPFC_SCZ_LIBD4100/processed-data/2_MockPNN/Training_tiles/Manual_annotations/Annotations/'
img_dir = '/dcs04/lieber/marmaypag/spatialDLPFC_SCZ_LIBD4100/raw-data/images/2_MockPNN/Training_tiles/'

df = pd.read_csv(os.path.join(csv_dir, os.listdir(csv_dir)[0]))
img_test = os.path.join(img_dir, os.listdir(img_dir)[23])
img_test = '/dcs04/lieber/marmaypag/spatialDLPFC_SCZ_LIBD4100/raw-data/images/2_MockPNN/Training_tiles/20220712_VIF_MockPNN_Strong_Scan1_[10087,51668]_component_data_01.tif'

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

out_img_gry = skimage.color.rgb2gray(out_img) # convert to gray to find contours and increase contrast

out_img255 = np.array(out_img_gry * 255, dtype = np.uint8) # change scale to 0-255 for find contours
out_img_clr = skimage.color.gray2rgb(np.array(out_img_gry * 255, dtype = np.uint8)) # convert to color to draw colored bb
hierachy1, img_threshold1 = cv2.threshold(np.array(out_img_gry * 255, dtype = np.uint8), 50, 255, cv2.THRESH_BINARY) # ADAPTIVE_THRESH_MEAN_C
contours1,_ = cv2.findContours(img_threshold1, cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_SIMPLE)
# cv2.drawContours(out_imgc, contours1, -1, (0, 255, 0), 2, cv2.LINE_AA) # color scheme: BGR
for cnt in contours1:
    coords = cv2.boundingRect(cnt) # x,y,w,h
    x1,y1,w1,h1 = cv2.boundingRect(cnt)
    if (w1*h1) >= 300:
        print(x1,y1,w1,h1)
        out_img1 = cv2.rectangle(out_img_clr, (x1-10,y1-10), (x1+w1+10, y1+h1+10), (0,255,0), 1) # change the color to black (0,0,0) if bb is not needed
        rect = cv2.minAreaRect(cnt)
        box = cv2.boxPoints(rect)
        box = np.int0(box)
        out_img1 = cv2.drawContours(out_img1,[box],0,(0,0,255),3) # comment out if contour box is not needed

fig,ax = plt.subplots(figsize = (20,20))
ax.imshow(out_img1)
fig.show()


# Initialize empty list
lst_intensities = []

# For each list of contour points...
for i in range(len(contours1)):
    # Create a mask image that contains the contour filled in
    cimg = np.zeros_like(out_img_gry)
    cv2.drawContours(cimg, contours1, i, color=255, thickness=-1)
    pts = np.where(cimg == 255)
    lst_intensities.append(out_img_gry[pts[0], pts[1]])



lst = []
for j in range(len(lst_intensities)):
    l =  lst_intensities[j] * 255
    lst.append(l)

high = []
for k in range(len(lst)):
    if lst[k].mean() >= 90:
        print(k, lst[k].mean())
        hi = lst[k].mean()
        high.append(hi)


fig,ax = plt.subplots(nrows = 1, ncols = 2, figsize = (20,20))
ax[0].hist(np.array(high).ravel(),256,[0,255])
ax[0].title.set_text('Histogram of mean pixel intensities of bright PNNs')
ax[0].set(xlabel = 'pixel intensities(grayscale)' , ylabel = 'frequency')
ax[1].hist(np.array(low).ravel(),256,[0,255])
ax[1].title.set_text('Histogram of mean pixel intensities of low intensity PNNs')
ax[1].set(xlabel = 'pixel intensities(grayscale)' , ylabel = 'frequency')
fig.show()
