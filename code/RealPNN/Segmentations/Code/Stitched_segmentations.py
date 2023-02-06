'''
For Stitched Visium-IF tissue sections from VistoSeg SplitSlide output
Channel0 = AF
Channel1 = Claudin - 5 (Alex 488),
Channel2 = DAPI,
Channel3 = NeuN,
Channel4 = WFA
'''



import numpy as np
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

###### DAPI segmentations
neun = Image.open(img_C1)
neun.seek(3)
im_neun = np.array(neun, dtype = 'uint8')
neun_clr = skimage.color.gray2rgb(im_neun)
fig,ax = plt.subplots(figsize = (20,20))
ax.imshow(neun, cmap = 'gray')
fig.show()
shifted, thresh, gray = morph_transform(neun_clr)
fig,ax = plt.subplots(figsize = (20,20))
ax.imshow(thresh, cmap = 'gray')
fig.show()

labels, localMax = find_labels(thresh)
dpx, dpy, dpw, dph, area, ws_img_bb = draw_rect_dapi(labels, gray, dapi_clr)
cv2.imwrite('/users/ukaipa/PNN/One_img/dapi_stitched_segmented_C1_run2.tif', ws_img_bb)




im_wfa, wfa = read_norm(img_C1, 4)
fig,ax = plt.subplots(figsize = (20,20))
ax.imshow(wfa, cmap = 'gray')
fig.show()
wfa255 = wfa * 255
ret,th1 = cv.threshold(wfa255,127,255,cv2.THRESH_BINARY)
th2 = cv2.adaptiveThreshold(wfa255,255,cv2.ADAPTIVE_THRESH_MEAN_C,cv2.THRESH_BINARY,11,2)
th3 = cv2.adaptiveThreshold(wfa255,255,cv2.ADAPTIVE_THRESH_GAUSSIAN_C,cv2.THRESH_BINARY,11,2)
# Otsu's thresholding
ret4,th4 = cv2.threshold(wfa255,0,255,cv2.THRESH_BINARY+cv2.THRESH_OTSU)
# Otsu's thresholding after Gaussian filtering
blur = cv2.GaussianBlur(wfa255,(5,5),0)
ret5,th5 = cv2.threshold(blur,0,255,cv2.THRESH_BINARY+cv2.THRESH_OTSU)

fig,ax = plt.subplots(figsize = (20,20))
ax.imshow(wfa255, cmap = 'gray')
fig.show()


im_dapi = np.array(dapi, dtype = 'uint8')

dapi, dapi_clr = read_norm(img_D1, 2)
fig,ax = plt.subplots(figsize = (20,20))
ax.imshow(dapi, cmap = 'gray')
fig.show()
shifted, thresh, gray = morph_transform(dapi)
fig,ax = plt.subplots(figsize = (20,20))
ax.imshow(gray, cmap = 'gray')
fig.show()
labels, localMax = find_labels(thresh)
dpx, dpy, dpw, dph, area, ws_img_bb = draw_rect_dapi(labels, gray, dapi_clr)
cv2.imwrite('/users/ukaipa/PNN/One_img/dapi_stitched_segmented_C1_run2.tif', ws_img_bb)

dapi = read_norm(img_C1, 2)
dpx, dpy, dpw, dph, area, ws_img_bb = draw_contours(dapi, 2, contours = None,  color = None, thickness = None, dapi_clr = dapi_clr)
fig,ax = plt.subplots(figsize = (20,20))
ax.imshow(ws_img_bb)
fig.show()


######### NeuN segmentations
neun = read_norm(img_C1, 3)
fig,ax = plt.subplots(figsize = (20,20))
ax.imshow(neun, cmap = 'gray')
fig.show()
shifted1, thresh1, gray1 = morph_transform(neun)
cv2.imwrite('/users/ukaipa/PNN/One_img/neun_shifted.tif', shifted1)
cv2.imwrite('/users/ukaipa/PNN/One_img/neun_thresh.tif', thresh1)
cv2.imwrite('/users/ukaipa/PNN/One_img/neun_gray.tif', gray1)
fig,ax = plt.subplots(figsize = (20,20))
ax.imshow(gray1, cmap = 'gray')
fig.show()
labels1 = find_labels(thresh1)
neunx, neuny, neunw, neunh, neun_area, ws_neun_bb = draw_rect_dapi(labels1, gray1, neun)

neun_contours = detect_contours(neun)
nx,ny,nw,nh,narea, seg_neun = draw_contours(neun, 2, neun_contours,(0,255,0), 2)
plot_img(neun, seg_neun)



########## Claudin segmentation
cla = read_norm(img_C1, 1)
fig,ax = plt.subplots(figsize = (20,20))
ax.imshow(cla, cmap = 'gray')
fig.show()
cla_contours = detect_contours(cla)
clx,cly,clw,clh, cl_area, seg_cla = draw_contours(im_cla, 1, cla_contours , (255,0,0), 2)
img_info_claudin = create_df(clx,cly,clw,clh, cl_area, img_test, 'claudin')


