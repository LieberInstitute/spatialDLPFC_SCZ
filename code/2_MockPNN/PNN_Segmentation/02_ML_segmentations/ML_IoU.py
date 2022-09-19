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

# define the training tiles directory
img_dir = '/dcs04/lieber/marmaypag/spatialDLPFC_SCZ_LIBD4100/raw-data/images/2_MockPNN/Training_tiles/'
img_dir_NTC = '/dcs04/lieber/marmaypag/spatialDLPFC_SCZ_LIBD4100/raw-data/images/2_MockPNN/20220712_VIF_MockPNN_Strong_NTC_C1_Br5182_MLtraining/'
img_test = '/dcs04/lieber/marmaypag/spatialDLPFC_SCZ_LIBD4100/raw-data/images/2_MockPNN/Training_tiles/20220712_VIF_MockPNN_Strong_Scan1_[6925,49106]_component_data_24.tif'
# define the csvs directory
csv_dir = '/dcs04/lieber/marmaypag/spatialDLPFC_SCZ_LIBD4100/processed-data/2_MockPNN/Training_tiles/Manual_annotations/Annotations/'
csv_test = '/dcs04/lieber/marmaypag/spatialDLPFC_SCZ_LIBD4100/processed-data/2_MockPNN/Training_tiles/Manual_annotations/Annotations/20220712_VIF_MockPNN_Strong_Scan1_[6925,49106]_component_data_24.csv'

# read the tile and the manual annotation csv
img_wfa = Image.open(img_test)
img_wfa.seek(3)
wfa = cv2.normalize(np.array(img_wfa, dtype = 'float32'), np.zeros(np.array(img_wfa, dtype = 'float32').shape, np.double), 1.0, 0.0, cv2.NORM_MINMAX)
img_claudin = Image.open(img_test)
img_claudin.seek(1)
claudin = cv2.normalize(np.array(img_claudin, dtype = 'float32'), np.zeros(np.array(img_claudin, dtype = 'float32').shape, np.double), 1.0, 0.0, cv2.NORM_MINMAX)

# read the manual annotations csv into dataframe
df_manual_test = pd.read_csv(csv_test)

# Calculating all 4 coordinates of the BB
df_manual_test = df_manual_test.rename(columns = {'BX': 'x1', 'BY': 'y1'})
df_manual_test['x2'] = df_manual_test['x1'] + df_manual_test['Width']
df_manual_test['y2'] = df_manual_test['y1']
df_manual_test['x3'] = df_manual_test['x1']
df_manual_test['y3'] = df_manual_test['y1'] + df_manual_test['Height']
df_manual_test['x4'] = df_manual_test['x2']
df_manual_test['y4'] = df_manual_test['y3']

# convert x,y,bx,by from floating point to integers (doing it after, reduces round off errors)
df_manual_test['X'], df_manual_test['Y'], df_manual_test['x1'], df_manual_test['y1'] = np.int0(np.ceil(df_manual_test['X'])), np.int0(np.ceil(df_manual_test['Y'])), np.int0(np.ceil(df_manual_test['x1'])), np.int0(np.ceil(df_manual_test['y1']))
df_manual_test['x2'], df_manual_test['y2'], df_manual_test['x3'], df_manual_test['y3'], df_manual_test['x4'], df_manual_test['y4'] = np.int0(np.ceil(df_manual_test['x2'])), np.int0(np.ceil(df_manual_test['y2'])), np.int0(np.ceil(df_manual_test['x3'])), np.int0(np.ceil(df_manual_test['y3'])), np.int0(np.ceil(df_manual_test['x4'])), np.int0(np.ceil(df_manual_test['y4']))

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
clx, cly, clw, clh = [], [], [], []
for cnt in contours:
    x,y,w,h = cv2.boundingRect(cnt)
    print("CLAUDIN CONTOURS", x,y,w,h)
    clx.append(x)
    cly.append(y)
    clw.append(w)
    clh.append(h)
    out_img = cv2.rectangle(wfac, (x-10,y-10), (x+w+10, y+h+10), (0,0,0), -1)
# out_img now has black colored on most of the blood vessels; out has changed contrast pixels where PNNs are highlighted
# so now we find contours for the PNNs using the out_img

col_names = ['img_file_name','type_of_object_str', 'X', 'Y', 'W', 'H', 'no_of_claudin']
object_name = 'Blood_vessels' # name of the objects stored in the dataframe
file_name = '20220712_VIF_MockPNN_Strong_Scan1_[6925,49106]_component_data_24.tif' # image file name

dict = {col_names[0]: file_name, col_names[1]: object_name, col_names[2]: clx, col_names[3]: cly, col_names[4]: clw, col_names[5]: clh, col_names[6]: len(clx)}
df_claudin_ml = pd.DataFrame(dict, columns = col_names)


fig,ax = plt.subplots(nrows = 1, ncols = 2,figsize = (20,20))
ax[0].imshow(claudin)
ax[1].imshow(out_img)
fig.show()


out_img_gry = skimage.color.rgb2gray(out_img) # convert to gray to find contours and increase contrast
# out_img_gry[out_img_gry <= 0.2] = 0.0 # decrease contrast of background
# out_img_gry[out_img_gry >= 0.3] = 1.0 # increase the contrast of PNNs
# fig,ax = plt.subplots(nrows = 1, ncols = 2, figsize = (20,20))
# ax[0].imshow(out_img)
# ax[1].imshow(out_img_gry)
# fig.show()

out_img255 = np.array(out_img_gry * 255, dtype = np.uint8) # change scale to 0-255 for find contours
out_img_clr = skimage.color.gray2rgb(np.array(out_img_gry * 255, dtype = np.uint8)) # convert to color to draw colored bb
hierachy1, img_threshold1 = cv2.threshold(np.array(out_img_gry * 255, dtype = np.uint8), 100, 255, cv2.THRESH_BINARY) # ADAPTIVE_THRESH_MEAN_C
contours1,_ = cv2.findContours(img_threshold1, cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_SIMPLE)
# cv2.drawContours(out_imgc, contours1, -1, (0, 255, 0), 2, cv2.LINE_AA) # color scheme: BGR
wfx, wfy, wfw, wfh = [], [], [], []
for cnt in contours1:
    coords = cv2.boundingRect(cnt) # x,y,w,h
    x1,y1,w1,h1 = cv2.boundingRect(cnt)
    if (w1*h1) >= 300:
        print(x1,y1,w1,h1)
        wfx.append(x1)
        wfy.append(y1)
        wfw.append(w1)
        wfh.append(h1)
        out_img1 = cv2.rectangle(out_img_clr, (x1-10,y1-10), (x1+w1+10, y1+h1+10), (0,255,0), 1) # change the color to black (0,0,0) if bb is not needed
        rect = cv2.minAreaRect(cnt)
        box = cv2.boxPoints(rect)
        box = np.int0(box)
        out_img1 = cv2.drawContours(out_img1,[box],0,(0,0,255),3) # comment out if contour box is not needed


fig,ax = plt.subplots(figsize = (20,20))
ax.imshow(out_img1)
fig.show()

# Populate the data in the dataframe
col_names = ['img_file_name','type_of_object_str', 'X', 'Y', 'W', 'H', 'no_of_pnns']
object_name = 'PNN' # name of the objects stored in the dataframe
file_name = '20220712_VIF_MockPNN_Strong_Scan1_[6925,49106]_component_data_24.tif' # image file name

dict = {col_names[0]: file_name, col_names[1]: object_name, col_names[2]: wfx, col_names[3]: wfy, col_names[4]: wfw, col_names[5]: wfh, col_names[6]: len(wfx)}
df_wfa_ml = pd.DataFrame(dict, columns = col_names)

# Works without the big condition -- FOUND boxes that may actually overlap
for i1, rows1 in df_manual_test['X'].iteritems():
    for i2, rows2 in df_wfa_ml['X'].iteritems():
         # print(i1, rows1, i2, rows2, rows1-rows2)
         if (rows1-rows2) >= 0 and (rows1-rows2) <= 50 and (df_manual_test['Y'][i1]-df_wfa_ml['Y'][i2]) >= 0 and  (df_manual_test['Y'][i1]-df_wfa_ml['Y'][i2]) <= 50:
            # print(i1, rows1, i2, rows2, rows1-rows2) # 20 718 35 716 2
            print("Manual", i1, df_manual_test['X'][i1], df_manual_test['Y'][i1], df_manual_test['BX'][i1], df_manual_test['BY'][i1])
            print("ML", i2, df_wfa_ml['X'][i2], df_wfa_ml['Y'][i2], df_wfa_ml['W'][i2], df_wfa_ml['H'][i2])

'''''
Manual 1 331 499 314 478
ML 39 293 473 29 23
Manual 7 849 514 840 502
ML 41 839 472 23 14
Manual 18 282 421 266 404
ML 46 272 421 23 30
'''''


# For visual verification, the next steps are:
# to compare the # (i1,i2) to manual v/s predicted bb
# look at fiji for the coordinates of predicted bb and the marked number in the screenshot of manual annotations

def bb_intersection_over_union(boxA, boxB):
    # determine the (x, y)-coordinates of the intersection rectangle
    xA = max(boxA[0], boxB[0])
    yA = max(boxA[1], boxB[1])
    xB = min(boxA[2], boxB[2])
    yB = min(boxA[3], boxB[3])

    # compute the area of intersection rectangle
    interArea = abs(max((xB - xA, 0)) * max((yB - yA), 0))
    if interArea == 0:
        return 0
    # compute the area of both the prediction and ground-truth
    # rectangles
    boxAArea = abs((boxA[2] - boxA[0]) * (boxA[3] - boxA[1]))
    boxBArea = abs((boxB[2] - boxB[0]) * (boxB[3] - boxB[1]))

    # compute the intersection over union by taking the intersection
    # area and dividing it by the sum of prediction + ground-truth
    # areas - the interesection area
    iou = interArea / float(boxAArea + boxBArea - interArea)

    # return the intersection over union value
    return iou
