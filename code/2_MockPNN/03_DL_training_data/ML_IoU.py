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
import pyhere
from pathlib import Path
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
from scipy.signal import convolve2d
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
img_dir = pyhere.here('raw-data', 'images', '2_MockPNN', 'Training_tiles')
img_dir_NTC = pyhere.here('raw-data', 'images', '2_MockPNN', '20220712_VIF_MockPNN_Strong_NTC_C1_Br5182_MLtraining')
img_NTC = pyhere.here('raw-data', 'images', '2_MockPNN', '20220712_VIF_MockPNN_Strong_NTC_C1_Br5182_MLtraining', '20220712_VIF_MockPNN_Strong_NTC_Scan1_[11013,50974]_component_data.tif')
img_SCZ = pyhere.here('raw-data', 'images', '2_MockPNN', '20220712_VIF_MockPNN_Strong_SCZ_C1_Br2039_MLtraining', '20220712_VIF_MockPNN_Strong_SCZ_Scan1_[10629,49106]_component_data.tif')
img_test = pyhere.here('raw-data', 'images', '2_MockPNN', 'Training_tiles', '20220712_VIF_MockPNN_Strong_Scan1_[6925,51188]_component_data_02.tif') #20220712_VIF_MockPNN_Strong_Scan1_[10087,51668]_component_data_01
csv_test = pyhere.here('processed-data', '2_MockPNN', 'Training_tiles', 'Manual_annotations', 'Annotations', '20220712_VIF_MockPNN_Strong_Scan1_[6925,51188]_component_data_02.csv') #20220712_VIF_MockPNN_Strong_Scan1_[10087,51668]_component_data_01

'20220712_VIF_MockPNN_Strong_Scan1_[12480,49800]_component_data_22'
'20220712_VIF_MockPNN_Strong_Scan1_[11939,48892]_component_data_07'
'20220712_VIF_MockPNN_Strong_Scan1_[6925,51188]_component_data_02'
# read the tile and the manual annotation csv
img_wfa = Image.open(img_test)
img_wfa.seek(3)
wfa = cv2.normalize(np.array(img_wfa, dtype = 'float32'), np.zeros(np.array(img_wfa, dtype = 'float32').shape, np.double), 1.0, 0.0, cv2.NORM_MINMAX)
img_claudin = Image.open(img_test)
img_claudin.seek(1)
claudin = cv2.normalize(np.array(img_claudin, dtype = 'float32'), np.zeros(np.array(img_claudin, dtype = 'float32').shape, np.double), 1.0, 0.0, cv2.NORM_MINMAX)

# Preprocessing the csv
conv_factor = 2.0112375738 # the fiji annotations are measured in microns and they need to be translated to pixels (1860/924.81)
df_manual_test = pd.read_csv(csv_test) # read the manual annotations csv into dataframe
df_manual_test = df_manual_test.rename(columns = {'X': 'xc', 'Y': 'yc', 'BX': 'x1', 'BY': 'y1', 'Perim.': 'Perimeter'}) # xc,yc are the centroids of the BB
df_manual_test.loc[:,['xc']], df_manual_test.loc[:,['yc']], df_manual_test.loc[:,['x1']], df_manual_test.loc[:,['y1']], df_manual_test['Width'], df_manual_test['Height'] = df_manual_test['xc']*conv_factor, df_manual_test['yc']*conv_factor, df_manual_test['x1']*conv_factor, df_manual_test['y1']*conv_factor, df_manual_test['Width']*conv_factor, df_manual_test['Height']*conv_factor
df_manual_test['x2'] = (df_manual_test['x1'] + df_manual_test['Width'])
df_manual_test['y2'], df_manual_test['x3'] = df_manual_test['y1'], df_manual_test['x1']
df_manual_test['y3'] = (df_manual_test['y1'] + df_manual_test['Height'])
df_manual_test['x4'], df_manual_test['y4']  = df_manual_test['x2'], df_manual_test['y3'] # Calculating all 4 coordinates of the BB
df_manual_test['xc'], df_manual_test['yc'], df_manual_test['x1'], df_manual_test['y1'] = np.int0(np.ceil(df_manual_test['xc'])), np.int0(np.ceil(df_manual_test['yc'])), np.int0(np.ceil(df_manual_test['x1'])), np.int0(np.ceil(df_manual_test['y1'])) # convert x,y,bx,by from floating point to integers (doing it after, reduces round off errors)
df_manual_test['x2'], df_manual_test['y2'], df_manual_test['x3'], df_manual_test['y3'], df_manual_test['x4'], df_manual_test['y4'] = np.int0(np.ceil(df_manual_test['x2'])), np.int0(np.ceil(df_manual_test['y2'])), np.int0(np.ceil(df_manual_test['x3'])), np.int0(np.ceil(df_manual_test['y3'])), np.int0(np.ceil(df_manual_test['x4'])), np.int0(np.ceil(df_manual_test['y4']))
df_manual_test = df_manual_test[['Area', 'Perimeter', 'Mean', 'Min', 'Max', 'xc', 'yc', 'x1', 'y1', 'x2', 'y2', 'x3', 'y3' , 'x4' , 'y4', 'Width', 'Height', 'Ch']] # rearranging the columns

# Increasing the contrast
claudin[claudin <= claudin.mean()] = 0.0
claudin[claudin >= claudin.mean()] = 1.0

# Plot the normalized/pre-processed image
fig,ax = plt.subplots(nrows = 1, ncols = 2,figsize = (20,20))
ax[0].imshow(claudin)
ax[1].imshow(wfa)
fig.show()

# Detecting contours
wfa_clr = skimage.color.gray2rgb(wfa)
claudin_clr = skimage.color.gray2rgb((np.array((claudin * 255), dtype = np.uint8))) # convert to color to draw colored bb
hierachy, img_threshold = cv2.threshold(np.array(claudin * 255, dtype = np.uint8), 10, 255, cv2.THRESH_BINARY)
contours,_ = cv2.findContours(img_threshold, cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_SIMPLE)
clx, cly, clw, clh = [], [], [], []
for cnt in contours:
    x,y,w,h = cv2.boundingRect(cnt)
    # print("CLAUDIN CONTOURS", x,y,w,h)
    clx.append(x)
    cly.append(y)
    clw.append(w)
    clh.append(h)
    out_img = cv2.rectangle(wfa_clr, (x-10,y-10), (x+w+10, y+h+10), (0,0,0), -2)
# out_img now has black colored on most of the blood vessels; out has changed contrast pixels where PNNs are highlighted
# so now we find contours for the PNNs using the out_img

col_names = ['img_file_name','type_of_object_str', 'x1', 'y1', 'Width', 'Height', 'total_number_claudin']
object_name = 'Blood_vessels' # name of the objects stored in the dataframe
file_name = os.path.basename(img_test) # image file name


dict = {col_names[0]: file_name, col_names[1]: object_name, col_names[2]: clx, col_names[3]: cly, col_names[4]:
    clw, col_names[5]: clh, col_names[6]: len(clx)}
df_claudin_ml = pd.DataFrame(dict, columns = col_names)

df_claudin_ml['x2'] = df_claudin_ml['x1'] + df_claudin_ml['Width']
df_claudin_ml['y2'], df_claudin_ml['x3'] = df_claudin_ml['y1'], df_claudin_ml['x1']
df_claudin_ml['y3'] = df_claudin_ml['y1'] + df_claudin_ml['Height']
df_claudin_ml['x4'], df_claudin_ml['y4'] = df_claudin_ml['x2'], df_claudin_ml['y3']
df_claudin_ml = df_claudin_ml[['img_file_name', 'type_of_object_str', 'x1', 'y1', 'x2', 'y2', 'x3', 'y3', 'x4', 'y4', 'Width', 'Height', 'total_number_claudin']]

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



# for 11: area = 268.973 xc = 840.782	yc = 80.548 X = 832.329	y = 72.593	w = 16.905	h = 15.911
out_img255 = np.array(out_img_gry * 255, dtype = np.uint8) # change scale to 0-255 for find contours
out_img_clr = skimage.color.gray2rgb(np.array(out_img_gry * 255, dtype = np.uint8)) # convert to color to draw colored bb
hierachy1, img_threshold1 = cv2.threshold(np.array(out_img_gry * 255, dtype = np.uint8), 100, 255, cv2.THRESH_BINARY) # ADAPTIVE_THRESH_MEAN_C
contours1,_ = cv2.findContours(img_threshold1, cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_SIMPLE)
# cv2.drawContours(out_imgc, contours1, -1, (0, 255, 0), 2, cv2.LINE_AA) # color scheme: BGR
wfx, wfy, wfw, wfh, pnn_area = [], [], [], [], []
for cnt in contours1:
    x1,y1,w1,h1 = cv2.boundingRect(cnt)
    area = cv2.contourArea(cnt)
    # print(area)
    if area >= 100 and area <= 20000:
        # print(x1,y1,w1,h1)
        wfx.append(x1)
        wfy.append(y1)
        wfw.append(w1)
        wfh.append(h1)
        # pnn_area.append(area)
        out_img1 = cv2.rectangle(out_img_clr, (x1-10,y1-10), (x1+w1+10, y1+h1+10), (0,0,0), 1) # change the color to black (0,0,0) if bb is not needed
        rect = cv2.minAreaRect(cnt)
        box = cv2.boxPoints(rect)
        box = np.int0(box)
        out_img1 = cv2.drawContours(out_img1,[box],0,(0,0,255),3) # comment out if contour box is not needed

fig,ax = plt.subplots(figsize = (20,20))
ax.imshow(out_img1)
fig.show()

# Populate the data in the dataframe
col_names = ['img_file_name','type_of_object_str', 'x1', 'y1', 'Width', 'Height', 'Area', 'total_number_pnns']
object_name = 'PNN' # name of the objects stored in the dataframe
file_name = os.path.basename(img_test) # image file name

dict = {col_names[0]: file_name, col_names[1]: object_name, col_names[2]: wfx, col_names[3]: wfy, col_names[4]: wfw, col_names[5]: wfh, col_names[6]: pnn_area, col_names[7]: len(wfx)}
df_wfa_ml = pd.DataFrame(dict, columns = col_names)

df_wfa_ml['x2'] = df_wfa_ml['x1'] + df_wfa_ml['Width']
df_wfa_ml['y2'], df_wfa_ml['x3'] = df_wfa_ml['y1'], df_wfa_ml['x1']
df_wfa_ml['y3'] = df_wfa_ml['y1'] + df_wfa_ml['Height']
df_wfa_ml['x4'], df_wfa_ml['y4'] = df_wfa_ml['x2'], df_wfa_ml['y3']
df_wfa_ml = df_wfa_ml[['img_file_name', 'type_of_object_str', 'x1', 'y1', 'x2', 'y2', 'x3', 'y3', 'x4', 'y4', 'Width', 'Height', 'Area', 'total_number_pnns']]


# Find boxes that may actually overlap
count = 0
for i1, rows1 in df_manual_test['x1'].iteritems():
    for i2, rows2 in df_wfa_ml['x1'].iteritems():
         # print(i1, rows1, i2, rows2, rows1-rows2)
         if (rows1-rows2) >= 0 and (rows1-rows2) <= 50:  # and (df_manual_test['y1'][i1]-df_wfa_ml['y1'][i2]) >= 0 and  (df_manual_test['y1'][i1]-df_wfa_ml['y1'][i2]) <= 50:
            # print(i1, rows1, i2, rows2, rows1-rows2) # 20 718 35 716 2
            print("Manual", i1, df_manual_test['x1'][i1], df_manual_test['y1'][i1], df_manual_test['x4'][i1], df_manual_test['y4'][i1])
            print("ML", i2, df_wfa_ml['x1'][i2], df_wfa_ml['y1'][i2], df_wfa_ml['x4'][i2], df_wfa_ml['y4'][i2])
            count += 1
            l = df_manual_test.iloc[i1]
            box1 = [[l['x1'], l['y1']], [l['x2'], l['y2']],
                    [l['x3'], l['y3']], [l['x4'], l['y4']]]
            ml = df_wfa_ml.iloc[i2]
            box2 = [[ml['x1'], ml['y1']], [ml['x2'], ml['y2']],
                    [ml['x3'], ml['y3']], [ml['x4'], ml['y4']]]
            poly_1 = Polygon(box1)
            poly_2 = Polygon(box2)
            print(poly_1, poly_2, poly_1.intersection(poly_2).area, poly_1.union(poly_2).area)
            # print(calculate_iou(box1, box2))




def calculate_iou(box_1, box_2):
    poly_1 = Polygon(box_1)
    poly_2 = Polygon(box_2)
    print(poly_1, poly_2, poly_1.intersection(poly_2).area, poly_1.union(poly_2).area)
    iou = poly_1.intersection(poly_2).area / poly_1.union(poly_2).area
    return iou


# Manual 1 314 478 349 520
# ML 19 312 380 356 428
l = df_manual_test.iloc[2]
box1 = [[l['x1'], l['y1']], [l['x2'], l['y2']],
        [l['x3'], l['y3']], [l['x4'], l['y4']]]
ml = df_wfa_ml.iloc[17]
box2 = [[ml['x1'], ml['y1']], [ml['x2'], ml['y2']],
        [ml['x3'], ml['y3']], [ml['x4'], ml['y4']]]



# draw a rectangle from the manual annotations csv on the contour detected image
def draw_rect(df_manual_test, contour_img):
    for box in range(len(df_manual_test['x1'])):
        print(box)
        rect = cv2.rectangle(contour_img, (df_manual_test['x1'][box], df_manual_test['y1'][box]), (df_manual_test['x4'][box], df_manual_test['y4'][box]), (255,0,0), 3)
    fig,ax = plt.subplots(figsize = (20,20))
    ax.imshow(out_img1)
    fig.show()
    return contour_img


