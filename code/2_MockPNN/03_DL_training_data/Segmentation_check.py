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
import helper_functions
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

img_dir = '/dcs04/lieber/marmaypag/spatialDLPFC_SCZ_LIBD4100/raw-data/images/2_MockPNN/Training_tiles/'
csv_dir = '/dcs04/lieber/marmaypag/spatialDLPFC_SCZ_LIBD4100/processed-data/2_MockPNN/Training_tiles/Manual_annotations/Annotations/'
img_test = pyhere.here('raw-data', 'images', '2_MockPNN', 'Training_tiles', '20220712_VIF_MockPNN_Strong_Scan1_[6384,53057]_component_data_11.tif')
csv_test = '/dcs04/lieber/marmaypag/spatialDLPFC_SCZ_LIBD4100/processed-data/2_MockPNN/Training_tiles/Manual_annotations/Annotations/20220712_VIF_MockPNN_Strong_Scan1_[6384,53057]_component_data_11.csv'
img_test = '/dcs04/lieber/marmaypag/spatialDLPFC_SCZ_LIBD4100/raw-data/images/2_MockPNN/Training_tiles/20220712_VIF_MockPNN_Strong_Scan1_[6384,53057]_component_data_11.tif'

####### All for 1 image ######
# segment all the channels and put the code here
# open and normalise
im_dapi = read_norm(img_test, 0)
im_cla = read_norm(img_test, 1)
im_neun = read_norm(img_test, 2)


### works well and shows the plots but without decimal places on the x-axis
n, bins, patches = plt.hist(im_wfa, 30, range = [0,0.02], facecolor='gray', align='mid')
pylab.rc("axes", linewidth=8.0)
pylab.rc("lines", markeredgewidth=2.0)
plt.xlabel('pix int', fontsize=14)
plt.ylabel('# of targets', fontsize=14)
pylab.xticks(fontsize=15, rotation = 'vertical')
pylab.yticks(fontsize=15)
plt.grid(True)
plt.show()


fig,ax = plt.subplots(nrows = 1, ncols = 2,figsize = (20,20))
ax[0].imshow(im_neun)
ax[1].imshow(im_wfa)
fig.show()

# dapi segmentation (functions to be run first in the dapi segmentations code)
dapi_clr = skimage.color.gray2rgb((np.array((im_dapi * 255), dtype = np.uint8))) # convert to color to draw colored bb
shifted, thresh, gray = morph_transform(dapi_clr)
labels = find_labels(thresh)
dpx, dpy, dpw, dph, area, seg_dapi = draw_rect_dapi(labels, gray, dapi_clr)
img_info_dapi = create_df(dpx, dpy, dpw, dph, area, img_test, 'DAPI')

plot_img(im_dapi, seg_dapi)

# claudin segmentation
claudin_clr = skimage.color.gray2rgb((np.array((im_cla * 255), dtype = np.uint8))) # convert to color to draw colored bb
hierachy, img_threshold = cv2.threshold((np.array((im_cla * 255), dtype = np.uint8)), 100, 255, cv2.THRESH_BINARY)
contours,_ = cv2.findContours(img_threshold, cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_SIMPLE)
# cv2.drawContours(out_imgc, contours1, -1, (0, 255, 0), 2, cv2.LINE_AA) # color scheme: BGR len(contours)
clx, cly, clw, clh = [],[],[],[]
for cnt in contours:
      x,y,w,h = cv2.boundingRect(cnt)
      if(w*h >= 100):
            clx.append(x)
            cly.append(y)
            clw.append(w)
            clh.append(h)
            out_img_bb = cv2.rectangle(claudin_clr, (x,y), (x+w+5, y+h+5), (0,0,255), 2)
            rect = cv2.minAreaRect(cnt)
            box = cv2.boxPoints(rect)
            box = np.int0(box)
            out_img_cnt = cv2.drawContours(out_img_bb,[box],0,(0,255,0),1)

# Put the BB details in the csv
col_names = ['img_file_name','type_of_object_str', 'x1', 'y1', 'Width', 'Height', 'total_number_claudin']
object_name = 'Blood_vessels' # name of the objects stored in the dataframe
file_name = os.path.basename(img_test) # image file name

dict = {col_names[0]: file_name, col_names[1]: object_name, col_names[2]: clx, col_names[3]: cly, col_names[4]: clw, col_names[5]: clh}
img_info_claudin = pd.DataFrame(dict, columns = col_names)
# compute the rest of the coordinates of the BB
img_info_claudin['x2'] = img_info_claudin['x1'] + img_info_claudin['Width']
img_info_claudin['y2'], img_info_claudin['x3'] = img_info_claudin['y1'], img_info_claudin['x1']
img_info_claudin['y3'] = img_info_claudin['y1'] + img_info_claudin['Height']
img_info_claudin['x4'], img_info_claudin['y4'] = img_info_claudin['x2'], img_info_claudin['y3']
img_info_claudin = img_info_claudin[['img_file_name', 'type_of_object_str', 'x1', 'y1', 'x2', 'y2', 'x3', 'y3', 'x4', 'y4', 'Width', 'Height']]


# segment pnns and get their coordinates in a df
# claudin + wfa segmentation
wfa_clr = skimage.color.gray2rgb(im_wfa)
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
out_img_gry = skimage.color.rgb2gray(out_img) # convert to gray to find contours and increase contrast
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
    if area >= 1000:
        out_img1 = cv2.rectangle(out_img_clr, (x1-10,y1-10), (x1+w1+10, y1+h1+10), (0,0,0), -1) # eliminating all the big objects
    elif 100 <= area < 2000: # size threshold
    # if area >= 100 and area <= 2000: # area >= 2000
        # print(x1,y1,w1,h1)
        wfx.append(x1)
        wfy.append(y1)
        wfw.append(w1)
        wfh.append(h1)
        pnn_area.append(area)
        out_img1 = cv2.rectangle(out_img_clr, (x1-10,y1-10), (x1+w1+10, y1+h1+10), (0,0,0), 2) # change the color to black (0,0,0) if bb is not needed
        rect = cv2.minAreaRect(cnt)
        box = cv2.boxPoints(rect)
        box = np.int0(box)
        out_img1 = cv2.drawContours(out_img1,[box],0,(0,0,0),1) # comment out if contour box is not needed



# find the coordinates in other channels that actually overlap with the pnns
# plot histogram and find the avg pixel intensities of other channels in within the coordinates that overlap with pnns
# put the pixel intensities of the overlaped coordinates from other channels into a df
# check the overlap of all channels for presence/absence of pixel intensities
