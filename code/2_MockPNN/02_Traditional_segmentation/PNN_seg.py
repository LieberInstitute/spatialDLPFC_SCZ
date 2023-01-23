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
csv_test = '/dcs04/lieber/marmaypag/spatialDLPFC_SCZ_LIBD4100/processed-data/2_MockPNN/Training_tiles/Manual_annotations/Annotations/20220712_VIF_MockPNN_Strong_Scan1_[10087,53057]_component_data_23.csv'
img_test = '/dcs04/lieber/marmaypag/spatialDLPFC_SCZ_LIBD4100/raw-data/images/2_MockPNN/Training_tiles/20220712_VIF_MockPNN_Strong_Scan1_[10087,53057]_component_data_23.tif'
#           '20220712_VIF_MockPNN_Strong_Scan1_[6384,53057]_component_data_11.tif'


# read and normalize the images
im_cla = read_norm(img_test, 1)
orig_wfa, im_wfa = read_norm(img_test, 3) #, im_wfa, img_arr_adj

# plot_img(orig_wfa, im_wfa)
# fig,ax = plt.subplots(figsize = (20,20))
# ax.imshow(orig_wfa, cmap = 'gray')
# fig.show()

# segment PNNs using Claudin-5 for a single image
cla_wfa_contour = detect_contours(im_cla)
clx,cly,clw,clh, cl_area, seg_cla_wfa = draw_contours(im_wfa, 1, cla_wfa_contour, (0,0,0), -1)
# plot_img(im_cla, seg_cla_wfa)
out_img_gry = skimage.color.rgb2gray(seg_cla_wfa) # convert to gray to find contours and increase contrast
wfa_contours = detect_contours(out_img_gry)
wfx, wfy, wfw, wfh, pnn_area, seg_wfa = draw_contours(out_img_gry, 3, wfa_contours, (0,0,255), 2)
# plot_img(im_wfa, seg_wfa)
img_info_wfa = create_df(wfx, wfy, wfw, wfh, pnn_area, img_test, 'PNN')
# wfa, img_info1_wfa = all_pix_pnns(img_info_wfa, seg_wfa, orig_wfa) # to get all the pixels and their mean pixel intensities
im_wfa_clr = skimage.color.gray2rgb(im_wfa)
approx, contours1, shapes, contour_img1 = detect_shape_pnns(seg_wfa, img_info_wfa, wfa_contours)
fig,ax = plt.subplots(figsize = (20,20))
ax.imshow(contour_img1, cmap = 'gray')
fig.show()


# segment PNNs from all images in a directory
dst_wfa_img = '/dcs04/lieber/marmaypag/spatialDLPFC_SCZ_LIBD4100/processed-data/2_MockPNN/Training_tiles/ML_annotations/Images/WFA_segmentations/'
for img_name in os.listdir(img_dir):
    if img_name.endswith('.tif'):
          im_cla = read_norm(os.path.join(img_dir, img_name), 1)
          im_wfa = read_norm(os.path.join(img_dir, img_name), 3)
          cla_wfa_contour = detect_contours(im_cla)
          clx,cly,clw,clh, cl_area, seg_cla_wfa = draw_contours(im_wfa, 1, cla_wfa_contour, (0,0,0), -1)
          out_img_gry = skimage. color.rgb2gray(seg_cla_wfa) # convert to gray to find contours and increase contrast
          wfa_contours = detect_contours(out_img_gry)
          wfx, wfy, wfw, wfh, pnn_area, seg_wfa = draw_contours(out_img_gry, 3, wfa_contours, (0,0,255), 2)
          # cv2.imwrite((dst_wfa_img + img_name.split('.')[0] + '.tif'), seg_wfa)
          img_info_wfa = create_df(wfx, wfy, wfw, wfh, pnn_area, img_test, 'PNN')
          print("img name:",img_name.split('.')[0],"PNNs:",len(img_info_wfa))


#### Draw the manual annotations contours to check how much of background is being cut off and also check FPs
manual_pnns = manual_annot(csv_test)
orig_wfa_clr = skimage.color.gray2rgb(orig_wfa)
manual = draw_rect(manual_pnns, orig_wfa_clr, (255,0,0))
img, df = all_pix_pnns(manual_pnns, manual,  orig_wfa)
x = [0.63221335, 0.5955076,0.6419167,0.67618203,0.7545275,0.75522727,0.7332946,0.7341727,0.72368777]

fig = plt.figure(figsize = (5, 5))
plt.bar(list(range(0,9)), x, color = 'blue', width = 0.2)
plt.xticks(np.arange(0,9, 1), labels = list(range(0,8)))
plt.xlabel("Number of segmented PNNs")
plt.ylabel("Mean pixel intensities")
plt.title("Mean pixel intensities plot for annotated PNNs")
plt.show()

fig,ax = plt.subplots(figsize = (20,20))
ax.imshow(img, cmap = 'gray')
fig.show()
