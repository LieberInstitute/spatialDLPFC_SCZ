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
img_test = '/dcs04/lieber/marmaypag/spatialDLPFC_SCZ_LIBD4100/raw-data/images/2_MockPNN/Training_tiles/20220712_VIF_MockPNN_Strong_Scan1_[10087,53057]_component_data_23.tif'
#           '20220712_VIF_MockPNN_Strong_Scan1_[6384,53057]_component_data_11.tif'


# read and normalize the images
im_cla = read_norm(img_test, 1)
orig_wfa, im_wfa = read_norm(img_test, 3) #, im_wfa, img_arr_adj
# plot_img(im_cla, im_wfa)
# fig,ax = plt.subplots(figsize = (20,20))
# ax.imshow(im_wfa, cmap = 'gray')
# fig.show()

# segment PNNs using Claudin-5 for a single image
approx_cla, cla_wfa_contour, shape_cla = detect_contours(im_cla)
clx,cly,clw,clh, cl_area, seg_cla_wfa = draw_contours(im_wfa, 1, cla_wfa_contour, (0,0,0), -1)
# plot_img(im_cla, seg_cla_wfa)
out_img_gry = skimage.color.rgb2gray(seg_cla_wfa) # convert to gray to find contours and increase contrast
approx_wfa,wfa_contours, shape_wfa = detect_contours(out_img_gry)
wfx, wfy, wfw, wfh, pnn_area, seg_wfa = draw_contours(out_img_gry, 3, wfa_contours, (0,0,255), 2)
# plot_img(im_wfa, seg_wfa)
img_info_wfa = create_df(wfx, wfy, wfw, wfh, pnn_area, img_test, 'PNN')

normalised_img, ch_num, contours = None,  color = None, thickness = None, dapi_clr = None

# out_img_gry = skimage.color.rgb2gray(im_wfa)
wfa, img_info1_wfa = all_pix_pnns(img_info_wfa, seg_wfa, orig_wfa) # to get all the pixels and their mean pixel intensities
# fig,ax = plt.subplots(figsize = (20,20))
# ax.imshow(contour_img, cmap = 'gray')
# fig.show()


# segment PNNs from all images in a directory
dst_wfa_img = '/dcs04/lieber/marmaypag/spatialDLPFC_SCZ_LIBD4100/processed-data/2_MockPNN/Training_tiles/ML_annotations/Images/WFA_segmentations/'
for img_name in os.listdir(img_dir):
    if img_name.endswith('.tif'):
          # print(img_name.split('.')[0])
          im_cla = read_norm(os.path.join(img_dir, img_name), 1)
          im_wfa = read_norm(os.path.join(img_dir, img_name), 3)
          # print("read cla and wfa")
          cla_wfa_contour = detect_contours(im_cla)
          clx,cly,clw,clh, cl_area, seg_cla_wfa = draw_contours(im_wfa, 1, cla_wfa_contour, (0,0,0), -1)
          # print("claudins segmented")
          out_img_gry = skimage.color.rgb2gray(seg_cla_wfa) # convert to gray to find contours and increase contrast
          wfa_contours = detect_contours(out_img_gry)
          wfx, wfy, wfw, wfh, pnn_area, seg_wfa = draw_contours(out_img_gry, 3, wfa_contours, (0,0,255), 2)
          # print("PNNs segmented")
          # cv2.imwrite((dst_wfa_img + img_name.split('.')[0] + '.tif'), seg_wfa)
          img_info_wfa = create_df(wfx, wfy, wfw, wfh, pnn_area, img_test, 'PNN')
          print("img name:",img_name.split('.')[0],"PNNs:",len(img_info_wfa))


#### Trying the shape detector
shape = "unidentified"
peri = cv2.arcLength(c, True) # c is the contour
approx = cv2.approxPolyDP(c, 0.04 * peri, True)
