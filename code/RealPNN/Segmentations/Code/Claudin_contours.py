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
from skimage.morphology import (erosion,dilation,opening,closing,white_tophat,black_tophat,skeletonize,convex_hull_image)
from skimage.draw import polygon_perimeter
from itertools import product
from collections import defaultdict
from shapely.geometry import Point
from shapely.geometry.polygon import Polygon
from stitched_functions import read_img, watershed_segmentation, save_coordinates
from stitched_functions import draw_contours, all_pixels


# directory path
Image.MAX_IMAGE_PIXELS = None # increase the max image pixels to avoid decompression error
source_dir = '/dcs04/lieber/marmaypag/spatialDLPFC_SCZ_LIBD4100/processed-data/VistoSeg/captureAreas/'
dst_dir_claudin = '/dcs04/lieber/marmaypag/spatialDLPFC_SCZ_LIBD4100/processed-data/RealPNN/capture_area_segmentations/Claudin/'


# image paths
img_A1 = pyhere.here('processed-data', 'VistoSeg', 'captureAreas','V12F14-053_A1.tif')
img_B1 = '/dcs04/lieber/marmaypag/spatialDLPFC_SCZ_LIBD4100/processed-data/VistoSeg/captureAreas/V12F14-053_B1.tif'
img_C2 = '/dcs04/lieber/marmaypag/spatialDLPFC_SCZ_LIBD4100/processed-data/VistoSeg/captureAreas/V12F14-057_C1.tif'
# img_C1 = pyhere.here('processed-data', 'VistoSeg', 'captureAreas','V12F14-057_C1.tif') # /dcs04/lieber/marmaypag/spatialDLPFC_SCZ_LIBD4100/processed-data/VistoSeg/captureAreas/V12F14-057_C1.tif
# img_D1 = pyhere.here('processed-data', 'VistoSeg', 'captureAreas','V12F14-057_D1.tif') # /dcs04/lieber/marmaypag/spatialDLPFC_SCZ_LIBD4100/processed-data/VistoSeg/captureAreas/V12F14-057_D1.tif
# img_dir = pyhere.here('processed-data', 'VistoSeg', 'captureAreas')

# csv paths
csv_A1 = '/dcs04/lieber/marmaypag/spatialDLPFC_SCZ_LIBD4100/processed-data/RealPNN/capture_area_segmentations/Claudin/Data_files/V12F14-053_A1_info.csv'

# claudin segmentations by detecting contours for 1 image
# im_claudin = read_img.read_and_preprocess(img_C1, 1)
# plot_im(im_claudin)
# cla_contours = detect_contours.return_contours(im_claudin)
# clx,cly,clw,clh, cl_area, seg_cla = draw_contours.draw_detected_contours(im_claudin, 1, cla_contours , (255,0,0), 2)
# claudin_df = save_coordinates.create_df(clx,cly,clw,clh, cl_area, im_claudin, 'claudin')
#
# # claudin segmentations by detecting contours for all images in the directory
# for img_path in os.listdir(img_dir):
#     if img_path.endswith(".tif"):
#         im_claudin = read_img.read_and_preprocess(img_path, 1)
#         print("read", os.path.basename(img_path))
#         # plot_im(im_claudin)
#         cla_contours = detect_contours.return_contours(im_claudin)
#         clx,cly,clw,clh, cl_area, seg_cla = draw_contours.draw_detected_contours(im_claudin, 1, cla_contours , (255,0,0), 2)
#         img_info_claudin = save_coordinates.create_df(clx,cly,clw,clh, cl_area, im_claudin, 'claudin')

# watershed segmentations claudin
# img_claudin, claudin_shifted, claudin_gray, claudin_thresh = read_img.read_and_preprocess(img_D1, 3)
# plot_im(img_claudin)
# fig,ax = plt.subplots(figsize = (20,20))
# ax.imshow(neun_thresh, cmap = 'gray')
# fig.show()
# claudin_labels, claudin_localmax = watershed_segmentation.find_labels(claudin_thresh)
# clx, cly, clw, clh, cl_area, claudin_segmented = watershed_segmentation.draw_rect_dapi(claudin_labels, claudin_gray, claudin_neun)
# cv2.imwrite('/users/ukaipa/PNN/One_img/claudin_stitched_segmented_D1_run1.tif', claudin_segmented)
# print("segmented image saved")
# claudin_df = save_coordinates.create_df(clx, cly, clw, clh, cl_area, img_claudin, 'Claudin')



# find contours for all images in the dir
# Image.MAX_IMAGE_PIXELS = None
# for img_path in os.listdir(source_dir):
#     if img_path.endswith(".tif"):
#         claudin_img = Image.open(os.path.join(source_dir, img_path))
#         claudin_img.seek(1)
#         claudin = np.array(claudin_img, dtype = 'uint8')
#         claudin_c = cv2.cvtColor(claudin,cv2.COLOR_BGR2RGB)
#         gray = cv2.cvtColor(claudin_c,cv2.COLOR_RGB2GRAY)
#         _,thresh = cv2.threshold(gray, np.mean(gray), 255, cv2.THRESH_BINARY | cv2.THRESH_OTSU) #_INV
#         claudin_contours,_ = cv2.findContours(thresh, cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_SIMPLE)
#         print("found", len(claudin_contours), "in", img_path)
#         # dp_cnt = cv2.drawContours(dapi_c, contours, -1, (0, 255, 0), 2)
#         clx,cly,clw,clh, cl_area, claudin_segmented = draw_contours.draw_all_contours(claudin_c, claudin_contours, (0,255,0), 2)
#         claudin_df = save_coordinates.create_df(clx,cly,clw,clh, cl_area, img_path.split('.')[0], 'Claudin-5')
#         claudin_df.to_csv(dst_dir_claudin + img_path.split('.')[0] + '_info.csv')
        # cv2.imwrite(dst_dir_neun + img_path + '_neun_contours_segmented.tif', neun_segmented)


# tested out deriving all pixels from within the contour
Image.MAX_IMAGE_PIXELS = None
claudin_img = Image.open(img_A1)
claudin_img.seek(2)
claudin = np.array(claudin_img, dtype = 'uint8') # (17799, 16740)
claudin_c = cv2.cvtColor(claudin,cv2.COLOR_BGR2RGB)
gray = cv2.cvtColor(claudin_c,cv2.COLOR_RGB2GRAY)
_,thresh = cv2.threshold(gray, np.mean(gray), 255, cv2.THRESH_BINARY | cv2.THRESH_OTSU) #_INV
claudin_contours,_ = cv2.findContours(thresh, cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_SIMPLE)
print(len(claudin_contours))
claudin_contoured_img = draw_contours.draw_all_contours(claudin_c, claudin_contours, (0, 255, 0), 2) #dp_cnt = cv2.drawContours(dapi_c, contours, -1, (0, 255, 0), 2)
# cv2.imwrite(dst_dir_dapi + os.basename(img_A1)[0] + '_claudin_contours_thresholded.tif', claudin_contoured_img)
# claudin_df = save_coordinates.create_df(clx,cly,clw,clh, cl_area, img_A1, 'Claudin-5')
A1 = pd.read_csv(csv_A1)
contour_img, claudin_df_all, mean_pix_int_list = all.pixels.all_pix_pnns(A1, claudin_contoured_img, claudin)
claudin_df_all.to_csv(dst_dir_claudin + img_path.split('.')[0] + '_pix_info.csv')
