
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
dst_dir_neun = '/dcs04/lieber/marmaypag/spatialDLPFC_SCZ_LIBD4100/processed-data/RealPNN/capture_area_segmentations/NeuN/'

# test images
img_A1 = pyhere.here('processed-data', 'VistoSeg', 'captureAreas','V12F14-053_A1.tif')
img_B1 = '/dcs04/lieber/marmaypag/spatialDLPFC_SCZ_LIBD4100/processed-data/VistoSeg/captureAreas/V12F14-053_B1.tif'
img_C2 = '/dcs04/lieber/marmaypag/spatialDLPFC_SCZ_LIBD4100/processed-data/VistoSeg/captureAreas/V12F14-057_C1.tif'

# csv paths
csv_A1 = '/dcs04/lieber/marmaypag/spatialDLPFC_SCZ_LIBD4100/processed-data/RealPNN/capture_area_segmentations/Claudin/Data_files/V12F14-053_A1_info.csv'


# img_neun, neun_shifted, neun_gray, neun_thresh = read_img.read_and_preprocess(img_D1, 3)
# plot_im(img_neun)
# fig,ax = plt.subplots(figsize = (20,20))
# ax.imshow(neun_thresh, cmap = 'gray')
# fig.show()

# neun_labels, neun_localmax = find_labels(neun_thresh) #watershed_segmentation.
# nnx, nny, nnw, nnh, nn_area, neun_segmented = draw_rect_dapi(neun_labels, neun_gray, img_neun) #watershed_segmentation.
# cv2.imwrite('/users/ukaipa/PNN/One_img/neun_stitched_segmented_D1_1629382.tif', neun_segmented)
# print("segmented image saved")

# dapi_df = save_coordinates.create_df(nnx, nny, nnw, nnh, nn_area, img_neun, 'NeuN')

# neun segmentations by detecting contours for all images in the directory
# for img_path in os.listdir(img_dir):
#     if img_path.endswith(".tif"):
#         im_neun = read_img.read_and_preprocess(img_path, 3)
#         print("read", os.path.basename(img_path))
#         # plot_im(im_claudin)
#         neun_contours = detect_contours.return_contours(im_neun)
#         nnx, nny, nnw, nnh, nn_area, neun_segmented = draw_contours.draw_detected_contours(im_neun, 3, neun_contours , (255,0,0), 2)
#         img_info_claudin = save_coordinates.create_df(nnx, nny, nnw, nnh, nn_area, im_neun, 'NeuN')


# find contours for all images in the dir
# Image.MAX_IMAGE_PIXELS = None
# for img_path in os.listdir(source_dir):
#     if img_path.endswith(".tif"):
#         neun_img = Image.open(os.path.join(source_dir, img_path))
#         neun_img.seek(3)
#         neun = np.array(neun_img, dtype = 'uint8')
#         neun_c = cv2.cvtColor(neun,cv2.COLOR_BGR2RGB)
#         gray = cv2.cvtColor(neun_c,cv2.COLOR_RGB2GRAY)
#         _,thresh = cv2.threshold(gray, np.mean(gray), 255, cv2.THRESH_BINARY | cv2.THRESH_OTSU) #_INV
#         neun_contours,_ = cv2.findContours(thresh, cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_SIMPLE)
#         print("found", len(neun_contours), "in", img_path)
#         # dp_cnt = cv2.drawContours(dapi_c, contours, -1, (0, 255, 0), 2)
#         nnx, nny, nnw, nnh, nn_area, neun_segmented = draw_contours.draw_all_contours(neun_c, neun_contours, (0,255,0), 2)
#         neun_df = save_coordinates.create_df(nnx, nny, nnw, nnh, nn_area, img_path.split('.')[0], 'NeuN')
#         neun_df.to_csv(dst_dir_neun + img_path.split('.')[0] + '_info.csv')
        # cv2.imwrite(dst_dir_neun + img_path + '_neun_contours_segmented.tif', neun_segmented)


# tested out deriving all pixels from within the contour
Image.MAX_IMAGE_PIXELS = None
neun_img = Image.open(img_A1)
neun_img.seek(2)
neun = np.array(neun_img, dtype = 'uint8') # (17799, 16740)
neun_c = cv2.cvtColor(neun,cv2.COLOR_BGR2RGB)
gray = cv2.cvtColor(neun_c,cv2.COLOR_RGB2GRAY)
_,thresh = cv2.threshold(gray, np.mean(gray), 255, cv2.THRESH_BINARY | cv2.THRESH_OTSU) #_INV
neun_contours,_ = cv2.findContours(thresh, cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_SIMPLE)
print(len(neun_contours))
neun_contoured_img = draw_contours.draw_all_contours(neun_c, neun_contours, (0, 255, 0), 2) #dp_cnt = cv2.drawContours(dapi_c, contours, -1, (0, 255, 0), 2)
# cv2.imwrite(dst_dir_dapi + os.basename(img_A1)[0] + '_claudin_contours_thresholded.tif', claudin_contoured_img)
# claudin_df = save_coordinates.create_df(clx,cly,clw,clh, cl_area, img_A1, 'Claudin-5')
A1 = pd.read_csv(csv_A1)
contour_img, neun_df_all, mean_pix_int_list = all.pixels.all_pix_pnns(A1, neun_contoured_img, neun)
neun_df_all.to_csv(dst_dir_neun + img_A1.split('.')[0] + '_pix_info.csv')
