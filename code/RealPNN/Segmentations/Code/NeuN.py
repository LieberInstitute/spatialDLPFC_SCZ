
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
img_dir = '/dcs04/lieber/marmaypag/spatialDLPFC_SCZ_LIBD4100/processed-data/VistoSeg/captureAreas/'
dst_dir_neun = '/dcs04/lieber/marmaypag/spatialDLPFC_SCZ_LIBD4100/processed-data/RealPNN/single_channels_segmented/NeuN/slide4/'

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
Image.MAX_IMAGE_PIXELS = None
img_A1 = '/dcs04/lieber/marmaypag/spatialDLPFC_SCZ_LIBD4100/processed-data/VistoSeg/captureAreas/V12F14-053_A1.tif'
for img_path in os.listdir(img_dir):
    if img_path.endswith(".tif") and ('V13M06-279') in img_path:
        neun_img = Image.open(os.path.join(img_dir, img_path))
        neun_img.seek(3)
        neun = np.array(neun_img, dtype = 'uint8')
        neun_c = cv2.cvtColor(neun,cv2.COLOR_BGR2RGB)
        gray = cv2.cvtColor(neun_c,cv2.COLOR_RGB2GRAY)
        _,thresh = cv2.threshold(gray, np.mean(gray), 255, cv2.THRESH_BINARY | cv2.THRESH_OTSU) #_INV
        neun_contours,_ = cv2.findContours(thresh, cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_SIMPLE)
        print(len(neun_contours), img_path)
        nn_cnt = cv2.drawContours(neun_c, neun_contours, -1, (0, 0, 0), 1)
        # fig,ax = plt.subplots(figsize = (20,20))
        # ax.imshow(nn_cnt) # , cmap = 'gray'
        # plt.title(img_path.split('.')[0])
        # fig.show()
        for cnt in neun_contours:
            x,y,w,h = cv2.boundingRect(cnt)
            neun_area = cv2.contourArea(cnt)
            if neun_area >=100:
                neun_cnt = cv2.rectangle(neun_c, (x,y), (x+w, y+h), (0,255,0), 2)
            elif neun_area <100:
                neun_cnt = cv2.rectangle(neun_c, (x,y), (x+w, y+h), (255,0,0), -1)
        # gray_segmented = cv2.cvtColor(neun_cnt,cv2.COLOR_RGB2GRAY)
        # thresh_segmented = cv2.threshold(gray_segmented, np.mean(gray_segmented), 255, cv2.THRESH_BINARY | cv2.THRESH_OTSU)[1] #_INV
        # binary_segmented = cv2.normalize(np.array(thresh_segmented, dtype = 'uint8'), np.zeros(np.array(thresh_segmented, dtype = 'uint8').shape, np.double), 1.0, 0.0, cv2.NORM_MINMAX)
        # fig,ax = plt.subplots(figsize = (20,20))
        # ax.imshow(neun_cnt) # , cmap = 'gray'
        # fig.show()
        # nnx, nny, nnw, nnh, nn_area, neun_segmented = draw_contours.draw_all_contours(neun_c, neun_contours, (0,255,0), 2)
        # neun_df = save_coordinates.create_df(nnx, nny, nnw, nnh, nn_area, img_path.split('.')[0], 'NeuN')
        # neun_df.to_csv(dst_dir_neun + img_path.split('.')[0] + '_info.csv')
        cv2.imwrite(dst_dir_neun + img_path.split('.')[0] + '_neun_clr_segmented.tif', neun_cnt)


# tested out deriving all pixels from within the contour
Image.MAX_IMAGE_PIXELS = None
neun_img = Image.open(img_A1)
neun_img.seek(3)
neun = np.array(neun_img, dtype = 'uint8') # (17799, 16740)
neun_c = cv2.cvtColor(neun,cv2.COLOR_BGR2RGB)
gray = cv2.cvtColor(neun_c,cv2.COLOR_RGB2GRAY)
_,thresh = cv2.threshold(gray, np.mean(gray), 255, cv2.THRESH_BINARY | cv2.THRESH_OTSU) #_INV
neun_contours,_ = cv2.findContours(thresh, cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_SIMPLE)
# create a blank mask image of the same size as the input image
# mask = np.zeros(neun.shape[:2], dtype=np.uint8)
# # draw the contour on the mask image
# cv2.drawContours(mask, neun_contours, -1, 255, cv2.FILLED)
# # extract the pixels within the contour region using the mask
# neun_binary_img = cv2.bitwise_and(neun_c, neun_c, mask=mask)
neun_contoured_img = cv2.drawContours(neun_c, neun_contours, -1, (0, 255, 0), 2) # green all of the pixels
# for row in range(neun_contoured_img.shape[0]):
#     for col in range(neun_contoured_img.shape[1]):
#         pixel_val = neun_contoured_img[row][col][0]  # assuming the image is in RGB format
#         if pixel_val >= 10 and pixel_val < 30:
#             neun_contoured_img = cv2.drawContours(neun_c, neun_contours, -1, (255, 255, 0), 2) # yellow only bright pixels



# create a blank mask image of the same size as the input image
mask = np.zeros(neun_contoured_img.shape[:2], dtype=np.uint8)

# loop through the image and set the pixels with values between 10 and 30 to white
for row in range(neun_contoured_img.shape[0]):
    for col in range(neun_contoured_img.shape[1]):
        pixel_val = neun_contoured_img[row][col][0]  # assuming the image is in RGB format
        if pixel_val >= 10 and pixel_val < 30:
            mask[row][col] = 255

# find contours on the mask image --works for size threshold
contours, _ = cv2.findContours(mask, cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_SIMPLE)
for cnt in neun_contours:
    x,y,w,h = cv2.boundingRect(cnt)
    neun_area = cv2.contourArea(cnt)
    if neun_area >=50:
        neun_cnt = cv2.rectangle(neun_c, (x,y), (x+w, y+h), (0,255,0), 2)
    elif neun_area <50:
        neun_cnt = cv2.rectangle(neun_c, (x,y), (x+w, y+h), (255,255,125), 2)

# draw the contours on the original image
neun_contoured_img = cv2.drawContours(neun_c, contours, -1, (255, 255, 0), 2)

cv2.imwrite(dst_dir_neun + os.path.basename(img_A1).split('.')[0] + '_neun_size.tif', neun_cnt)

fig,ax = plt.subplots(figsize = (20,20))
ax.imshow(neun_cnt) #, cmap = 'gray'
fig.show()


# for pix intensity threshold
#find the contours,
#get all pixles within the contour, get mean pixels intensity of the contour,
#plot histogram, based on that, get the threshold value, set that to be the threhold for neun segmentation


# for row in range(neun.shape[0]):
#     for col in range(neun.shape[1]):
#         if neun[row][col] >= 10 and neun[row][col] < 30:
#             # print(len(neun_contours))
#             neun_contoured_img = cv2.drawContours(neun_c, neun_contours, -1, (0, 255, 0), 2)


# fig,ax = plt.subplots(figsize = (20,20))
# ax.imshow(neun_contoured_img) #, cmap = 'gray'
# fig.show()

# cv2.imwrite(dst_dir_dapi + os.basename(img_A1)[0] + '_claudin_contours_thresholded.tif', claudin_contoured_img)
# claudin_df = save_coordinates.create_df(clx,cly,clw,clh, cl_area, img_A1, 'Claudin-5')
# A1 = pd.read_csv(csv_A1)
# contour_img, neun_df_all, mean_pix_int_list = all_pixels.all_pix_pnns(A1, neun_contoured_img, neun)
# neun_df_all.to_csv(dst_dir_neun + img_A1.split('.')[0] + '_pix_info.csv')


# plot histogram
def hist_plot(img, bins=255):
    range = (img.min(), img.max()) # 50
    histogram, bin_edges = np.histogram(img.ravel(), bins=bins, range=range)
    plt.figure()
    plt.title("Grayscale Histogram")
    plt.xlabel("grayscale value")
    plt.ylabel("pixel count")
    plt.xlim([img.min(), img.max()]) # 50
    plt.plot(bin_edges[0:-1], histogram)
    plt.show()

hist_plot(neun)

