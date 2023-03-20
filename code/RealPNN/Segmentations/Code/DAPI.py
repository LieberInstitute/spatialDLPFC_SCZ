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
import matplotlib.pyplot as plt
import cv2
import scipy
from scipy.spatial.distance import *
import skimage
import tifffile
from scipy.io import savemat
from skimage import feature, segmentation, draw, measure, morphology
from stitched_functions import read_img, watershed_segmentation, draw_contours
from stitched_functions import draw_contours, save_coordinates

# directory path
Image.MAX_IMAGE_PIXELS = None # increase the max image pixels to avoid decompression error
source_dir = '/dcs04/lieber/marmaypag/spatialDLPFC_SCZ_LIBD4100/processed-data/VistoSeg/captureAreas/'
dst_dir_dapi = '/dcs04/lieber/marmaypag/spatialDLPFC_SCZ_LIBD4100/processed-data/RealPNN/capture_area_segmentations/DAPI/DAPI_binarized/'

# test file paths
img_A1 = pyhere.here('processed-data', 'VistoSeg', 'captureAreas','V12F14-053_A1.tif')
img_B1 = '/dcs04/lieber/marmaypag/spatialDLPFC_SCZ_LIBD4100/processed-data/VistoSeg/captureAreas/V12F14-053_B1.tif'

# watershed segmentation on test image
# img_dapi, dapi_shifted, dapi_gray, dapi_thresh = read_img.read_and_preprocess(img_A1, 2)
# dapi_labels, dapi_localmax = find_labels(dapi_thresh) #watershed_segmentation.
# dpx, dpy, dpw, dph, dp_area, dapi_segmented = draw_rect_from_labels(dapi_labels, dapi_gray, img_dapi)
# cv2.imwrite(dst_dir_dapi + os.basename(img_A1)[0] + '_dapi_watershed_segmented.tif', dapi_segmented)
# dapi_df = save_coordinates.create_df(dpx, dpy, dpw, dph, dp_area, img_dapi, 'DAPI')


# contour detection on test image
Image.MAX_IMAGE_PIXELS = None
dapi_img = Image.open(img_A1)
dapi_img.seek(2)
dapi = np.array(dapi_img, dtype = 'uint8') # (17799, 16740)
dapi_c = cv2.cvtColor(dapi,cv2.COLOR_BGR2RGB)
gray = cv2.cvtColor(dapi_c,cv2.COLOR_RGB2GRAY)
_,thresh = cv2.threshold(gray, np.mean(gray), 255, cv2.THRESH_BINARY | cv2.THRESH_OTSU) #_INV
contours,_ = cv2.findContours(thresh, cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_SIMPLE)
print(len(contours))
# dapi_cnt = cv2.drawContours(dapi_c, contours, -1, (0, 0, 0), 1) #black
for cnt in contours:
    x,y,w,h = cv2.boundingRect(cnt)
    area = cv2.contourArea(cnt)
    if area >=50:
        # dapi_cnt = cv2.drawContours(dapi_c, contours, -1, (0, 0, 0), 1)
        dp_cnt = cv2.rectangle(dapi_c, (x,y), (x+w, y+h), (0,0,0), 1)
    elif area<50:
        # dapi_cnt = cv2.drawContours(dapi_c, contours, -1, (0, 0, 0), -1)
        dp_cnt = cv2.rectangle(dapi_c, (x,y), (x+w, y+h), (0,0,0), -1)
gray_segmented = cv2.cvtColor(dp_cnt,cv2.COLOR_RGB2GRAY)
thresh_segmented = cv2.threshold(gray_segmented, np.mean(gray_segmented), 255, cv2.THRESH_BINARY | cv2.THRESH_OTSU)[1] #_INV
binary_segmented = cv2.normalize(np.array(thresh_segmented, dtype = 'uint8'), np.zeros(np.array(thresh_segmented, dtype = 'uint8').shape, np.double), 1.0, 0.0, cv2.NORM_MINMAX)
cv2.imwrite(dst_dir_dapi + os.path.basename(img_A1).split('.')[0] + '_dapi.tif', thresh_segmented)
thresh_segmented.save(dst_dir_dapi + os.path.basename(img_A1).split('.')[0] + '_dapi.tif', 'tif')

_, binary_segmented_thresh = cv2.threshold(binary_segmented, 0.5, 1, cv2.THRESH_BINARY)
cv2.imwrite(filename, binary_segmented_thresh)
io.imsave(dst_dir_dapi + os.path.basename(img_A1).split('.')[0] + '_dapi.png', thresh_segmented, plugin='pil', check_contrast=False, compress=9)



# cv2.imwrite(dst_dir_dapi + os.path.basename(img_A1).split('.')[0] + '_dapi_segmented.tif', dp_cnt)
fig,ax = plt.subplots(figsize = (20,20))
ax.imshow(thresh_segmented, cmap = 'gray') #
fig.show()

for cont in contours:
    if cv2.contourArea(cont) >= 20:
        dapi_contoured_img = draw_contours.draw_dapi_contours(dapi_c, contours, (0, 255, 0), 2) #dp_cnt = cv2.drawContours(dapi_c, contours, -1, (0, 255, 0), 2)
cv2.imwrite(dst_dir_dapi + os.basename(img_A1)[0] + '_dapi_contours_thresholded.tif', dapi_contoured_img)
# dapi_df = save_coordinates.create_df(dpx, dpy, dpw, dph, dp_area, img_dapi, 'DAPI')

# cv2.imwrite(dst_dir_dapi + os.basename(img_A1)[0] + '_dapi_contours_detected.tif', dapi_contoured_img)
# dapi_df = save_coordinates.create_df(dpx, dpy, dpw, dph, dp_area, img_dapi, 'DAPI')
a = '/dcs04/lieber/marmaypag/spatialDLPFC_SCZ_LIBD4100/processed-data/RealPNN/capture_area_segmentations/DAPI/Segmented_images/V12F14-053_A1_dapi_size_thresholded.tif'
aim = Image.open(a)
fig,ax = plt.subplots(figsize = (20,20))
ax.imshow(binary_segmented, cmap = 'gray')
fig.show()



# find contours for all images in the dir
Image.MAX_IMAGE_PIXELS = None
for img_path in os.listdir(source_dir):
    if img_path.endswith(".tif"):
        dapi_img = Image.open(os.path.join(source_dir, img_path))
        dapi_img.seek(2)
        dapi = np.array(dapi_img, dtype = 'uint8')
        dapi_c = cv2.cvtColor(dapi,cv2.COLOR_BGR2RGB)
        gray = cv2.cvtColor(dapi_c,cv2.COLOR_RGB2GRAY)
        _,thresh = cv2.threshold(gray, np.mean(gray), 255, cv2.THRESH_BINARY | cv2.THRESH_OTSU) #_INV
        dapi_contours,_ = cv2.findContours(thresh, cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_SIMPLE)
        print("found", len(dapi_contours), "in", img_path)
        # dp_cnt = cv2.drawContours(dapi_c, contours, -1, (0, 255, 0), 2)
        for cnt in dapi_contours:
            x,y,w,h = cv2.boundingRect(cnt)
            area = cv2.contourArea(cnt)
            if area >=50:
                dp_cnt = cv2.rectangle(dapi_c, (x,y), (x+w, y+h), (0,0,0), 1)
            elif area<50:
                dp_cnt = cv2.rectangle(dapi_c, (x,y), (x+w, y+h), (0,0,0), -1)
        gray_segmented = cv2.cvtColor(dp_cnt,cv2.COLOR_RGB2GRAY)
        thresh_segmented = cv2.threshold(gray_segmented, np.mean(gray_segmented), 255, cv2.THRESH_BINARY | cv2.THRESH_OTSU)[1] #_INV
        # binary_segmented = cv2.normalize(np.array(thresh_segmented, dtype = 'uint8'), np.zeros(np.array(thresh_segmented, dtype = 'uint8').shape, np.double), 1.0, 0.0, cv2.NORM_MINMAX)
        cv2.imwrite(dst_dir_dapi + img_path.split('.')[0] + '_dapi_binarized.tif', thresh_segmented)

#             if cv2.contourArea(cnt) >=20:
#                 dpx, dpy, dpw, dph, dp_area, dp_segmented = draw_contours.draw_all_contours(dapi_c, dapi_contours, (255,125,155), 2)
#                 dapi_df = save_coordinates.create_df(dpx, dpy, dpw, dph, dp_area, img_path.split('.')[0], 'DAPI')
#                 dapi_df.to_csv(dst_dir_dapi + img_path.split('.')[0] + '_info.csv')
                # cv2.imwrite(dst_dir + img_path + '_dapi_contours_segmented.tif', dp_cnt)



# dapi segmentations by detecting contours for all images in the directory
# for img_path in os.listdir(img_dir):
#     if img_path.endswith(".tif"):
#         im_dapi, dapi_clr = read_img.read_and_preprocess(img_path, 2)
#         print("read", os.path.basename(img_path))
#         # plot_im(im_claudin)
#         dapi_contours = detect_contours.return_contours(im_dapi)
#         dpx, dpy, dpw, dph, dp_area, dapi_segmented = draw_contours.draw_detected_contours(im_dapi, 2, dapi_contours , (255,0,0), 2)
#         img_info_dapi = save_coordinates.create_df(dpx, dpy, dpw, dph, dp_area, dapi_segmented, im_dapi, 'DAPI')

# dapi edge detection for all images in the directory
# edges_all_images(source_dir, dst_dir)

# get all pixels info and append it to the existing csv
# csv_info_053_B1 = '/dcs04/lieber/marmaypag/spatialDLPFC_SCZ_LIBD4100/processed-data/RealPNN/capture_area_segmentations/DAPI/Data_files/V12F14-053_B1_info.csv'


# Image.MAX_IMAGE_PIXELS = None
# dapi_img = Image.open(img_B1)
# dapi_img.seek(2)
# dapi = np.array(dapi_img, dtype = 'uint8') # (17799, 16740)
# dapi_c = cv2.cvtColor(dapi,cv2.COLOR_BGR2RGB)
# gray = cv2.cvtColor(dapi_c,cv2.COLOR_RGB2GRAY)
# _,thresh = cv2.threshold(gray, np.mean(gray), 255, cv2.THRESH_BINARY | cv2.THRESH_OTSU) #_INV
# contours,_ = cv2.findContours(thresh, cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_SIMPLE)
# print(len(contours))

# B1_053 = pd.read_csv(csv_info_053_B1)
# contour_img, img_info_df, mean_pix_int_list = all_pix_pnns(img_info_df, contour_img, original_img)
