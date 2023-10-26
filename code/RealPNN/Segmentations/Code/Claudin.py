'''
For Stitched Visium-IF tissue sections from VistoSeg SplitSlide output
Channel0 = AF
Channel1 = Claudin - 5 (Alex 488),
Channel2 = DAPI,
Channel3 = NeuN,
Channel4 = WFA
'''


from scipy import ndimage
import numpy as np
import pyhere
import pandas as pd
import PIL
from PIL import Image
import os
import matplotlib
import matplotlib.pyplot as plt
# import sys
import cv2
import math
import scipy
# from scipy.spatial.distance import *
import skimage
# from skimage import *
from itertools import product
from collections import defaultdict
# from stitched_functions import read_img, watershed_segmentation, save_coordinates
# from stitched_functions import draw_contours, all_pixels


# directory path
Image.MAX_IMAGE_PIXELS = None # increase the max image pixels to avoid decompression error
img_dir = pyhere.here('processed-data', 'VistoSeg', 'captureAreas')
img_dir = '/dcs04/lieber/marmaypag/spatialDLPFC_SCZ_LIBD4100/processed-data/VistoSeg/captureAreas/'
# dst_dir_claudin = pyhere.here('processed-data', 'RealPNN', 'capture_area_segmentations', 'Claudin', 'claudin_binarized')
# dst_dir_claudin = '/dcs04/lieber/marmaypag/spatialDLPFC_SCZ_LIBD4100/processed-data/RealPNN/capture_area_segmentations/Claudin/claudin_binarized/'
dst_dir_claudin = '/dcs04/lieber/marmaypag/spatialDLPFC_SCZ_LIBD4100/processed-data/RealPNN/single_channels_segmented/Claudin/slide4/'

# image paths
img_A1_r1 = pyhere.here('processed-data', 'VistoSeg', 'captureAreas','V12F14-053_A1.tif')
img_D1_r2 = pyhere.here('processed-data', 'VistoSeg', 'captureAreas', 'V12D07-334_D1.tif')
img_B1_r2 = pyhere.here('processed-data', 'VistoSeg', 'captureAreas', 'V12D07-334_B1.tif')

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


# find contours for all images in the dir
Image.MAX_IMAGE_PIXELS = None
for img_path in os.listdir(img_dir):
    if img_path.endswith(".tif") and ('V13M06-279') in img_path:
        claudin_img = Image.open(os.path.join(img_dir, img_path))
        claudin_img.seek(1)
        claudin = np.array(claudin_img, dtype = 'uint8')
        print("Raw image stats-", "\nMean:\t", claudin.mean(), "\nMax:\t", claudin.max(), "\nMin:\t", claudin.min())
        claudin_c = cv2.cvtColor(claudin,cv2.COLOR_BGR2RGB)
        gray = cv2.cvtColor(claudin_c,cv2.COLOR_RGB2GRAY)
        _,thresh = cv2.threshold(gray, np.mean(gray), 255, cv2.THRESH_BINARY | cv2.THRESH_OTSU) #_INV
        claudin_contours,_ = cv2.findContours(thresh, cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_SIMPLE)
        print("found", len(claudin_contours), "in", img_path)
        cl_cnt = cv2.drawContours(claudin_c, claudin_contours, -1, (0, 255, 0), 2)
        # cv2.imwrite(dst_dir_claudin + img_path.split('.')[0] + '_claudin_segmented_red.tif', cl_cnt)
        for cnt in claudin_contours:
            x,y,w,h = cv2.boundingRect(cnt)
            area = cv2.contourArea(cnt)
            if area<5000:
                cla_rect = cv2.rectangle(claudin_c, (x,y), (x+w, y+h), (0,0,0), 1) #green
            elif area>=5000:
                cla_rect = cv2.rectangle(claudin_c, (x,y), (x+50+w+100, y+50+h+100), (0,0,0), -1) #yellow
        gray_segmented = cv2.cvtColor(cl_cnt,cv2.COLOR_RGB2GRAY)
        thresh_segmented = cv2.threshold(gray_segmented, np.mean(gray_segmented), 255, cv2.THRESH_BINARY | cv2.THRESH_OTSU)[1] #_INV
        # binary_segmented = cv2.normalize(np.array(thresh_segmented, dtype = 'uint8'), np.zeros(np.array(thresh_segmented, dtype = 'uint8').shape, np.double), 1.0, 0.0, cv2.NORM_MINMAX)
        # fig,ax = plt.subplots(figsize = (20,20))
        # ax.imshow(thresh_segmented, cmap = 'gray')
        # fig.show()
        # clx,cly,clw,clh, cl_area, claudin_segmented = draw_contours.draw_all_contours(claudin_c, claudin_contours, (0,255,0), 2)
        # claudin_df = save_coordinates.create_df(clx,cly,clw,clh, cl_area, img_path.split('.')[0], 'Claudin-5')
        # claudin_df.to_csv(dst_dir_claudin + img_path.split('.')[0] + '_info.csv')
        cv2.imwrite(dst_dir_claudin + img_path.split('.')[0] + '_claudin_binarized.tif', thresh_segmented)


fig,ax = plt.subplots(figsize = (20,20))
ax.imshow(binary_segmented, cmap = 'gray')

fig.show()


# detecting claudin contours -- taken from WFA code
claudin_img = Image.open(img_B1_r2)
claudin_img.seek(1)
claudin = np.array(claudin_img, dtype = 'uint8')
print("Raw image stats-", "\nMean:\t", claudin.mean(), "\nMax:\t", claudin.max(), "\nMin:\t", claudin.min())
fig,ax = plt.subplots(figsize = (20,20))
ax.imshow(claudin, cmap = 'gray')
fig.show()

claudin_c = cv2.cvtColor(claudin,cv2.COLOR_BGR2RGB)
gray = cv2.cvtColor(claudin_c,cv2.COLOR_RGB2GRAY)
# thresh = cv2.adaptiveThreshold(gray, 255, cv2.ADAPTIVE_THRESH_MEAN_C, cv2.THRESH_BINARY, 21, 25)
_,thresh = cv2.threshold(gray, 50, 100, cv2.THRESH_BINARY | cv2.THRESH_OTSU) #_INV,
claudin_contours,_ = cv2.findContours(thresh, cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_SIMPLE)
cla_cnt = cv2.drawContours(claudin_c, claudin_contours, -1, (255, 153, 255), 2) #pink
print("found", len(claudin_contours))
area_ = []
for cnt in claudin_contours:
    x,y,w,h = cv2.boundingRect(cnt)
    area = cv2.contourArea(cnt)
    area_.append(area)
    if area<1000:
        cla_rect = cv2.rectangle(claudin_c, (x,y), (x+w, y+h), (0,0,0), -1) #green
    elif area>=5000:
        cla_rect = cv2.rectangle(claudin_c, (x,y), (x+50+w+100, y+50+h+100), (255,255,0), 2) #yellow

cv2.imwrite(dst_dir_claudin + 'V12D07-334_B1' + '_claudin_size_thresholded.tif', cla_rect)

fig,ax = plt.subplots(figsize = (20,20))
ax.imshow(cla_cnt, cmap = 'gray') #
fig.show()


# trying other thresholds to remove background noise from folds and artifacts
se=cv2.getStructuringElement(cv2.MORPH_RECT , (8,8))
bg=cv2.morphologyEx(claudin, cv2.MORPH_DILATE, se)
out_gray=cv2.divide(claudin, bg, scale=255)
out_binary=cv2.threshold(out_gray, 0, 255, cv2.THRESH_OTSU )[1]


# blur
blur = cv2.GaussianBlur(gray, (0,0), sigmaX=33, sigmaY=33)

# divide
divide = cv2.divide(gray, blur, scale=255)

# otsu threshold
thresh = cv2.threshold(divide, 0, 255, cv2.THRESH_BINARY+cv2.THRESH_OTSU)[1]

# apply morphology
kernel = cv2.getStructuringElement(cv2.MORPH_RECT, (3,3))
morph = cv2.morphologyEx(thresh, cv2.MORPH_CLOSE, kernel)

############# find the noisy area and try to maks out the noise only from that area ###########
import cv2
import numpy as np

# Load the image
claudin_img = Image.open(img_B1_r2)
claudin_img.seek(1)
claudin = np.array(claudin_img, dtype = 'uint8')

# convert to color first, to avoid the error
claudin_c = cv2.cvtColor(claudin,cv2.COLOR_BGR2RGB)

# Convert the image to grayscale
gray = cv2.cvtColor(claudin_c, cv2.COLOR_BGR2GRAY)

# Calculate the average pixel intensity of the entire image
average_intensity = np.mean(gray)

# Threshold the image to create a binary mask of the noisy region
_, mask = cv2.threshold(gray, average_intensity, 255, cv2.THRESH_BINARY)

# Perform morphological operations to remove noise from the mask
kernel = cv2.getStructuringElement(cv2.MORPH_RECT, (5, 5))
opening = cv2.morphologyEx(mask, cv2.MORPH_OPEN, kernel, iterations=2)

# Find contours in the mask
contours, _ = cv2.findContours(opening, cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_SIMPLE)

# Draw the contours on the original image
output = claudin_c.copy()
cv2.drawContours(output, contours, -1, (0, 255, 0), 2)

# Display the result
fig,ax = plt.subplots(figsize = (20,20))
ax.imshow(output) # , cmap = 'gray'
fig.show()
