'''
For Stitched Visium-IF tissue sections from VistoSeg SplitSlide output
Channel0 = AF
Channel1 = Claudin - 5 (Alex 488),
Channel2 = DAPI,
Channel3 = NeuN,
Channel4 = WFA
'''


import PIL
from PIL import Image, ImageFont, ImageDraw, ImageEnhance, ImageSequence
import os
import numpy as np
import pandas as pd
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
from skimage.filters import threshold_otsu
from skimage.morphology import (erosion,dilation,opening,closing,white_tophat,black_tophat,skeletonize,convex_hull_image)
from skimage.draw import polygon_perimeter
from itertools import product
from collections import defaultdict
# from shapely.geometry import Point
# from shapely.geometry.polygon import Polygon # shapely module error

# directory path
Image.MAX_IMAGE_PIXELS = None
source_dir = '/dcs04/lieber/marmaypag/spatialDLPFC_SCZ_LIBD4100/processed-data/VistoSeg/captureAreas/'
img_dir = '/dcs04/lieber/marmaypag/spatialDLPFC_SCZ_LIBD4100/processed-data/VistoSeg/captureAreas/'
dst_dir_wfa = '/dcs04/lieber/marmaypag/spatialDLPFC_SCZ_LIBD4100/processed-data/RealPNN/single_channels_segmented/WFA/slide10/'

# WFA threshold
Image.MAX_IMAGE_PIXELS = None
for img_path in os.listdir(source_dir):
    if img_path.endswith(".tif") and ('V13F27-295') in img_path:
        wfa_img = Image.open(os.path.join(source_dir, img_path))
        wfa_img.seek(4)
        wfa = np.array(wfa_img, dtype = 'uint8')
        wfa_c = cv2.cvtColor(wfa,cv2.COLOR_BGR2RGB)
        # adjusted = cv2.convertScaleAbs(wfa, alpha=0.3, beta=10) # decreased the contrast of the original image for better segmentation
        hierachy, img_threshold = cv2.threshold(wfa,  80, 255, cv2.THRESH_BINARY) # 150;; 80
        img_th_c = cv2.cvtColor(img_threshold,cv2.COLOR_BGR2RGB)
        # fig,ax = plt.subplots(figsize = (20,20))
        # ax.imshow(img_threshold, cmap = 'gray')
        # fig.show()
        # cv2.imwrite(dst_dir_wfa + img_path.split('.')[0] + '_thresholded.tif', img_threshold)
        wfa_contours,_ = cv2.findContours(img_threshold, cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_SIMPLE)
        print("found", len(wfa_contours), "in", img_path)
        # wfa_cnt = cv2.drawContours(img_th_c, wfa_contours, -1, (0, 0, 0), 1) # yellow all contours
        area_ = []
        count_outside_range = 0
        for cnt in wfa_contours:
            x,y,w,h = cv2.boundingRect(cnt)
            area = cv2.contourArea(cnt)
            area_.append(area) # max(area_) = 1809854.0
            # print(area)
            # max, avg = max(area_), (sum(area_)/len(area_)).astype('uint8')
            if area<=1000: # area<150 and area>=10000 and area>=100000
                # print(f"Area {area} satisfies condition 1")
                wfa_cnt = cv2.rectangle(img_th_c, (x,y), (x+w, y+h), (0,0,0), -1) # rectangle
            else:
                # print(f"Area {area} satisfies condition 2")
                wfa_cnt = cv2.rectangle(img_th_c, (x,y), (x+w, y+h), (0,0,0), 1) # rectangle
        gray_segmented_wfa = cv2.cvtColor(img_th_c,cv2.COLOR_RGB2GRAY)
        thresh_segmented_wfa = cv2.threshold(gray_segmented_wfa, 80, 255, cv2.THRESH_BINARY)[1] #_INV # | cv2.THRESH_OTSU
        cv2.imwrite(dst_dir_wfa + img_path.split('.')[0] + '_wfa_segmented.tif', thresh_segmented_wfa)

        # fig,ax = plt.subplots(figsize = (20,20))
        # ax.imshow(wfa_cnt)
        # fig.show()
            # else:
            #     # print(f"Area {area} satisfies the else condition")
            #     wfa_cnt = cv2.rectangle(img_th_c, (x,y), (x+w, y+h), (0,0,0), -1) # rectangle
                # if area >= 100000:
                #     text = f"{area:.2f}"
                #     text_size = cv2.getTextSize(text, cv2.FONT_HERSHEY_SIMPLEX, 5.0, 2)[0]
                #     text_x = x + (w - text_size[0]) // 2
                #     text_y = y + (h + text_size[1]) // 2
                #     cv2.putText(wfa_cnt, text, (text_x, text_y), cv2.FONT_HERSHEY_SIMPLEX, 5.0, (0, 255, 0), 2)
                #     count_outside_range += 1
        gray_segmented_wfa = cv2.cvtColor(wfa_cnt,cv2.COLOR_RGB2GRAY)
        thresh_segmented_wfa = cv2.threshold(gray_segmented_wfa, 80, 255, cv2.THRESH_BINARY)[1] #_INV # | cv2.THRESH_OTSU
        # print("Contours with area greater than 100000:", count_outside_range)
        cv2.imwrite(dst_dir_wfa + img_path.split('.')[0] + '_wfa_seg.tif', thresh_segmented_wfa)
        # approx, contours, shape, contour_img = detect_shape_pnns(img_th_c, wfa_contours)

# wfa contours for 1 image
Image.MAX_IMAGE_PIXELS = None
wfa_img = Image.open(img_A1)
wfa_img.seek(4)
wfa = np.array(wfa_img, dtype = 'uint8')
wfa_c = cv2.cvtColor(wfa,cv2.COLOR_BGR2RGB)
hierachy, img_threshold = cv2.threshold(wfa,  100, 150, cv2.THRESH_BINARY)
img_th_c = cv2.cvtColor(img_threshold,cv2.COLOR_BGR2RGB)
# fig,ax = plt.subplots(figsize = (20,20))
# ax.imshow(img_th_c)
# fig.show()
wfa_contours,_ = cv2.findContours(img_threshold, cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_SIMPLE)
print("found", len(wfa_contours), "in", os.path.basename(img_A1))
wfa_cnt = cv2.drawContours(img_th_c, wfa_contours, -1, (0, 0, 0), 1) # (255, 255, 0) yellow all contours
fig,ax = plt.subplots(figsize = (20,20))
ax.imshow(wfa_cnt)
fig.show()
for cnt in wfa_contours:
    x,y,w,h = cv2.boundingRect(cnt)
    area = cv2.contourArea(cnt)
    if area >=100:
        wfa_cnt = cv2.rectangle(wfa_c, (x,y), (x+w, y+h), (0,255,0), 2)
    elif area<100:
        wfa_cnt = cv2.rectangle(wfa_c, (x,y), (x+w, y+h), (0,0,0), -1)
    # if area<1000:
    #     cla_rect = cv2.rectangle(img_th_c, (x,y), (x+w, y+h), (0,0,255), 2) #green
    # elif area>=10000:
    #     cla_rect = cv2.rectangle(img_th_c, (x,y), (x+50+w+100, y+50+h+100), (255,255,0), 2) #yellow
gray_segmented_wfa = cv2.cvtColor(wfa_cnt,cv2.COLOR_RGB2GRAY)
thresh_segmented_wfa = cv2.threshold(gray_segmented_wfa, np.mean(gray_segmented_wfa), 255, cv2.THRESH_BINARY | cv2.THRESH_OTSU)[1] #_INV
cv2.imwrite(dst_dir_wfa + 'A1_bounding_box_bin.tif', thresh_segmented_wfa)
 #(255, 153, 255)pink
fig,ax = plt.subplots(figsize = (20,20))
ax.imshow(wfa_cnt)
fig.show()
gray_segmented_wfa = cv2.cvtColor(wfa_cnt,cv2.COLOR_RGB2GRAY)
thresh_segmented_wfa = cv2.threshold(gray_segmented_wfa, np.mean(gray_segmented_wfa), 255, cv2.THRESH_BINARY | cv2.THRESH_OTSU)[1] #_INV
# fig,ax = plt.subplots(figsize = (20,20))
# ax.imshow(thresh_segmented_wfa, cmap = 'gray')
# fig.show()
# detecting claudin contours
claudin_img = Image.open(img_C1)
claudin_img.seek(1)
claudin = np.array(claudin_img, dtype = 'uint8')
claudin_c = cv2.cvtColor(claudin,cv2.COLOR_BGR2RGB)
gray = cv2.cvtColor(claudin_c,cv2.COLOR_RGB2GRAY)
_,thresh = cv2.threshold(gray, np.mean(gray), 255, cv2.THRESH_BINARY | cv2.THRESH_OTSU) #_INV
claudin_contours,_ = cv2.findContours(thresh, cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_SIMPLE)
# cla_cnt = cv2.drawContours(img_th_c, claudin_contours, -1, (255, 153, 255), 2) #pink
# print("found", len(claudin_contours))
area_ = []
for cnt in claudin_contours:
    x,y,w,h = cv2.boundingRect(cnt)
    area = cv2.contourArea(cnt)
    area_.append(area)
    if area<1000:
        cla_rect = cv2.rectangle(img_th_c, (x,y), (x+w, y+h), (0,0,255), 2) #green
    elif area>=10000:
        cla_rect = cv2.rectangle(img_th_c, (x,y), (x+50+w+100, y+50+h+100), (255,255,0), 2) #yellow
fig,ax = plt.subplots(figsize = (20,20))
ax.imshow(cla_rect)
fig.show()
cv2.imwrite('/dcs04/lieber/marmaypag/spatialDLPFC_SCZ_LIBD4100/processed-data/RealPNN/capture_area_segmentations/WFA/B1_Claudin_bounding_box_test.tif', cla_rect)
cla_area = np.array(area_)
# drawing claudin contours on the color thresholded image of WFA
# cla_cnt = cv2.drawContours(img_th_c, claudin_contours, -1, (0, 255, 0), 2)
# fig,ax = plt.subplots(figsize = (20,20))
# ax.imshow(cla_cnt)
# fig.show()
# drawing DAPI on the color thresholded image of WFA
dapi_img = Image.open(img_A1)
dapi_img.seek(2)
dapi = np.array(dapi_img, dtype = 'uint8')
dapi_c = cv2.cvtColor(dapi,cv2.COLOR_BGR2RGB)
thresh_c = cv2.cvtColor(thresh_segmented_wfa, cv2.COLOR_BGR2RGB)
gray = cv2.cvtColor(dapi_c,cv2.COLOR_RGB2GRAY)
_,thresh = cv2.threshold(gray, np.mean(gray), 255, cv2.THRESH_BINARY | cv2.THRESH_OTSU) #_INV
dapi_contours,_ = cv2.findContours(thresh, cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_SIMPLE)
print("found", len(dapi_contours), "DAPI")
dapi_cnt = cv2.drawContours(wfa_cnt, dapi_contours, -1, (255, 0, 0), 2) #pink = (255, 153, 255)
for cnt in dapi_contours:
    if wfa_contours == dapi_contours:
        wfa_dapi_masked = cv2.drawContours(img_th_c, wfa_cnt, -1, (0, 0, 0), -1) #pink
fig,ax = plt.subplots(figsize = (20,20))
ax.imshow(dapi_cnt)
fig.show()
cv2.imwrite('/dcs04/lieber/marmaypag/spatialDLPFC_SCZ_LIBD4100/processed-data/RealPNN/capture_area_segmentations/Trials/dapi_wfa_clr.tif', dapi_cnt)
cv2.imwrite('/dcs04/lieber/marmaypag/spatialDLPFC_SCZ_LIBD4100/processed-data/RealPNN/capture_area_segmentations/Claudin/A1_Claudin_bounding_box_test.tif', cont_claudin)
gray_segmented_dapi = cv2.cvtColor(dapi_cnt,cv2.COLOR_RGB2GRAY)
thresh_segmented_dapi = cv2.threshold(gray_segmented_dapi, np.mean(gray_segmented_dapi), 255, cv2.THRESH_BINARY | cv2.THRESH_OTSU)[1] #_INV
fig,ax = plt.subplots(figsize = (20,20))
ax.imshow(thresh_segmented_dapi, cmap = 'gray')
fig.show()
# drawing neun contours
neun_img = Image.open(img_B2)
neun_img.seek(3)
neun = np.array(neun_img, dtype = 'uint8')
neun_c = cv2.cvtColor(neun,cv2.COLOR_BGR2RGB)
gray = cv2.cvtColor(neun_c,cv2.COLOR_RGB2GRAY)
_,thresh = cv2.threshold(gray, np.mean(gray), 255, cv2.THRESH_BINARY | cv2.THRESH_OTSU) #_INV
neun_contours,_ = cv2.findContours(thresh, cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_SIMPLE)
print("found", len(neun_contours), "in", img_path)
final = cv2.drawContours(dapi_cnt, neun_contours, -1, (255, 255, 0), 2) #yellow
fig,ax = plt.subplots(figsize = (20,20))
ax.imshow(final)
fig.show()
cv2.imwrite('/dcs04/lieber/marmaypag/spatialDLPFC_SCZ_LIBD4100/processed-data/RealPNN/capture_area_segmentations/Claudin/claudin_C1.tif', claudin)
cv2.imwrite('/dcs04/lieber/marmaypag/spatialDLPFC_SCZ_LIBD4100/processed-data/RealPNN/capture_area_segmentations/NeuN/B1_NeuN.tif', neun)
# creating a composite of DAPI, NeuN and WFA
import numpy as np
import matplotlib.pyplot as plt
# Load the image channels
channels = []
with Image.open(img_A1) as image:
    num_frames = image.n_frames
    for i in range(num_frames):
        image.seek(i)
        channel = np.array(image, dtype='uint8')
        channels.append(channel)
# Assign color channels
blue_channel = channels[2]  # DAPI assigned to blue
red_channel = channels[3]  # NeuN assigned to red
green_channel = thresh_segmented_wfa  # channels[4] - WFA assigned to green
# Create the composite image
composite = np.zeros(channels[2].shape + (3,), dtype='uint8')
composite[..., 0] = blue_channel  # Blue channel = DAPI
composite[..., 1] = red_channel  # Red channel = NeuN
composite[..., 2] = green_channel  # Green channel = WFA
# save the composite image
cv2.imwrite('/dcs04/lieber/marmaypag/spatialDLPFC_SCZ_LIBD4100/processed-data/RealPNN/capture_area_segmentations/DAPI_NeuN_segmentedWFA_composite.tif', composite)
# Plot the composite image
fig, ax = plt.subplots(figsize=(20, 20))
ax.imshow(composite)
ax.set_axis_off()
plt.show()



##### V12D07-334_C1
Image.MAX_IMAGE_PIXELS = None
for img_path in os.listdir(source_dir):
    if img_path.endswith(".tif") and ('V12D07-334_D1') in img_path:
        wfa_img = Image.open(os.path.join(source_dir, img_path))
        wfa_img.seek(4)
        wfa = np.array(wfa_img, dtype = 'uint8')
        wfa_c = cv2.cvtColor(wfa,cv2.COLOR_BGR2RGB)
        hierachy, img_threshold = cv2.threshold(wfa,  80, 255, cv2.THRESH_BINARY) # 150
        img_th_c = cv2.cvtColor(img_threshold,cv2.COLOR_BGR2RGB)
        wfa_contours,_ = cv2.findContours(img_threshold, cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_SIMPLE)
        print("found", len(wfa_contours), "in", img_path)
        # wfa_cnt = cv2.drawContours(img_th_c, wfa_contours, -1, (0, 0, 0), 1) # yellow all contours
        area_ = []
        count_outside_range = 0
        for cnt in wfa_contours:
            x,y,w,h = cv2.boundingRect(cnt)
            area = cv2.contourArea(cnt)
            area_.append(area) # max(area_) = 1809854.0
            # max, avg = max(area_), (sum(area_)/len(area_)).astype('uint8')
            # if area in [535359.00, 1397823.00, 148415.00]: # area<150 and area>=10000
            #     # print(f"Area {area} satisfies condition 1")
            #     wfa_cnt = cv2.rectangle(img_th_c, (x,y), (x+w, y+h), (0,0,0), -1) # rectangle
            if area >= 100 and area<=100000: #1809854.00
                # print(f"Area {area} satisfies condition 2")
                wfa_cnt = cv2.rectangle(img_th_c, (x,y), (x+w, y+h), (0,0,0), 1) # rectangle
            else:
                # print(f"Area {area} satisfies the else condition")
                wfa_cnt = cv2.rectangle(img_th_c, (x,y), (x+w, y+h), (0,0,0), -1) # rectangle
                # if area >= 100000:
                #     text = f"{area:.2f}"
                #     text_size = cv2.getTextSize(text, cv2.FONT_HERSHEY_SIMPLEX, 5.0, 2)[0]
                #     text_x = x + (w - text_size[0]) // 2
                #     text_y = y + (h + text_size[1]) // 2
                #     cv2.putText(wfa_cnt, text, (text_x, text_y), cv2.FONT_HERSHEY_SIMPLEX, 5.0, (0, 255, 0), 2)
                #     count_outside_range += 1
        gray_segmented_wfa = cv2.cvtColor(wfa_cnt,cv2.COLOR_RGB2GRAY)
        thresh_segmented_wfa = cv2.threshold(gray_segmented_wfa, 50, 255, cv2.THRESH_BINARY)[1] #_INV # | cv2.THRESH_OTSU
        cv2.imwrite(dst_dir_wfa + img_path.split('.')[0] + '__wfa.tif', thresh_segmented_wfa)

fig,ax = plt.subplots(figsize = (20,20))
ax.imshow(wfa_cnt)
fig.show()
