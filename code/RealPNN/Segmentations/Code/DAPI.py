'''
For Stitched Visium-IF tissue sections from VistoSeg SplitSlide output
Channel0 = AF
Channel1 = Claudin - 5 (Alex 488),
Channel2 = DAPI,
Channel3 = NeuN,
Channel4 = WFA
'''


from scipy import ndimage
import imutils
import numpy as np
import argparse
from argparse import ArgumentParser
import pyhere
from pyhere import here
from pathlib import Path
import pandas as pd
import PIL
from PIL import Image
import os
import matplotlib.pyplot as plt
import cv2
import scipy
from scipy.spatial.distance import *
import skimage
from skimage import feature, segmentation, draw, measure, morphology
from stitched_functions import read_img, watershed_segmentation, draw_contours
from stitched_functions import draw_contours, save_coordinates

# directory path
Image.MAX_IMAGE_PIXELS = None # increase the max image pixels to avoid decompression error
source_dir = '/dcs04/lieber/marmaypag/spatialDLPFC_SCZ_LIBD4100/processed-data/VistoSeg/captureAreas/'
dst_dir_dapi = '/dcs04/lieber/marmaypag/spatialDLPFC_SCZ_LIBD4100/processed-data/RealPNN/capture_area_segmentations/DAPI/DAPI_binarized/'

# test file paths
img_A1 = pyhere.here('processed-data', 'VistoSeg', 'captureAreas','V12D07-334_A1.tif')
img_B1 = '/dcs04/lieber/marmaypag/spatialDLPFC_SCZ_LIBD4100/processed-data/VistoSeg/captureAreas/V12F14-053_B1.tif'


# find contours for all images in the dir
Image.MAX_IMAGE_PIXELS = None
for img_path in os.listdir(source_dir):
    if img_path.endswith(".tif") and ('V12D07') in img_path:
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
            if area >=100:
                dp_cnt = cv2.rectangle(dapi_c, (x,y), (x+w, y+h), (255,0,0), 1)
            elif area <50:
                dp_cnt = cv2.rectangle(dapi_c, (x,y), (x+w, y+h), (0,0,0), -1)
        gray_segmented = cv2.cvtColor(dp_cnt,cv2.COLOR_RGB2GRAY)
        thresh_segmented = cv2.threshold(gray_segmented, np.mean(gray_segmented), 255, cv2.THRESH_BINARY | cv2.THRESH_OTSU)[1] #_INV
        cv2.imwrite(dst_dir_dapi + img_path.split('.')[0] + '_dapi_binarized.tif', thresh_segmented)



# For one test image
Image.MAX_IMAGE_PIXELS = None
dapi_img = Image.open(img_A1)
dapi_img.seek(2)
dapi = np.array(dapi_img, dtype = 'uint8')
dapi_c = cv2.cvtColor(dapi,cv2.COLOR_BGR2RGB)
gray = cv2.cvtColor(dapi_c,cv2.COLOR_RGB2GRAY)
_,thresh = cv2.threshold(gray, np.mean(gray), 255, cv2.THRESH_BINARY | cv2.THRESH_OTSU) #_INV
dapi_contours,_ = cv2.findContours(thresh, cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_SIMPLE)
print("found", len(dapi_contours))
# dp_cnt = cv2.drawContours(dapi_c, dapi_contours, -1, (0, 255, 0), 2)
for cnt in dapi_contours:
   x,y,w,h = cv2.boundingRect(cnt)
   area = cv2.contourArea(cnt)
   if area >=50:
      dp_cnt = cv2.rectangle(dapi_c, (x,y), (x+w, y+h), (255,0,0), 1)
   elif area<50:
      dp_cnt = cv2.rectangle(dapi_c, (x,y), (x+w, y+h), (0,0,0), -1)

gray_segmented = cv2.cvtColor(dp_cnt,cv2.COLOR_RGB2GRAY)
thresh_segmented = cv2.threshold(gray_segmented, np.mean(gray_segmented), 255, cv2.THRESH_BINARY | cv2.THRESH_OTSU)[1] #_INV
fig,ax = plt.subplots(figsize = (20,20))
ax.imshow(thresh_segmented, cmap='gray')
fig.show()

