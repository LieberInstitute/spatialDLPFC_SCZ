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
from skimage import feature, segmentation, draw, measure, morphology
from stitched_functions import read_img, watershed_segmentation
from stitched_functions import *

# directory path
Image.MAX_IMAGE_PIXELS = None # increase the max image pixels to avoid decompression error
source_dir = '/dcs04/lieber/marmaypag/spatialDLPFC_SCZ_LIBD4100/processed-data/VistoSeg/captureAreas/'
dst_dir = '/dcs04/lieber/marmaypag/spatialDLPFC_SCZ_LIBD4100/processed-data/RealPNN/capture_area_segmentations/DAPI/'

# test file paths
img_A1 = pyhere.here('processed-data', 'VistoSeg', 'captureAreas','V12F14-053_A1.tif')

# test image watershed
img_dapi, dapi_shifted, dapi_gray, dapi_thresh = read_img.read_and_preprocess(img_D1, 2)
dapi_labels, dapi_localmax = find_labels(dapi_thresh) #watershed_segmentation.
dpx, dpy, dpw, dph, dp_area, dapi_segmented = draw_rect_from_labels(dapi_labels, dapi_gray, img_dapi)
cv2.imwrite('/users/ukaipa/PNN/One_img/dapi_stitched_segmented_D1_1148569.tif', dapi_segmented)
dapi_df = save_coordinates.create_df(dpx, dpy, dpw, dph, dp_area, img_dapi, 'dapi')


# detect contours on test image

dapi_img = Image.open(img_C2)
dapi_img.seek(2)
dapi = np.array(dapi_img, dtype = 'uint8') # (17799, 16740)
dapi_c = cv2.cvtColor(dapi,cv2.COLOR_BGR2RGB)
gray = cv2.cvtColor(dapi_c,cv2.COLOR_RGB2GRAY)
_,thresh = cv2.threshold(gray, np.mean(gray), 255, cv2.THRESH_BINARY | cv2.THRESH_OTSU) #_INV
contours,_ = cv2.findContours(thresh, cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_SIMPLE)
print(len(contours))
dp_cnt = cv2.drawContours(dapi_c, contours, -1, (0, 255, 0), 2)
fig,ax = plt.subplots(figsize = (20,20))
ax.imshow(thresh, cmap = 'gray') #
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
        contours,_ = cv2.findContours(thresh, cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_SIMPLE)
        print("found", len(contours), "in", img_path)
        dp_cnt = cv2.drawContours(dapi_c, contours, -1, (0, 255, 0), 2)
        cv2.imwrite(dst_dir + img_path + '_dapi_contours_segmented.tif', dp_cnt)






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
edges_all_images(source_dir, dst_dir)
