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
from PIL import Image, ImageFont, ImageDraw, ImageEnhance, ImageSequence
import os
import matplotlib.pyplot as plt
import cv2
import scipy
from scipy.spatial.distance import *
from skimage import feature, segmentation, draw, measure, morphology
from skimage.draw import polygon_perimeter
from stitched_functions import read_img, watershed_segmentation
from stitched_functions import *

# directory path
source_dir = pyhere.here('processed-data', 'VistoSeg','captureAreas') #'/dcs04/lieber/marmaypag/spatialDLPFC_SCZ_LIBD4100/processed-data/VistoSeg/captureAreas/'
dst_dir = pyhere.here('processed-data', 'RealPNN', 'capture_area_segmentations', 'AF', 'AF_segmented_binary') #'/dcs04/lieber/marmaypag/spatialDLPFC_SCZ_LIBD4100/processed-data/RealPNN/capture_area_segmentations/AF/AF_segmented_binary/'

# file paths for test images
img_C1 = pyhere.here('processed-data', 'VistoSeg', 'captureAreas','V12F14-057_C1.tif') # /dcs04/lieber/marmaypag/spatialDLPFC_SCZ_LIBD4100/processed-data/VistoSeg/captureAreas/V12F14-057_C1.tif
img_D1 = pyhere.here('processed-data', 'VistoSeg', 'captureAreas','V12F14-057_D1.tif') # /dcs04/lieber/marmaypag/spatialDLPFC_SCZ_LIBD4100/processed-data/VistoSeg/captureAreas/V12F14-057_D1.tif

# directory path
img_dir = pyhere.here('processed-data', 'VistoSeg', 'captureAreas')

img_af, af_shifted, af_gray, af_thresh = read_img.read_and_preprocess(img_D1, 0)
plot_im(af_thresh)
af_labels, af_localmax = watershed_segmentation.find_labels(af_thresh)
afx, afy, afw, afh, af_area, af_segmented = watershed_segmentation.draw_rect_dapi(af_labels, af_gray, img_af)
af_df = save_coordinates.create_df(afx, afy, afw, afh, af_area, img_af, 'autofluorescence')

# segmentations by detecting contours for all images in the directory
for img_path in os.listdir(img_dir):
    if img_path.endswith(".tif"):
        im_af = read_img.read_and_preprocess(img_path, 0)
        print("read", os.path.basename(img_path))
        # af_labels, af_localmax = watershed_segmentation.find_labels(af_thresh)
        af_contours = detect_contours.return_contours(im_af)
        afx, afy, afw, afh, af_area, af_segmented = draw_contours.draw_detected_contours(im_af, 0, af_contours, (255,125,155), 2)
        af_df = save_coordinates.create_df(afx, afy, afw, afh, af_area, img_af, 'autofluorescence')


# find contours for all images in the dir
Image.MAX_IMAGE_PIXELS = None
for img_path in os.listdir(source_dir):
    if img_path.endswith(".tif"):
        af_img = Image.open(os.path.join(source_dir, img_path))
        af_img.seek(0)
        af = np.array(af_img, dtype = 'uint8')
        af_c = cv2.cvtColor(af,cv2.COLOR_BGR2RGB)
        gray = cv2.cvtColor(af_c,cv2.COLOR_RGB2GRAY)
        _,thresh = cv2.threshold(gray, np.mean(gray), 255, cv2.THRESH_BINARY | cv2.THRESH_OTSU) #_INV
        af_contours,_ = cv2.findContours(thresh, cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_SIMPLE)
        print("found", len(af_contours), "in", img_path)
        af_cnt = cv2.drawContours(af_c, af_contours, -1, (255, 0, 0), 2)
        gray_segmented = cv2.cvtColor(af_cnt,cv2.COLOR_RGB2GRAY)
        thresh_segmented = cv2.threshold(gray_segmented, np.mean(gray_segmented), 255, cv2.THRESH_BINARY | cv2.THRESH_OTSU)[1] #_INV
        # afx, afy, afw, afh, af_area, af_segmented = draw_contours.draw_all_contours(af_c, af_contours, (255,125,155), 2)
        # af_df = save_coordinates.create_df(afx, afy, afw, afh, af_area, img_path.split('.')[0], 'autofluorescence')
        # af_df.to_csv(dst_dir + img_path.split('.')[0] + '_info.csv')
        cv2.imwrite(dst_dir + img_path.split('.')[0] + '_af_binarized.tif', thresh_segmented)
