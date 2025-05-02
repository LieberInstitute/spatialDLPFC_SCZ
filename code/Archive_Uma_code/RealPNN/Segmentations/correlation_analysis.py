import numpy as np
import pandas as pd
import PIL
import pyhere
from pathlib import Path
from PIL import Image, ImageFont, ImageDraw, ImageEnhance, ImageSequence
import os
import matplotlib
import matplotlib.pyplot as plt
import sys
import cv2
import math
import scipy
from scipy.spatial.distance import *
from scipy.signal import convolve2d
from scipy.optimize import linear_sum_assignment
import skimage
from skimage import feature, segmentation, draw, measure, morphology
from skimage.morphology import (erosion,dilation,opening,closing,white_tophat,black_tophat,skeletonize,convex_hull_image)
from skimage.draw import polygon_perimeter
import itertools
from itertools import product
from collections import defaultdict
from shapely.geometry import Point
from shapely.geometry.polygon import Polygon


#1 - get the manual annotations image <-- run this
csv_test = '/dcs04/lieber/marmaypag/spatialDLPFC_SCZ_LIBD4100/processed-data/Experimentation_archive/2_MockPNN/Training_tiles/Manual_annotations/Annotations/20220712_VIF_MockPNN_Strong_Scan1_[7851,52577]_component_data_20.csv'
img_test = '/dcs04/lieber/marmaypag/spatialDLPFC_SCZ_LIBD4100/raw-data/images/2_MockPNN/Training_tiles_raw_no_annotations/20220712_VIF_MockPNN_Strong_Scan1_[7851,52577]_component_data_20.tif'

#2a - segment all the single tiles from both ntc and scz samples and save them in a folder 
# (then later, all the manual annoations can be overlaid on the segmented tiles) <-- run this
Image.MAX_IMAGE_PIXELS = None
source_dir = '/dcs04/lieber/marmaypag/spatialDLPFC_SCZ_LIBD4100/raw-data/images/2_MockPNN/Training_tiles_raw_no_annotations/' # folder where all the selected tiles exist
dst_dir = '/dcs04/lieber/marmaypag/spatialDLPFC_SCZ_LIBD4100/raw-data/images/2_MockPNN/Training_tiles_segmented_CV_IoU_test/' # where all the segmented tiles are located
dst_dir_csv = '/dcs04/lieber/marmaypag/spatialDLPFC_SCZ_LIBD4100/raw-data/images/2_MockPNN/Training_tiles_cv_segmentation_csv_files/' # where all the segmented tiles info is located


##### testing on file 2 of the manual annotations, # of annotations = 8
csv_file2 = '/dcs04/lieber/marmaypag/spatialDLPFC_SCZ_LIBD4100/processed-data/Experimentation_archive/2_MockPNN/Training_tiles/Manual_annotations/Annotations_in_pixels/20220712_VIF_MockPNN_Strong_Scan1_[6925,51188]_component_data_02_wfa__seg_manual_annotations.csv'
img_file2_ma = '/dcs04/lieber/marmaypag/spatialDLPFC_SCZ_LIBD4100/processed-data/Experimentation_archive/2_MockPNN/Training_tiles/Manual_annotations/Images/20220712_VIF_MockPNN_Strong_Scan1_[6925,51188]_component_data_02.tif'
img_file2_raw = '/dcs04/lieber/marmaypag/spatialDLPFC_SCZ_LIBD4100/raw-data/images/2_MockPNN/Training_tiles_raw_no_annotations/20220712_VIF_MockPNN_Strong_Scan1_[6925,51188]_component_data_02.tif'
csv_file2_cv = '/dcs04/lieber/marmaypag/spatialDLPFC_SCZ_LIBD4100/raw-data/images/2_MockPNN/Training_tiles_cv_segmentation_csv_files/20220712_VIF_MockPNN_Strong_Scan1_[6925,51188]_component_data_02wfa_seg.csv'

#### second test on file 21, # of annotations = 42
csv_file21_cv = '/dcs04/lieber/marmaypag/spatialDLPFC_SCZ_LIBD4100/raw-data/images/2_MockPNN/Training_tiles_cv_segmentation_csv_files/20220712_VIF_MockPNN_Strong_Scan1_[8235,50974]_component_data_21wfa_seg.csv'
img_file21_raw = '/dcs04/lieber/marmaypag/spatialDLPFC_SCZ_LIBD4100/raw-data/images/2_MockPNN/Training_tiles_raw_no_annotations/20220712_VIF_MockPNN_Strong_Scan1_[8235,50974]_component_data_21.tif'
csv_file21_ma = '/dcs04/lieber/marmaypag/spatialDLPFC_SCZ_LIBD4100/processed-data/Experimentation_archive/2_MockPNN/Training_tiles/Manual_annotations/Annotations_in_pixels/20220712_VIF_MockPNN_Strong_Scan1_[8235,50974]_component_data_21_wfa__seg_manual_annotations.csv'

##### third test on file 20, # of annotations = 9
csv_file20_cv = '/dcs04/lieber/marmaypag/spatialDLPFC_SCZ_LIBD4100/raw-data/images/2_MockPNN/Training_tiles_cv_segmentation_csv_files/20220712_VIF_MockPNN_Strong_Scan1_[7851,52577]_component_data_20wfa_seg.csv'
img_file20_raw = '/dcs04/lieber/marmaypag/spatialDLPFC_SCZ_LIBD4100/raw-data/images/2_MockPNN/Training_tiles_raw_no_annotations/20220712_VIF_MockPNN_Strong_Scan1_[7851,52577]_component_data_20.tif'
csv_file20_ma = '/dcs04/lieber/marmaypag/spatialDLPFC_SCZ_LIBD4100/processed-data/Experimentation_archive/2_MockPNN/Training_tiles/Manual_annotations/Annotations_in_pixels/20220712_VIF_MockPNN_Strong_Scan1_[7851,52577]_component_data_20_wfa__seg_manual_annotations.csv'


#### test on all files
csv_cv = 
img_raw = 
csv_ma = 


# read the files
file2_df = pd.read_csv(csv_file20)
file2_cv = pd.read_csv(csv_file20_cv)
image = Image.open(img_file20_raw)
image.seek(3)
image = cv2.normalize(np.array(image, dtype = 'float32'), np.zeros(np.array(image, dtype = 'float32').shape, np.double), 1.0, 0.0, cv2.NORM_MINMAX)
wfa_c = cv2.cvtColor(np.array(image*255, dtype = np.uint8),cv2.COLOR_BGR2RGB)

# Compute mean pixel intensities for each bounding box ---> doing this in a grayscale image
for index, row in file2_df.iterrows():
    x1, y1, x4, y4 = int(row['x1']), int(row['y1']), int(row['x4']), int(row['y4'])
    # Extract region of interest (ROI)
    roi = image[y1:y4, x1:x4]
    # Draw rectangle on the image
    color = (255, 255, 255)  # white color for the rectangle
    thickness = 2
    cv2.rectangle(image, (x1, y1), (x4, y4), color, thickness)
    # Compute mean pixel intensities
    mean_intensity = cv2.mean(roi)
    print(f"Bounding Box {index + 1}: Mean Intensity - {mean_intensity}")
    # Add text with mean intensity on top of the rectangle
    text = f"{mean_intensity[0]:.2f}"  # Assuming grayscale image, adjust if using RGB
    font = cv2.FONT_HERSHEY_SIMPLEX
    font_scale = 0.8
    font_thickness = 2
    text_size = cv2.getTextSize(text, font, font_scale, font_thickness)[0]
    # text_position = ((x1 + x4 - text_size[0]) // 2, (y1 + y4 + text_size[1]) // 2)
    text_position = ((x1 + x4 - text_size[0]) // 2, y1 - 5)
    cv2.putText(image, text, text_position, font, font_scale, color, font_thickness)

cv2.imwrite('/dcs04/lieber/marmaypag/spatialDLPFC_SCZ_LIBD4100/raw-data/images/2_MockPNN/Test/box_file2.tif', image)

# read the files
file2_df = pd.read_csv(csv_file20_ma)
file2_cv = pd.read_csv(csv_file20_cv)
image = Image.open(img_file20_raw)
image.seek(3)
image = cv2.normalize(np.array(image, dtype = 'float32'), np.zeros(np.array(image, dtype = 'float32').shape, np.double), 1.0, 0.0, cv2.NORM_MINMAX)
wfa_c = cv2.cvtColor(np.array(image*255, dtype = np.uint8),cv2.COLOR_BGR2RGB)


### doing the same thing on RGB images
for index, row in file2_df.iterrows():
    x1, y1, x4, y4 = int(row['x1']), int(row['y1']), int(row['x4']), int(row['y4'])
    # Extract region of interest (ROI)
    roi = image[y1:y4, x1:x4]
    # Compute mean pixel intensities
    mean_intensity = cv2.mean(roi)[0]
    print(f"Bounding Box {index + 1}: Mean Intensity - {mean_intensity}")
    # Draw rectangle on the image based on mean intensity
    if mean_intensity < 0.3:
        color = (0, 255, 255)  #### change this and make it draw a black box over the lower mean intensities
        color_name = "Yellow" # Yellow color for mean intensity <= 0.3
    else:
        color = (0, 255, 0)  # Green color for other mean intensity values
        color_name = "green"
    wfa_c = cv2.rectangle(wfa_c, (x1, y1), (x4, y4), color, 2)
    # Add text with mean intensity on top of the rectangle
    text_intensity = f"{mean_intensity:.2f}"
    text_index = f"Box {index + 1}"   
    text_combined = f"{text_intensity}, {text_index}" 
    font = cv2.FONT_HERSHEY_SIMPLEX
    font_scale = 0.8
    font_thickness = 2
    color_txt = (255, 255, 0) # cyan text
    text_size = cv2.getTextSize(text_combined, font, font_scale, font_thickness)[0]
    # text_position = ((x1 + x4 - text_size[0]) // 2, (y1 + y4 + text_size[1]) // 2)
    text_position = ((x1 + x4 - text_size[0]) // 2, y1 - 5)
    wfa_c = cv2.putText(wfa_c, text_combined, text_position, font, font_scale, color_txt, font_thickness)
    print(f"Bounding Box {index + 1}: Mean Intensity - {mean_intensity:.2f}, Color - {color_name}")


# Draw rectangles for CV bounding boxes
for _, row in file2_cv.iterrows():
    x1, y1, x4, y4 = int(row['x1']), int(row['y1']), int(row['x4']), int(row['y4'])    
    # Draw rectangle on the image
    color = (0, 0, 255)  # Red color for the CV BB
    thickness = 2
    wfa_c = cv2.rectangle(wfa_c, (x1, y1), (x4, y4), color, thickness)


cv2.imwrite('/dcs04/lieber/marmaypag/spatialDLPFC_SCZ_LIBD4100/raw-data/images/2_MockPNN/Test/box_box_yellow_file20.tif', wfa_c)
