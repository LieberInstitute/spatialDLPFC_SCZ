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
import skimage
from skimage import feature, segmentation, draw, measure, morphology
from skimage.morphology import (erosion,dilation,opening,closing,white_tophat,black_tophat,skeletonize,convex_hull_image)
from skimage.draw import polygon_perimeter
import itertools
from itertools import product
from collections import defaultdict
from shapely.geometry import Point
from shapely.geometry.polygon import Polygon

#1 - get the manual annotations image
csv_test = '/dcs04/lieber/marmaypag/spatialDLPFC_SCZ_LIBD4100/processed-data/Experimentation_archive/2_MockPNN/Training_tiles/Manual_annotations/Annotations/20220712_VIF_MockPNN_Strong_Scan1_[7851,52577]_component_data_20.csv'
img_test = '/dcs04/lieber/marmaypag/spatialDLPFC_SCZ_LIBD4100/raw-data/images/2_MockPNN/Training_tiles_raw_no_annotations/20220712_VIF_MockPNN_Strong_Scan1_[7851,52577]_component_data_20.tif'

#2a - segment all the single tiles from both ntc and scz samples and save them in a folder (then later, all the manual annoations can be overlaid on the segmented tiles)
Image.MAX_IMAGE_PIXELS = None
source_dir = '/dcs04/lieber/marmaypag/spatialDLPFC_SCZ_LIBD4100/raw-data/images/2_MockPNN/Training_tiles_raw_no_annotations/' # folder where all the selected tiles exist
dst_dir = '/dcs04/lieber/marmaypag/spatialDLPFC_SCZ_LIBD4100/raw-data/images/2_MockPNN/Training_tiles_segmented_CV_IoU_test/' # where all the segmented tiles are located

# segment all the individual tiles
for img_path in os.listdir(source_dir):
    if img_path.endswith(".tif"):
        # print(img_path)
        wfa_img = Image.open(os.path.join(source_dir, img_path))
        wfa_img.seek(3)
        # wfa = np.array(wfa_img, dtype = 'uint8')
        wfa = cv2.normalize(np.array(wfa_img, dtype = 'float32'), np.zeros(np.array(wfa_img, dtype = 'float32').shape, np.double), 1.0, 0.0, cv2.NORM_MINMAX)
        # cv2.imwrite(dst_dir + img_path.split('.')[0] + '_raw.tif', wfa)
        wfa_c = cv2.cvtColor(np.array(wfa*255, dtype = np.uint8),cv2.COLOR_BGR2RGB)
        # cv2.imwrite(dst_dir + img_path.split('.')[0] + '_colored.tif', wfa_c)
        # adjusted = cv2.convertScaleAbs(wfa, alpha=0.3, beta=10) # decreased the contrast of the original image for better segmentation
        wfa_gry = skimage.color.rgb2gray(wfa_c)
        hierachy, img_threshold = cv2.threshold(np.array(wfa_gry*255, dtype = np.uint8),  100, 255, cv2.THRESH_BINARY) # 150
        img_th_c = cv2.cvtColor(img_threshold,cv2.COLOR_BGR2RGB)
        img_threshold_gry = cv2.cvtColor(img_th_c, cv2.COLOR_BGR2GRAY)
        # fig,ax = plt.subplots(figsize = (20,20))
        # ax.imshow(img_threshold, cmap = 'gray')
        # fig.show()
        # cv2.imwrite(dst_dir + img_path.split('.')[0] + '_thresholded.tif', img_threshold)
        wfa_contours,_ = cv2.findContours(img_threshold_gry, cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_SIMPLE)
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
            if area<500: # ; area>=3000 and area>=3000
                # print(f"Area {area} satisfies condition 1")
                wfa_cnt = cv2.rectangle(img_th_c, (x,y), (x+w, y+h), (0,0,0), -1) # rectangle
            if area>=3000:
                wfa_cnt = cv2.rectangle(img_th_c, (x,y), (x+w, y+h), (0,0,0), -1) # rectangle
            elif 300 <= area <= 3000:
                # print(f"Area {area} satisfies condition 2")
                wfa_cnt = cv2.rectangle(img_th_c, (x,y), (x+w, y+h), (0,0,255), 2) # red rectangle
        # gray_segmented_wfa = cv2.cvtColor(img_th_c,cv2.COLOR_RGB2GRAY)
        # thresh_segmented_wfa = cv2.threshold(gray_segmented_wfa, 80, 255, cv2.THRESH_BINARY)[1] #_INV # | cv2.THRESH_OTSU
        cv2.imwrite(dst_dir + img_path.split('.')[0] + '_wfa__seg.tif', wfa_cnt)

# match the tile numbers from annotation images to segmented images and overlay the boxes
# draw a rectangle from the manual annotations csv on the contour detected image --> works to overlay annotations over segmentations
def draw_rect(df_manual_test, contour_img):
    for box in range(len(df_manual_test['x1'])):
        # print(box)
        rect = cv2.rectangle(contour_img, (df_manual_test['x1'][box], df_manual_test['y1'][box]), (df_manual_test['x4'][box], df_manual_test['y4'][box]), (0,255,0), 3)
    # fig,ax = plt.subplots(figsize = (20,20))
    # ax.imshow(out_img1)
    # fig.show()
    return contour_img

# Get a list of filenames in each folder
folder1_path = '/dcs04/lieber/marmaypag/spatialDLPFC_SCZ_LIBD4100/raw-data/images/2_MockPNN/Training_tiles_segmented_CV_IoU_test/'
folder2_path = '/dcs04/lieber/marmaypag/spatialDLPFC_SCZ_LIBD4100/processed-data/Experimentation_archive/2_MockPNN/Training_tiles/Manual_annotations/Annotations/'


# Create a dictionary to store the filenames with matching parts
matching_files = {}

# Iterate through the filenames in both folders
for filename1 in os.listdir(folder1_path):
    for filename2 in os.listdir(folder2_path):
        # Extract the part enclosed in square brackets from both filenames
        part1 = filename1.split('[')[-1].split(']')[0]
        part2 = filename2.split('[')[-1].split(']')[0]
        # If the extracted parts match, add the filenames to the dictionary
        if part1 == part2:
            matching_files[part1] = (filename1, filename2)

# folder to save all the overlaid images
dst_dir_overlay = '/dcs04/lieber/marmaypag/spatialDLPFC_SCZ_LIBD4100/raw-data/images/2_MockPNN/Training_tiles_segmented_overlaid_annotations/'

# Process the matched files
for part, (file1, file2) in matching_files.items():
    # Perform operations on the matched files here
    print(f"Matching part: {part}")
    print(f"Matching file in folder 1: {file1}")
    print(f"Matching file in folder 2: {file2}")
    # Preprocessing the manual annotations csv
    conv_factor = 2.0112375738 # the fiji annotations are measured in microns and they need to be translated to pixels (1860/924.81)
    df_manual_test = pd.read_csv(os.path.join(folder2_path,file2)) # read the manual annotations csv into dataframe
    df_manual_test = df_manual_test.rename(columns = {'X': 'xc', 'Y': 'yc', 'BX': 'x1', 'BY': 'y1', 'Perim.': 'Perimeter'}) # xc,yc are the centroids of the BB
    df_manual_test.loc[:,['xc']], df_manual_test.loc[:,['yc']], df_manual_test.loc[:,['x1']], df_manual_test.loc[:,['y1']], df_manual_test['Width'], df_manual_test['Height'] = df_manual_test['xc']*conv_factor, df_manual_test['yc']*conv_factor, df_manual_test['x1']*conv_factor, df_manual_test['y1']*conv_factor, df_manual_test['Width']*conv_factor, df_manual_test['Height']*conv_factor
    df_manual_test['x2'] = (df_manual_test['x1'] + df_manual_test['Width'])
    df_manual_test['y2'], df_manual_test['x3'] = df_manual_test['y1'], df_manual_test['x1']
    df_manual_test['y3'] = (df_manual_test['y1'] + df_manual_test['Height'])
    df_manual_test['x4'], df_manual_test['y4']  = df_manual_test['x2'], df_manual_test['y3'] # Calculating all 4 coordinates of the BB
    df_manual_test['xc'], df_manual_test['yc'], df_manual_test['x1'], df_manual_test['y1'] = np.int0(np.ceil(df_manual_test['xc'])), np.int0(np.ceil(df_manual_test['yc'])), np.int0(np.ceil(df_manual_test['x1'])), np.int0(np.ceil(df_manual_test['y1'])) # convert x,y,bx,by from floating point to integers (doing it after, reduces round off errors)
    df_manual_test['x2'], df_manual_test['y2'], df_manual_test['x3'], df_manual_test['y3'], df_manual_test['x4'], df_manual_test['y4'] = np.int0(np.ceil(df_manual_test['x2'])), np.int0(np.ceil(df_manual_test['y2'])), np.int0(np.ceil(df_manual_test['x3'])), np.int0(np.ceil(df_manual_test['y3'])), np.int0(np.ceil(df_manual_test['x4'])), np.int0(np.ceil(df_manual_test['y4']))
    df_manual_test = df_manual_test[['Area', 'Perimeter', 'Mean', 'Min', 'Max', 'xc', 'yc', 'x1', 'y1', 'x2', 'y2', 'x3', 'y3' , 'x4' , 'y4', 'Width', 'Height', 'Ch']] # rearranging the columns
    segmented_img = Image.open(os.path.join(folder1_path,file1))
    overlay_img = draw_rect(df_manual_test, (np.array(segmented_img, dtype = 'uint8')))
    cv2.imwrite(dst_dir_overlay + file1.split('.')[0] + '_overlay_rect.tif', overlay_img) # actual = green box, predicted = red box


    

#2b - segment br5182 = ntc and br2039 = scz using the CV algorithm (the whole tissue section) and save it
#3 - overlay the manual annotation boxes on the segmented images, just on the wfa channel
#4 - this gives figure a,b,c
#5 - for the math,
#6 - find the number of overlaps in all 13 ntc and 12 scz tiles
#7 - find the number of overlaps in the zero PNN tiles as well
#8 - find the IoUs for all of the overlaps
#9 - find the avg of IoUs 
#10 - find how many are aboe 80% or high/low? ==>do this only based on the previous step (avgs)
#11- visulize the low, medium and high IoU ones
