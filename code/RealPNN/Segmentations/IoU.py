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

#1 - get the manual annotations image
csv_test = '/dcs04/lieber/marmaypag/spatialDLPFC_SCZ_LIBD4100/processed-data/Experimentation_archive/2_MockPNN/Training_tiles/Manual_annotations/Annotations/20220712_VIF_MockPNN_Strong_Scan1_[7851,52577]_component_data_20.csv'
img_test = '/dcs04/lieber/marmaypag/spatialDLPFC_SCZ_LIBD4100/raw-data/images/2_MockPNN/Training_tiles_raw_no_annotations/20220712_VIF_MockPNN_Strong_Scan1_[7851,52577]_component_data_20.tif'

#2a - segment all the single tiles from both ntc and scz samples and save them in a folder (then later, all the manual annoations can be overlaid on the segmented tiles)
Image.MAX_IMAGE_PIXELS = None
source_dir = '/dcs04/lieber/marmaypag/spatialDLPFC_SCZ_LIBD4100/raw-data/images/2_MockPNN/Training_tiles_raw_no_annotations/' # folder where all the selected tiles exist
dst_dir = '/dcs04/lieber/marmaypag/spatialDLPFC_SCZ_LIBD4100/raw-data/images/2_MockPNN/Training_tiles_segmented_CV_IoU_test/' # where all the segmented tiles are located
dst_dir_csv = '/dcs04/lieber/marmaypag/spatialDLPFC_SCZ_LIBD4100/raw-data/images/2_MockPNN/Training_tiles_cv_segmentation_csv_files/' # where all the segmented tiles info is located

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
        wfx, wfy, wfw, wfh, pnn_area = [], [], [], [], []
        for cnt in wfa_contours:
            x,y,w,h = cv2.boundingRect(cnt)
            area = cv2.contourArea(cnt)
            area_.append(area) # max(area_) = 1809854.0
            # print(area)
            # max, avg = max(area_), (sum(area_)/len(area_)).astype('uint8')
            if area<500: # ; area>=3000 and area>=3000
                # print(f"Area {area} satisfies condition 1")
                wfa_cnt = cv2.rectangle(img_th_c, (x,y), (x+w, y+h), (0,0,0), -1) # black rectangle
            if area>=3000:
                wfa_cnt = cv2.rectangle(img_th_c, (x,y), (x+w, y+h), (0,0,0), -1) # black rectangle
            elif 300 <= area <= 3000:
                # print(f"Area {area} satisfies condition 2")
                wfa_cnt = cv2.rectangle(img_th_c, (x,y), (x+w, y+h), (0,0,255), 2) # red rectangle
                wfx.append(x)
                wfy.append(y)
                wfw.append(w)
                wfh.append(h)
                pnn_area.append(area)
        # gray_segmented_wfa = cv2.cvtColor(img_th_c,cv2.COLOR_RGB2GRAY)
        # thresh_segmented_wfa = cv2.threshold(gray_segmented_wfa, 80, 255, cv2.THRESH_BINARY)[1] #_INV # | cv2.THRESH_OTSU
        cv2.imwrite(dst_dir + img_path.split('.')[0] + '_wfa__seg.tif', wfa_cnt)
        # Populate the data in the dataframe
        col_names = ['img_file_name','type_of_object', 'x1', 'y1', 'Width', 'Height', 'Area', 'total_number_pnns']
        object_name = 'PNN' # name of the objects stored in the dataframe
        file_name = os.path.basename(img_test) # image file name
        dict = {col_names[0]: file_name, col_names[1]: object_name, col_names[2]: wfx, col_names[3]: wfy, col_names[4]: wfw, col_names[5]: wfh, col_names[6]: pnn_area, col_names[7]: len(wfx)}
        df_wfa_cv = pd.DataFrame(dict, columns = col_names) 
        df_wfa_cv['x2'] = df_wfa_cv['x1'] + df_wfa_cv['Width']
        df_wfa_cv['y2'], df_wfa_cv['x3'] = df_wfa_cv['y1'], df_wfa_cv['x1']
        df_wfa_cv['y3'] = df_wfa_cv['y1'] + df_wfa_cv['Height']
        df_wfa_cv['x4'], df_wfa_cv['y4'] = df_wfa_cv['x2'], df_wfa_cv['y3']
        df_wfa_cv = df_wfa_cv[['img_file_name', 'type_of_object', 'x1', 'y1', 'x2', 'y2', 'x3', 'y3', 'x4', 'y4', 'Width', 'Height', 'Area', 'total_number_pnns']]
        df_wfa_cv.to_csv(dst_dir_csv + img_path.split('.')[0] + '_wfa_seg.csv')


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
annotations_pixels = '/dcs04/lieber/marmaypag/spatialDLPFC_SCZ_LIBD4100/processed-data/Experimentation_archive/2_MockPNN/Training_tiles/Manual_annotations/Annotations_in_pixels/'

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
    df_manual_test.to_csv(annotations_pixels + file1.split('.')[0] + '_manual_annotations.csv')
    # segmented_img = Image.open(os.path.join(folder1_path,file1))
    # overlay_img = draw_rect(df_manual_test, (np.array(segmented_img, dtype = 'uint8')))
    # cv2.imwrite(dst_dir_overlay + file1.split('.')[0] + '_overlay_rect.tif', overlay_img) # actual = green box, predicted = red box


#2b - segment br5182 = ntc and br2039 = scz using the CV algorithm (the whole tissue section) and save it
#3 - overlay the manual annotation boxes on the segmented images, just on the wfa channel
#4 - this gives figure a,b,c
#5 - for the math,
#6 - find the number of overlaps in all 13 ntc and 12 scz tiles
def calculate_iou(box1, box2):
    # Calculate the coordinates of the intersection rectangle
    x1 = max(box1[0], box2[0])
    y1 = max(box1[1], box2[1])
    x2 = min(box1[2], box2[2])
    y2 = min(box1[3], box2[3])
    # Calculate the area of the intersection rectangle
    intersection_area = max(0, x2 - x1) * max(0, y2 - y1)
    # Calculate the areas of both bounding boxes
    area_box1 = (box1[2] - box1[0]) * (box1[3] - box1[1])
    area_box2 = (box2[2] - box2[0]) * (box2[3] - box2[1])
    # Calculate the IoU
    iou = intersection_area / float(area_box1 + area_box2 - intersection_area)
    return iou

# Define the list of ground truth bounding boxes and predicted bounding boxes
ground_truth_boxes = [(x1, y1, x2, y2), (x1, y1, x2, y2), ...]  # Replace with your actual values
predicted_boxes = [(x1, y1, x2, y2), (x1, y1, x2, y2), ...]  # Replace with your actual values

# match the file names from the ground truth and predicted boxes 
ground_truth_folder = '/dcs04/lieber/marmaypag/spatialDLPFC_SCZ_LIBD4100/processed-data/Experimentation_archive/2_MockPNN/Training_tiles/Manual_annotations/Annotations/'
predicted_folder = '/dcs04/lieber/marmaypag/spatialDLPFC_SCZ_LIBD4100/raw-data/images/2_MockPNN/Training_tiles_cv_segmentation_csv_files/'


matching_files_boxes = {}

# Iterate through the filenames in both folders
for filename1 in os.listdir(annotations_pixels):
    for filename2 in os.listdir(predicted_folder):
        # Extract the part enclosed in square brackets from both filenames
        part1 = filename1.split('[')[-1].split(']')[0]
        part2 = filename2.split('[')[-1].split(']')[0]
        # If the extracted parts match, add the filenames to the dictionary
        if part1 == part2:
            matching_files[part1] = (filename1, filename2)

# Process the matched files
for part, (file1, file2) in matching_files_boxes.items():
    # Perform operations on the matched files here
    print(f"Matching part: {part}")
    print(f"Matching file in folder 1: {file1}")
    print(f"Matching file in folder 2: {file2}")
    df_gt = pd.read_csv(os.path.join(annotations_pixels,file1.split('.')[0] + '_manual_annotations.csv'))
    df_pred = pd.read_csv(os.path.join(predicted_folder, file2.split('.')[0] + 'wfa_seg.csv'))
    ground_truth_boxes = df_gt[['x1', 'y1', 'x2', 'y2']].values.tolist()
    predicted_boxes = df_pred[['x1', 'y1', 'x2', 'y2']].values.tolist()
    # Initialize a matrix to store the IoU values
    iou_matrix = np.zeros((len(ground_truth_boxes), len(predicted_boxes)))
    # Calculate IoU for all pairs of ground truth and predicted bounding boxes
    for i, gt_box in enumerate(ground_truth_boxes):
        for j, pred_box in enumerate(predicted_boxes):
            iou_matrix[i, j] = calculate_iou(gt_box, pred_box)


# Initialize a matrix to store the IoU values
iou_matrix = np.zeros((len(ground_truth_boxes), len(predicted_boxes)))

# Calculate IoU for all pairs of ground truth and predicted bounding boxes
for i, gt_box in enumerate(ground_truth_boxes):
    for j, pred_box in enumerate(predicted_boxes):
        iou_matrix[i, j] = calculate_iou(gt_box, pred_box)

# Count how many ground truth boxes overlap with each predicted box
gt_matches = np.sum(iou_matrix > threshold, axis=0)

# Count how many predicted boxes overlap with each ground truth box
pred_matches = np.sum(iou_matrix > threshold, axis=1)

# You can then access the counts of matches for each predicted and ground truth box.
# gt_matches[i] contains the number of ground truth boxes that overlap with predicted box i.
# pred_matches[j] contains the number of predicted boxes that overlap with ground truth box j.

#7 - find the number of overlaps in the zero PNN tiles as well
#8 - find the IoUs for all of the overlaps
#9 - find the avg of IoUs 
#10 - find how many are aboe 80% or high/low? ==>do this only based on the previous step (avgs)
#11- visulize the low, medium and high IoU ones
