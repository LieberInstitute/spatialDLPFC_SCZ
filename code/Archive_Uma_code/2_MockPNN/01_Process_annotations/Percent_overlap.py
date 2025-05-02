'''
For Visium-IF
Channel0 = DAPI, DAPI
Channel1 = Claudin5 (Alex 488),
Channel2 = NeuN (Alexa 555),
Channel3 = WFA (Alexa 647),
Channel4 = AF (Autofluorescence), sample AF
Channel5 = Thumbnail
'''


import numpy as np
import helper_functions
import pyhere
import re
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
import itertools
import collections
from itertools import product
from collections import defaultdict, Counter
from shapely.geometry import Point
from shapely.geometry.polygon import Polygon


img_dir = '/dcs04/lieber/marmaypag/spatialDLPFC_SCZ_LIBD4100/raw-data/images/2_MockPNN/Training_tiles/'
# img_dir = pyhere.here('raw_data', 'images', '2_Mock_PNN', 'Training_tiles')
csv_dir = '/dcs04/lieber/marmaypag/spatialDLPFC_SCZ_LIBD4100/processed-data/2_MockPNN/Training_tiles/Manual_annotations/Annotations/'
# csv_dir = pyhere.here('processed-data', '2_Mock_PNN', 'Training_tiles', 'Manual_annotations', 'Annotations')
img_test = pyhere.here('raw-data', 'images', '2_MockPNN', 'Training_tiles', '20220712_VIF_MockPNN_Strong_Scan1_[6384,53057]_component_data_11.tif')
csv_test = '/dcs04/lieber/marmaypag/spatialDLPFC_SCZ_LIBD4100/processed-data/2_MockPNN/Training_tiles/Manual_annotations/Annotations/20220712_VIF_MockPNN_Strong_Scan1_[6384,53057]_component_data_11.csv'


# read manual annotation files as well for the corresponding image file
for img_name in os.listdir(img_dir):
    for csv_name in os.listdir(csv_dir):
        if int(img_name.split('_')[8].split('.')[0]) == int(csv_name.split('_')[8].split('.')[0]):
            # print(int(img_name.split('_')[8].split('.')[0]), int(csv_name.split('_')[8].split('.')[0]))
            print(img_name, csv_name)
            images = Image.open(os.path.join(img_dir, img_name))
            csv = manual_annot(os.path.join(csv_dir, csv_name))
            print("Manual annotated PNNs:", len(csv))
            images.seek(1) # seek the frame of interest #1=Claudin
            cla = cv2.normalize(np.array(images, dtype = 'float32'), np.zeros(np.array(images, dtype = 'float32').shape, np.double), 1.0, 0.0, cv2.NORM_MINMAX) # normalisation using local min_max
            img_cla = Image.fromarray(cla) # reconstruct the image
            images.seek(3) # seek the frame of interest #3=WFA
            wfa = cv2.normalize(np.array(images, dtype = 'float32'), np.zeros(np.array(images, dtype = 'float32').shape, np.double), 1.0, 0.0, cv2.NORM_MINMAX) # normalisation using local min_max
            img_wfa = Image.fromarray(wfa) # reconstruct the image
            wfa255 = np.array(wfa * 255, dtype = np.uint8)
            wfac = skimage.color.gray2rgb(wfa)
            hierachy, img_threshold = cv2.threshold(np.array(cla * 255, dtype = np.uint8), 10, 255, cv2.THRESH_BINARY)
            contours,_ = cv2.findContours(img_threshold, cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_SIMPLE)
            print("Claudin contours:",len(contours))
            for cnt in contours:
                x,y,w,h = cv2.boundingRect(cnt)
                # print("CLAUDIN CONTOURS", x,y,w,h)
                out_img = cv2.rectangle(wfac, (x-10,y-10), (x+w+10, y+h+10), (0,0,0), -1)
                # out_img now has black colored on most of the blood vessels; out has changed contrast pixels where PNNs are highlighted
                # so now we find contours for the PNNs using the out_img
            out_img_gry = skimage.color.rgb2gray(out_img) # convert to gray to find contours and increase contrast
            # fig,ax = plt.subplots(figsize = (20,20))
            # ax.imshow(out_img_gry)
            # fig.show()
            # hist_plot(out_img)
            # for row in range(out_img_gry.shape[0]):
            #     for col in range(out_img_gry.shape[1]):
            #         if out_img_gry[row][col] > 0.02:
            #             out_img_gry[row][col] = 0.0
            #         else:
            #             out_img_gry[row][col] = 1.0
            # fig,ax = plt.subplots(figsize = (20,20))
            # ax.imshow(out_img_gry)
            # fig.show()
            # out_img_gry[out_img_gry > 0.2] = 0.0 # decrease contrast of background
            # out_img_gry[out_img_gry <= 0.2] = 1.0 # increase the contrast of PNNs
            # fig,ax = plt.subplots(figsize = (20,20))
            # ax.imshow(out_img_gry)
            # fig.show()
            out_img255 = np.array(out_img_gry * 255, dtype = np.uint8) # change scale to 0-255 for find contours
            out_img_clr = skimage.color.gray2rgb(np.array(out_img_gry * 255, dtype = np.uint8)) # convert to color to draw colored bb
            hierachy1, img_threshold1 = cv2.threshold(np.array(out_img_gry * 255, dtype = np.uint8), 100, 255, cv2.ADAPTIVE_THRESH_MEAN_C)
            contours1,_ = cv2.findContours(img_threshold1, cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_SIMPLE)
            # cv2.drawContours(out_imgc, contours1, -1, (0, 255, 0), 2, cv2.LINE_AA) # color scheme: BGR
            x,y,w,h,area = [], [], [], [], []
            for cnt in contours1:
                coords = cv2.boundingRect(cnt) # x,y,w,h
                x1,y1,w1,h1 = cv2.boundingRect(cnt)
                area1 = cv2.contourArea(cnt)
                if area1 >= 1000:
                    out_img1 = cv2.rectangle(out_img_clr, (x1-10,y1-10), (x1+w1+10, y1+h1+10), (0,0,0), -1) # eliminating all the big objects
                elif area1 >= 100 and area1 < 2000:
                    # if (w1*h1) >= 300:
                    x.append(x1)
                    y.append(y1)
                    w.append(w1)
                    h.append(h1)
                    area.append(area1)
                    out_img1 = cv2.rectangle(out_img_clr, (x1-10,y1-10), (x1+w1+10, y1+h1+10), (255,0,0), 2) # change the color to black (0,0,0) if bb is not needed
                    rect = cv2.minAreaRect(cnt)
                    box = cv2.boxPoints(rect)
                    box = np.int0(box)
                    out_img1 = cv2.drawContours(out_img1,[box],0,(0,0,0),1) # comment out if contour box is not needed
            print("total number of PNNs in ", img_name, "is", len(x))
            # fig,ax = plt.subplots(figsize = (20,20))
            # ax.imshow(out_img1)
            # fig.show()
            draw_rect(csv, out_img1, (0,255,0)) # draw manual annotations BB on PNN segmented image GREEN
            # fig,ax = plt.subplots(figsize = (20,20))
            # ax.imshow(out_img1)
            # fig.show()
            df_wfa_ml = create_df(x,y,w,h, area, os.path.join(img_dir, img_name), 'PNN')
            df_wfa_ml.to_csv('/dcs04/lieber/marmaypag/spatialDLPFC_SCZ_LIBD4100/processed-data/2_MockPNN/Training_tiles/ML_annotations/Images/Original/'+ img_name.split('.')[0] + '.csv')

            overlap_ml_ann = []
            for pnn in range(len(df_wfa_ml)):
                for ann in range(len(csv)):
                    xmin1, xmin2, xmax1, xmax2, ymin1, ymin2, ymax1, ymax2 = df_wfa_ml['x1'][pnn],csv['x1'][ann],df_wfa_ml['x4'][pnn],csv['x4'][ann],df_wfa_ml['y1'][pnn],csv['y1'][ann],df_wfa_ml['y4'][pnn],csv['y4'][ann]
                    if xmax1 >= xmin2 and xmax2 >= xmin1 and ymax1 >= ymin2 and ymax2 >= ymin1:
                        overlap_ml_ann.append(ann)
            print("Number of manual and segmented box overlap",len(Counter(overlap_ml_ann).keys()))
            # fig,ax = plt.subplots(figsize = (20,20))
            # ax.imshow(out_img1)
            # fig.show()
            # DAPI SEGMENTATION
            dapi, dapi_clr = read_norm(os.path.join(img_dir, img_name), 0)
            shifted, thresh, gray = morph_transform(dapi_clr)
            labels = find_labels(thresh)
            dpx, dpy, dpw, dph, area, segmented_dapi = draw_rect_dapi(labels, gray, dapi_clr)
            img_info_dapi = create_df(dpx, dpy, dpw, dph, area, os.path.join(img_dir, img_name), 'DAPI')
            draw_rect(img_info_dapi, out_img1, (0,0,255)) # draw dapi BB on the PNN segmented image
            overlap_ml_dapi = []
            for pnn1 in range(len(df_wfa_ml)):
                for dap in range(len(img_info_dapi)):
                    xmin1, xmin2, xmax1, xmax2, ymin1, ymin2, ymax1, ymax2 = df_wfa_ml['x1'][pnn1], img_info_dapi['x1'][dap], df_wfa_ml['x4'][pnn1], img_info_dapi['x4'][dap], df_wfa_ml['y1'][pnn1], img_info_dapi['y1'][dap],df_wfa_ml['y4'][pnn1], img_info_dapi['y4'][dap]
                    if xmax1 >= xmin2 and xmax2 >= xmin1 and ymax1 >= ymin2 and ymax2 >= ymin1:
                        overlap_ml_dapi.append(pnn1)
            print("Number of DAPI and segmented box overlap",len(Counter(overlap_ml_dapi).keys()))
            # NEUN SEGMENTATION
            neun= read_norm(os.path.join(img_dir, img_name), 2)
            neun_contours = detect_contours(neun)
            nx,ny,nw,nh,narea, neun_segmented = draw_contours(neun_contours, neun, (0,255,0), 2)
            df_ml_neun = create_df(nx,ny,nw,nh, narea,os.path.join(img_dir, img_name) , 'NeuN')
            draw_rect(df_ml_neun, out_img1, (255,215,0)) # draw neun boxes on the same PNN and dapi segmented image
            overlap_ml_neun = []
            for pnn2 in range(len(df_wfa_ml)):
                for nn in range(len(df_ml_neun)):
                    xmin1, xmin2, xmax1, xmax2, ymin1, ymin2, ymax1, ymax2 = df_wfa_ml['x1'][pnn2], df_ml_neun['x1'][nn], df_wfa_ml['x4'][pnn2], df_ml_neun['x4'][nn], df_wfa_ml['y1'][pnn2], df_ml_neun['y1'][nn],df_wfa_ml['y4'][pnn2], df_ml_neun['y4'][nn]
                    if xmax1 >= xmin2 and xmax2 >= xmin1 and ymax1 >= ymin2 and ymax2 >= ymin1:
                        overlap_ml_neun.append(pnn2)
            print("Number of NeuN and segmented box overlap",len(Counter(overlap_ml_neun).keys()))
            fig,ax = plt.subplots(figsize = (20,20))
            ax.imshow(out_img1)
            fig.show()


            # for i in range(len(df_wfa_ml)): # PNN
            #     for k in range(len(csv)):
            #         for j in range(len(img_info_dapi)): # DAPI
            #             xmin1, xmax1, xmin2, xmax2 = df_wfa_ml['x1'][i], df_wfa_ml['x4'][i], csv['x1'][k], csv['x4'][k]
            #             ymin1, ymax1, ymin2, ymax2 = df_wfa_ml['y1'][i], df_wfa_ml['y4'][i], csv['y1'][k], csv['y4'][k]
            #             # xmin1, xmax1, xmin2, xmax2 = df_wfa_ml['x1'][i], df_wfa_ml['x4'][i], img_info_dapi['x1'][k], img_info_dapi['x4'][k]
            #             # ymin1, ymax1, ymin2, ymax2 = df_wfa_ml['y1'][i], df_wfa_ml['y4'][i], img_info_dapi['y1'][k], img_info_dapi['y4'][k]
            #             if xmax1 >= xmin2 and xmax2 >= xmin1 and ymax1 >= ymin2 and ymax2 >= ymin1:
            #                 print(xmin1, xmax1, xmin2, xmax2, i, k)
            fig,ax = plt.subplots(figsize = (20,20))
            ax.imshow(out_img1)
            fig.show()


# i = manual annotations
# j = dapi
# k = NeuN
# l = segmentations

for i in range(len(csv)):
    for j in range(len(df_wfa_ml)):
        for k in range(len(img_info_dapi)):
            xmin1, xmin2, xmin3 = csv['x1'][i], df_wfa_ml['x1'][j], img_info_dapi['x1'][k]
            ymin1, ymin2, ymin3 = csv['y1'][i], df_wfa_ml['y1'][j], img_info_dapi['y1'][k]
            xmax1, xmax2, xmax3 = csv['x4'][i], df_wfa_ml['x4'][j], img_info_dapi['x4'][k]
            ymax1, ymax2, ymax3 = csv['y4'][i], df_wfa_ml['y4'][j], img_info_dapi['y4'][k]


# loop through the whole directory to find the overlaps between WFA and DAPI
img_dir = '/dcs04/lieber/marmaypag/spatialDLPFC_SCZ_LIBD4100/raw-data/images/2_MockPNN/Training_tiles/'
csv_dir = '/dcs04/lieber/marmaypag/spatialDLPFC_SCZ_LIBD4100/processed-data/2_MockPNN/Training_tiles/Manual_annotations/Annotations/'

for img_name in os.listdir(img_dir):
    for csv_name in os.listdir(csv_dir):
        if int(img_name.split('_')[8].split('.')[0]) == int(csv_name.split('_')[8].split('.')[0]):
            # print(int(img_name.split('_')[8].split('.')[0]), int(csv_name.split('_')[8].split('.')[0]))
            print(img_name, csv_name)
            dapi, dapi_clr = read_norm(os.path.join(img_dir, img_name), 0)
            print(img_name)
            csv = manual_annot(os.path.join(csv_dir, csv_name))
            print(len(csv))
            shifted, thresh, gray = morph_transform(dapi_clr)
            labels = find_labels(thresh)
            dpx, dpy, dpw, dph, area, segmented_dapi = draw_rect_dapi(labels, gray, dapi_clr)
            img_info_dapi = create_df(dpx, dpy, dpw, dph, area, os.path.join(img_dir, img_name), 'DAPI')
            draw_rect(csv, segmented_dapi)
            # df_wfa_ml = create_df(x,y,w,h, area, os.path.join(img_dir, img_name), 'PNN')
            box_lists_dapi = []
            for i in range(len(img_info_dapi)): # PNN
                for k in range(len(csv)):
                    # for j in range(len(img_info_dapi)): # DAPI
                        xmin1, xmax1, xmin2, xmax2 = img_info_dapi['x1'][i], img_info_dapi['x4'][i], csv['x1'][k], csv['x4'][k]
                        ymin1, ymax1, ymin2, ymax2 = img_info_dapi['y1'][i], img_info_dapi['y4'][i], csv['y1'][k], csv['y4'][k]
                        # xmin1, xmax1, xmin2, xmax2 = df_wfa_ml['x1'][i], df_wfa_ml['x4'][i], img_info_dapi['x1'][k], img_info_dapi['x4'][k]
                        # ymin1, ymax1, ymin2, ymax2 = df_wfa_ml['y1'][i], df_wfa_ml['y4'][i], img_info_dapi['y1'][k], img_info_dapi['y4'][k]
                        # xmin1, xmax1, xmin2, xmax2 = df_wfa_ml['x1'][i], df_wfa_ml['x4'][i], img_info_dapi['x1'][k], img_info_dapi['x4'][k]
                        # ymin1, ymax1, ymin2, ymax2 = df_wfa_ml['y1'][i], df_wfa_ml['y4'][i], img_info_dapi['y1'][k], img_info_dapi['y4'][k]
                        if xmax1 >= xmin2 and xmax2 >= xmin1 and ymax1 >= ymin2 and ymax2 >= ymin1:
                            # print(xmin1, xmax1, xmin2, xmax2, i, k)
                            box_lists_dapi.append(k)
            print(Counter(box_lists_dapi))

            fig,ax = plt.subplots(figsize = (20,20))
            ax.imshow(segmented_dapi)
            fig.show()



# loop through the whole directory for NeuN and WFA overlaps
img_dir = '/dcs04/lieber/marmaypag/spatialDLPFC_SCZ_LIBD4100/raw-data/images/2_MockPNN/Training_tiles/'
csv_dir = '/dcs04/lieber/marmaypag/spatialDLPFC_SCZ_LIBD4100/processed-data/2_MockPNN/Training_tiles/Manual_annotations/Annotations/'

for img_name in os.listdir(img_dir):
    for csv_name in os.listdir(csv_dir):
        if int(img_name.split('_')[8].split('.')[0]) == int(csv_name.split('_')[8].split('.')[0]):
            # print(int(img_name.split('_')[8].split('.')[0]), int(csv_name.split('_')[8].split('.')[0]))
            print(img_name, csv_name)
            csv = manual_annot(os.path.join(csv_dir, csv_name))
            print(len(csv))
            neun= read_norm(os.path.join(img_dir, img_name), 2)
            neun_contours = detect_contours(neun)
            nx,ny,nw,nh,narea, neun_segmented = draw_contours(neun_contours, neun, (0,255,0), 2)
            df_ml_neun = create_df(nx,ny,nw,nh, narea,os.path.join(img_dir, img_name) , 'NeuN')
            draw_rect(csv, neun_segmented)
            box_lists_neun = []
            for i in range(len(df_ml_neun)): # PNN
                for k in range(len(csv)):
                    # for j in range(len(img_info_dapi)): # DAPI
                        xmin1, xmax1, xmin2, xmax2 = df_ml_neun['x1'][i], df_ml_neun['x4'][i], csv['x1'][k], csv['x4'][k]
                        ymin1, ymax1, ymin2, ymax2 = df_ml_neun['y1'][i], df_ml_neun['y4'][i], csv['y1'][k], csv['y4'][k]
                        # xmin1, xmax1, xmin2, xmax2 = df_wfa_ml['x1'][i], df_wfa_ml['x4'][i], img_info_dapi['x1'][k], img_info_dapi['x4'][k]
                        # ymin1, ymax1, ymin2, ymax2 = df_wfa_ml['y1'][i], df_wfa_ml['y4'][i], img_info_dapi['y1'][k], img_info_dapi['y4'][k]
                        # xmin1, xmax1, xmin2, xmax2 = df_wfa_ml['x1'][i], df_wfa_ml['x4'][i], img_info_dapi['x1'][k], img_info_dapi['x4'][k]
                        # ymin1, ymax1, ymin2, ymax2 = df_wfa_ml['y1'][i], df_wfa_ml['y4'][i], img_info_dapi['y1'][k], img_info_dapi['y4'][k]
                        if xmax1 >= xmin2 and xmax2 >= xmin1 and ymax1 >= ymin2 and ymax2 >= ymin1:
                            # print(xmin1, xmax1, xmin2, xmax2, i, k)
                            box_lists_neun.append(k)
            print(Counter(box_lists_neun))

            fig,ax = plt.subplots(figsize = (20,20))
            ax.imshow(neun_segmented)
            fig.show()

