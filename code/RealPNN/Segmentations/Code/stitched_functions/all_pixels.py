env'''
For Visium-IF
Channel0 = DAPI, DAPI
Channel1 = Claudin5 (Alex 488),
Channel2 = NeuN (Alexa 555),
Channel3 = WFA (Alexa 647),
Channel4 = AF (Autofluorescence), sample AF
Channel5 = Thumbnail
'''

from __future__ import print_function
from skimage.feature import peak_local_max
from skimage.segmentation import find_boundaries, watershed
from scipy import ndimage
import argparse
from argparse import ArgumentParser
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
from itertools import product
from collections import defaultdict
from shapely.geometry import Point
from shapely.geometry.polygon import Polygon

# draw a white rectangle filled using the coordinates from the csv
from collections import Counter
pix_list, locs_list, mean_pix_int_list  = [], [], []
def all_pix_pnns(img_info_df, contour_img, original_img):
    for box in range(len(img_info_df['x1'])):
        print(box)
        gray_image = np.full((17799, 16740), 0, dtype=np.uint8) # new blank image so the original pix intensities are retained
        rect = cv2.rectangle(gray_image, (img_info_df['x1'][box], img_info_df['y1'][box]), (img_info_df['x4'][box], img_info_df['y4'][box]), (255,255,255), -1) # draw white filled rect on the copy of the image
        cv2.putText(contour_img, ('%d'%box), (img_info_df['x1'][box],img_info_df['y1'][box]), cv2.FONT_HERSHEY_SIMPLEX, 2, (125, 246, 55), 3)
        # gray_seg_wfa = skimage.color.rgb2gray(contour_img)
        # plot_img(gray_image, contour_img)
        locs = np.argwhere(gray_image == 255)
        print(locs.shape)
        for i in range(locs.shape[0]):
            for j in range(locs.shape[1] -1):
                if gray_image[locs[i,j],locs[i,j+1]] == 255: # gray image has white filled boxes
                    # print(original_img[locs[i,j],locs[i,j+1]])
                    pix_list.append(original_img[locs[i,j],locs[i,j+1]]) # append all pix intensities of coordinates inside the PNN box
        print("pix mean:", (np.array(pix_list)).mean()) # convert list to array and find the mean pix intensities
        locs_list.append(locs) # append the all pixels of all PNNs detected
        mean_pix_int_list.append((np.array(pix_list)).mean()) # and their mean intensities
    print("Number of PNNs segmented:", len(locs_list))
    print("lengths", len(locs_list), len(mean_pix_int_list))
    img_info_df['pixels'] = locs_list
    img_info_df['mean_pixel_int'] = mean_pix_int_list
    img_info_df = img_info_df[['img_file_name', 'type_of_object_str', 'x1', 'y1', 'x2', 'y2', 'x3', 'y3', 'x4', 'y4', 'xc', 'yc', 'Width', 'Height', 'area', 'mean_pixel_int', 'pixels']]
    fig = plt.figure(figsize = (5, 5))
    plt.bar(list(range(len(img_info_df))), img_info_df['mean_pixel_int'], color = 'blue', width = 0.2)
    plt.xticks(np.arange(0,len(img_info_df), 1), labels = list(range(len(img_info_df))))
    plt.xlabel("Number of segmented PNNs")
    plt.ylabel("Mean pixel intensities")
    plt.title("Mean pixel intensities plot for segmented PNNs")
    plt.show()
    return contour_img, img_info_df, mean_pix_int_list # this returns a color image with PNN contours marked along with numbers
