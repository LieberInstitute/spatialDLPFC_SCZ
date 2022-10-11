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
import numpy.ma as ma
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
import imageio
from skimage import feature, segmentation, draw, measure, morphology
from skimage.morphology import (erosion,dilation,opening,closing,white_tophat,black_tophat,skeletonize,convex_hull_image)
from skimage.draw import polygon_perimeter
import tifffile
import imagecodecs
import itertools
from itertools import product
from collections import defaultdict
from shapely.geometry import Point
from shapely.geometry.polygon import Polygon

# Read all images
img_dir = '/dcs04/lieber/marmaypag/spatialDLPFC_SCZ_LIBD4100/raw-data/images/2_MockPNN/Training_tiles/'
# src_folder = '/dcs04/lieber/marmaypag/spatialDLPFC_SCZ_LIBD4100/raw-data/images/Visium_IF_Inform_no_overlap/'
dst_folder = '/dcs04/lieber/marmaypag/spatialDLPFC_SCZ_LIBD4100/processed-data/2_MockPNN/Training_tiles/ML_annotations/Images/Original/'
# dst_folder_WFA = '/dcs04/lieber/marmaypag/spatialDLPFC_SCZ_LIBD4100/processed-data/WFA_normalized/'

for file_name in os.listdir(img_dir):
    images = Image.open(os.path.join(img_dir, file_name))
    print("Filename", os.path.basename(file_name))
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
    print(len(contours))
    for cnt in contours:
        x,y,w,h = cv2.boundingRect(cnt)
        # print("CLAUDIN CONTOURS", x,y,w,h)
        out_img = cv2.rectangle(wfac, (x-10,y-10), (x+w+10, y+h+10), (0,0,0), -1)
    # out_img now has black colored on most of the blood vessels; out has changed contrast pixels where PNNs are highlighted
    # so now we find contours for the PNNs using the out_img
    out_img_gry = skimage.color.rgb2gray(out_img) # convert to gray to find contours and increase contrast
    out_img_gry[out_img_gry <= 0.2] = 0.0 # decrease contrast of background
    out_img_gry[out_img_gry >= 0.3] = 1.0 # increase the contrast of PNNs
    out_img255 = np.array(out_img_gry * 255, dtype = np.uint8) # change scale to 0-255 for find contours
    out_img_clr = skimage.color.gray2rgb(np.array(out_img_gry * 255, dtype = np.uint8)) # convert to color to draw colored bb
    hierachy1, img_threshold1 = cv2.threshold(np.array(out_img_gry * 255, dtype = np.uint8), 100, 255, cv2.ADAPTIVE_THRESH_MEAN_C)
    contours1,_ = cv2.findContours(img_threshold1, cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_SIMPLE)
    # cv2.drawContours(out_imgc, contours1, -1, (0, 255, 0), 2, cv2.LINE_AA) # color scheme: BGR
    for cnt in contours1:
      coords = cv2.boundingRect(cnt) # x,y,w,h
      x1,y1,w1,h1 = cv2.boundingRect(cnt)
      if (w1*h1) >= 300:
            out_img1 = cv2.rectangle(out_img_clr, (x1-10,y1-10), (x1+w1+10, y1+h1+10), (0,0,0), 1) # change the color to black (0,0,0) if bb is not needed
            rect = cv2.minAreaRect(cnt)
            box = cv2.boxPoints(rect)
            box = np.int0(box)
            out_img1 = cv2.drawContours(out_img1,[box],0,(0,0,0),1) # comment out if contour box is not needed
    fig,ax = plt.subplots(figsize = (20,20))
    ax.imshow(out_img1)
    fig.show()
    print(os.path.basename(file_name).split('.')[0].split('_')[5])
    cv2.imwrite(dst_folder + os.path.basename(file_name).split('.')[0] +'_PNN_segmented' + '.tif', out_img1)

