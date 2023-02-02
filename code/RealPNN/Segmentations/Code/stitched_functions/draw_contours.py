'''
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

# extract the watershed algorithm labels
def draw_rect_dapi(labels, gray, dapi):
    dpx, dpy, dpw, dph, area = [], [], [], [], []
    # print("1) entering the label loop")
    for counter, label in enumerate(np.unique(labels)):
        if label == 0: # label marked 0 are background
            continue
        mask = np.zeros(gray.shape, dtype="uint8") # create masks that only have the detected labels as foreground and 0 as background
        mask[labels == label] = 255
        # print("2) found the masks", counter)
        # detect contours in the mask and grab the largest one
        cnts = cv2.findContours(mask.copy(), cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_SIMPLE) # detect the watershed contours
        # print("3) found contours", counter)
        cnts = imutils.grab_contours(cnts) # extract only the contours
        # print("4) grabbed contours",counter)
        c = max(cnts, key=cv2.contourArea) # get the area
        x,y,w,h = cv2.boundingRect(c) # BB coordinates
        area1 = cv2.contourArea(c)
        # print("5) found BB coordinates and area", counter)
        if 2 >= area1 <= 500: # adding shape filter to filter out the smaller contours and the biggest cluter contours
            # print("6) appending BB coordinates", counter)
            dpx.append(x)
            dpy.append(y)
            dpw.append(w)
            dph.append(h)
            area.append(cv2.contourArea(c))
            # print("7) appended")
            ws_img_bb = cv2.rectangle(skimage.color.gray2rgb(dapi), (x,y), (x+w, y+h), (0,255,0), 1) # if a colored BB is not required then, change color to (0,0,0) and thickness to 1
            # print("8) drawing rectangles")
            print("drawing rectangle number",counter, "with area", area1)
    return dpx, dpy, dpw, dph, area, ws_img_bb

# draw the extracted contours onto the image
def draw_contours(normalised_img, ch_num, contours = None,  color = None, thickness = None, dapi_clr = None):
    if ch_num == 1: #Claudin
        color_img = skimage.color.gray2rgb(normalised_img, dtype = np.uint8)
        x, y, w, h, area = [],[],[],[],[]
        for cnt in contours:
            x_, y_, w_, h_ = cv2.boundingRect(cnt)
            area_ = cv2.contourArea(cnt)
            if area_ >= 1000:
                contour_img = cv2.rectangle(color_img, (x_-10,y_-10), (x_+w_+10, y_+h_+10), (0,0,0), -1) # eliminating all the big objects
            elif 90 <= area_ < 2500: # size threshold
                # area_ = cv2.contourArea(cnt)
                # print(ax,ay,aw,ah)
                x.append(x_)
                y.append(y_)
                w.append(w_)
                h.append(h_)
                area.append(area_)
                bb_img = cv2.rectangle(color_img, (x_,y_), (x_+w_+10, y_+h_+10), color, thickness) #(255,0,0), 2-- to draw colored boxes
                box = np.int0(cv2.boxPoints(cv2.minAreaRect(cnt)))
                contour_img = cv2.drawContours(bb_img,[box],0,(0,0,0),1) # change the color and thickness here if contours need to be visible
        return x, y, w, h, area, contour_img
    elif ch_num == 3: #NeuN
        color_img = skimage.color.gray2rgb(normalised_img, dtype = np.uint8)
        x, y, w, h, area = [],[],[],[],[]
        for cnt in contours:
            x_, y_, w_, h_ = cv2.boundingRect(cnt)
            area_ = cv2.contourArea(cnt)
            if area_ >= 100:
                # area_ = cv2.contourArea(cnt)
                # print(ax,ay,aw,ah)
                x.append(x_)
                y.append(y_)
                w.append(w_)
                h.append(h_)
                area.append(area_)
                bb_img = cv2.rectangle(color_img, (x_,y_), (x_+w_+10, y_+h_+10), color, thickness) #(255,0,0), 2-- to draw colored boxes
                box = np.int0(cv2.boxPoints(cv2.minAreaRect(cnt)))
                contour_img = cv2.drawContours(bb_img,[box],0,(0,0,0),1) # change the color and thickness here if contours need to be visible
        return x, y, w, h, area, contour_img
    else: #WFA
        color_img = skimage.color.gray2rgb((np.array((normalised_img * 255), dtype = np.uint8)))
        x, y, w, h, area = [],[],[],[],[]
        for cnt in contours:
            x_, y_, w_, h_ = cv2.boundingRect(cnt)
            area_ = cv2.contourArea(cnt)
            if area_ >= 1000:
                contour_img = cv2.rectangle(color_img, (x_-10,y_-10), (x_+w_+10, y_+h_+10), (0,0,0), -1) # eliminating all the big objects
            elif 90 <= area_ < 2500: # size threshold
                # print(ax,ay,aw,ah)
                x.append(x_)
                y.append(y_)
                w.append(w_)
                h.append(h_)
                area.append(area_)
                bb_img = cv2.rectangle(color_img, (x_,y_), (x_+w_+10, y_+h_+10), color, thickness) #(255,0,0), 2-- to draw colored boxes
                box = np.int0(cv2.boxPoints(cv2.minAreaRect(cnt)))
                contour_img = cv2.drawContours(bb_img,[box],0,(0,0,0),1) # change the color and thickness here if contours need to be visible
                # cv2.putText(contour_img, label, (x_,y_), cv2.FONT_HERSHEY_SIMPLEX, 0.7, (125, 246, 55), 3)
        return x, y, w, h, area, contour_img
