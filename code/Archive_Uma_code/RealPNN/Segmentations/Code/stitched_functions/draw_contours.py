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

# draw the extracted contours onto the image
def draw_detected_contours(normalised_img, ch_num, contours = None,  color = None, thickness = None, dapi_clr = None):
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


def draw_all_contours(color_img, contours, box_color = (255, 0, 0), box_thickness = 2):
    x_,y_,width_,height_,area_ = [],[],[],[],[] #x,y,width,height,area
    print(len(contours))
    contour_img = cv2.drawContours(color_img, contours, -1, box_color, box_thickness)
    for number_cnt, cnt in enumerate(contours):
        # x, y, w, h = cv2.boundingRect(cnt)
        # area = cv2.contourArea(cnt)
        # x_.append(x)
        # y_.append(y)
        # width_.append(w)
        # height_.append(h)
        # area_.append(area)
        # print("appended", number_cnt, "row to csv")
        contour_img = cv2.drawContours(color_img, contours, -1, box_color, box_thickness)
    return contour_img #x_,y_,width_,height_,area_,


def draw_dapi_contours(color_img, contours, box_color = (255, 0, 0), box_thickness = 2):
    x_,y_,width_,height_,area_ = [],[],[],[],[] #x,y,width,height,area
    print(len(contours))
    for number_cnt, cnt in enumerate(contours):
        if cv2.contourArea(cnt) >= 30:
            x, y, w, h = cv2.boundingRect(cnt)
            area = cv2.contourArea(cnt)
            x_.append(x)
            y_.append(y)
            width_.append(w)
            height_.append(h)
            area_.append(area)
            print("appended", number_cnt, "row to csv")
            contour_img = cv2.drawContours(color_img, contours, -1, box_color, box_thickness)
    return x_,y_,width_,height_,area_,contour_img

