
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

def detect_shape_pnns(contour_img, contours): #img_info_df,
    # color_img = skimage.color.gray2rgb(normalised_img)
    shape = "unidentified"
    for cnts in contours:
        peri = cv2.arcLength(cnts, True) # c is the contour
        approx = cv2.approxPolyDP(cnts, 0.04 * peri, True)# gray_seg_wfa = skimage.color.rgb2gray(contour_img)
        c = max(cnts, key=cv2.contourArea) # get the area
        x,y,w,h = cv2.boundingRect(c)
        area_ = cv2.contourArea(cnts)
        # print("peri, approx", peri, approx)
        if 90 <= area_ <= 2500:
            if len(approx) == 3:
                shape = "triangle"
                cv2.putText(contour_img, shape, (x,y), cv2.FONT_HERSHEY_SIMPLEX, 1, (125, 246, 55), 3)
            elif len(approx) == 4:
                (x, y, w, h) = cv2.boundingRect(approx)
                ar = w / float(h)
                shape = "square" if ar >= 0.95 and ar <= 1.05 else "rectangle"
                cv2.putText(contour_img, shape, (x,y), cv2.FONT_HERSHEY_SIMPLEX, 1, (125, 246, 55), 3)
            elif len(approx) == 5:
                shape = "pentagon"
                cv2.putText(contour_img, shape, (x,y), cv2.FONT_HERSHEY_SIMPLEX, 1, (125, 246, 55), 3)
            else:
                shape = "circle"
                cv2.putText(contour_img, shape, (x,y), cv2.FONT_HERSHEY_SIMPLEX, 1, (125, 246, 55), 3)
        print("shape", shape)
    # rect = cv2.rectangle(contour_img, (img_info_df['x1'][box], img_info_df['y1'][box]), (img_info_df['x4'][box], img_info_df['y4'][box]), (255,255,255), 3) # draw white filled rect on the copy of the image
    # cv2.putText(contour_img, shape, (img_info_df['x1'][box], img_info_df['y1'][box]), cv2.FONT_HERSHEY_SIMPLEX, 3, (125, 246, 55), 3)
    return approx, contours, shape, contour_img
