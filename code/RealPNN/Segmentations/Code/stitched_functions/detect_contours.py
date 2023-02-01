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

# detect contours in the normalised_img
def contours(original_img): ### create a separate function for shape detection and run it through this loop for contour detection
    hierachy, img_threshold = cv2.threshold(original_img, 100, 255, cv2.THRESH_BINARY)
    contours,_ = cv2.findContours(img_threshold, cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_SIMPLE)
    print("contours detected", len(contours))
    return contours


def find_labels(threshold):
    D = ndimage.distance_transform_edt(threshold) # Euclidean distance from binary to nearest 0-pixel
    print("distance measured")
    localMax = peak_local_max(D, indices=False, min_distance=5, labels=threshold) # find the local maxima for all the individual objects
    print("local max found")
    markers = ndimage.label(localMax, structure=np.ones((3, 3)))[0] # 8-connectivity connected component analysis
    print("local max markers found")
    labels = watershed(-D, markers, mask=threshold)
    print("{} unique segments found".format(len(np.unique(labels)) - 1))
    return labels
