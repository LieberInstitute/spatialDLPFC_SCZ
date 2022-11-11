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
from itertools import product
from collections import defaultdict
from shapely.geometry import Point
from shapely.geometry.polygon import Polygon

img_test = pyhere.here('raw-data', 'images', '2_MockPNN', 'Training_tiles', '20220712_VIF_MockPNN_Strong_Scan1_[6384,53057]_component_data_11.tif')
csv_test = '/dcs04/lieber/marmaypag/spatialDLPFC_SCZ_LIBD4100/processed-data/2_MockPNN/Training_tiles/Manual_annotations/Annotations/20220712_VIF_MockPNN_Strong_Scan1_[6384,53057]_component_data_11.csv'

# read and normalise images
dapi, dapi_clr = read_norm(img_test, 0)
claudin = read_norm(img_test, 1)
neun = read_norm(img_test, 2)
wfa = read_norm(img_test, 3)

# read the csv got manual annotations
file_11 = manual_annot(csv_test)

# create a new image
new_im = np.zeros(dapi.shape, np.double)

# draw a white rectangle filled using the coordinates from the csv
def draw_rect(df_manual_test, contour_img):
    for box in range(len(df_manual_test['x1'])):
        print(box)
        rect = cv2.rectangle(contour_img, (df_manual_test['x1'][box], df_manual_test['y1'][box]), (df_manual_test['x4'][box], df_manual_test['y4'][box]), (255,255,255), -1)
    return contour_img

rect_img = draw_rect(file_11, new_im)

fig,ax = plt.subplots(figsize = (20,20))
ax.imshow(rect_img)
fig.show()

# get those coordinates in an array
locs = np.where(rect_img == 255)
pixels = new_im[locs]
print(np.mean(pixels), len(locs[0]))

# get the pixel values of those coordinates to find if dapi signal is present
dapi_box_mean # from the dapi segmentation code

