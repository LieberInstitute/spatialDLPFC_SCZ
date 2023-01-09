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
import pylab
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
import itertools
from itertools import product
from collections import defaultdict
from shapely.geometry import Point
from shapely.geometry.polygon import Polygon

img_dir = '/dcs04/lieber/marmaypag/spatialDLPFC_SCZ_LIBD4100/raw-data/images/2_MockPNN/Training_tiles/'
csv_dir = '/dcs04/lieber/marmaypag/spatialDLPFC_SCZ_LIBD4100/processed-data/2_MockPNN/Training_tiles/Manual_annotations/Annotations/'
img_test = pyhere.here('raw-data', 'images', '2_MockPNN', 'Training_tiles', '20220712_VIF_MockPNN_Strong_Scan1_[6384,53057]_component_data_11.tif')
csv_test = '/dcs04/lieber/marmaypag/spatialDLPFC_SCZ_LIBD4100/processed-data/2_MockPNN/Training_tiles/Manual_annotations/Annotations/20220712_VIF_MockPNN_Strong_Scan1_[6384,53057]_component_data_11.csv'
img_test = '/dcs04/lieber/marmaypag/spatialDLPFC_SCZ_LIBD4100/raw-data/images/2_MockPNN/Training_tiles/20220712_VIF_MockPNN_Strong_Scan1_[6384,53057]_component_data_11.tif'


im_cla = read_norm(img_test, 1)

# segment pnns and get their coordinates in a df
# claudin + wfa segmentation
cla_wfa_contour = detect_contours(im_cla)
clx,cly,clw,clh, cl_area, seg_cla_wfa = draw_contours(cla_wfa_contour, im_wfa, (0,0,0), -1)
plot_img(seg_cla, seg_cla_wfa)
out_img_gry = skimage.color.rgb2gray(seg_cla_wfa) # convert to gray to find contours and increase contrast
wfa_contours = detect_contours(out_img_gry)
wfx, wfy, wfw, wfh, pnn_area, seg_wfa = draw_contours(wfa_contours, out_img_gry, (0,0,255), 2)
plot_img(im_wfa, seg_wfa)
img_info_wfa = create_df(wfx, wfy, wfw, wfh, pnn_area, img_test, 'PNN')
