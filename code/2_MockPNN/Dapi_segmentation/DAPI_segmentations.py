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
import mpldatacursor
import glob
import sys
import cv2
import math
import scipy
from scipy.spatial.distance import *
# from scipy import ndimage, distance
import imageio
import skimage
from skimage import *
from skimage import feature, segmentation, draw, measure, morphology
from skimage.morphology import (erosion,dilation,opening,closing,white_tophat,black_tophat,skeletonize,convex_hull_image)
from skimage.draw import polygon_perimeter
import tifffile as tif
import imagecodecs
import itertools
from itertools import product
from collections import defaultdict
from shapely.geometry import Point
from shapely.geometry.polygon import Polygon


img_dir = '/dcs04/lieber/marmaypag/spatialDLPFC_SCZ_LIBD4100/raw-data/images/2_MockPNN/Training_tiles/'
img_dir_NTC = '/dcs04/lieber/marmaypag/spatialDLPFC_SCZ_LIBD4100/raw-data/images/2_MockPNN/20220712_VIF_MockPNN_Strong_NTC_C1_Br5182_MLtraining/'
img_NTC = '/dcs04/lieber/marmaypag/spatialDLPFC_SCZ_LIBD4100/raw-data/images/2_MockPNN/20220712_VIF_MockPNN_Strong_NTC_C1_Br5182_MLtraining/20220712_VIF_MockPNN_Strong_NTC_Scan1_[11013,50974]_component_data.tif'
img_SCZ = '/dcs04/lieber/marmaypag/spatialDLPFC_SCZ_LIBD4100/raw-data/images/2_MockPNN/20220712_VIF_MockPNN_Strong_SCZ_C1_Br2039_MLtraining/20220712_VIF_MockPNN_Strong_SCZ_Scan1_[10629,49106]_component_data.tif'

# Read and normalize the image
img_dapi = Image.open(img_NTC)
img_dapi.seek(0) # channel 0 = DAPI
dapi = cv2.normalize(np.array(img_dapi, dtype = 'float32'), np.zeros(np.array(img_dapi, dtype = 'float32').shape, np.double), 1.0, 0.0, cv2.NORM_MINMAX)


# Increasing the contrast (merges the dapi--this step might not be necessary)
dapi[dapi <= dapi.mean()] = 0.0
dapi[dapi >= 0.001] = 1.0

# Plot the normalized/pre-processed image
fig,ax = plt.subplots(figsize = (20,20))
ax.imshow(dapi)
fig.show()


