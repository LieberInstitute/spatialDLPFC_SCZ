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
img_neun = Image.open(img_NTC)
img_neun.seek(2) # channel 0 = DAPI
neun = cv2.normalize(np.array(img_neun, dtype = 'float32'), np.zeros(np.array(img_neun, dtype = 'float32').shape, np.double), 1.0, 0.0, cv2.NORM_MINMAX)


# Increasing the contrast --this step might not be necessary
# neun[neun <= neun.mean()] = 0.0
# neun[neun >= neun.mean()] = 1.0

# Plot the normalized/pre-processed image
fig,ax = plt.subplots(figsize = (20,20))
ax.imshow(neun)
fig.show()

# Find contours/segmentations for NeuN layer
neun_gray = skimage.color.rgb2gray(neun) # convert to gray to find contours
neun_clr = skimage.color.gray2rgb(np.array((neun * 255), dtype = np.uint8)) # convert to color to draw colored bb
hierachy, img_threshold = cv2.threshold((np.array((neun * 255), dtype = np.uint8)), 100, 255, cv2.THRESH_BINARY)
contours,_ = cv2.findContours(img_threshold, cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_SIMPLE)
# cv2.drawContours(out_imgc, contours1, -1, (0, 255, 0), 2, cv2.LINE_AA) # color scheme: BGR
nux, nuy, nuw, nuh = [],[],[],[]
for cnt in contours:
      x,y,w,h = cv2.boundingRect(cnt)
      if(w*h >= 100):
            print(x,y,w,h)
            nux.append(x)
            nuy.append(y)
            nuw.append(w)
            nuh.append(h)
            print(x,y,w,h)
            out_img1_neun = cv2.rectangle(neun_clr, (x,y), (x+w+5, y+h+5), (255,0,0), 2)
            rect = cv2.minAreaRect(cnt)
            box = cv2.boxPoints(rect)
            box = np.int0(box)
            out_img2_neun = cv2.drawContours(out_img1_neun,[box],0,(0,0,255),1)

fig,ax = plt.subplots(figsize = (20,20))
fig,ax = plt.subplots(nrows = 1, ncols = 2, figsize = (20,20))
ax[0].imshow(claudin)
ax[0].title.set_text('Original')
ax[1].imshow(out_img2_neun)
ax[1].title.set_text('Segemented')
fig.show()

# Populate the data in the dataframe
col_names = ['img_file_name','type_of_object_str', 'x1', 'y1', 'Width', 'Height', 'total_number_neun']
object_name = 'NeuN' # name of the objects stored in the dataframe
file_name = os.path.basename(img_test) # image file name

dict = {col_names[0]: file_name, col_names[1]: object_name, col_names[2]: nux, col_names[3]: nuy, col_names[4]: nuw, col_names[5]: nuh, col_names[6]: len(nux)}
img_info_neun = pd.DataFrame(dict, columns = col_names)
img_info_neun['x2'] = img_info_neun['x1'] + img_info_neun['Width']
img_info_neun['y2'], img_info_neun['x3'] = img_info_neun['y1'], img_info_neun['x1']
img_info_neun['y3'] = img_info_neun['y1'] + img_info_neun['Height']
img_info_neun['x4'], img_info_neun['y4'] = img_info_neun['x2'], img_info_neun['y3']
img_info_neun = img_info_neun[['img_file_name', 'type_of_object_str', 'x1', 'y1', 'x2', 'y2', 'x3', 'y3', 'x4', 'y4', 'Width', 'Height', 'total_number_neun']]




