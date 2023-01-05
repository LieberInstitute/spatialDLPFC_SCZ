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
import pyhere
from pathlib import Path
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


img_dir = pyhere.here('raw-data', 'images', '2_MockPNN', 'Training_tiles')
img_dir_NTC = pyhere.here('raw-data', 'images', '2_MockPNN', '20220712_VIF_MockPNN_Strong_NTC_C1_Br5182_MLtraining')
img_NTC = pyhere.here('raw-data', 'images', '2_MockPNN', '20220712_VIF_MockPNN_Strong_NTC_C1_Br5182_MLtraining', '20220712_VIF_MockPNN_Strong_NTC_Scan1_[11013,50974]_component_data.tif')
img_SCZ = pyhere.here('raw-data', 'images', '2_MockPNN', '20220712_VIF_MockPNN_Strong_SCZ_C1_Br2039_MLtraining', '20220712_VIF_MockPNN_Strong_SCZ_Scan1_[10629,49106]_component_data.tif')
img_test = pyhere.here('raw-data', 'images', '2_MockPNN', 'Training_tiles', '20220712_VIF_MockPNN_Strong_SCZ_Scan1_[6925,49106]_component_data_24.tif')
img_test = '/dcs04/lieber/marmaypag/spatialDLPFC_SCZ_LIBD4100/raw-data/images/2_MockPNN/Training_tiles/20220712_VIF_MockPNN_Strong_Scan1_[12864,50280]_component_data_17.tif'


# Read and normalize the image
img_claudin = Image.open(img_test)
img_claudin.seek(1) # channel 1 = Claudin 5
claudin = cv2.normalize(np.array(img_claudin, dtype = 'float32'), np.zeros(np.array(img_claudin, dtype = 'float32').shape, np.double), 1.0, 0.0, cv2.NORM_MINMAX)


# Increasing the contrast
claudin[claudin <= claudin.mean()] = 0.0
claudin[claudin >= claudin.mean()] = 1.0

# Plot the normalized/pre-processed image
fig,ax = plt.subplots(figsize = (20,20))
ax.imshow(claudin)
fig.show()

# Find contours/segmentations for Claudin5 layer
# claudin_gray = skimage.color.rgb2gray(claudin) # convert to gray to find contours
claudin_clr = skimage.color.gray2rgb((np.array((claudin * 255), dtype = np.uint8))) # convert to color to draw colored bb
hierachy, img_threshold = cv2.threshold((np.array((claudin * 255), dtype = np.uint8)), 100, 255, cv2.THRESH_BINARY)
contours,_ = cv2.findContours(img_threshold, cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_SIMPLE)
# cv2.drawContours(out_imgc, contours1, -1, (0, 255, 0), 2, cv2.LINE_AA) # color scheme: BGR len(contours)
clx, cly, clw, clh = [],[],[],[]
for cnt in contours:
      x,y,w,h = cv2.boundingRect(cnt)
      if(w*h >= 100):
            clx.append(x)
            cly.append(y)
            clw.append(w)
            clh.append(h)
            out_img_bb = cv2.rectangle(claudin_clr, (x,y), (x+w+5, y+h+5), (0,0,255), 2)
            rect = cv2.minAreaRect(cnt)
            box = cv2.boxPoints(rect)
            box = np.int0(box)
            out_img_cnt = cv2.drawContours(out_img_bb,[box],0,(0,255,0),1)

cv2.imwrite('/users/ukaipa/PNN/One_img/Images/Cla_seg.tif', out_img_cnt) # save the image if needed

fig,ax = plt.subplots(nrows = 1, ncols = 2, figsize = (20,20))
ax[0].imshow(claudin)
ax[0].title.set_text('Original')
ax[1].imshow(out_img_cnt)
ax[1].title.set_text('Segemented')
fig.show()

# Put the BB details in the csv
col_names = ['img_file_name','type_of_object_str', 'x1', 'y1', 'Width', 'Height', 'total_number_claudin']
object_name = 'Blood_vessels' # name of the objects stored in the dataframe
file_name = os.path.basename(img_test) # image file name

dict = {col_names[0]: file_name, col_names[1]: object_name, col_names[2]: clx, col_names[3]: cly, col_names[4]: clw, col_names[5]: clh}
img_info_claudin = pd.DataFrame(dict, columns = col_names)
# compute the rest of the coordinates of the BB
img_info_claudin['x2'] = img_info_claudin['x1'] + img_info_claudin['Width']
img_info_claudin['y2'], img_info_claudin['x3'] = img_info_claudin['y1'], img_info_claudin['x1']
img_info_claudin['y3'] = img_info_claudin['y1'] + img_info_claudin['Height']
img_info_claudin['x4'], img_info_claudin['y4'] = img_info_claudin['x2'], img_info_claudin['y3']
img_info_claudin = img_info_claudin[['img_file_name', 'type_of_object_str', 'x1', 'y1', 'x2', 'y2', 'x3', 'y3', 'x4', 'y4', 'Width', 'Height']]

# export the dataframe to csv
img_info_claudin.to_csv("/dcs04/lieber/marmaypag/spatialDLPFC_SCZ_LIBD4100/processed-data/2_MockPNN/Training_tiles/ML_annotations/Annotations/20220712_VIF_MockPNN_Strong_Scan1_[12864,50280]_component_data_17_claudin.csv")

###### segmentation using the helper_functions
im_cla = read_norm(img_test, 1)
cla_contours = detect_contours(im_cla)
clx,cly,clw,clh, cl_area, contour_cla = draw_contours(cla_contours, im_cla, 1, (255,0,0), 2)
img_info_claudin = create_df(clx,cly,clw,clh, cl_area, img_test, 'claudin')
