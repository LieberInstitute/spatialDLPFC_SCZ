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
import collections
from itertools import product
from collections import defaultdict, Counter
from shapely.geometry import Point
from shapely.geometry.polygon import Polygon


img_dir = pyhere.here('raw-data', 'images', '2_MockPNN', 'Training_tiles')
img_dir_NTC = pyhere.here('raw-data', 'images', '2_MockPNN', '20220712_VIF_MockPNN_Strong_NTC_C1_Br5182_MLtraining')
img_NTC = pyhere.here('raw-data', 'images', '2_MockPNN', '20220712_VIF_MockPNN_Strong_NTC_C1_Br5182_MLtraining', '20220712_VIF_MockPNN_Strong_NTC_Scan1_[11013,50974]_component_data.tif')
img_SCZ = pyhere.here('raw-data', 'images', '2_MockPNN', '20220712_VIF_MockPNN_Strong_SCZ_C1_Br2039_MLtraining', '20220712_VIF_MockPNN_Strong_SCZ_Scan1_[10629,49106]_component_data.tif')
img_test = pyhere.here('raw-data', 'images', '2_MockPNN', 'Training_tiles', '20220712_VIF_MockPNN_Strong_SCZ_Scan1_[6925,49106]_component_data_24.tif')

img_test = '/dcs04/lieber/marmaypag/spatialDLPFC_SCZ_LIBD4100/raw-data/images/2_MockPNN/Training_tiles/20220712_VIF_MockPNN_Strong_Scan1_[6925,49106]_component_data_24.tif'

# loop through the whole directory to segment only NeuN
img_dir = '/dcs04/lieber/marmaypag/spatialDLPFC_SCZ_LIBD4100/raw-data/images/RealPNN/round1/20220814_VIF_PNN_S2_SCZ/'
csv_dst = '/dcs04/lieber/marmaypag/spatialDLPFC_SCZ_LIBD4100/processed-data/RealPNN/NeuN_segmentations/Image_csvs/'
img_dst = '/dcs04/lieber/marmaypag/spatialDLPFC_SCZ_LIBD4100/processed-data/RealPNN/NeuN_segmentations/Image_segmentations'

for img_name in os.listdir(img_dir):
    if img_name.endswith('.tif'):
        print(img_name)
        neun= read_norm(os.path.join(img_dir, img_name), 2)
        neun_contours = detect_contours(neun)
        nx,ny,nw,nh,narea, neun_segmented = draw_contours(neun, 2, neun_contours, (0,255,0), 2)
        df_ml_neun = create_df(nx,ny,nw,nh, narea,os.path.join(img_dir, img_name) , 'NeuN')
        df_ml_neun.to_csv(path_or_buf = (csv_dst + img_name.split('.')[0] + '.csv')) # df to csv and save it in the csv_dst folder
        cv2.imwrite((img_dst + img_name.split('.')[0] + '.tif'), neun_segmented)



###### using helper_functions for 1 image
neun = read_norm(img_test, 2)
fig,ax = plt.subplots(figsize = (20,20))
ax.imshow(neun, cmap = 'gray')
fig.show()
neun_contours = detect_contours(neun)
nx,ny,nw,nh,narea, seg_neun = draw_contours(neun, 2, neun_contours,(0,255,0), 2)
plot_img(neun, seg_neun)
img_info_neun = create_df(nx,ny,nw,nh, narea, img_test, 'NeuN')
