import numpy as np
import numpy.ma as ma
import pandas as pd
import csv
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
import re
from itertools import product
from collections import defaultdict
from shapely.geometry import Point
from shapely.geometry.polygon import Polygon

csv_dir = '/dcs04/lieber/marmaypag/spatialDLPFC_SCZ_LIBD4100/processed-data/2_MockPNN/Training_tiles/Manual_annotations/Annotations/'
img_dir = '/dcs04/lieber/marmaypag/spatialDLPFC_SCZ_LIBD4100/raw-data/images/2_MockPNN/Training_tiles/'
img_test = '/dcs04/lieber/marmaypag/spatialDLPFC_SCZ_LIBD4100/raw-data/images/2_MockPNN/Training_tiles/20220712_VIF_MockPNN_Strong_Scan1_[10087,51668]_component_data_01.tif'

# function to extract all the coordinates from the csv
def extract_coords(csv_path, object_name):
    img_info_df = pd.read_csv(csv_path)
    col_names = ['img_file_name','type_of_object_str', 'x1', 'y1', 'Width', 'Height']
    object_name = object_name # name of the objects stored in the dataframe
    file_name = os.path.basename(csv_path).split('.')[0] # image file name
    dict = {col_names[0]: file_name, col_names[1]: object_name, col_names[2]: np.int0(np.ceil(img_info_df['X'])), col_names[3]: np.int0(np.ceil(img_info_df['Y'])), col_names[4]: np.int0(np.ceil(img_info_df['Width'])), col_names[5]: np.int0(np.ceil(img_info_df['Height']))}
    img_info_df = pd.DataFrame(dict, columns = col_names)
    img_info_df['x2'] = img_info_df['x1'] + img_info_df['Width']
    img_info_df['y2'], img_info_df['x3'] = img_info_df['y1'], img_info_df['x1']
    img_info_df['y3'] = img_info_df['y1'] + img_info_df['Height']
    img_info_df['x4'], img_info_df['y4'] = img_info_df['x2'], img_info_df['y3']
    img_info_df = img_info_df[['img_file_name', 'type_of_object_str', 'x1', 'y1', 'x2', 'y2', 'x3', 'y3', 'x4', 'y4', 'Width', 'Height']]
    return(img_info_df)


# Initialize empty list
lst_intensities = []

# For each list of contour points...
for i in range(len(contours1)):
    # Create a mask image that contains the contour filled in
    cimg = np.zeros_like(out_img_gry)
    cv2.drawContours(cimg, contours1, i, color=255, thickness=-1)
    pts = np.where(cimg == 255)
    lst_intensities.append(out_img_gry[pts[0], pts[1]])



lst = []
for j in range(len(lst_intensities)):
    l =  lst_intensities[j] * 255
    lst.append(l)

high = []
for k in range(len(lst)):
    if lst[k].mean() >= 90:
        print(k, lst[k].mean())
        hi = lst[k].mean()
        high.append(hi)


fig,ax = plt.subplots(nrows = 1, ncols = 2, figsize = (20,20))
ax[0].hist(np.array(high).ravel(),256,[0,255])
ax[0].title.set_text('Histogram of mean pixel intensities of bright PNNs')
ax[0].set(xlabel = 'pixel intensities(grayscale)' , ylabel = 'frequency')
ax[1].hist(np.array(low).ravel(),256,[0,255])
ax[1].title.set_text('Histogram of mean pixel intensities of low intensity PNNs')
ax[1].set(xlabel = 'pixel intensities(grayscale)' , ylabel = 'frequency')
fig.show()


# plot histograms of Sang Ho's annotations for bright and dim PNNs according to their diagnosis
img_dir = '/dcs04/lieber/marmaypag/spatialDLPFC_SCZ_LIBD4100/raw-data/images/2_MockPNN/Training_tiles/'
csv_dir = '/dcs04/lieber/marmaypag/spatialDLPFC_SCZ_LIBD4100/processed-data/2_MockPNN/Training_tiles/Manual_annotations/Annotations/'
cont = []
for tile_name in os.listdir(img_dir):
    for ann_name in os.listdir(csv_dir):
        if (int(re.split(r'[_.]',os.path.basename(tile_name))[8])) % 2 == 0:
            if int(re.split(r'[_.]',os.path.basename(tile_name))[8]) == int(re.split(r'[_.]',os.path.basename(ann_name))[8]):
                # print(int(re.split(r'[_.]',os.path.basename(tile_name))[8]), int(re.split(r'[_.]',os.path.basename(ann_name))[8]))
                img_pnn = Image.open(os.path.join(img_dir, tile_name))
                img_pnn.seek(3)
                print(tile_name, ann_name)
                pnn = cv2.normalize(np.array(img_pnn, dtype = 'float32'), np.zeros(np.array(img_pnn, dtype = 'float32').shape, np.double), 1.0, 0.0, cv2.NORM_MINMAX)
                pnn_clr = skimage.color.gray2rgb(np.array(pnn * 255, dtype = np.uint8))
                print(tile_name, ann_name)
                csv_pnn = extract_coords(os.path.join(csv_dir, ann_name), 'PNN')
                x1,y1,x2,y2,x3,y3,x4,y4 = csv_pnn['x1'], csv_pnn['y1'],csv_pnn['x2'], csv_pnn['y2'],csv_pnn['x3'], csv_pnn['y3'],csv_pnn['x4'], csv_pnn['y4']
                print(tile_name, ann_name)
                for i,x1_,y1_,x2_,y2_,x3_,y3_,x4_,y4_ in zip(range(len(csv_pnn)),x1,y1,x2,y2,x3,y3,x4,y4):
                    cont = np.array([x1_,y1_,x2_,y2_,x3_,y3_,x4_,y4_])
                    print(cont)
                    out_ = cv2.drawContours(pnn_clr,cont,0,(0,0,255),3)







