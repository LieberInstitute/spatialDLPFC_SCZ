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

training_tiles_path = '/dcs04/lieber/marmaypag/spatialDLPFC_SCZ_LIBD4100/raw-data/images/2_MockPNN/Training_tiles/'
tile_annotations = '/dcs04/lieber/marmaypag/spatialDLPFC_SCZ_LIBD4100/processed-data/2_MockPNN/Training_tiles/Manual_annotations/Annotations/'

# Read the corresponding tiles to its annotations file (needs work to work for all images)
for tile_name in os.listdir(training_tiles_path):
    for ann_name in os.listdir(tile_annotations):
        if int(re.split(r'[_.]',os.path.basename(tile_name))[8]) == int(re.split(r'[_.]',os.path.basename(ann_name))[8]):
            print(tile_name,ann_name)
            img_pnn = Image.open(os.path.join(training_tiles_path, tile_name))
            img_pnn.seek(3)
            pnn = cv2.normalize(np.array(img_pnn, dtype = 'float32'), np.zeros(np.array(img_pnn, dtype = 'float32').shape, np.double), 1.0, 0.0, cv2.NORM_MINMAX)
            csv_pnn = pd.read_csv(os.path.join(tile_annotations, ann_name))
            i = 0
            x,y,w,h = np.int0(np.ceil(csv_pnn['BX'])), np.int0(np.ceil(csv_pnn['BY'])), np.int0(np.ceil(csv_pnn['Width'])), np.int0(np.ceil(csv_pnn['Height']))
            print(x,y,w,h)
            cv2.imwrite('/users/ukaipa/PNN/One_img/Cropped_annotations/img_pnn{}.tif'.format(i), pnn[y:y+h,x:x+w]) # if you want the annotations to be cropped out and saved
            i += 1


# This works for one image
im, ann = '20220712_VIF_MockPNN_Strong_Scan1_[12864,50280]_component_data_17.tif', '20220712_VIF_MockPNN_Strong_Scan1_[12864,50280]_component_data_17.csv'#'20220712_VIF_MockPNN_Strong_Scan1_[6925,49106]_component_data_24.csv'
img = Image.open(os.path.join(training_tiles_path, im))
img.seek(3)
pnn = cv2.normalize(np.array(img, dtype = 'float32'), np.zeros(np.array(img, dtype = 'float32').shape, np.double), 1.0, 0.0, cv2.NORM_MINMAX)
annot = pd.read_csv(os.path.join(tile_annotations, ann))
x,y,w,h = np.int0(np.ceil(annot['BX'])), np.int0(np.ceil(annot['BY'])), np.int0(np.ceil(annot['Width'])), np.int0(np.ceil(annot['Height']))
print(x,y,w,h)
i = 0
for ct in range(len(x)):
    print(x[ct], y[ct], w[ct], h[ct], w[ct]*h[ct])
    big_arr = np.zeros((100,100))
    small_arr = pnn[y[ct]:y[ct]+h[ct],x[ct]:x[ct]+w[ct]]
    big_arr[5:small_arr.shape[0] + 5, 5:small_arr.shape[1] + 5] = small_arr
    cv2.imwrite('/users/ukaipa/PNN/One_img/Cropped_annotations/img_pnn_big_{}.tif'.format(i), big_arr) # if you want the annotations to be cropped out and saved
    i += 1


