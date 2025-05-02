import numpy as np
import numpy.ma as ma
import pandas as pd
import PIL
import pyhere
from pathlib import Path
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
from scipy.signal import convolve2d
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
import matplotlib.patches as patches


# define the training tiles directory
img_dir = pyhere.here('raw-data', 'images', '2_MockPNN', 'Training_tiles')
img_dir_NTC = pyhere.here('raw-data', 'images', '2_MockPNN', '20220712_VIF_MockPNN_Strong_NTC_C1_Br5182_MLtraining')
img_NTC = pyhere.here('raw-data', 'images', '2_MockPNN', '20220712_VIF_MockPNN_Strong_NTC_C1_Br5182_MLtraining', '20220712_VIF_MockPNN_Strong_NTC_Scan1_[11013,50974]_component_data.tif')
img_SCZ = pyhere.here('raw-data', 'images', '2_MockPNN', '20220712_VIF_MockPNN_Strong_SCZ_C1_Br2039_MLtraining', '20220712_VIF_MockPNN_Strong_SCZ_Scan1_[10629,49106]_component_data.tif')
img_test = pyhere.here('raw-data', 'images', '2_MockPNN', 'Training_tiles', '20220712_VIF_MockPNN_Strong_Scan1_[6384,53057]_component_data_11.tif') #20220712_VIF_MockPNN_Strong_Scan1_[10087,51668]_component_data_01, 20220712_VIF_MockPNN_Strong_Scan1_[6925,49106]_component_data_24.tif
csv_test = pyhere.here('processed-data', '2_MockPNN', 'Training_tiles', 'Manual_annotations', 'Annotations', '20220712_VIF_MockPNN_Strong_Scan1_[6384,53057]_component_data_11.csv') #20220712_VIF_MockPNN_Strong_Scan1_[10087,51668]_component_data_01

df_11 = pd.read_csv(csv_test)

