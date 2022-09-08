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

# Read the corresponding tiles to its annotations file
for tile_name in os.listdir(training_tiles_path):
    for ann_name in os.listdir(tile_annotations):
        if int(re.split(r'[_.]',os.path.basename(t))[8]) == int(re.split(r'[_.]',os.path.basename(a))[8]):
            print(tile_name,ann_name)

# read the WFA channel
# crop the PNNs from WFA
# resize them
# crop the BV from WFA
# detect contours for both
# get coordinates for both
# put them all in 1 csv

