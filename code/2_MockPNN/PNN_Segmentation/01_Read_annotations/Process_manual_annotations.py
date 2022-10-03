from __future__ import print_function
from skimage.feature import peak_local_max
from skimage.segmentation import find_boundaries, watershed
from scipy import ndimage
import numpy as np
from PIL import Image, ImageFont, ImageDraw, ImageEnhance, ImageSequence
import os
import skimage
import matplotlib
import matplotlib.pyplot as plt
import argparse
import imutils
import cv2
import tensorflow as tf
import keras

# import one of the manual annotation csvs--pref which has lesser rows
img_test = '/dcs04/lieber/marmaypag/spatialDLPFC_SCZ_LIBD4100/raw-data/images/2_MockPNN/Training_tiles/20220712_VIF_MockPNN_Strong_Scan1_[12864,50280]_component_data_17.tif'
csv_test = '/dcs04/lieber/marmaypag/spatialDLPFC_SCZ_LIBD4100/processed-data/2_MockPNN/Training_tiles/Manual_annotations/Annotations/20220712_VIF_MockPNN_Strong_Scan1_[12864,50280]_component_data_17.csv'
test = pd.read_csv(csv_test)

# create the tf train object for each individual row
# write a loop for this process
# create a function for one file
# call the function for all the csvs
# check if this works
