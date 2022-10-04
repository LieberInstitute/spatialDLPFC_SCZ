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

# function to extract all the coordinates from the csv
def extracting_coords(csv_path, object_name):
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

# run this for a single file
file17 = extracting_coords(csv_test, 'PNN')

# create the tf train object for each individual row
img_name_list = tf.train.BytesList(value = [b'20220712_VIF_MockPNN_Strong_Scan1_[12864,50280]_component_data_17'])
x1_list = tf.train.Int64List(value = np.int0(np.ceil([x for x in file17['x1']])))
y1_list = tf.train.Int64List(value = np.int0(np.ceil([y for y in file17['y1']])))
x2_list = tf.train.Int64List(value = np.int0(np.ceil([x for x in file17['x2']])))
y2_list = tf.train.Int64List(value = np.int0(np.ceil([y for y in file17['y2']])))
x3_list = tf.train.Int64List(value = np.int0(np.ceil([x for x in file17['x3']])))
y3_list = tf.train.Int64List(value = np.int0(np.ceil([y for y in file17['y3']])))
x4_list = tf.train.Int64List(value = np.int0(np.ceil([x for x in file17['x4']])))
y4_list = tf.train.Int64List(value = np.int0(np.ceil([y for y in file17['y4']])))

# convert list to features
img_name = tf.train.Feature(bytes_list = img_name_list)
x1 = tf.train.Feature(int64_list = x1_list)
y1 = tf.train.Feature(int64_list = y1_list)
x2 = tf.train.Feature(int64_list = x2_list)
y2 = tf.train.Feature(int64_list = y2_list)
x3 = tf.train.Feature(int64_list = x3_list)
y3 = tf.train.Feature(int64_list = y3_list)
x4 = tf.train.Feature(int64_list = x4_list)
y4 = tf.train.Feature(int64_list = y4_list)

# write a loop for this process
# create a function for one file
# call the function for all the csvs
# check if this works
