from __future__ import print_function
from skimage.feature import peak_local_max
from skimage.segmentation import find_boundaries, watershed
from scipy import ndimage
import numpy as np
import pandas as pd
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
img_test = pyhere.here('raw-data', 'images', '2_MockPNN', '20220712_VIF_MockPNN_Strong_NTC_C1_Br5182_MLtraining', '20220712_VIF_MockPNN_Strong_Scan1_[12864,50280]_component_data_17.tif')
csv_test = pyhere.here('processed-data', '2_MockPNN', 'Training_tiles', 'Manual_annotations', 'Annotations', '20220712_VIF_MockPNN_Strong_Scan1_[12864,50280]_component_data_17.csv')

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

# writing into a tf record
with tf.io.TFRecordWriter('x1.tfrecord') as writer:
  writer.write(ex.SerializeToString())


# write a loop for this process
obj_name = b'PNN'
example = tf.train.Example(features=tf.train.Features(feature={'img_name': tf.train.Feature(bytes_list = tf.train.BytesList(value = [m.encode('utf-8') for m in file17['img_file_name']])),
                                                               'label': tf.train.Feature(bytes_list = tf.train.BytesList(value = [obj_name])),
                                                               'x1':tf.train.Feature(int64_list = tf.train.Int64List(value = np.int0(np.ceil([x for x in file17['x1']])))),
                                                               'y1':tf.train.Feature(int64_list = tf.train.Int64List(value = np.int0(np.ceil([y for y in file17['y1']])))),
                                                               'w':tf.train.Feature(int64_list = tf.train.Int64List(value = np.int0(np.ceil([y for y in file17['Width']])))),
                                                               'h':tf.train.Feature(int64_list = tf.train.Int64List(value = np.int0(np.ceil([y for y in file17['Height']]))))}))

# for loop to create the format of label,x,y,w,h
for x,y,w,h in zip(file17['x1'], file17['y1'], file17['w'], file17['h']):
    print(x,y,w,h)
# import the df for blood vessels
# create the feature objects for blood vessels
# then create 1 example object for 1 tile
# loop through the rest of the manual annotations slides
# sanity check
# convert to Pascal-VOC
