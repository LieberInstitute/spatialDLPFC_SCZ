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
import tensorflow.compat.v1 as tfc
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
pnn_df = extracting_coords(csv_test, 'PNN')

# tf.compat.v1.disable_eager_execution() # only if eager execution is not needed (eager execution is enabled by default in tf2)
# tf.compat.v1.disable_v2_behavior() # if using a tf1 function
# write a loop for this process
# tf.compat.v1.enable_eager_execution()
object_label = b'PNN'
def create_tf_example(object_label, data_frame):
    example = tf.train.Example(features=tf.train.Features(feature={'img_name': tf.train.Feature(bytes_list = tf.train.BytesList(value = [m.encode('utf-8') for m in data_frame['img_file_name']])),
                                                               'label': tf.train.Feature(bytes_list = tf.train.BytesList(value = [object_label])),
                                                               'x1':tf.train.Feature(int64_list = tf.train.Int64List(value = np.int0(np.ceil([x for x in data_frame['x1']])))),
                                                               'y1':tf.train.Feature(int64_list = tf.train.Int64List(value = np.int0(np.ceil([y for y in data_frame['y1']])))),
                                                               'w':tf.train.Feature(int64_list = tf.train.Int64List(value = np.int0(np.ceil([w for w in data_frame['Width']])))),
                                                               'h':tf.train.Feature(int64_list = tf.train.Int64List(value = np.int0(np.ceil([h for h in data_frame['Height']]))))}))
    return example

# create a tfexample for blood vessels
claudin_csv = '/dcs04/lieber/marmaypag/spatialDLPFC_SCZ_LIBD4100/processed-data/2_MockPNN/Training_tiles/ML_annotations/Annotations/20220712_VIF_MockPNN_Strong_Scan1_[12864,50280]_component_data_17_claudin.csv'
claudin_df = pd.read_csv(claudin_csv)
cla_ex = create_tf_example(b'blood_vessels', claudin_df)


# writing into a tf record
with tf.io.TFRecordWriter('csv1.tfrecord') as writer:
  writer.write(example.SerializeToString())
  # print(tf.train.Example.FromString(example))

# start a tf sess
sess = tf.compat.v1.InteractiveSession()

# read the tf record file
tf.compat.v1.enable_eager_execution()
reader = tf.compat.v1.TFRecordReader()
filename_queue = tf.data.TFRecordDataset(['csv1.tfrecord'])

# displays the tfrecords
for raw_record in filename_queue.take(1):
    exampl = tf.train.Example()
    exampl.ParseFromString(raw_record.numpy())
    print(exampl)

# fix this!
tf.compat.v1.disable_v2_behavior()
_, serialized_example = reader.read(filename_queue)

# nested dict for tf features in the right order (find a feature that accepts a nested dict as an example object)
dict_list=[]
object_label = 'PNN'
for i in range(len(pnn_df['x1'])):
    dict = {'boxes': {'label': tf.train.Feature(bytes_list = tf.train.BytesList(value = [b'PNN'])),
                      'x': tf.train.Feature(int64_list = tf.train.Int64List(value = np.int0(np.ceil([pnn_df['x1'][i]])))),
                      'y': tf.train.Feature(int64_list = tf.train.Int64List(value = np.int0(np.ceil([pnn_df['y1'][i]])))),
                      'w':  tf.train.Feature(int64_list = tf.train.Int64List(value = np.int0(np.ceil([pnn_df['Width'][i]])))),
                      'h': tf.train.Feature(int64_list = tf.train.Int64List(value = np.int0(np.ceil([pnn_df['Height'][i]]))))}}
    dict_list.append(dict)

exm = tf.train.Example(features = tf.train.Feature(bytes_list = tf.train.BytesList(value = [dict_list])))




# works, but in the wrong order
example = tf.train.Example(features=tf.train.Features(feature={'img_name': tf.train.Feature(bytes_list = tf.train.BytesList(value = [m.encode('utf-8') for m in data_frame['img_file_name']])),
                                                               'label': tf.train.Feature(bytes_list = tf.train.BytesList(value = [object_label])),
                                                               'x1':tf.train.Feature(int64_list = tf.train.Int64List(value = np.int0(np.ceil([x for x in data_frame['x1']])))),
                                                               'y1':tf.train.Feature(int64_list = tf.train.Int64List(value = np.int0(np.ceil([y for y in data_frame['y1']])))),
                                                               'w':tf.train.Feature(int64_list = tf.train.Int64List(value = np.int0(np.ceil([w for w in data_frame['Width']])))),
                                                               'h':tf.train.Feature(int64_list = tf.train.Int64List(value = np.int0(np.ceil([h for h in data_frame['Height']]))))}))



# loop through the rest of the manual annotations slides
# sanity check
# convert to Pascal-VOC
