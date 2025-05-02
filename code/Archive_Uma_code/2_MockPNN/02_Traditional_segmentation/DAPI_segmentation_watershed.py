# import the necessary packages
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
import pyhere

# construct the argument parse and parse the arguments
# ap = argparse.ArgumentParser()
# ap.add_argument("-i", "--image", required=True,
# 	help="path to input image")
# args = vars(ap.parse_args())

# define the image path
img_test = pyhere.here('raw-data', 'images', '2_MockPNN', 'Training_tiles', '20220712_VIF_MockPNN_Strong_Scan1_[6384,53057]_component_data_11.tif')
csv_test = '/dcs04/lieber/marmaypag/spatialDLPFC_SCZ_LIBD4100/processed-data/2_MockPNN/Training_tiles/Manual_annotations/Annotations/20220712_VIF_MockPNN_Strong_Scan1_[6384,53057]_component_data_11.csv'
# loop through the whole directory to segment only DAPI
img_dir = '/dcs04/lieber/marmaypag/spatialDLPFC_SCZ_LIBD4100/raw-data/images/RealPNN/round1/20220814_VIF_PNN_S2_SCZ/'
csv_dst = '/dcs04/lieber/marmaypag/spatialDLPFC_SCZ_LIBD4100/processed-data/RealPNN/DAPI_segmentations/Image_csvs/'
img_dst = '/dcs04/lieber/marmaypag/spatialDLPFC_SCZ_LIBD4100/processed-data/RealPNN/DAPI_segmentations/Image_segmentations/'

for img_name in os.listdir(img_dir):
    if img_name.endswith('.tif'):
        # print(int(img_name.split('_')[8].split('.')[0]), int(csv_name.split('_')[8].split('.')[0]))
        print(img_name)
        dapi, dapi_clr = read_norm(os.path.join(img_dir, img_name), 0)
        # csv = manual_annot(os.path.join(csv_dir, csv_name))
        # print(len(csv))
        shifted, thresh, gray = morph_transform(dapi_clr)
        labels = find_labels(thresh)
        dpx, dpy, dpw, dph, area, segmented_dapi = draw_rect_dapi(labels, gray, dapi_clr)
        img_info_dapi = create_df(dpx, dpy, dpw, dph, area, os.path.join(img_dir, img_name), 'DAPI')
        img_info_dapi.to_csv(path_or_buf = (csv_dst + img_name.split('.')[0] + '.csv')) # df to csv and save it in the csv_dst folder
        cv2.imwrite((img_dst + img_name.split('.')[0] + '.tif'), segmented_dapi) # save the segmented images in the img_dst folder


fig,ax = plt.subplots(figsize = (20,20))
ax.imshow(segmented_dapi)
fig.show()



# dapi segmentations for a single image
dapi, dapi_clr = read_norm(img_test, 0)
# fig,ax = plt.subplots(figsize = (20,20))
# ax.imshow(dapi, cmap = 'gray')
# fig.show()
dpx, dpy, dpw, dph, area, ws_img_bb = draw_contours(dapi, 0, contours = None,  color = None, thickness = None, dapi_clr = dapi_clr)
# fig,ax = plt.subplots(figsize = (20,20))
# ax.imshow(ws_img_bb)
# fig.show()
plot_img(dapi, ws_img_bb)
img_info_dapi = create_df(dpx, dpy, dpw, dph, area, img_test, 'DAPI') # populate the data in the dataframe


# for loop for looping through the total num of BB/len of the csv
dapi_box_means = []
for bb in range(len(img_info_dapi)-191):
    # print("entered the loop", bb)
    new_im = np.zeros(dapi.shape, np.double)
    # print("created a new image", bb)
    rect_img = draw_single_rect(img_info_dapi.iloc[bb], new_im)
    # print("drew a rectangle", bb)
    locs = np.where(rect_img == 255)
    # print("found pix == 255", bb)
    # print("Number of DAPI",len(locs[0]))
    dapi_box, x_, y_ = [], [], []
    # print("entering 2nd loop")
    for x,y in zip(locs[0], locs[1]):
        if rect_img[x, y] == 255:
            # print(x,y, dapi[x,y])
            x_.append(x)
            y_.append(y)
            print("appending", bb)
            coords = np.array((x_,y_))
            dapi_box.append(dapi[x,y])
    dapi_box = np.array(dapi_box)
    print("final appending", bb)
    coords.append(coords)
    # print(dapi_box.mean())
    dapi_box_means.append(dapi_box.mean())

# convert dapi means from list to array
dapi_box_mean = np.array(dapi_box_means)

