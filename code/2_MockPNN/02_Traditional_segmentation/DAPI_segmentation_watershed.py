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

# load and preprocess the image
# img_dapi = Image.open(img_test)
# img_dapi.seek(0) # channel 0 = DAPI
# dapi = cv2.normalize(np.array(img_dapi, dtype = 'float32'), np.zeros(np.array(img_dapi, dtype = 'float32').shape, np.double), 1.0, 0.0, cv2.NORM_MINMAX)
# dapi_clr = skimage.color.gray2rgb((np.array((dapi * 255), dtype = np.uint8))) # convert to color to draw colored bb

dapi, dapi_clr = read_norm(img_test, 0)
fig,ax = plt.subplots(figsize = (20,20))
ax.imshow(dapi)
fig.show()


# perform pyramid mean shifting
def morph_transform(image_clr):
    shifted = cv2.pyrMeanShiftFiltering(image_clr, 21, 51) #dapi_clr
    gray = cv2.cvtColor(shifted, cv2.COLOR_BGR2GRAY)
    thresh = cv2.threshold(gray, 0, 255, cv2.THRESH_BINARY | cv2.THRESH_OTSU)[1]
    fig,ax = plt.subplots(figsize = (20,20))
    ax.imshow(image_clr)
    fig.show()
    return shifted, thresh, gray

shifted, thresh, gray = morph_transform(dapi_clr)



# find labels in the image
def find_labels(threshold):
    D = ndimage.distance_transform_edt(threshold) # Euclidean distance from binary to nearest 0-pixel
    localMax = peak_local_max(D, indices=False, min_distance=5, labels=threshold) # find the local maxima for all the individual objects
    markers = ndimage.label(localMax, structure=np.ones((3, 3)))[0] # 8-connectivity connected component analysis
    labels = watershed(-D, markers, mask=threshold)
    print("{} unique segments found".format(len(np.unique(labels)) - 1))
    return labels

labels = find_labels(thresh)

# extract the watershed algorithm labels
def draw_rect_dapi(labels, gray, dapi): # add area
    dpx, dpy, dpw, dph, area = [], [], [], [], []
    for label in np.unique(labels):
        if label == 0: # label marked 0 are background
            continue
        mask = np.zeros(gray.shape, dtype="uint8") # create masks that only have the detected labels as foreground and 0 as background
        mask[labels == label] = 255
        # detect contours in the mask and grab the largest one
        cnts = cv2.findContours(mask.copy(), cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_SIMPLE) # detect the watershed contours
        cnts = imutils.grab_contours(cnts) # extract only the contours
        c = max(cnts, key=cv2.contourArea) # get the area
        x,y,w,h = cv2.boundingRect(c) # BB coordinates
        area.append(cv2.contourArea(c))
        dpx.append(x)
        dpy.append(y)
        dpw.append(w)
        dph.append(h)
        ws_img_bb = cv2.rectangle(dapi, (x,y), (x+w, y+h), (255,0,0), 2) # draw BB
    return dpx, dpy, dpw, dph, area, ws_img_bb

dpx, dpy, dpw, dph, area, segmented_dapi = draw_rect_dapi(labels, gray, dapi_clr)

# Plot the segmentation result
fig,ax = plt.subplots(figsize = (20,20))
ax.imshow(segmented_dapi)
ax.title.set_text('Watershed Segmentation')
fig.show()


# Populate the data in the dataframe
img_info_dapi = create_df(dpx, dpy, dpw, dph, area, img_test, 'DAPI')

# for loop for looping through the total num of BB/len of the csv
dapi_box_means = [] # not looping through all of the BBs
for bb in range(len(img_info_dapi)):
    print("entered the loop", bb)
    new_im = np.zeros(dapi.shape, np.double)
    print("created a new image", bb)
    rect_img = draw_single_rect(img_info_dapi.iloc[bb], new_im)
    print("drew a rectangle", bb)
    locs = np.where(rect_img == 255)
    print("found pix == 255", bb)
    print("Number of DAPI",len(locs[0]))
    dapi_box = []
    # print("entering 2nd loop")
    for x,y in zip(locs[0], locs[1]):
        if rect_img[x, y] == 255:
            # print(x,y, dapi[x,y])
            dapi_box.append(dapi[x,y])
    dapi_box = np.array(dapi_box)
    print(dapi_box.mean())
    dapi_box_means.append(dapi_box.mean())
