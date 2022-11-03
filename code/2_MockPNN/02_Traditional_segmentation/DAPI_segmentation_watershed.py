# import the necessary packages
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


# construct the argument parse and parse the arguments
# ap = argparse.ArgumentParser()
# ap.add_argument("-i", "--image", required=True,
# 	help="path to input image")
# args = vars(ap.parse_args())

# define the image path
img_test = pyhere.here('raw-data', 'images', '2_MockPNN', 'Training_tiles', '20220712_VIF_MockPNN_Strong_Scan1_[6925,49106]_component_data_24.tif')

# load and preprocess the image
# img_dapi = Image.open(img_test)
# img_dapi.seek(0) # channel 0 = DAPI
# dapi = cv2.normalize(np.array(img_dapi, dtype = 'float32'), np.zeros(np.array(img_dapi, dtype = 'float32').shape, np.double), 1.0, 0.0, cv2.NORM_MINMAX)
# dapi_clr = skimage.color.gray2rgb((np.array((dapi * 255), dtype = np.uint8))) # convert to color to draw colored bb

dapi = read_norm(img_test, 0)

# perform pyramid mean shifting
def morph_transform(image):
	shifted = cv2.pyrMeanShiftFiltering(image, 21, 51) #dapi_clr
	gray = cv2.cvtColor(shifted, cv2.COLOR_BGR2GRAY)
	thresh = cv2.threshold(gray, 0, 255, cv2.THRESH_BINARY | cv2.THRESH_OTSU)[1]
	fig,ax = plt.subplots(figsize = (20,20))
	ax.imshow(dapi_clr)
	fig.show()
	return shifted, thresh

# Otsu's thresholding
# gray = cv2.cvtColor(shifted, cv2.COLOR_BGR2GRAY)
# thresh = cv2.threshold(gray, 0, 255, cv2.THRESH_BINARY | cv2.THRESH_OTSU)[1]
# fig,ax = plt.subplots(figsize = (20,20))
# ax.imshow(thresh)
# fig.show()


# find labels in the image
D = ndimage.distance_transform_edt(thresh) # Euclidean distance from binary to nearest 0-pixel
localMax = peak_local_max(D, indices=False, min_distance=5, labels=thresh) # find the local maxima for all the individual objects
markers = ndimage.label(localMax, structure=np.ones((3, 3)))[0] # 8-connectivity connected component analysis
labels = watershed(-D, markers, mask=thresh)
print("[INFO] {} unique segments found".format(len(np.unique(labels)) - 1))

# extract the watershed algorithm labels
dpx, dpy, dpw, dph = [], [], [], []
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
	dpx.append(x)
	dpy.append(y)
	dpw.append(w)
	dph.append(h)
	ws_img_bb = cv2.rectangle(dapi_clr, (x,y), (x+w, y+h), (255,0,0), 2) # draw BB

# Plot the segmentation result
fig,ax = plt.subplots(figsize = (20,20))
ax.imshow(ws_img_bb)
ax.title.set_text('Watershed Segmentation')
fig.show()


# Populate the data in the dataframe
col_names = ['img_file_name','type_of_object_str', 'x1', 'y1', 'Width', 'Height', 'total_number_dapi']
object_name = 'DAPI' # name of the objects stored in the dataframe
file_name = os.path.basename(img_test) # image file name

dict = {col_names[0]: file_name, col_names[1]: object_name, col_names[2]: dpx, col_names[3]: dpy, col_names[4]: dpw, col_names[5]: dph, col_names[6]: len(dpx)}
img_info_dapi = pd.DataFrame(dict, columns = col_names)
img_info_dapi['x2'] = img_info_dapi['x1'] + img_info_dapi['Width']
img_info_dapi['y2'], img_info_dapi['x3'] = img_info_dapi['y1'], img_info_dapi['x1']
img_info_dapi['y3'] = img_info_dapi['y1'] + img_info_dapi['Height']
img_info_dapi['x4'], img_info_dapi['y4'] = img_info_dapi['x2'], img_info_dapi['y3']
img_info_dapi = img_info_dapi[['img_file_name', 'type_of_object_str', 'x1', 'y1', 'x2', 'y2', 'x3', 'y3', 'x4', 'y4', 'Width', 'Height', 'total_number_dapi']]
