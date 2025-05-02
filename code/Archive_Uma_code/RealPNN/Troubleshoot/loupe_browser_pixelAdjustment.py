from scipy import ndimage
import imutils
import numpy as np
import argparse
from argparse import ArgumentParser
import pyhere
from pyhere import here
from pathlib import Path
import pandas as pd
import PIL
from PIL import Image
import os
import matplotlib.pyplot as plt
import cv2
import scipy
from scipy.spatial.distance import *
import skimage
from skimage import feature, segmentation, draw, measure, morphology

im_A1 = '/dcs04/lieber/marmaypag/spatialDLPFC_SCZ_LIBD4100/processed-data/VistoSeg/captureAreas/V12D07-334_A1.tif'

Image.MAX_IMAGE_PIXELS = None
im = Image.open(im_A1)
im.seek(3)
im_arr = np.array(im, dtype = 'uint8')
fig,ax = plt.subplots(figsize = (20,20))
ax.imshow(im_arr, cmap = 'gray') #
fig.show()

# define the alpha and beta
alpha = 1.5 # Contrast control
beta = 10 # Brightness control

# call convertScaleAbs function
adjusted = cv2.convertScaleAbs(im_arr, alpha=alpha, beta=beta)

# define the contrast and brightness value
contrast = 100. # Contrast control ( 0 to 127)
brightness = 70. # Brightness control (0-100)

# call addWeighted function. use beta = 0 to effectively only operate on one image
out = cv2.addWeighted( im_arr, contrast, im_arr, 0, brightness)

cv2.imwrite('/dcs04/lieber/marmaypag/spatialDLPFC_SCZ_LIBD4100/processed-data/RealPNN/single_channels_segmented/V12D07-334_A1_neun_pix_adj.tif', out)

fig,ax = plt.subplots(figsize = (20,20))
ax.imshow(out, cmap = 'gray')
fig.show()

