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
# define the contrast and brightness value
contrast = 100. # Contrast control ( 0 to 127)
brightness = 70. # Brightness control (0-100)
im.seek(3) # 3 = neun
im_arr = np.array(im, dtype = 'uint8')
for i in range(5):
    im.seek(i)
    im_arr = np.array(im, dtype = 'uint8')
    out = cv2.addWeighted( im_arr, contrast, im_arr, 0, brightness)
    fig,ax = plt.subplots(figsize = (20,20))
    ax.imshow(out, cmap = 'gray')
    fig.show()

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

import matplotlib
import matplotlib.pyplot as plt
import numpy as np

from skimage import data, img_as_float
from skimage import exposure

matplotlib.rcParams['font.size'] = 8


def plot_img_and_hist(image, axes, bins=256):
    """Plot an image along with its histogram and cumulative histogram.

    """
    image = img_as_float(image)
    ax_img, ax_hist = axes
    ax_cdf = ax_hist.twinx()
    # Display image
    ax_img.imshow(image, cmap=plt.cm.gray)
    ax_img.set_axis_off()
    # Display histogram
    ax_hist.hist(image.ravel(), bins=bins, histtype='step', color='black')
    ax_hist.ticklabel_format(axis='y', style='scientific', scilimits=(0, 0))
    ax_hist.set_xlabel('Pixel intensity')
    ax_hist.set_xlim(0, 1)
    ax_hist.set_yticks([])
    # Display cumulative distribution
    img_cdf, bins = exposure.cumulative_distribution(image, bins)
    ax_cdf.plot(bins, img_cdf, 'r')
    ax_cdf.set_yticks([])
    return ax_img, ax_hist, ax_cdf


# Load an example image
Image.MAX_IMAGE_PIXELS = None
im = Image.open(im_A1)
im.seek(3) # 3 = neun
im_arr = np.array(im, dtype = 'uint8')
img = im_arr

# Gamma
gamma_corrected = exposure.adjust_gamma(img, 2)

# Logarithmic
logarithmic_corrected = exposure.adjust_log(img, 1)

# Display results
fig = plt.figure(figsize=(8, 5))
axes = np.zeros((2, 3), dtype=object)
axes[0, 0] = plt.subplot(2, 3, 1)
axes[0, 1] = plt.subplot(2, 3, 2, sharex=axes[0, 0], sharey=axes[0, 0])
axes[0, 2] = plt.subplot(2, 3, 3, sharex=axes[0, 0], sharey=axes[0, 0])
axes[1, 0] = plt.subplot(2, 3, 4)
axes[1, 1] = plt.subplot(2, 3, 5)
axes[1, 2] = plt.subplot(2, 3, 6)

ax_img, ax_hist, ax_cdf = plot_img_and_hist(img, axes[:, 0])
ax_img.set_title('Low contrast image')

y_min, y_max = ax_hist.get_ylim()
ax_hist.set_ylabel('Number of pixels')
ax_hist.set_yticks(np.linspace(0, y_max, 5))

ax_img, ax_hist, ax_cdf = plot_img_and_hist(gamma_corrected, axes[:, 1])
ax_img.set_title('Gamma correction')

ax_img, ax_hist, ax_cdf = plot_img_and_hist(logarithmic_corrected, axes[:, 2])
ax_img.set_title('Logarithmic correction')

ax_cdf.set_ylabel('Fraction of total intensity')
ax_cdf.set_yticks(np.linspace(0, 1, 5))

# prevent overlap of y-axis labels
fig.tight_layout()
plt.show()


fig,ax = plt.subplots(figsize = (20,20))
ax.imshow(gamma_corrected, cmap = 'gray')
fig.show()
