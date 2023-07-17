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

