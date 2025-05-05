'''
For Stitched Visium-IF tissue sections from VistoSeg SplitSlide output
Channel0 = AF
Channel1 = Claudin - 5,
Channel2 = DAPI,
Channel3 = NeuN,
Channel4 = WFA
'''

'''
Slide 1 = V12F14-053
Slide 2 = V12F14-057
Slide 3 = V12D07-334
Slide 4 = V13M06_279
Slide 5 = V13M06_280
Slide 6 = 
Slide 7 = 
Slide 8 = 
Slide 9 = 
Slide 10 = 
Slide 11 = 
Slide 12 = 
Slide 13 = 
Slide 14 = 
Slide 15 = 
Slide 16 = 

'''


import numpy as np
import pyhere
from pyhere import here
from pylab import xticks
from pathlib import Path
import pandas as pd
import PIL
from PIL import Image, ImageFont, ImageDraw, ImageEnhance, ImageSequence
import os
import matplotlib
import matplotlib.pyplot as plt
import sys
import cv2
import math
import scipy


#### splislide troubleshoot
img = '/dcs04/lieber/marmaypag/spatialDLPFC_SCZ_LIBD4100/raw-data/images/RealPNN/Inform_0%_overlap/V12D07-334_channel_1.tif'
Image.MAX_IMAGE_PIXELS = None
im = Image.open(img)
arr = np.array(im, dtype = 'uint8')
arr_c = cv2.cvtColor(arr,cv2.COLOR_BGR2RGB)
gray = cv2.cvtColor(arr_c,cv2.COLOR_RGB2GRAY)
_,thresh = cv2.threshold(gray, np.mean(gray), 255, cv2.THRESH_BINARY | cv2.THRESH_OTSU) #_INV
fig,ax = plt.subplots(figsize = (20,20))
ax.imshow(thresh, cmap = 'gray')
fig.show()

