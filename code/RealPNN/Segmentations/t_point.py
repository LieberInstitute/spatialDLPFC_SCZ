import PIL
from PIL import Image, ImageFont, ImageDraw, ImageEnhance, ImageSequence
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import rayleigh, kstest
from scipy.ndimage import gaussian_gradient_magnitude
from scipy.optimize import curve_fit
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
import sys
import cv2
import math
import scipy
import os

# directory path
Image.MAX_IMAGE_PIXELS = None
source_dir = '/dcs04/lieber/marmaypag/spatialDLPFC_SCZ_LIBD4100/processed-data/VistoSeg/captureAreas/'
img_dir = '/dcs04/lieber/marmaypag/spatialDLPFC_SCZ_LIBD4100/processed-data/VistoSeg/captureAreas/'
dst_dir_wfa = '/dcs04/lieber/marmaypag/spatialDLPFC_SCZ_LIBD4100/processed-data/RealPNN/single_channels_segmented/WFA/test/'
dst_dir_hist_all_values = '/dcs04/lieber/marmaypag/spatialDLPFC_SCZ_LIBD4100/processed-data/RealPNN/single_channels_segmented/WFA/histograms/'
dst_dir_samui = '/dcs04/lieber/marmaypag/spatialDLPFC_SCZ_LIBD4100/processed-data/RealPNN/single_channels_segmented/WFA/samui_samples/'
dst_dir_hist_cutoff_values = '/dcs04/lieber/marmaypag/spatialDLPFC_SCZ_LIBD4100/processed-data/RealPNN/single_channels_segmented/WFA/histograms_for_all_tissue_sections_only_cutoff_values/'


# importing one image
'''
This test needs to be done to check if the image data follows the Rayleigh distribution and the Kolmogorov-Smirnov (KS) test
The KS statistic being close to 1 and the P-value being 0 indicates that we can reject the null hypothesis that the data comes from a Rayleigh distribution
The shape of histogram suggests that the distribution of pixel intensities in the images is not well-modeled by a Rayleigh distribution.
'''
Image.MAX_IMAGE_PIXELS = None
source_dir = '/dcs04/lieber/marmaypag/spatialDLPFC_SCZ_LIBD4100/processed-data/VistoSeg/captureAreas/V12D07-334_A1.tif'
wfa_img = Image.open(source_dir)
print("sample num: ", os.path.splitext(os.path.basename(source_dir))[0])
wfa_img.seek(4)
wfa = np.array(wfa_img, dtype = 'uint8')
plt.clf()  # Clear the previous figure
# Compute the histogram from the second peak onwards
second_peak = 25.90  # Your identified second peak value
# hist = cv2.calcHist([wfa], [0], None, [256], [second_peak, 256])

hist, bin_edges = np.histogram(wfa, bins=256, range=(0, 256), density=True)
bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2

# Fit a Rayleigh distribution to your data to find the scale parameter sigma
# For demonstration purposes, let's assume we have an estimated sigma
sigma_est = np.std(wfa) / np.sqrt((4 - np.pi) / 2)

# Generate a Rayleigh distribution curve
rayleigh_curve = rayleigh.pdf(bin_centers, scale=sigma_est)

# Plot the histogram and the Rayleigh distribution curve
plt.figure(figsize=(10, 6))
plt.bar(bin_centers, hist, width=bin_edges[1] - bin_edges[0], color='blue', alpha=0.6, label='Image Histogram')
plt.plot(bin_centers, rayleigh_curve, 'r-', label=f'Rayleigh Distribution (σ={sigma_est:.2f})')
plt.title('Image Histogram vs. Rayleigh Distribution')
plt.xlabel('Pixel Intensity')
plt.ylabel('Frequency')
plt.legend()
plt.savefig(os.path.join(dst_dir_wfa, os.path.splitext(os.path.basename(source_dir))[0] + '_test_for_rayleigh_Ks.png')) # save the figures

# Perform a Kolmogorov-Smirnov test
statistic, p_value = kstest(hist, 'rayleigh', args=(0, sigma_est))

print(f'KS Statistic: {statistic}, P-value: {p_value}') # KS Statistic: 0.9998436557747595, P-value: 0.0
# A high P-value (close to 1) indicates a good fit.
