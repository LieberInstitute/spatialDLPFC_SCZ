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

# importing one image
# Load the image
# directory path
Image.MAX_IMAGE_PIXELS = None
source_dir = '/dcs04/lieber/marmaypag/spatialDLPFC_SCZ_LIBD4100/processed-data/VistoSeg/captureAreas/V12D07-334_A1.tif'
wfa_img = Image.open(os.path.join(source_dir, img_path))
print("sample num: ", os.path.splitext(os.path.basename(source_dir))[0])
wfa_img.seek(4)
wfa = np.array(wfa_img, dtype = 'uint8')
plt.clf()  # Clear the previous figure
# Compute the histogram from the second peak onwards
second_peak = 25.90  # Your identified second peak value
hist = cv2.calcHist([wfa], [0], None, [256], [second_peak, 256])
hist, bin_edges = np.histogram(image, bins=256, range=(0, 256), density=True)
bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2

# Fit a Rayleigh distribution to your data to find the scale parameter sigma
# For demonstration purposes, let's assume we have an estimated sigma
sigma_est = np.std(image) / np.sqrt((4 - np.pi) / 2)

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
plt.show()

# Perform a Kolmogorov-Smirnov test
statistic, p_value = kstest(hist, 'rayleigh', args=(0, sigma_est))

print(f'KS Statistic: {statistic}, P-value: {p_value}')
# A high P-value (close to 1) indicates a good fit.
