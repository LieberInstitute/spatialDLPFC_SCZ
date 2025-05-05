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
dst_dir_wfa = '/dcs04/lieber/marmaypag/spatialDLPFC_SCZ_LIBD4100/processed-data/RealPNN/single_channels_segmented/WFA/different_thresh/'
dst_dir_hist_all_values = '/dcs04/lieber/marmaypag/spatialDLPFC_SCZ_LIBD4100/processed-data/RealPNN/single_channels_segmented/WFA/histograms/'
dst_dir_samui = '/dcs04/lieber/marmaypag/spatialDLPFC_SCZ_LIBD4100/processed-data/RealPNN/single_channels_segmented/WFA/samui_samples/'
dst_dir_hist_cutoff_values = '/dcs04/lieber/marmaypag/spatialDLPFC_SCZ_LIBD4100/processed-data/RealPNN/single_channels_segmented/WFA/histograms_for_all_tissue_sections_only_cutoff_values/'


# importing one image
'''
This test needs to be done to check if the image data follows the Rayleigh distribution and the Kolmogorov-Smirnov (KS) test
The KS statistic being close to 1 and the P-value being 0 indicates that we can reject the null hypothesis that the data comes from a Rayleigh distribution
The shape of histogram suggests that the distribution of pixel intensities in the images is not well-modeled by a Rayleigh distribution.
Calculating the statistic with the whole histogram where the first peak starts at 0, leads to the failure of the test and acceptance of the null hypothesis
The data fails the test even after removing the first peak at 0
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
first_peak_index = np.argmax(hist) # find the first peak
first_peak_value = bin_edges[first_peak_index] # the first peak is usually 0
pixel_values = np.arange(len(hist)) # Dynamically calculate exclusion_range based on standard deviation
std_dev = np.std(pixel_values)
exclusion_range = int(0 * std_dev)  # Adjust the multiplier as needed
# exclusion_range = 20 # adjust this value to exclude a region around the first peak
hist[range(max(0, first_peak_index - exclusion_range), min(len(hist), first_peak_index + exclusion_range + 1))] = 0 
second_peak_index = np.argmax(hist) # Find the index of the second peak
second_peak_value = bin_edges[second_peak_index] # Get the value of the second peak
second_peak_y_value = hist[second_peak_index]
bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2
# Now create a modified histogram that starts from the second peak
second_peak_bin = second_peak_index  # This should be the bin index, not the bin value
mod_hist = hist[second_peak_bin:]
mod_bins = bin_edges[second_peak_bin:]
# Estimate the Rayleigh scale parameter sigma from the modified histogram
# Assuming the mode (peak) of Rayleigh distribution is at the second peak
sigma_est = mod_bins[np.argmax(mod_hist)] / np.sqrt(2 * np.log(2))
# Generate the Rayleigh distribution curve for the modified histogram
rayleigh_curve = rayleigh.pdf(mod_bins, scale=sigma_est)
# Plot the modified histogram and the Rayleigh distribution curve for comparison
plt.figure(figsize=(10, 6))
plt.bar(mod_bins[:-1], mod_hist / np.sum(mod_hist), width=np.diff(mod_bins), color='blue', alpha=0.7, label='Modified Image Histogram')
plt.plot(mod_bins, rayleigh_curve, 'r-', label=f'Rayleigh Distribution (σ={sigma_est:.2f})')
plt.title('Modified Image Histogram vs. Rayleigh Distribution')
plt.xlabel('Pixel Intensity')
plt.ylabel('Frequency')
plt.legend()
plt.savefig(os.path.join(dst_dir_wfa, os.path.splitext(os.path.basename(source_dir))[0] + '_test_for_rayleigh_Ks_from_second_peak.png')) # save the figures

# Normalize the histogram to form a probability distribution
norm_mod_hist = mod_hist / np.sum(mod_hist)
cumulative_data = np.cumsum(norm_mod_hist)

# The KS test compares the CDFs, so we need to use the cumulative data
ks_statistic, p_value = kstest(cumulative_data, 'rayleigh', args=(0, sigma_est))

print(f'KS Statistic: {ks_statistic}, P-value: {p_value}')



############ actual work using the second peak on V12D07-334 A1 image ######## THIs WORKs!!!!
import cv2
import numpy as np
import matplotlib.pyplot as plt

# Load the image
# directory path
Image.MAX_IMAGE_PIXELS = None
source_dir = '/dcs04/lieber/marmaypag/spatialDLPFC_SCZ_LIBD4100/processed-data/VistoSeg/captureAreas/V12D07-334_A1.tif'
wfa_img = Image.open(source_dir)
print("sample num: ", os.path.splitext(os.path.basename(source_dir))[0])
wfa_img.seek(4)
wfa = np.array(wfa_img, dtype = 'uint8')
plt.clf()  # Clear the previous figure
# Compute the histogram from the second peak onwards
second_peak = 25.90  # Your identified second peak value
hist = cv2.calcHist([wfa], [0], None, [256], [second_peak, 256])

# Initialize variables for the T-point algorithm
M = int(np.ceil(second_peak))
L = len(hist) - 1
best_error = np.inf
best_k = None

# Iterate over each possible breakpoint k
for k in range(M+1, L):
    # Perform linear regression on both sides of k
    left_slope, left_intercept = np.polyfit(bins[M:k], hist[M:k], 1)
    right_slope, right_intercept = np.polyfit(bins[k+1:L+1], hist[k+1:L+1], 1)   
    # Calculate predicted values
    left_predicted = left_slope * bins[M:k] + left_intercept
    right_predicted = right_slope * bins[k+1:L+1] + right_intercept    
    # Calculate errors
    left_errors = (hist[M:k] - left_predicted) ** 2
    right_errors = (hist[k+1:L+1] - right_predicted) ** 2    
    # Sum the errors to get the total error for this breakpoint k
    total_error = np.sum(left_errors) + np.sum(right_errors)    
    # Check if this is the best breakpoint so far
    if total_error < best_error:
        best_error = total_error
        best_k = k

# The best breakpoint is our threshold T
T = best_k

# Apply the threshold to segment the image
_, thresholded_image = cv2.threshold(wfa, T, 255, cv2.THRESH_BINARY)

# Display the segmented image
plt.imshow(thresholded_image, cmap='gray')
plt.show()
