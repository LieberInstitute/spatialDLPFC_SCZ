import os
import re
import cv2
import numpy as np
import pandas as pd
import PIL
import matplotlib
import matplotlib.pyplot as plt
from PIL import Image, ImageFont, ImageDraw, ImageEnhance, ImageSequence

# Define the directories
csv_cv_dir = '/dcs04/lieber/marmaypag/spatialDLPFC_SCZ_LIBD4100/raw-data/images/2_MockPNN/Training_tiles_cv_segmentation_csv_files/'
csv_ma_dir = '/dcs04/lieber/marmaypag/spatialDLPFC_SCZ_LIBD4100/processed-data/Experimentation_archive/2_MockPNN/Training_tiles/Manual_annotations/Annotations_in_pixels/'
img_raw_dir = '/dcs04/lieber/marmaypag/spatialDLPFC_SCZ_LIBD4100/raw-data/images/2_MockPNN/Training_tiles_raw_no_annotations/'

# Get a list of filenames in each directory
csv_cv_files = os.listdir(csv_cv_dir)
csv_ma_files = os.listdir(csv_ma_dir)
img_raw_files = os.listdir(img_raw_dir)

# Function to extract numbers within square brackets from a filename
def extract_numbers(filename):
    match = re.search(r'\[(\d+),(\d+)\]', filename)
    if match:
        return tuple(map(int, match.groups()))
    else:
        return None

# Dictionary to store matching filenames
matching_files = {}

# Iterate over files in each directory and find matches
for csv_cv_file in csv_cv_files:
    numbers_cv = extract_numbers(csv_cv_file)
    if numbers_cv:
        for csv_ma_file in csv_ma_files:
            numbers_ma = extract_numbers(csv_ma_file)
            if numbers_ma and numbers_cv == numbers_ma:
                for img_raw_file in img_raw_files:
                    numbers_img_raw = extract_numbers(img_raw_file)
                    if numbers_img_raw and numbers_cv == numbers_img_raw:
                        matching_files[numbers_cv] = {
                            'csv_cv': os.path.join(csv_cv_dir, csv_cv_file),
                            'csv_ma': os.path.join(csv_ma_dir, csv_ma_file),
                            'img_raw': os.path.join(img_raw_dir, img_raw_file)
                        }
                        break


# function to find mean pixel intensity of pixels with abox
def calculate_mean_intensity(row):
    x1, y1, x4, y4 = int(row['x1']), int(row['y1']), int(row['x4']), int(row['y4'])
    roi = image[y1:y4, x1:x4]
    mean_intensity = cv2.mean(roi)[0]
    return mean_intensity


# List to store all mean intensity values from all files
all_mean_intensity_values = []
# Now, matching_files dictionary contains matching filenames grouped by the numbers within square brackets
# You can perform further operations on the matched files
for numbers, files in matching_files.items():
    # Read the CSV files
    file_cv_df = pd.read_csv(files['csv_cv'])
    file_ma_df = pd.read_csv(files['csv_ma'])
    # Read the image
    image = Image.open(files['img_raw'])
    image.seek(3)
    image = cv2.normalize(np.array(image, dtype='float32'), np.zeros(np.array(image, dtype='float32').shape, np.double), 1.0, 0.0, cv2.NORM_MINMAX)
    wfa_c = cv2.cvtColor(np.array(image * 255, dtype=np.uint8), cv2.COLOR_BGR2RGB)
    # Create an array to store mean intensities
    mean_intensity_values = []
   # Create a list to store normalized mean intensities for the current file
    normalized_mean_intensity_values = []
    # Iterate over bounding boxes and compute mean pixel intensities
    # for index, row in file_ma_df.iterrows():
    #     x1, y1, x4, y4 = int(row['x1']), int(row['y1']), int(row['x4']), int(row['y4'])
    #     # Extract region of interest (ROI)
    #     roi = image[y1:y4, x1:x4]
    #     # Compute mean pixel intensities
    #     mean_intensity = cv2.mean(roi)[0]
    #     mean_intensity_values.append(mean_intensity)
    # Extend the list with mean intensity values for the current file
    all_mean_intensity_values.extend(mean_intensity_values)
    # Filter rows with mean_pixel_intensity >= 0.25
    # Add a new column 'mean_pixel_intensity' to the DataFrame
    file_ma_df['mean_pixel_intensity'] = file_ma_df.apply(calculate_mean_intensity, axis=1)
    filtered_df = file_ma_df[file_ma_df['mean_pixel_intensity'] >= 0.25]
    # Save the filtered DataFrame to a new CSV file
    output_csv_path = '/dcs04/lieber/marmaypag/spatialDLPFC_SCZ_LIBD4100/raw-data/images/2_MockPNN/pixel_intensities_filtered_csvs/filtered_data_{numbers}.csv'
    filtered_df.to_csv(output_csv_path, index=False)
    # Plot histogram
    # plt.hist(mean_intensity_values, bins=50, color='blue', alpha=0.7)
    # plt.title(f'Histogram of Mean Pixel Intensities - Files with numbers {numbers}')
    # plt.xlabel('Mean Pixel Intensity')
    # plt.ylabel('Frequency')
    # # Save histogram as an image file
    # output_histogram_path = f'/dcs04/lieber/marmaypag/spatialDLPFC_SCZ_LIBD4100/raw-data/images/2_MockPNN/Test/histogram_{numbers}.png'
    # plt.savefig(output_histogram_path)
    # plt.clf()  # Clear the current figure to start a new one for the next iteration
    # print(f"Files with numbers {numbers}:")
    # print(f"CSV (CV): {files['csv_cv']}")
    # print(f"CSV (MA): {files['csv_ma']}")
    # print(f"Image (Raw): {files['img_raw']}")
    # print(f"Histogram saved to: {output_histogram_path}")
    # print("\n")
# Plot a single histogram for all mean intensity values
bins = np.linspace(0, 1, 51)
hist, bins = np.histogram(all_mean_intensity_values, bins=bins, range=(0, 1))
# Plot the histogram using plt.bar
plt.bar(bins[:-1], hist, width=np.diff(bins), color='blue', alpha=0.7, density=True)
plt.title('Histogram of All Mean Pixel Intensities')
plt.xlabel('Mean Pixel Intensity')
plt.ylabel('Frequency')
# Save histogram as an image file
output_histogram_path = '/dcs04/lieber/marmaypag/spatialDLPFC_SCZ_LIBD4100/raw-data/images/2_MockPNN/Test/histogram_all_normalized__.png'
plt.savefig(output_histogram_path)




### removing the rows with mean pixle intensities <=0.25 and creating a new filtered df
# function to find mean pixel intensity of pixels with abox
def calculate_mean_intensity(row):
    x1, y1, x4, y4 = int(row['x1']), int(row['y1']), int(row['x4']), int(row['y4'])
    roi = image[y1:y4, x1:x4]
    mean_intensity = cv2.mean(roi)[0]
    return mean_intensity



# Now, matching_files dictionary contains matching filenames grouped by the numbers within square brackets
# You can perform further operations on the matched files
for numbers, files in matching_files.items():
    # Read the CSV files
    file_cv_df = pd.read_csv(files['csv_cv'])
    file_ma_df = pd.read_csv(files['csv_ma'])
    # Read the image
    image = Image.open(files['img_raw'])
    image.seek(3)
    image = cv2.normalize(np.array(image, dtype='float32'), np.zeros(np.array(image, dtype='float32').shape, np.double), 1.0, 0.0, cv2.NORM_MINMAX)
    wfa_c = cv2.cvtColor(np.array(image * 255, dtype=np.uint8), cv2.COLOR_BGR2RGB)
    # Filter rows with mean_pixel_intensity >= 0.25
    # Add a new column 'mean_pixel_intensity' to the DataFrame
    file_ma_df['mean_pixel_intensity'] = file_ma_df.apply(calculate_mean_intensity, axis=1)
    filtered_df = file_ma_df[file_ma_df['mean_pixel_intensity'] >= 0.25]
    # Save the filtered DataFrame to a new CSV file
    output_csv_path = f'/dcs04/lieber/marmaypag/spatialDLPFC_SCZ_LIBD4100/raw-data/images/2_MockPNN/pixel_intensities_filtered_csvs/filtered_data_{numbers}.csv'
    filtered_df.to_csv(output_csv_path, index=False)
    