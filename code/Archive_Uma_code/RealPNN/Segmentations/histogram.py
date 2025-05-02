import os
import re
import cv2
import numpy as np
import pandas as pd
import PIL
import matplotlib
import matplotlib.pyplot as plt
from PIL import Image, ImageFont, ImageDraw, ImageEnhance, ImageSequence
import glob

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
# def calculate_mean_intensity(row):
#     x1, y1, x4, y4 = int(row['x1']), int(row['y1']), int(row['x4']), int(row['y4'])
#     roi = image[y1:y4, x1:x4]
#     mean_intensity = cv2.mean(roi)[0]
#     return mean_intensity

## find mean, sum and max pixel intensities
def calculate_intensity_metrics(row, image):
    x1, y1, x4, y4 = int(row['x1']), int(row['y1']), int(row['x4']), int(row['y4'])
    roi = image[y1:y4, x1:x4]    
    mean_intensity = cv2.mean(roi)[0]
    max_intensity = roi.max()
    sum_intensity = roi.sum()
    return  mean_intensity, max_intensity, sum_intensity #


# List to store all mean intensity values from all files
all_mean_intensity_values, all_sum_values, all_max_values = [],[], []
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
    mean_intensity_values, sum_values, max_values = [], [], []
    # Iterate over bounding boxes and compute mean pixel intensities
    for index, row in file_ma_df.iterrows():
        x1, y1, x4, y4 = int(row['x1']), int(row['y1']), int(row['x4']), int(row['y4'])
        # Extract region of interest (ROI)
        roi = wfa_c[y1:y4, x1:x4]
        # Compute mean pixel intensities
        mean_intensity = cv2.mean(roi)[0]
        max_intensity = roi.max()
        sum_intensity = roi.sum()
        mean_intensity_values.append(mean_intensity)
        max_values.append(max_intensity)
        sum_values.append(sum_intensity)
    # Extend the list with mean intensity values for the current file
    all_mean_intensity_values.extend(mean_intensity_values)
    all_sum_values.extend(sum_values)
    all_max_values.extend(max_values)
    # Filter rows with mean_pixel_intensity >= 0.25
    # Add a new column 'mean_pixel_intensity' to the DataFrame
    # file_ma_df['mean_pixel_intensity'], file_ma_df['max_pixel_intensity'], file_ma_df['sum_pixel_intensity'] = zip(*file_ma_df.apply(lambda row: calculate_intensity_metrics(row, image), axis=1)) #  
    # output_csv_path_ = f'/dcs04/lieber/marmaypag/spatialDLPFC_SCZ_LIBD4100/raw-data/images/2_MockPNN/pixel_intensities_filtered_csvs/Max_Mean_sum_{numbers}.csv'
    # file_ma_df.to_csv(output_csv_path_, index=False)

    filtered_df = file_ma_df[file_ma_df['mean_pixel_intensity'] >= 0.25]
    # Save the filtered DataFrame to a new CSV file
    output_csv_path = f'/dcs04/lieber/marmaypag/spatialDLPFC_SCZ_LIBD4100/raw-data/images/2_MockPNN/pixel_intensities_filtered_csvs/filtered_data_{numbers}.csv'
    filtered_df.to_csv(output_csv_path, index=False)
    # Plot histogram for each tile
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
output_histogram_path = '/dcs04/lieber/marmaypag/spatialDLPFC_SCZ_LIBD4100/raw-data/images/2_MockPNN/Test/histogram_all_mean_values.png'
plt.savefig(output_histogram_path)


## plotting mean, max and sum pixels
# Set up subplots
fig, axs = plt.subplots(3, 1, figsize=(20, 20))

# Plot histogram for all_mean_intensity_values
axs[0].hist(all_mean_intensity_values, bins=50, alpha=0.7, color='blue')
axs[0].set_title('Histogram of Mean Intensity Values')
axs[0].set_xlabel('Mean Intensity')
axs[0].set_ylabel('Frequency')

# Plot histogram for all_sum_values
axs[1].hist(all_sum_values, bins=50, alpha=0.7, color='green')
axs[1].set_title('Histogram of Sum Values')
axs[1].set_xlabel('Sum Value')
axs[1].set_ylabel('Frequency')

# Plot histogram for all_max_values
axs[2].hist(all_max_values, bins=50, alpha=0.7, color='red')
axs[2].set_title('Histogram of Max Values')
axs[2].set_xlabel('Max Value')
axs[2].set_ylabel('Frequency')

# Adjust layout to prevent overlap of titles and labels
plt.tight_layout()

# Save the combined histogram as an image file
output_histogram_path = '/dcs04/lieber/marmaypag/spatialDLPFC_SCZ_LIBD4100/raw-data/images/2_MockPNN/Test/histogram_all_mean_sum_max_values.png'
plt.savefig(output_histogram_path)


## individual plots# Plot histogram for all_mean_intensity_values
plt.figure(figsize=(20, 20))
plt.hist(all_mean_intensity_values, bins=50, alpha=0.7, color='blue')
plt.title('Histogram of Mean Intensity Values')
plt.xlabel('Mean Intensity', fontsize=14)
plt.ylabel('Frequency', fontsize=14)
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)
plt.savefig('/dcs04/lieber/marmaypag/spatialDLPFC_SCZ_LIBD4100/raw-data/images/2_MockPNN/Test/histogram_all_mean_values.png')


# Plot histogram for all_sum_values
plt.figure(figsize=(20, 20))
plt.hist(all_sum_values, bins=50, alpha=0.7, color='green')
plt.title('Histogram of Sum Values')
plt.xlabel('Sum Value', fontsize=14)
plt.ylabel('Frequency', fontsize=14)
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)
plt.savefig('/dcs04/lieber/marmaypag/spatialDLPFC_SCZ_LIBD4100/raw-data/images/2_MockPNN/Test/histogram_all_sum_values.png')


# Plot histogram for all_max_values
plt.figure(figsize=(20, 20))
plt.hist(all_max_values, bins=50, alpha=0.7, color='red')
plt.title('Histogram of Max Values')
plt.xlabel('Max Value', fontsize=14)
plt.ylabel('Frequency', fontsize=14)
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)
plt.savefig('/dcs04/lieber/marmaypag/spatialDLPFC_SCZ_LIBD4100/raw-data/images/2_MockPNN/Test/histogram_all_max_values.png')




#### plotting the mean, sum and max pixel intensities for all tiles
# combining the dfs
# Combine all dataframes into one
# Specify the directory containing CSV files
csv_directory = '/dcs04/lieber/marmaypag/spatialDLPFC_SCZ_LIBD4100/raw-data/images/2_MockPNN/pixel_intensities_filtered_csvs/'

# Use glob to find all CSV files in the directory
csv_files = glob.glob(csv_directory + '*.csv')

# Initialize an empty list to store dataframes
all_dfs = []

# Iterate through CSV files and read them into dataframes
for csv_file in csv_files:
    print(csv_file)
    df = pd.read_csv(csv_file)
    print(len(df))
    all_dfs.append(df)

# Combine all dataframes into one
# combined_df = pd.concat(all_dfs, ignore_index=True)

# Exclude empty or all-NA columns from each dataframe <-- do this to avoid warnings
all_dfs = [df.dropna(axis=1, how='all') for df in all_dfs]
# Drop the "Unnamed: 0" column
combined_df = combined_df.drop("Unnamed: 0", axis=1, errors="ignore")


# plot all histograms
def plot_and_save_histogram(data, column, color, output_path):
    plt.hist(data[column], bins=50, alpha=0.7, color=color) # , label=column
    plt.title(f'Histogram of {column}')
    plt.xlabel('Intensity')
    plt.ylabel('Frequency')
    plt.legend()
    plt.savefig(output_path)

# Plot histogram for mean_pixel_intensity
plt.hist(combined_df['mean_pixel_intensity'], bins=50, alpha=0.7, color='blue') #, label='Mean Pixel Intensity'
plt.title('Histogram of Mean Pixel Intensity')
plt.xlabel('Intensity')
plt.ylabel('Frequency')
# plt.legend()
output_histogram_path_mean = '/dcs04/lieber/marmaypag/spatialDLPFC_SCZ_LIBD4100/raw-data/images/2_MockPNN/Test/histogram_mean_.png'
plt.savefig(output_histogram_path_mean)
# plt.show()

# Plot histogram for max_pixel_intensity
plt.hist(combined_df['max_pixel_intensity'], bins=50, alpha=0.7, color='green') # , label='Max Pixel Intensity'
plt.title('Histogram of Max Pixel Intensity')
plt.xlabel('Intensity')
plt.ylabel('Frequency')
# plt.legend()
output_histogram_path_max = '/dcs04/lieber/marmaypag/spatialDLPFC_SCZ_LIBD4100/raw-data/images/2_MockPNN/Test/histogram_max_.png'
plt.savefig(output_histogram_path_max)
# plt.show()

# Plot histogram for sum_pixel_intensity
plt.hist(combined_df['sum_pixel_intensity'], bins=50, alpha=0.7, color='red') # , label='Sum Pixel Intensity'
plt.title('Histogram of Sum Pixel Intensity')
plt.xlabel('Intensity')
plt.ylabel('Frequency')
# plt.legend()
output_histogram_path_sum = '/dcs04/lieber/marmaypag/spatialDLPFC_SCZ_LIBD4100/raw-data/images/2_MockPNN/Test/histogram_sum_.png'
plt.savefig(output_histogram_path_sum)
# plt.show()



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
    