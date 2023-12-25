import os
import re
import cv2
import numpy as np
import pandas as pd
import PIL
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from PIL import Image, ImageFont, ImageDraw, ImageEnhance, ImageSequence
import glob



# read the full tissue image for NTC and sCZ --> Raw
ntc_path = '/dcs04/lieber/marmaypag/spatialDLPFC_SCZ_LIBD4100/raw-data/images/2_MockPNN/20220712_VIF_MockPNN_Strong_NTC_C1_Br5182_MLtraining/V12F14-053.tif'
scz_path = '/dcs04/lieber/marmaypag/spatialDLPFC_SCZ_LIBD4100/raw-data/images/2_MockPNN/20220712_VIF_MockPNN_Strong_SCZ_C1_Br2039_MLtraining/V12F14-057.tif'

# read the segmented full tissue images of ntc and scz
ntc_seg = '/dcs04/lieber/marmaypag/spatialDLPFC_SCZ_LIBD4100/raw-data/images/2_MockPNN/Test_images/V12F14-053_wfa_ntc.tif'
scz_seg = '/dcs04/lieber/marmaypag/spatialDLPFC_SCZ_LIBD4100/raw-data/images/2_MockPNN/Test_images/V12F14-057_wfa_scz.tif'

# Each of the individual tiles for NTC 
ntc_tile_shape = [1860, 1396]
scz_tile_shape = [1860, 1396]
ntc_tissue_shape = [16740, 13960]
scz_tissue_shape = [14880, 16752]


#### plot the gridlines --> using the same code for ntc and scz images
# Load the image
Image.MAX_IMAGE_PIXELS = None
image = Image.open(scz_seg)
# image.seek(4)
wfa_scz = np.array(image, dtype = 'uint8')
wfa_c_scz = cv2.cvtColor(wfa_scz,cv2.COLOR_BGR2RGB)
wfa_ntc = np.array(image, dtype = 'uint8')
wfa_c_ntc = cv2.cvtColor(wfa_ntc,cv2.COLOR_BGR2RGB)

cv2.imwrite('/dcs04/lieber/marmaypag/spatialDLPFC_SCZ_LIBD4100/raw-data/images/2_MockPNN/Test_images/wfa_raw_scz.tif', wfa_scz)


# Get image dimensions
image_width, image_height = image.size
# Box dimensions
box_width = 1860
box_height = 1396
# Calculate the number of boxes in both dimensions
num_boxes_horizontal = image_width // box_width
num_boxes_vertical = image_height // box_height
# Create a list to store the starting coordinates of each box
box_coordinates = []
# Create a list to store the position of each box
box_positions = []
# Create a figure and axis
fig, ax = plt.subplots(1)
# Plot the image
ax.imshow(image)


## draw both gridlines and also save the box coordinates
# Draw grid lines and store box coordinates
for i in range(num_boxes_vertical):
    for j in range(num_boxes_horizontal):
        x = j * box_width
        y = i * box_height
        # ax.add_patch(patches.Rectangle((x, y), box_width, box_height, linewidth=1, edgecolor='white', facecolor='none'))
        box_coordinates.append((x, y)) # list of all the top left coordiantes of all the grid boxes
        box_positions.append((i, j)) # stores the positions of the boxes like (0,0), (0,1) etc


# Set axis limits
ax.set_xlim([0, image_width])
ax.set_ylim([image_height, 0])

# Remove axis ticks
ax.set_xticks([])
ax.set_yticks([])

# Save the plotted image with the grid
plt.savefig('/dcs04/lieber/marmaypag/spatialDLPFC_SCZ_LIBD4100/raw-data/images/2_MockPNN/Test_images/grid_image_scz.png', bbox_inches='tight', pad_inches=0.1)


## choosing specific boxes of the grid to compare to the annotations
# Specify the target box position
# target_positions_tiles_to_tissue = [(0,3), (0,7), (1,7), (2,4), (2,8), (3,2), (4,2), (4,3), (4,7), (7,2), (7,7), (6,1)] # wrong x,y ==> but this is what is to be used for the x,y coordinates for scalling up the pixels from tiles to tissue
target_positions = [(3,0), (7,0), (7,1), (4,2), (8,2), (2,3), (2,4), (3,4), (7,4), (2,7), (7,7), (9,7)] # correct ones

# Draw red boxes around target positions
for pos in target_positions:
    i, j = pos
    x = j * box_width
    y = i * box_height
    rect = patches.Rectangle((x, y), box_width, box_height, linewidth=1, edgecolor='red', facecolor='none')
    ax.add_patch(rect)
    # Add text with box position
    ax.text(x + 0.5 * box_width, y + 0.5 * box_height, f'({i},{j})', color='red',
            fontsize=8, ha='center', va='center')

# Set axis limits
ax.set_xlim([0, image_width])
ax.set_ylim([image_height, 0])

# Remove axis ticks
ax.set_xticks([])
ax.set_yticks([])

# Save the plotted image with the grid and red boxes
plt.savefig('/dcs04/lieber/marmaypag/spatialDLPFC_SCZ_LIBD4100/raw-data/images/2_MockPNN/Test_images/grid_with_red_boxes_corrected.png', bbox_inches='tight', pad_inches=0.1)


# Ensure the target position is within the grid
# if target_position[0] <= num_boxes_vertical and target_position[1] <= num_boxes_horizontal:
#     # Calculate box coordinates for the target position
#     target_box_x = target_position[1] * box_width
#     target_box_y = target_position[0] * box_height
#     # Crop the image based on the target box coordinates
#     cropped_image = image.crop((target_box_x, target_box_y, target_box_x + box_width, target_box_y + box_height))
#     cropped_image_arr = (np.array(cropped_image, dtype = 'uint8').shape)
#     # Save the cropped image
#     save_path = '/dcs04/lieber/marmaypag/spatialDLPFC_SCZ_LIBD4100/raw-data/images/2_MockPNN/Test_images/cropped_32_ntc_.png'
#     cropped_image.save(save_path)
#     print(f"Cropped image saved at {save_path}")


# Loop through target positions, crop, and save
for target_position in target_positions:
    i, j = target_position
    if i < num_boxes_vertical and j < num_boxes_horizontal:
        # Calculate box coordinates for the target position
        target_box_x = j * box_width
        target_box_y = i * box_height
        # Crop the image based on the target box coordinates
        cropped_image = image.crop((target_box_x, target_box_y, target_box_x + box_width, target_box_y + box_height))
        cropped_image_arr = np.array(cropped_image, dtype='uint8').shape
        # Save the cropped image
        save_path = f'/dcs04/lieber/marmaypag/spatialDLPFC_SCZ_LIBD4100/raw-data/images/2_MockPNN/Test_images/cropped_image_{i}_{j}.png'
        cropped_image.save(save_path)
        print(f"Cropped image saved at {save_path}")


## convert the pixels from tiles to tissue section
## formula n4 = ((n5 + (x*b1)), (n6 + (y*b2))
## n4 = new coords; n5, n6 = x,y coordinates that needs to be converted; x,y = box position in the fulll tissue section (eg., (3,2), etc);
## b1,b2 = width, height of each tile
target_positions_tiles_to_tissue = [(0,3), (0,7), (1,7), (2,4), (2,8), (3,2), (4,2), (4,3), (4,7), (7,2), (7,7), (6,1)] # wrong x,y ==> but this is what is to be used for the x,y coordinates for scalling up the pixels from tiles to tissue
# test_csv_path_11 = '/dcs04/lieber/marmaypag/spatialDLPFC_SCZ_LIBD4100/processed-data/Experimentation_archive/2_MockPNN/Training_tiles/Manual_annotations/Annotations_in_pixels/20220712_VIF_MockPNN_Strong_Scan1_[6384,53057]_component_data_11_wfa__seg_manual_annotations.csv'
# csv_11 = pd.read_csv(test_csv_path_11)
# csv_11.head(5) # box position for csv_11 = (7,0)
# csv_11 old values = xc to x4 = 1692  163  1675  147  1709  147  1675  179  1709  179
######### this formula works if you flip the x and y axis
csv_dir_manual_annotations = '/dcs04/lieber/marmaypag/spatialDLPFC_SCZ_LIBD4100/processed-data/Experimentation_archive/2_MockPNN/Training_tiles/Manual_annotations/Annotations_in_pixels/'
# Box dimensions
box_width = 1860
box_height = 1396

for csv_path in os.listdir(csv_dir_manual_annotations):
    if "NTC" in csv_path:
        print(csv_path)
        csv_file = pd.read_csv(os.path.join(csv_dir_manual_annotations,csv_path))
        for csv_row in range(len(csv_file)):
            # print(csv_row)
            for list_len in range(len(target_positions_tiles_to_tissue)):
                csv_file.loc[csv_row, 'x1'] += target_positions_tiles_to_tissue[list_len][0] * box_width
                csv_file.loc[csv_row, 'y1'] += target_positions_tiles_to_tissue[list_len][1] * box_height
                csv_file.loc[csv_row, 'x2'] += target_positions_tiles_to_tissue[list_len][0] * box_width
                csv_file.loc[csv_row, 'y2'] += target_positions_tiles_to_tissue[list_len][1] * box_height
                csv_file.loc[csv_row, 'x3'] += target_positions_tiles_to_tissue[list_len][0] * box_width
                csv_file.loc[csv_row, 'y3'] += target_positions_tiles_to_tissue[list_len][1] * box_height
                csv_file.loc[csv_row, 'x4'] += target_positions_tiles_to_tissue[list_len][0] * box_width
                csv_file.loc[csv_row, 'y4'] += target_positions_tiles_to_tissue[list_len][1] * box_height
                wfa_cnt = cv2.rectangle(wfa_c_ntc, (csv_file['x1'].iloc[0], csv_file['y1'].iloc[0]), (csv_file['x1'].iloc[0] + (int(csv_file['Width'].iloc[0])), csv_file['y1'].iloc[0] + (int(csv_file['Height'].iloc[0]))), (0,0,255), 2) # red rectangle



                # erro with this code
                # csv_file.loc[csv_row, 'x1'] = csv_file.loc[csv_row, 'x1'] + (target_positions_tiles_to_tissue[list_len][0]*box_width)  
                # csv_file.loc[csv_row, 'y1'] = csv_file.loc[csv_row, 'y1'] + (target_positions_tiles_to_tissue[list_len][1]*box_height) 
                # csv_file.loc[csv_row, 'x2'] = csv_file.loc[csv_row, 'x2'] + (target_positions_tiles_to_tissue[list_len][0]*box_width)
                # csv_file.loc['y2'][csv_row] = csv_file.loc['y2'][csv_row] + (target_positions_tiles_to_tissue[list_len][1]*box_height)
                # csv_file.loc['x3'][csv_row] = csv_file.loc['x3'][csv_row] + (target_positions_tiles_to_tissue[list_len][0]*box_width)
                # csv_file.loc['y3'][csv_row] = csv_file.loc['y3'][csv_row] + (target_positions_tiles_to_tissue[list_len][1]*box_height)
                # csv_file.loc['x4'][csv_row] = csv_file.loc['x4'][csv_row] + (target_positions_tiles_to_tissue[list_len][0]*box_width)
                # csv_file.loc['y4'][csv_row] = csv_file.loc['y4'][csv_row] + (target_positions_tiles_to_tissue[list_len][1]*box_height)
                

cv2.imwrite('/dcs04/lieber/marmaypag/spatialDLPFC_SCZ_LIBD4100/raw-data/images/2_MockPNN/tile_tissue_bb_match/all_files_red_rect.tif', wfa_c_ntc)



### doing this manually
target_positions_tiles_to_tissue_ntc = [(0,3), (0,7), (1,7), (2,4), (2,8), (3,2), (4,2), (4,3), (4,7), (7,2), (7,7), (6,1)] # wrong x,y ==> but this is what is to be used for the x,y coordinates for scalling up the pixels from tiles to tissue
target_positions_tiles_to_tissue_scz = [(0,2), (0,5), (1,0), (1,2), (1,7), (1,10), (3,1), (3,11), (4,4), (5,4), (5,7), (6,3)]

test_csv_path_11 = '/dcs04/lieber/marmaypag/spatialDLPFC_SCZ_LIBD4100/processed-data/Experimentation_archive/2_MockPNN/Training_tiles/Manual_annotations/Annotations_in_pixels/20220712_VIF_NTC_MockPNN_Strong_Scan1_[12864,50280]_component_data_17_wfa__seg_manual_annotations.csv'
csv_11 = pd.read_csv(test_csv_path_11)
csv_11.head(5) # box position for csv_11 = (7,0)
print(len(csv_11))

for csv_row in range(len(csv_11)):
    print("BEFORE: row number, x1, y1", csv_row, csv_11.loc[csv_row, 'x1' ], csv_11.loc[csv_row, 'y1'])
    csv_11.loc[csv_row, 'x1' ] = csv_11.loc[csv_row, 'x1'] + (4*1860)  
    csv_11.loc[csv_row, 'y1'] = csv_11.loc[csv_row, 'y1'] + (3*1396) 
    x = csv_11['x1'].iloc[csv_row]
    y = csv_11['y1'].iloc[csv_row]
    width = int(csv_11['Width'].iloc[csv_row])
    height = int(csv_11['Height'].iloc[csv_row])
    print("AFTER", x,y,width,height)
    wfa_cnt = cv2.rectangle(wfa_c_ntc, (x, y), (x + width, y + height), (0,0,255), 2) # red rectangle
    print("box drawn")

    
cv2.imwrite('/dcs04/lieber/marmaypag/spatialDLPFC_SCZ_LIBD4100/raw-data/images/2_MockPNN/tile_tissue_bb_match/overlay_tiles_on_tissue_scz.tif', wfa_c_scz)



    csv_11['x2'][csv_row] = csv_11['x2'][csv_row] + (0*1860)
    csv_11['y2'][csv_row] = csv_11['y2'][csv_row] + (3*1396)
    csv_11['x3'][csv_row] = csv_11['x3'][csv_row] + (0*1860)
    csv_11['y3'][csv_row] = csv_11['y3'][csv_row] + (3*1396)
    csv_11['x4'][csv_row] = csv_11['x4'][csv_row] + (0*1860)
    csv_11['y4'][csv_row] = csv_11['y4'][csv_row] + (3*1396)

# csv_11['x1'][csv_row] = csv_11['x1'][csv_row] + (0*1860)  
#     csv_11['y1'][csv_row] = csv_11['y1'][csv_row] + (3*1396) 
#     csv_11['x2'][csv_row] = csv_11['x2'][csv_row] + (0*1860)
#     csv_11['y2'][csv_row] = csv_11['y2'][csv_row] + (3*1396)
#     csv_11['x3'][csv_row] = csv_11['x3'][csv_row] + (0*1860)
#     csv_11['y3'][csv_row] = csv_11['y3'][csv_row] + (3*1396)
#     csv_11['x4'][csv_row] = csv_11['x4'][csv_row] + (0*1860)
#     csv_11['y4'][csv_row] = csv_11['y4'][csv_row] + (3*1396)



# Extract relevant information from the DataFrame
x = csv_11['x1'].iloc[0]
y = csv_11['y1'].iloc[0]
width = int(csv_11['Width'].iloc[0])
height = int(csv_11['Height'].iloc[0])
print(x,y,width,height)

wfa_cnt = cv2.rectangle(wfa_c, (x, y), (x + width, y + height), (0,0,255), 2) # red rectangle

# Add text above the rectangle
# text = "RED BOX HERE RIGHT PLACE"
# font = cv2.FONT_HERSHEY_SIMPLEX
# font_scale = 5
# font_thickness = 2
# text_size = cv2.getTextSize(text, font, font_scale, font_thickness)[0]
# text_x = x + (width - text_size[0]) // 2
# text_y = y - 5  # Adjust this value for the desired vertical distance from the rectangle
# cv2.putText(wfa_c, text, (text_x, text_y), font, font_scale, (255, 255, 255), font_thickness, cv2.LINE_AA)



cv2.imwrite('/dcs04/lieber/marmaypag/spatialDLPFC_SCZ_LIBD4100/raw-data/images/2_MockPNN/tile_tissue_bb_match/all_files_red_rect.tif', wfa_c_ntc)


