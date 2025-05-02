OMJsGJsRJH
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
ntc_raw = '/dcs04/lieber/marmaypag/spatialDLPFC_SCZ_LIBD4100/raw-data/images/2_MockPNN/20220712_VIF_MockPNN_Strong_NTC_C1_Br5182_MLtraining/V12F14-053.tif'
scz_raw = '/dcs04/lieber/marmaypag/spatialDLPFC_SCZ_LIBD4100/raw-data/images/2_MockPNN/20220712_VIF_MockPNN_Strong_SCZ_C1_Br2039_MLtraining/V12F14-057.tif'

# read the segmented full tissue images of ntc and scz
ntc_seg = '/dcs04/lieber/marmaypag/spatialDLPFC_SCZ_LIBD4100/raw-data/images/2_MockPNN/Test_images/V12F14-053_wfa_ntc.tif'
scz_seg = '/dcs04/lieber/marmaypag/spatialDLPFC_SCZ_LIBD4100/raw-data/images/2_MockPNN/Test_images/V12F14-057_scz_ntc.tif'

#### plot the gridlines --> using the same code for ntc and scz images
# Load the image
Image.MAX_IMAGE_PIXELS = None
image = Image.open(scz_seg)
wfa_scz = np.array(image, dtype = 'uint8')
wfa_c_scz = cv2.cvtColor(wfa_scz,cv2.COLOR_BGR2RGB)
# cv2.imwrite('/dcs04/lieber/marmaypag/spatialDLPFC_SCZ_LIBD4100/raw-data/images/2_MockPNN/Test_images/wfa_raw_scz.tif', wfa_scz)

# draw a grid on the segmented scz image
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
## highlight target tiles on the grid image
target_tiles_list_ntc = [(3,0), (7,0), (7,1), (4,2), (8,2), (2,3), (2,4), (3,4), (7,4), (2,7), (7,7), (9,7)]
target_tiles_list_scz = [(2,0), (5,0), (0,1), (2,1), (7,1), (10, 1), (1,3), (11,3), (4,4), (4,5), (7,5), (3,6)]

## draw both gridlines and also save the box coordinates
## draw the target boxes and their box positions
for i in range(num_boxes_vertical + 1):
    y = i * box_height
    cv2.line(wfa_c_scz, (0, y), (wfa_scz.shape[1], y), ((0, 255, 255)), 2)
    for j in range(num_boxes_horizontal + 1):
        x = j * box_width
        cv2.line(wfa_c_scz, (x, 0), (x, wfa_scz.shape[0]), ((0, 255, 255)), 2)
        # Append box coordinates and positions
        box_coordinates.append((x, y))  # Top-left coordinates
        box_positions.append((i, j))     # Box positions
        # Highlight target boxes in red
        if (i, j) in target_tiles_list_scz:
            x1 = x
            y1 = y
            x2 = x + box_width
            y2 = y + box_height
            cv2.rectangle(wfa_c_scz, (x1, y1), (x2, y2), (255, 0, 255), 10)
            # Put text on the image with box positions
            position_text = f"({i},{j})"
            text_position = (x1 + 5, y1) ##((x1 + x2) // 2, (y1 + y2) // 2)
            cv2.putText(wfa_c_scz, position_text, text_position, cv2.FONT_HERSHEY_SIMPLEX, 5, (255, 0, 255), 6)


# Save the NumPy array as a TIFF image 
cv2.imwrite('/dcs04/lieber/marmaypag/spatialDLPFC_SCZ_LIBD4100/raw-data/images/2_MockPNN/tile_tissue_bb_match/grid_image_with_target_tiles_scz.tif', wfa_c_scz)

### drawing the bouding boxes from tiles onto the tissue section
target_positions_tiles_to_tissue_ntc = [(0,3), (0,7), (1,7), (2,4), (2,8), (3,2), (4,2), (4,3), (4,7), (7,2), (7,7), (6,1)] # wrong x,y ==> but this is what is to be used for the x,y coordinates for scalling up the pixels from tiles to tissue
target_positions_tiles_to_tissue_scz = [(2,0), (0,5), (1,2), (1,7), (1,10), (3,1), (3,11), (4,4), (5,4), (5,7), (6,3)] ## fully wrong, need to chnage x, y to y,x

## match the csv to the tile position 
test_csv_path = '/dcs04/lieber/marmaypag/spatialDLPFC_SCZ_LIBD4100/processed-data/Experimentation_archive/2_MockPNN/Training_tiles/Manual_annotations/Annotations_in_pixels/20220712_VIF_SCZ_MockPNN_Strong_Scan1_[12480,49800]_component_data_22_wfa__seg_manual_annotations.csv'
csv_ = pd.read_csv(test_csv_path)
csv_.head(5) # box position for csv_11 = (7,0)
print(len(csv_))

# read the raw image to find mean pixel intensity from there
Image.MAX_IMAGE_PIXELS = None
wfa_img = Image.open(scz_raw)
wfa_img.seek(4)
wfa = np.array(wfa_img, dtype = 'uint8')
wfa_c = cv2.cvtColor(wfa,cv2.COLOR_BGR2RGB)

## read the saved grid image with target tiles highlighted
# img_path = '/dcs04/lieber/marmaypag/spatialDLPFC_SCZ_LIBD4100/raw-data/images/2_MockPNN/tile_tissue_bb_match/grid_image_with_target_tiles_ntc.tif'
# img = Image.open(img_path)
# img_arr = np.array(img, dtype = 'uint8')
# img_c = cv2.cvtColor(img_arr,cv2.COLOR_BGR2RGB)

## loop through the rows to overlay the tiles BB and also the mean pixel intensities
for index, csv_row in csv_.iterrows():
    print("BEFORE: row number, x1, y1", index, csv_row['x1'], csv_row['y1'])
    csv_.loc[index, 'x1'] = csv_row['x1'] + (6 * 1860)
    csv_.loc[index, 'y1'] = csv_row['y1'] + (3 * 1396)
    csv_.loc[index, 'x4'] = csv_row['x4'] + (6 * 1860)
    csv_.loc[index, 'y4'] = csv_row['y4'] + (3 * 1396)
    # x = csv_row['x1']
    # y = csv_row['y1']
    # width = int(csv_row['Width'])
    # height = int(csv_row['Height'])
    print("AFTER", csv_.loc[index, 'x1'], csv_.loc[index, 'y1'], csv_.loc[index, 'x4'], csv_.loc[index, 'y4'])
    # Draw rectangles on the image
    x1, y1, x4, y4 = csv_.loc[index, 'x1'], csv_.loc[index, 'y1'], csv_.loc[index, 'x4'], csv_.loc[index, 'y4']
    cRec = cv2.rectangle(wfa_c_scz, (int(csv_.loc[index, 'x1']), int(csv_.loc[index, 'y1'])), 
                  (int(csv_.loc[index, 'x4']), int(csv_.loc[index, 'y4'])), (0, 0, 255), 2)
    # Calculate mean pixel intensity within the bounding box
    roi = wfa[y1:y4, x1:x4]  # Changed from wfa to wfa_c_scz
    mean_intensity = int(cv2.mean(roi)[0])
    # Put text above the bounding box with mean pixel intensity
    text = f"{mean_intensity}"
    text_position = (x1, y1 - 5)  # Adjust for text placement
    cv2.putText(wfa_c_scz, text, text_position, cv2.FONT_HERSHEY_SIMPLEX, 1, (0, 0, 255), 2)


cv2.imwrite('/dcs04/lieber/marmaypag/spatialDLPFC_SCZ_LIBD4100/raw-data/images/2_MockPNN/tile_tissue_bb_match/overlay_tiles_on_tissue_scz_with_grid.tif', wfa_c_scz)


