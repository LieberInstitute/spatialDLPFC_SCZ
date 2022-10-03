import numpy as np
import pandas as pd
from PIL import Image, ImageFont, ImageDraw, ImageEnhance, ImageSequence
import os
import matplotlib
import argparse
import matplotlib.pyplot as plt
import mpldatacursor
import sys
import cv2
import math
import scipy

# segmentation function
def segment(normalised_img, color_img):
    hierachy, img_threshold = cv2.threshold((np.array((normalised_img * 255), dtype = np.uint8)), 100, 255, cv2.THRESH_BINARY)
    contours,_ = cv2.findContours(img_threshold, cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_SIMPLE)
    x, y, w, h, area = [],[],[],[],[]
    for cnt in contours:
        ax,ay,aw,ah = cv2.boundingRect(cnt)
        if(w*h >= 100):
            area1 = cv2.contourArea(cnt)
            print(x,y,w,h)
            x.append(ax)
            y.append(ay)
            w.append(aw)
            h.append(ah)
            area.append(area1)
            print(x,y,w,h)
            bb_img = cv2.rectangle(color_img, (x,y), (x+w+5, y+h+5), (255,0,0), 2)
            rect = cv2.minAreaRect(cnt)
            box = cv2.boxPoints(rect)
            box = np.int0(box)
            cnt_img = cv2.drawContours(out_img1_neun,[box],0,(0,0,255),1)

    col_names = ['img_file_name','type_of_object_str', 'x1', 'y1', 'Width', 'Height', 'total_number']
    object_name = obj_name # name of the objects stored in the dataframe
    file_name = os.path.basename(img_test) # image file name
    dict = {col_names[0]: file_name, col_names[1]: object_name, col_names[2]: x, col_names[3]: y, col_names[4]: w, col_names[5]: h, col_names[6]: len(x)}
    img_info_df = pd.DataFrame(dict, columns = col_names)
    img_info_df['x2'] = img_info_df['x1'] + img_info_df['Width']
    img_info_df['y2'], img_info_df['x3'] = img_info_df['y1'], img_info_df['x1']
    img_info_df['y3'] = img_info_df['y1'] + img_info_df['Height']
    img_info_df['x4'], img_info_df['y4'] = img_info_df['x2'], img_info_df['y3']
    img_info_df = img_info_df[['img_file_name', 'type_of_object_str', 'x1', 'y1', 'x2', 'y2', 'x3', 'y3', 'x4', 'y4', 'Width', 'Height', 'total_number']]
    return(cnt_img, img_info_df)



def plot_seg_img(original_img, segmented_img):
    fig,ax = plt.subplots(nrows = 1, ncols = 2, figsize = (20,20))
    ax[0].imshow(original_img)
    ax[0].title.set_text('Original')
    ax[1].imshow(segmented_img)
    ax[1].title.set_text('Segemented')
    fig.show()

