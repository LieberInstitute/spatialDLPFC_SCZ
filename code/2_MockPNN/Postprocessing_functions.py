'''
For Visium-IF
Channel0 = DAPI, DAPI
Channel1 = Claudin5 (Alex 488),
Channel2 = NeuN (Alexa 555),
Channel3 = WFA (Alexa 647),
Channel4 = AF (Autofluorescence), sample AF
Channel5 = Thumbnail
'''


# populate a dataframe with the coordinates info
def create_df(x,y,w,h, area, img_test, label):
    print("Segmented {0} {1}".format(len(x), label))
    col_names = ['img_file_name','type_of_object_str', 'x1', 'y1', 'Width', 'Height', 'area']
    file_name = os.path.basename(img_test) # image file name
    dict = {col_names[0]: file_name, col_names[1]: label, col_names[2]: x, col_names[3]: y, col_names[4]: w, col_names[5]: h, col_names[6]: area}
    img_info_df = pd.DataFrame(dict, columns = col_names)
    img_info_df['x2'] = img_info_df['x1'] + img_info_df['Width']
    img_info_df['y2'], img_info_df['x3'] = img_info_df['y1'], img_info_df['x1']
    img_info_df['y3'] = img_info_df['y1'] + img_info_df['Height']
    img_info_df['x4'], img_info_df['y4'] = img_info_df['x2'], img_info_df['y3']
    img_info_df['xc'] = np.int0(np.ceil((img_info_df['x1'] + img_info_df['x4'])/2))
    img_info_df['yc'] = np.int0(np.ceil((img_info_df['y1'] + img_info_df['y4'])/2))
    img_info_df = img_info_df[['img_file_name', 'type_of_object_str', 'x1', 'y1', 'x2', 'y2', 'x3', 'y3', 'x4', 'y4', 'xc', 'yc', 'Width', 'Height', 'area']]
    return img_info_df

# draw a white rectangle filled using the coordinates from the csv
from collections import Counter
pix_list, locs_list, mean_pix_int_list  = [], [], []
def all_pix_pnns(img_info_df, contour_img, original_img):
    for box in range(len(img_info_df['x1'])):
        print(box)
        gray_image = np.full((1396,1860), 0, dtype=np.uint8) # new blank image so the original pix intensities are retained
        rect = cv2.rectangle(gray_image, (img_info_df['x1'][box], img_info_df['y1'][box]), (img_info_df['x4'][box], img_info_df['y4'][box]), (255,255,255), -1) # draw white filled rect on the copy of the image
        cv2.putText(contour_img, ('%d'%box), (img_info_df['x1'][box],img_info_df['y1'][box]), cv2.FONT_HERSHEY_SIMPLEX, 2, (125, 246, 55), 3)
        # gray_seg_wfa = skimage.color.rgb2gray(contour_img)
        # plot_img(gray_image, contour_img)
        locs = np.argwhere(gray_image == 255)
        print(locs.shape)
        for i in range(locs.shape[0]):
            for j in range(locs.shape[1] -1):
                if gray_image[locs[i,j],locs[i,j+1]] == 255: # gray image has white filled boxes
                    # print(original_img[locs[i,j],locs[i,j+1]])
                    pix_list.append(original_img[locs[i,j],locs[i,j+1]]) # append all pix intensities of coordinates inside the PNN box
        print("pix mean:", (np.array(pix_list)).mean()) # convert list to array and find the mean pix intensities
        locs_list.append(locs) # append the all pixels of all PNNs detected
        mean_pix_int_list.append((np.array(pix_list)).mean()) # and their mean intensities
    print("Number of PNNs segmented:", len(locs_list))
    print("lengths", len(locs_list), len(mean_pix_int_list))
    img_info_df['pixels'] = locs_list
    img_info_df['mean_pixel_int'] = mean_pix_int_list
    img_info_df = img_info_df[['img_file_name', 'type_of_object_str', 'x1', 'y1', 'x2', 'y2', 'x3', 'y3', 'x4', 'y4', 'xc', 'yc', 'Width', 'Height', 'area', 'mean_pixel_int', 'pixels']]
    fig = plt.figure(figsize = (5, 5))
    plt.bar(list(range(len(img_info_df))), img_info_df['mean_pixel_int'], color = 'blue', width = 0.2)
    plt.xticks(np.arange(0,len(img_info_df), 1), labels = list(range(len(img_info_df))))
    plt.xlabel("Number of segmented PNNs")
    plt.ylabel("Mean pixel intensities")
    plt.title("Mean pixel intensities plot for segmented PNNs")
    plt.show()
    return contour_img, img_info_df, mean_pix_int_list # this returns a color image with PNN contours marked along with numbers


def detect_shape_pnns(contour_img, img_info_df, contours):
    # color_img = skimage.color.gray2rgb(normalised_img)
    shape = "unidentified"
    for cnts in contours:
        peri = cv2.arcLength(cnts, True) # c is the contour
        approx = cv2.approxPolyDP(cnts, 0.04 * peri, True)# gray_seg_wfa = skimage.color.rgb2gray(contour_img)
        c = max(cnts, key=cv2.contourArea) # get the area
        x,y,w,h = cv2.boundingRect(c)
        area_ = cv2.contourArea(cnts)
        # print("peri, approx", peri, approx)
        if 90 <= area_ <= 2500:
            if len(approx) == 3:
                shape = "triangle"
                cv2.putText(contour_img, shape, (x,y), cv2.FONT_HERSHEY_SIMPLEX, 1, (125, 246, 55), 3)
            elif len(approx) == 4:
                (x, y, w, h) = cv2.boundingRect(approx)
                ar = w / float(h)
                shape = "square" if ar >= 0.95 and ar <= 1.05 else "rectangle"
                cv2.putText(contour_img, shape, (x,y), cv2.FONT_HERSHEY_SIMPLEX, 1, (125, 246, 55), 3)
            elif len(approx) == 5:
                shape = "pentagon"
                cv2.putText(contour_img, shape, (x,y), cv2.FONT_HERSHEY_SIMPLEX, 1, (125, 246, 55), 3)
            else:
                shape = "circle"
                cv2.putText(contour_img, shape, (x,y), cv2.FONT_HERSHEY_SIMPLEX, 1, (125, 246, 55), 3)
        print("shape", shape)
    # rect = cv2.rectangle(contour_img, (img_info_df['x1'][box], img_info_df['y1'][box]), (img_info_df['x4'][box], img_info_df['y4'][box]), (255,255,255), 3) # draw white filled rect on the copy of the image
    # cv2.putText(contour_img, shape, (img_info_df['x1'][box], img_info_df['y1'][box]), cv2.FONT_HERSHEY_SIMPLEX, 3, (125, 246, 55), 3)
    return approx, contours, shape, contour_img


