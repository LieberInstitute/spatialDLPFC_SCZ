# segmentation function
# def segment_gen_dataframe(normalised_img, label):
#     color_img = skimage.color.gray2rgb((np.array((normalised_img * 255), dtype = np.uint8)))
#     hierachy, img_threshold = cv2.threshold((np.array((normalised_img * 255), dtype = np.uint8)), 100, 255, cv2.THRESH_BINARY)
#     contours,_ = cv2.findContours(img_threshold, cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_SIMPLE)
#     x, y, w, h, area = [],[],[],[],[]
#     for cnt in contours:
#         x_, y_, w_, h_ = cv2.boundingRect(cnt)
#         if(w_*h_ >= 100):
#             area_ = cv2.contourArea(cnt)
#             # print(ax,ay,aw,ah)
#             x.append(x_)
#             y.append(y_)
#             w.append(w_)
#             h.append(h_)
#             area.append(area_)
#             bb_img = cv2.rectangle(color_img, (x_,y_), (x_+w_+5, y_+h_+5), (255,0,0), 2)
#             box = np.int0(cv2.boxPoints(cv2.minAreaRect(cnt)))
#             contour_img = cv2.drawContours(bb_img,[box],0,(0,0,255),1)
#     col_names = ['img_file_name','type_of_object_str', 'x1', 'y1', 'Width', 'Height', 'area']
#     file_name = os.path.basename(img_test) # image file name
#     dict = {col_names[0]: file_name, col_names[1]: label, col_names[2]: x, col_names[3]: y, col_names[4]: w, col_names[5]: h, col_names[6]: area}
#     img_info_df = pd.DataFrame(dict, columns = col_names)
#     img_info_df['x2'] = img_info_df['x1'] + img_info_df['Width']
#     img_info_df['y2'], img_info_df['x3'] = img_info_df['y1'], img_info_df['x1']
#     img_info_df['y3'] = img_info_df['y1'] + img_info_df['Height']
#     img_info_df['x4'], img_info_df['y4'] = img_info_df['x2'], img_info_df['y3']
#     img_info_df = img_info_df[['img_file_name', 'type_of_object_str', 'x1', 'y1', 'x2', 'y2', 'x3', 'y3', 'x4', 'y4', 'Width', 'Height', 'area']]
#     return(contour_img, img_info_df)


# read and normalise the image
def read_norm(filepath, ch_num):
    img = Image.open(filepath)
    img.seek(ch_num)
    if ch_num == 0: # DAPI
        dapi = cv2.normalize(np.array(img, dtype = 'float32'), np.zeros(np.array(img, dtype = 'float32').shape, np.double), 1.0, 0.0, cv2.NORM_MINMAX)
        dapi_clr = skimage.color.gray2rgb((np.array((dapi * 255), dtype = np.uint8))) # convert to color to draw colored bb
        return dapi, dapi_clr
    if ch_num == 1: # claudin
        img_claudin = cv2.normalize(np.array(img, dtype = 'float32'), np.zeros(np.array(img, dtype = 'float32').shape, np.double), 1.0, 0.0, cv2.NORM_MINMAX)
        img_claudin[img_claudin <= img_claudin.mean()] = 0.0
        img_claudin[img_claudin >= img_claudin.mean()] = 1.0
        return img_claudin
    if ch_num == 2: #NeuN
        img_neun = cv2.normalize(np.array(img, dtype = 'float32'), np.zeros(np.array(img, dtype = 'float32').shape, np.double), 1.0, 0.0, cv2.NORM_MINMAX)
        # neun_gray = skimage.color.rgb2gray(img_neun) # convert to gray to find contours
        return img_neun
    else: # wfa
        img_arr = np.array(img, dtype = 'float32')
        img_arr[img_arr <= img_arr.mean()] = 0.0
        img_arr[img_arr >= 1.0] = img_arr.max()
        img_wfa = cv2.normalize(img_arr, np.zeros(img_arr.shape, np.double), 1.0, 0.0, cv2.NORM_MINMAX)
        return img_wfa


def plot_img(original_img, segmented_img):
    fig,ax = plt.subplots(nrows = 1, ncols = 2, figsize = (20,20))
    ax[0].imshow(original_img)
    ax[0].title.set_text('Original')
    ax[1].imshow(segmented_img)
    ax[1].title.set_text('Segemented')
    fig.show()

