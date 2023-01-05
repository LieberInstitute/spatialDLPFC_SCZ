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



# detect contours in the normalised_img
def detect_contours(normalised_img):
    hierachy, img_threshold = cv2.threshold((np.array((normalised_img * 255), dtype = np.uint8)), 100, 255, cv2.THRESH_BINARY)
    contours,_ = cv2.findContours(img_threshold, cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_SIMPLE)
    return contours



# draw the extracted contours onto the image
from __future__ import print_function
from skimage.feature import peak_local_max
from skimage.segmentation import find_boundaries, watershed
from scipy import ndimage
import imutils
def draw_contours(normalised_img, ch_num, contours = None,  color = None, thickness = None, dapi_clr = None):
    if ch_num == 0:
        shifted = cv2.pyrMeanShiftFiltering(dapi_clr, 21, 51) #dapi_clr
        gray = cv2.cvtColor(shifted, cv2.COLOR_BGR2GRAY)
        thresh = cv2.threshold(gray, 0, 255, cv2.THRESH_BINARY | cv2.THRESH_OTSU)[1]
        D = ndimage.distance_transform_edt(thresh) # Euclidean distance from binary to nearest 0-pixel
        localMax = peak_local_max(D, indices=False, min_distance=5, labels=thresh) # find the local maxima for all the individual objects
        markers = ndimage.label(localMax, structure=np.ones((3, 3)))[0] # 8-connectivity connected component analysis
        labels = watershed(-D, markers, mask=thresh)
        dpx, dpy, dpw, dph, area = [], [], [], [], []
        for label in np.unique(labels):
            if label == 0: # label marked 0 are background
                continue
            mask = np.zeros(gray.shape, dtype="uint8") # create masks that only have the detected labels as foreground and 0 as background
            mask[labels == label] = 255
            # detect contours in the mask and grab the largest one
            cnts = cv2.findContours(mask.copy(), cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_SIMPLE) # detect the watershed contours
            cnts = imutils.grab_contours(cnts) # extract only the contours
            c = max(cnts, key=cv2.contourArea) # get the area
            x,y,w,h = cv2.boundingRect(c) # BB coordinates
            area1 = cv2.contourArea(c)
            if area1 <= 100:
                dpx.append(x)
                dpy.append(y)
                dpw.append(w)
                dph.append(h)
                area.append(cv2.contourArea(c))
                ws_img_bb = cv2.rectangle(dapi_clr, (x,y), (x+w, y+h), (0,255,0), 1) # if a colored BB is not required then, change color to (0,0,0) and thickness to 1
        return dpx, dpy, dpw, dph, area, ws_img_bb
    elif ch_num == 3:
        color_img = skimage.color.gray2rgb((np.array((normalised_img * 255), dtype = np.uint8)))
        x, y, w, h, area = [],[],[],[],[]
        for cnt in contours:
            x_, y_, w_, h_ = cv2.boundingRect(cnt)
            area_ = cv2.contourArea(cnt)
            if area >= 1000:
                contour_img = cv2.rectangle(color_img, (x1-10,y1-10), (x1+w1+10, y1+h1+10), (0,0,0), -1) # eliminating all the big objects
            elif 100 <= area < 2000: # size threshold
                # print(ax,ay,aw,ah)
                x.append(x_)
                y.append(y_)
                w.append(w_)
                h.append(h_)
                area.append(area_)
                bb_img = cv2.rectangle(color_img, (x_,y_), (x_+w_+10, y_+h_+10), color, thickness) #(255,0,0), 2-- to draw colored boxes
                box = np.int0(cv2.boxPoints(cv2.minAreaRect(cnt)))
                contour_img = cv2.drawContours(bb_img,[box],0,(0,0,0),1) # change the color and thickness here if contours need to be visible
        return (x,y,w,h, area, contour_img)
    else:
        color_img = skimage.color.gray2rgb((np.array((normalised_img * 255), dtype = np.uint8)))
        x, y, w, h, area = [],[],[],[],[]
        for cnt in contours:
            x_, y_, w_, h_ = cv2.boundingRect(cnt)
            area_ = cv2.contourArea(cnt)
            if area_ >= 100:
                # area_ = cv2.contourArea(cnt)
                # print(ax,ay,aw,ah)
                x.append(x_)
                y.append(y_)
                w.append(w_)
                h.append(h_)
                area.append(area_)
                bb_img = cv2.rectangle(color_img, (x_,y_), (x_+w_+10, y_+h_+10), color, thickness) #(255,0,0), 2-- to draw colored boxes
                box = np.int0(cv2.boxPoints(cv2.minAreaRect(cnt)))
                contour_img = cv2.drawContours(bb_img,[box],0,(0,0,0),1) # change the color and thickness here if contours need to be visible
        return (x,y,w,h, area, contour_img)








# def draw_contours(contours, normalised_img, color, thickness):
#     color_img = skimage.color.gray2rgb((np.array((normalised_img * 255), dtype = np.uint8)))
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
#             bb_img = cv2.rectangle(color_img, (x_,y_), (x_+w_+10, y_+h_+10), color, thickness) #(255,0,0), 2-- to draw colored boxes
#             box = np.int0(cv2.boxPoints(cv2.minAreaRect(cnt)))
#             contour_img = cv2.drawContours(bb_img,[box],0,(0,0,0),1) # change the color and thickness here if contours need to be visible
#     return (x,y,w,h, area, contour_img)


######### DAPI segmentation functions
from __future__ import print_function
from skimage.feature import peak_local_max
from skimage.segmentation import find_boundaries, watershed
from scipy import ndimage
import imutils
def morph_transform(image_clr):
    shifted = cv2.pyrMeanShiftFiltering(image_clr, 21, 51) #dapi_clr
    gray = cv2.cvtColor(shifted, cv2.COLOR_BGR2GRAY)
    thresh = cv2.threshold(gray, 0, 255, cv2.THRESH_BINARY | cv2.THRESH_OTSU)[1]
    # fig,ax = plt.subplots(figsize = (20,20))
    # ax.imshow(image_clr)
    # fig.show()
    return shifted, thresh, gray



def find_labels(threshold):
    D = ndimage.distance_transform_edt(threshold) # Euclidean distance from binary to nearest 0-pixel
    localMax = peak_local_max(D, indices=False, min_distance=5, labels=threshold) # find the local maxima for all the individual objects
    markers = ndimage.label(localMax, structure=np.ones((3, 3)))[0] # 8-connectivity connected component analysis
    labels = watershed(-D, markers, mask=threshold)
    # print("{} unique segments found".format(len(np.unique(labels)) - 1))
    return labels



# extract the watershed algorithm labels
def draw_rect_dapi(labels, gray, dapi):
    dpx, dpy, dpw, dph, area = [], [], [], [], []
    for label in np.unique(labels):
        if label == 0: # label marked 0 are background
            continue
        mask = np.zeros(gray.shape, dtype="uint8") # create masks that only have the detected labels as foreground and 0 as background
        mask[labels == label] = 255
        # detect contours in the mask and grab the largest one
        cnts = cv2.findContours(mask.copy(), cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_SIMPLE) # detect the watershed contours
        cnts = imutils.grab_contours(cnts) # extract only the contours
        c = max(cnts, key=cv2.contourArea) # get the area
        x,y,w,h = cv2.boundingRect(c) # BB coordinates
        area1 = cv2.contourArea(c)
        if area1 <= 100:
            dpx.append(x)
            dpy.append(y)
            dpw.append(w)
            dph.append(h)
            area.append(cv2.contourArea(c))
            ws_img_bb = cv2.rectangle(dapi, (x,y), (x+w, y+h), (0,255,0), 1) # if a colored BB is not required then, change color to (0,0,0) and thickness to 1
    return dpx, dpy, dpw, dph, area, ws_img_bb



def draw_rect(df_manual_test, contour_img, color):
    for box in range(len(df_manual_test['x1'])):
        # print(box)
        rect = cv2.rectangle(contour_img, (df_manual_test['x1'][box], df_manual_test['y1'][box]), (df_manual_test['x4'][box], df_manual_test['y4'][box]), color, 2)
    # return contour_img



def draw_single_rect(df_manual_test, contour_img):
    cv2.rectangle(contour_img, (df_manual_test['x1'], df_manual_test['y1']), (df_manual_test['x4'], df_manual_test['y4']), (255,211,155), 2)
    # return contour_img



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



def manual_annot(filepath):
    conv_factor = 2.0112375738 # the fiji annotations are measured in microns which need to be translated to pixels (1860/924.81)
    df_manual_test = pd.read_csv(filepath) # read the manual annotations csv into dataframe
    df_manual_test = df_manual_test.rename(columns = {'X': 'xc', 'Y': 'yc', 'BX': 'x1', 'BY': 'y1', 'Perim.': 'Perimeter'}) # xc,yc are the centroids of the BB
    df_manual_test.loc[:,['xc']], df_manual_test.loc[:,['yc']], df_manual_test.loc[:,['x1']], df_manual_test.loc[:,['y1']], df_manual_test['Width'], df_manual_test['Height'] = df_manual_test['xc']*conv_factor, df_manual_test['yc']*conv_factor, df_manual_test['x1']*conv_factor, df_manual_test['y1']*conv_factor, df_manual_test['Width']*conv_factor, df_manual_test['Height']*conv_factor
    df_manual_test['x2'] = (df_manual_test['x1'] + df_manual_test['Width'])
    df_manual_test['y2'], df_manual_test['x3'] = df_manual_test['y1'], df_manual_test['x1']
    df_manual_test['y3'] = (df_manual_test['y1'] + df_manual_test['Height'])
    df_manual_test['x4'], df_manual_test['y4']  = df_manual_test['x2'], df_manual_test['y3'] # Calculating all 4 coordinates of the BB
    df_manual_test['xc'], df_manual_test['yc'], df_manual_test['x1'], df_manual_test['y1'] = np.int0(np.ceil(df_manual_test['xc'])), np.int0(np.ceil(df_manual_test['yc'])), np.int0(np.ceil(df_manual_test['x1'])), np.int0(np.ceil(df_manual_test['y1'])) # convert x,y,bx,by from floating point to integers (doing it after, reduces round off errors)
    df_manual_test['x2'], df_manual_test['y2'], df_manual_test['x3'], df_manual_test['y3'], df_manual_test['x4'], df_manual_test['y4'] = np.int0(np.ceil(df_manual_test['x2'])), np.int0(np.ceil(df_manual_test['y2'])), np.int0(np.ceil(df_manual_test['x3'])), np.int0(np.ceil(df_manual_test['y3'])), np.int0(np.ceil(df_manual_test['x4'])), np.int0(np.ceil(df_manual_test['y4']))
    df_manual_test = df_manual_test[['Area', 'Perimeter', 'Mean', 'Min', 'Max', 'xc', 'yc', 'x1', 'y1', 'x2', 'y2', 'x3', 'y3' , 'x4' , 'y4', 'Width', 'Height', 'Ch']] # rearranging the columns
    return df_manual_test



# plot the segmented and original image
def plot_img(original_img, segmented_img):
    fig,ax = plt.subplots(nrows = 1, ncols = 2, figsize = (20,20))
    ax[0].imshow(original_img)
    ax[0].title.set_text('Original')
    ax[1].imshow(segmented_img)
    ax[1].title.set_text('Segmented')
    fig.show()



# plot histogram
def hist_plot(img):
    range = (img.min(), img.max())
    histogram, bin_edges = np.histogram(img, bins=256, range=range)
    plt.figure()
    plt.title("Grayscale Histogram")
    plt.xlabel("grayscale value")
    plt.ylabel("pixel count")
    plt.xlim([img.min(), img.max()])
    plt.plot(bin_edges[0:-1], histogram)
    plt.show()




# plot histogram improved
import pylab
from pylab import xticks
def histo(img,range = [0,1]):
    n, bins, patches = plt.hist(img, 30, range = range, facecolor='gray', align='mid')
    pylab.rc("axes", linewidth=8.0)
    pylab.rc("lines", markeredgewidth=2.0)
    plt.xlabel('pix int', fontsize=14)
    plt.ylabel('# of targets', fontsize=14)
    pylab.xticks(fontsize=15, rotation = 'vertical')
    pylab.yticks(fontsize=15)
    plt.grid(True)
    plt.show()


