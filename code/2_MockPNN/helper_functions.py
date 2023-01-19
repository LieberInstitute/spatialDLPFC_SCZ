
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
        # img_arr_adj = img_arr
        histo(img_arr,range = [img_arr.min(),img_arr.max()])
        # img_arr[img_arr < 1.7598] = 0.0
        # img_arr[img_arr >= 0.6513] = 0.6513
        # img_arr[img_arr <= img_arr.mean()] = 0.0
        # img_arr[img_arr >= 1.0] = img_arr.max()
        img_wfa = cv2.normalize(img_arr, np.zeros(img_arr.shape, np.double), 1.0, 0.0, cv2.NORM_MINMAX)
        return img_arr, img_wfa



# detect contours in the normalised_img
def detect_contours(normalised_img): ### create a separate function for shape detection and run it through this loop for contour detection
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
            if area1 <= 500:
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
            if area_ >= 1000:
                contour_img = cv2.rectangle(color_img, (x_-10,y_-10), (x_+w_+10, y_+h_+10), (0,0,0), -1) # eliminating all the big objects
            elif 100 <= area_ < 2000: # size threshold
                # print(ax,ay,aw,ah)
                x.append(x_)
                y.append(y_)
                w.append(w_)
                h.append(h_)
                area.append(area_)
                bb_img = cv2.rectangle(color_img, (x_,y_), (x_+w_+10, y_+h_+10), color, thickness) #(255,0,0), 2-- to draw colored boxes
                box = np.int0(cv2.boxPoints(cv2.minAreaRect(cnt)))
                contour_img = cv2.drawContours(bb_img,[box],0,(0,0,0),1) # change the color and thickness here if contours need to be visible
                # cv2.putText(contour_img, label, (x_,y_), cv2.FONT_HERSHEY_SIMPLEX, 0.7, (125, 246, 55), 3)
        return x, y, w, h, area, contour_img
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
        return x, y, w, h, area, contour_img





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
        print(locs.shape, locs.mean())
        for i in range(locs.shape[0]):
            for j in range(locs.shape[1] -1):
                if gray_image[locs[i,j],locs[i,j+1]] == 255: # gray image has white filled boxes
                    # print(gray_seg_wfa[locs[i,j],locs[i,j+1]])
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
    return contour_img, img_info_df # this returns a color image with PNN contours marked along with numbers


def detect_shape_pnns(contour_img, img_info_df, contours):
    # color_img = skimage.color.gray2rgb(normalised_img)
    shape = "unidentified"
    for c in contours:
        peri = cv2.arcLength(c, True) # c is the contour
        approx = cv2.approxPolyDP(c, 0.04 * peri, True)# gray_seg_wfa = skimage.color.rgb2gray(contour_img)
        # print("peri, approx", peri, approx)
        if len(approx) == 3:
            shape = "triangle"
        elif len(approx) == 4:
            (x, y, w, h) = cv2.boundingRect(approx)
            ar = w / float(h)
            shape = "square" if ar >= 0.95 and ar <= 1.05 else "rectangle"
        elif len(approx) == 5:
            shape = "pentagon"
        else:
            shape = "circle"
        print("shape", shape)
    rect = cv2.rectangle(contour_img, (img_info_df['x1'][box], img_info_df['y1'][box]), (img_info_df['x4'][box], img_info_df['y4'][box]), (255,255,255), -1) # draw white filled rect on the copy of the image
    cv2.putText(contour_img, ('%d'%box), (img_info_df['x1'][box],img_info_df['y1'][box]), cv2.FONT_HERSHEY_SIMPLEX, 2, (125, 246, 55), 3)
    return approx, contours, shape



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
    ax[0].imshow(original_img, cmap='gray')
    ax[0].title.set_text('Original')
    ax[1].imshow(segmented_img, cmap='gray')
    ax[1].title.set_text('Segmented')
    fig.show()



# plot histogram
def hist_plot(img, bins=256):
    range = (img.min(), img.max())
    histogram, bin_edges = np.histogram(img.ravel(), bins=bins, range=range)
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
def histo(img, bins = 30, range = [0,1]):
    n, bins, patches = plt.hist(img.ravel(), bins = bins, range = range, facecolor='gray', align='mid') # (y, x, _)
    order = np.argsort(n)[::-1]
    # print(" highest bins:", n[order][:10])
    # print("  their ranges:", [ (bins[i+1])   for i in order[:10]]) #bins[i],
    img[img <= ([(bins[i+1])   for i in order[7:8]])] = 0.0 # select the bin, below which the pix intensities will be blackened
    img[img >= ([(bins[i+1])   for i in order[8:9]])] = 1.0
    print("the order to be used",[(bins[i+1])   for i in order[7:8]]) # print the chosen value below which all pix intensities are considered to be noise
    pylab.rc("axes", linewidth=8.0)
    pylab.rc("lines", markeredgewidth=2.0)
    xticks = [(bins[idx+1] + value)/2 for idx, value in enumerate(bins[:-1])]
    xticks_labels = [ "{:.1f}\nto\n{:.1f}".format(value, bins[idx+1]) for idx, value in enumerate(bins[:-1])]
    # plt.xticks(xticks, labels = xticks_labels)
    plt.xlabel('Pixel intensities', fontsize=14)
    plt.ylabel('Frequency', fontsize=14)
    plt.title('Pixel intensitiy plot')
    pylab.xticks(fontsize=15, rotation = 'vertical')
    pylab.yticks(fontsize=15)
    # for idx, value in enumerate(bins[:-1]):
    #     # print(idx, value)
    #     if 0.0 < n[idx] <= 0.25*pow(10, 4):
    #         print("in loop",idx, value)
    # plt.grid(True)
    # plt.show()


