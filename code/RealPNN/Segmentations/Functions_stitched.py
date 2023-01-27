'''
For Stitched Visium-IF tissue sections from VistoSeg SplitSlide output
Channel0 = AF
Channel1 = Claudin - 5 (Alex 488),
Channel2 = DAPI,
Channel3 = NeuN,
Channel4 = WFA
'''

# read and normalise the image
Image.MAX_IMAGE_PIXELS = None
def read_norm(filepath, ch_num):
    img = Image.open(filepath)
    img.seek(ch_num)
    if ch_num == 1: # claudin
        img_claudin = cv2.normalize(np.array(img, dtype = 'uint8'), np.zeros(np.array(img, dtype = 'uint8').shape, np.double), 1.0, 0.0, cv2.NORM_MINMAX)
        img_claudin[img_claudin <= img_claudin.mean()] = 0.0
        img_claudin[img_claudin >= img_claudin.mean()] = 1.0
        return img_claudin
    if ch_num == 2: # DAPI
        dapi_clr = skimage.color.gray2rgb((np.array(img, dtype = np.uint8))) # convert to color to draw colored bb
        dapi = cv2.normalize(np.array(img, dtype = 'uint8'), np.zeros(np.array(img, dtype = 'uint8').shape, np.double), 1.0, 0.0, cv2.NORM_MINMAX)
        return dapi, dapi_clr
    if ch_num == 3: #NeuN
        img_neun = cv2.normalize(np.array(img, dtype = 'uint8'), np.zeros(np.array(img, dtype = 'uint8').shape, np.double), 1.0, 0.0, cv2.NORM_MINMAX)
        # neun_gray = skimage.color.rgb2gray(img_neun) # convert to gray to find contours
        return img_neun
    else: # wfa
        img_arr = np.array(img, dtype = 'uint8')
        # histo(img_arr,range = [img_arr.min(),img_arr.max()])
        img_wfa = cv2.normalize(img_arr, np.zeros(img_arr.shape, np.double), 1.0, 0.0, cv2.NORM_MINMAX)
        return img_arr, img_wfa


# detect contours in the normalised_img
def detect_contours(normalised_img): ### create a separate function for shape detection and run it through this loop for contour detection
    hierachy, img_threshold = cv2.threshold((np.array((normalised_img * 255), dtype = np.uint8)), 100, 255, cv2.THRESH_BINARY)
    contours,_ = cv2.findContours(img_threshold, cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_SIMPLE)
    return contours

######### DAPI segmentation functions
from __future__ import print_function
from skimage.feature import peak_local_max
from skimage.segmentation import find_boundaries, watershed
from scipy import ndimage
import imutils
def morph_transform(image_clr):
    shifted = cv2.pyrMeanShiftFiltering(image_clr, 21, 51) #dapi_clr
    print("shifted")
    gray = cv2.cvtColor(shifted, cv2.COLOR_BGR2GRAY)
    print("grayed")
    thresh = cv2.threshold(gray, 0, 255, cv2.THRESH_BINARY | cv2.THRESH_OTSU)[1]
    print("thresholded")
    # fig,ax = plt.subplots(figsize = (20,20))
    # ax.imshow(image_clr)
    # fig.show()
    return shifted, thresh, gray



def find_labels(threshold):
    D = ndimage.distance_transform_edt(threshold) # Euclidean distance from binary to nearest 0-pixel
    print("distance measured")
    localMax = peak_local_max(D, indices=False, min_distance=5, labels=threshold) # find the local maxima for all the individual objects
    print("local max found")
    markers = ndimage.label(localMax, structure=np.ones((3, 3)))[0] # 8-connectivity connected component analysis
    print("local max markers found")
    labels = watershed(-D, markers, mask=threshold)
    print("{} unique segments found".format(len(np.unique(labels)) - 1))
    return labels



# extract the watershed algorithm labels
def draw_rect_dapi(labels, gray, dapi_clr):
    dpx, dpy, dpw, dph, area = [], [], [], [], []
    print("1) entering the label loop")
    for label in np.unique(labels):
        if label == 0: # label marked 0 are background
            continue
        mask = np.zeros(gray.shape, dtype="uint8") # create masks that only have the detected labels as foreground and 0 as background
        mask[labels == label] = 255
        print("2) found the masks")
        # detect contours in the mask and grab the largest one
        cnts = cv2.findContours(mask.copy(), cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_SIMPLE) # detect the watershed contours
        print("3) found contours")
        cnts = imutils.grab_contours(cnts) # extract only the contours
        print("4) grabbed contours")
        c = max(cnts, key=cv2.contourArea) # get the area
        x,y,w,h = cv2.boundingRect(c) # BB coordinates
        area1 = cv2.contourArea(c)
        print("5) found BB coordinates and area")
        if area1 <= 100:
            print("6) appending BB coordinates")
            dpx.append(x)
            dpy.append(y)
            dpw.append(w)
            dph.append(h)
            area.append(cv2.contourArea(c))
            print("7) appended")
            ws_img_bb = cv2.rectangle(dapi_clr, (x,y), (x+w, y+h), (0,255,0), 1) # if a colored BB is not required then, change color to (0,0,0) and thickness to 1
            print("8) drawing rectangles")
    return dpx, dpy, dpw, dph, area, ws_img_bb


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
            elif 90 <= area_ < 2500: # size threshold
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
        color_img = skimage.color.gray2rgb(neun, dtype = np.uint8)
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




# plot the segmented and original image
def plot_img(original_img, segmented_img):
    fig,ax = plt.subplots(nrows = 1, ncols = 2, figsize = (20,20))
    ax[0].imshow(original_img, cmap='gray')
    ax[0].title.set_text('Original')
    ax[1].imshow(segmented_img, cmap='gray')
    ax[1].title.set_text('Segmented')
    fig.show()
