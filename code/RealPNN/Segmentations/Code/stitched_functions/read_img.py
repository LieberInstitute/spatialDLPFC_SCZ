'''
For Visium-IF
Channel0 = DAPI, DAPI
Channel1 = Claudin5 (Alex 488),
Channel2 = NeuN (Alexa 555),
Channel3 = WFA (Alexa 647),
Channel4 = AF (Autofluorescence), sample AF
Channel5 = Thumbnail
'''

# read and normalise the image
Image.MAX_IMAGE_PIXELS = None
def preprocess(filepath, ch_num):
    img = Image.open(filepath)
    img.seek(ch_num)
    if ch_num == 1: # claudin
        # img_claudin = cv2.normalize(np.array(img, dtype = 'uint8'), np.zeros(np.array(img, dtype = 'uint8').shape, np.double), 1.0, 0.0, cv2.NORM_MINMAX)
        img_claudin = np.array(img, dtype = 'uint8')
        img_claudin[img_claudin <= img_claudin.mean()] = 0
        img_claudin[img_claudin >= img_claudin.mean()] = 255
        return img_claudin
    if ch_num == 2: # DAPI
        dapi_clr = skimage.color.gray2rgb((np.array(img, dtype = np.uint8))) # convert to color to draw colored bb
        dapi = cv2.normalize(np.array(img, dtype = 'uint8'), np.zeros(np.array(img, dtype = 'uint8').shape, np.double), 1.0, 0.0, cv2.NORM_MINMAX)
        image_clr = skimage.color.gray2rgb((np.array((original_img * 255), dtype = np.uint8)))
        shifted = cv2.pyrMeanShiftFiltering(image_clr, 21, 51) #dapi_clr
        print("shifted")
        gray = cv2.cvtColor(shifted, cv2.COLOR_BGR2GRAY)
        print("grayed")
        thresh = cv2.threshold(gray, 0, 255, cv2.THRESH_BINARY | cv2.THRESH_OTSU)[1]
        print("thresholded")
        return dapi, dapi_clr
    if ch_num == 3: #NeuN
        img_neun = np.array(img, dtype = 'uint8')
        # img_neun[img_neun <= img_neun.mean()] = 0
        # img_neun[img_neun >= img_neun.mean()] = 255
        # img_neun = cv2.normalize(np.array(img, dtype = 'uint8'), np.zeros(np.array(img, dtype = 'uint8').shape, np.double), 1.0, 0.0, cv2.NORM_MINMAX)
        # neun_gray = skimage.color.rgb2gray(img_neun) # convert to gray to find contours
        return img_neun
    else: # wfa
        img_arr = np.array(img, dtype = 'uint8')
        # histo(img_arr,range = [img_arr.min(),img_arr.max()])
        img_wfa = cv2.normalize(img_arr, np.zeros(img_arr.shape, np.double), 1.0, 0.0, cv2.NORM_MINMAX)
        return img_arr, img_wfa
