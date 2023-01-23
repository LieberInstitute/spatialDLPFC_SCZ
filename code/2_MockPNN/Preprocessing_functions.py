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
        img_wfa = cv2.normalize(img_arr, np.zeros(img_arr.shape, np.double), 1.0, 0.0, cv2.NORM_MINMAX)
        return img_arr,img_wfa


# plot histogram improved
import pylab
from pylab import xticks
def histo(img, bins = 30, range = [0,1]):
    n, bins, patches = plt.hist(img.ravel(), bins = bins, range = range, facecolor='gray', align='mid') # (y, x, _)
    order = np.argsort(n)[::-1]
    # print(" highest bins:", n[order][:10])
    print("  their ranges:", [ (bins[i+1])   for i in order[:10]]) #bins[i],
    # change the contrast such that the order[3:8] are only visible and rest are all masked
    img[img <= ([(bins[i+1])   for i in order[7:8]])] = 0.0 # select the bin, below which the pix intensities will be blackened
    img[img <= ([(bins[i+1])   for i in order[2:3]])] = 0.0
    # img[img <= ([(bins[i+1])   for i in order[2:3]])] = 0.0
    # img[img <= ([(bins[i+1])   for i in order[7:]])] = 0.0
    print("0 init done!")
    # img[img >= ([(bins[i+1])   for i in order[3:4]])] = 1.0
    # img[img >= ([(bins[i+1])   for i in order[4:5]])] = 1.0
    # img[img >= ([(bins[i+1])   for i in order[5:6]])] = 1.0
    # img[img >= ([(bins[i+1])   for i in order[6:7]])] = 1.0
    img[img >= ([(bins[i+1])   for i in order[4:5]])] = 1.0
    img[img >= ([(bins[i+1])   for i in order[3:4]])] = 1.0
    img[img >= ([(bins[i+1])   for i in order[5:6]])] = 1.0
    print("1 init done!")
    # for i in order:
    #     if i>=8:
    #         print("3",i, bins[i+1], img[img >= ([(bins[i+1])])])
    #         img[img >= ([(bins[i+1])])] = 0.0
    #         print("after 3", img[img >= ([(bins[i+1])])])
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




img_arr = Image.open(img_test)
img_arr.seek(3)
img = np.array(img_arr, dtype = 'float32')
range = [img.min(),img.max()]
n, bins, patches = plt.hist(img.ravel(), bins = 30, range = range, facecolor='gray', align='mid') # (y, x, _)
order = np.argsort(n)[::-1]
print(order)
# print(" highest bins:", n[order][:10])
print("  their ranges:", [ (bins[i+1])   for i in order[:]]) #bins[i],
for i in order:
    if i >=0 and i<3:
        print("1",i, bins[i+1], img[img <= ([(bins[i+1])])])
        img[img <= ([(bins[i+1])])] = 0.0
        print("after 1",img[img <= ([(bins[i+1])])])
    elif i>=8:
        print("3",i, bins[i+1], img[img >= ([(bins[i+1])])])
        img[img >= ([(bins[i+1])])] = 0.0
        print("after 3", img[img >= ([(bins[i+1])])])
    elif i>=3 and i<8:
        print("2",i, bins[i+1], img[img >= ([(bins[i+1])])])
        # img[img >= ([(bins[i+1])])] = 1.0
#
# fig,ax = plt.subplots(figsize = (20,20))
# ax.imshow(img, cmap = 'gray')
# fig.show()
#
#
#
# # change the contrast such that the order[3:8] are only visible and rest are all masked
#     img[img <= ([(bins[i+1])   for i in order[0:1]])] = 0.0 # select the bin, below which the pix intensities will be blackened
#     img[img <= ([(bins[i+1])   for i in order[1:2]])] = 0.0
#     img[img <= ([(bins[i+1])   for i in order[2:3]])] = 0.0
#     # img[img <= ([(bins[i+1])   for i in order[7:]])] = 0.0
#     print("0 init done!")
#     img[img >= ([(bins[i+1])   for i in order[3:4]])] = 1.0
#     img[img >= ([(bins[i+1])   for i in order[4:5]])] = 1.0
#     img[img >= ([(bins[i+1])   for i in order[5:6]])] = 1.0
#     img[img >= ([(bins[i+1])   for i in order[6:7]])] = 1.0
#     img[img >= ([(bins[i+1])   for i in order[7:8]])] = 1.0
#     print("the order to be used",[(bins[i+1])   for i in order[3:8]]) # print the chosen value below which all pix intensities are considered to be noise
#     pylab.rc("axes", linewidth=8.0)
#     pylab.rc("lines", markeredgewidth=2.0)
#     xticks = [(bins[idx+1] + value)/2 for idx, value in enumerate(bins[:-1])]
#     xticks_labels = [ "{:.1f}\nto\n{:.1f}".format(value, bins[idx+1]) for idx, value in enumerate(bins[:-1])]
#     # plt.xticks(xticks, labels = xticks_labels)
#     plt.xlabel('Pixel intensities', fontsize=14)


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
