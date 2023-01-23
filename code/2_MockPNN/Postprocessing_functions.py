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

