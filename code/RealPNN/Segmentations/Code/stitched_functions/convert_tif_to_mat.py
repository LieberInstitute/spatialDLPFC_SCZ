
import os
import numpy as np
import tifffile
from scipy.io import savemat

def convert_tiff_to_mat(directory):
    for filename in os.listdir(directory):
        if filename.endswith(".tif") or filename.endswith(".tiff"):
            tiff_path = os.path.join(directory, filename)
            mat_path = os.path.splitext(tiff_path)[0] + ".mat"
            tiff_data = tifffile.imread(tiff_path)
            mat_data = {"data": np.array(tiff_data)}
            savemat(mat_path, mat_data)


convert_tiff_to_mat('/dcs04/lieber/marmaypag/spatialDLPFC_SCZ_LIBD4100/processed-data/RealPNN/capture_area_segmentations/DAPI/Segmented_images/')
convert_tiff_to_mat('/dcs04/lieber/marmaypag/spatialDLPFC_SCZ_LIBD4100/processed-data/RealPNN/capture_area_segmentations/NeuN/NeuN_segmented_binary/')
convert_tiff_to_mat('/dcs04/lieber/marmaypag/spatialDLPFC_SCZ_LIBD4100/processed-data/RealPNN/capture_area_segmentations/Claudin/Segmented_images_binary/')
convert_tiff_to_mat('/dcs04/lieber/marmaypag/spatialDLPFC_SCZ_LIBD4100/processed-data/RealPNN/capture_area_segmentations/WFA/Segmented_images_binary/')
convert_tiff_to_mat('/dcs04/lieber/marmaypag/spatialDLPFC_SCZ_LIBD4100/processed-data/RealPNN/capture_area_segmentations/DAPI/DAPI_binarized/')
convert_tiff_to_mat('/dcs04/lieber/marmaypag/spatialDLPFC_SCZ_LIBD4100/processed-data/RealPNN/capture_area_segmentations/Claudin/claudin_binarized/')
