fname = 'path and name to your tif file';
DAPI = imread(fname,1);
Claudin5 = imread(fname,2);
NeuN = imread(fname,3);
WFA = imread(fname,4);
AF = imread(fname,5);
save(fullfile('/dcs04/lieber/marmaypag/spatialDLPFC_SCZ_LIBD4100/processed-data/RealPNN/capture_area_segmentations/all_channels_segemented/Test/', 'V12F14-053_A1_thresholded.mat'), 'DAPI', 'Claudin5','NeuN','WFA','AF', '-v7.3')