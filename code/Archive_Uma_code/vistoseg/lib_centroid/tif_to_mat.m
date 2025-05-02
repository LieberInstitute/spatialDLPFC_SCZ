fname = '/dcs04/lieber/marmaypag/spatialDLPFC_SCZ_LIBD4100/processed-data/RealPNN/all_channels_segemented/Test2/V12F14-053_D1.tif';
DAPI = imread(fname,1);
Claudin5 = imread(fname,2);
NeuN = imread(fname,3);
WFA = imread(fname,4);
AF = imread(fname,5);
save(fullfile('/dcs04/lieber/marmaypag/spatialDLPFC_SCZ_LIBD4100/processed-data/RealPNN/all_channels_segemented/Test2/', 'V12F14-053_D1.mat'), 'DAPI', 'Claudin5','NeuN','WFA','AF', '-v7.3')
