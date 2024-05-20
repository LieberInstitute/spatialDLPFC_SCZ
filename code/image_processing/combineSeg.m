% neunFiles = dir('/dcs04/lieber/marmaypag/spatialDLPFC_SCZ_LIBD4100/processed-data/RealPNN/single_channels_segmented/NeuN/slide3_final/*.tif');
% DAPIfiles = dir('/dcs04/lieber/marmaypag/spatialDLPFC_SCZ_LIBD4100/processed-data/RealPNN/single_channels_segmented/DAPI/slide3_final/*.tif');
D = '/dcs04/lieber/marmaypag/spatialDLPFC_SCZ_LIBD4100/processed-data/RealPNN/segmented_channels_stitched/slide3/';
% for i=1:4
%     fileP = fullfile(neunFiles(i).folder, neunFiles(i).name);
%     NeuN = imread(fileP);
%     DAPI = imread(strrep(strrep(fileP, 'NeuN', 'DAPI'), 'neun', 'dapi'));
%     save(fullfile(D,[neunFiles(i).name(1:end-19), '_thresholded.mat']), 'NeuN', 'DAPI')
% end

WFAseg = '/dcs04/lieber/marmaypag/spatialDLPFC_SCZ_LIBD4100/processed-data/image_processing/WFAseg';
myfiles = dir(fullfile(strrep(D,'slide3','*'),'*thresholded.mat'));
segs = '/dcs04/lieber/marmaypag/spatialDLPFC_SCZ_LIBD4100/processed-data/image_processing/DAPI_NeuN_WFA_Segs';
for i = 1:64
    currentSample = myfiles(i).name(1:end-16);
    fileP = fullfile(myfiles(i).folder, myfiles(i).name);
    load(fileP)
    load(fullfile(WFAseg, [currentSample, '_WFAseg3.mat']))
    WFA = pnn;
    tb = [tb; size(WFA,1), size(DAPI,1), size(NeuN, 1), size(WFA, 2), size(DAPI,2), size(NeuN, 2)];
    disp(tb)
    %save(fullfile(segs,[currentSample, '_segs.mat']))
end