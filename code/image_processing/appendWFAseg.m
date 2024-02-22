dt = '/dcs04/lieber/marmaypag/spatialDLPFC_SCZ_LIBD4100/processed-data/VistoSeg/captureAreas/';
dt1 = '/dcs04/lieber/marmaypag/spatialDLPFC_SCZ_LIBD4100/processed-data/image_processing/WFAseg/';
ot = '/dcs04/lieber/marmaypag/spatialDLPFC_SCZ_LIBD4100/processed-data/image_processing/samui/';


%copyfile('/dcs04/lieber/marmaypag/spatialDLPFC_SCZ_LIBD4100/processed-data/VistoSeg/captureAreas/*1.mat', ot)
myfiles1 = dir([dt,'*.mat']);
myfiles = dir([dt1,'*WFAseg.mat']);

for i = 1:numel(myfiles)
    fname=fullfile(myfiles(i).folder, myfiles(i).name);
    currentName = [myfiles(i).name(1:end-11),'.mat'];
    idx = find(strcmp({myfiles1.name}', currentName));
    fname1=fullfile(myfiles1(idx).folder, myfiles1(idx).name);
    load(fname)
    BW = uint8(BW);
    BW(BW==1) = 255;
    
    load(fname1)
    imwrite(mat2gray(DAPI),fullfile(ot,[currentName(1:end-4),'.tif']))
    imwrite(mat2gray(NeuN),fullfile(ot,[currentName(1:end-4),'.tif']),'writemode', 'append')
    imwrite(mat2gray(Claudin5),fullfile(ot,[currentName(1:end-4),'.tif']),'writemode', 'append')
    imwrite(mat2gray(WFA),fullfile(ot,[currentName(1:end-4),'.tif']),'writemode', 'append')
    imwrite(mat2gray(AF),fullfile(ot,[currentName(1:end-4),'.tif']),'writemode', 'append')
    imwrite(BW,fullfile(ot,[currentName(1:end-4),'.tif']),'writemode', 'append')
disp(currentName)
end

