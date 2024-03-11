dt = '/dcs04/lieber/marmaypag/spatialDLPFC_SCZ_LIBD4100/processed-data/VistoSeg/captureAreas/';
dt1 = '/dcs04/lieber/marmaypag/spatialDLPFC_SCZ_LIBD4100/processed-data/image_processing/WFAseg/';
ot = '/dcs04/lieber/marmaypag/spatialDLPFC_SCZ_LIBD4100/processed-data/image_processing/samui/';
at = '/dcs04/lieber/marmaypag/spatialDLPFC_SCZ_LIBD4100/processed-data/image_processing/AFseg/';

%copyfile('/dcs04/lieber/marmaypag/spatialDLPFC_SCZ_LIBD4100/processed-data/VistoSeg/captureAreas/*1.mat', ot)
WFAfiles = dir([dt,'*.mat']);

for i = 1:numel(WFAfiles)
    fname1=fullfile(WFAfiles(i).folder, WFAfiles(i).name);
    load(fname1)
    
    currentName = WFAfiles(i).name(1:end-4);

    % fname=fullfile(dt1, [currentName, '_WFAseg.mat']);
    % load(fname)
    % seg = uint8(BW);
    % seg(seg==1) = 255;
    
    % fname=fullfile(dt1, [currentName, '_WFAseg1.mat']);
    % load(fname)
    % seg1 = uint8(BW);
    % seg1(seg1==1) = 255;

    % fname=fullfile(dt1, [currentName, '_WFAseg2.mat']);
    % load(fname)
    % seg2 = uint8(pnn);
    % seg2(seg2==1) = 255;

    fname=fullfile(dt1, [currentName, '_WFAseg3.mat']);
    load(fname)
    seg3 = uint8(pnn);
    seg3(seg3==1) = 255;
    
    fname=fullfile(at, [currentName, '_AFseg.mat']);
    load(fname)
    AFseg = uint8(BW);
    AFseg(AFseg==1) = 255;

    imwrite(mat2gray(DAPI),fullfile(ot,[currentName,'.tif']))
    imwrite(mat2gray(NeuN),fullfile(ot,[currentName,'.tif']),'writemode', 'append')
    imwrite(mat2gray(Claudin5),fullfile(ot,[currentName,'.tif']),'writemode', 'append')
    imwrite(mat2gray(WFA),fullfile(ot,[currentName,'.tif']),'writemode', 'append')
    imwrite(mat2gray(AF),fullfile(ot,[currentName,'.tif']),'writemode', 'append')
    %imwrite(seg,fullfile(ot,[currentName,'.tif']),'writemode', 'append')
    %imwrite(seg1,fullfile(ot,[currentName,'.tif']),'writemode', 'append')
    %imwrite(seg2,fullfile(ot,[currentName,'.tif']),'writemode', 'append')
    imwrite(seg3,fullfile(ot,[currentName,'.tif']),'writemode', 'append')
    imwrite(AFseg,fullfile(ot,[currentName,'.tif']),'writemode', 'append')
disp(currentName)
end

