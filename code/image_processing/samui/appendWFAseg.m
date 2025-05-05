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

     fname=fullfile(dt1, [currentName, '_WFAseg.mat']);
     load(fname)
     tpoint = uint8(BW);
     tpoint(tpoint==1) = 255;
    
     fname=fullfile(dt1, [currentName, '_WFAseg1.mat']);
     load(fname)
     ostus = uint8(BW);
     ostus(ostus==1) = 255;

    % fname=fullfile(dt1, [currentName, '_WFAseg2.mat']);
    % load(fname)
    % seg2 = uint8(pnn);
    % seg2(seg2==1) = 255;

    fname=fullfile(dt1, [currentName, '_WFAseg3.mat']);
    load(fname)
    adapT = uint8(pnn);
    adapT(adapT==1) = 255;
    
    fname=fullfile(at, [currentName, '_AFseg.mat']);
    load(fname)
    AFmask = pnn;
    AFmask(BW) = 0;
    AFmask = uint8(AFmask);
    AFmask(AFmask==1) = 255;
    
    fname=fullfile(dt1, [currentName, '_AFseg_mask.mat']);
    load(fname)
    filt = uint8(filteredPnnM);
    filt(filt==1) = 255;

    imwrite(mat2gray(DAPI),fullfile(ot,[currentName,'.tif']))
    imwrite(mat2gray(NeuN),fullfile(ot,[currentName,'.tif']),'writemode', 'append')
    imwrite(mat2gray(Claudin5),fullfile(ot,[currentName,'.tif']),'writemode', 'append')
    imwrite(mat2gray(WFA),fullfile(ot,[currentName,'.tif']),'writemode', 'append')
    imwrite(mat2gray(AF),fullfile(ot,[currentName,'.tif']),'writemode', 'append')
    imwrite(tpoint,fullfile(ot,[currentName,'.tif']),'writemode', 'append')
    imwrite(ostus,fullfile(ot,[currentName,'.tif']),'writemode', 'append')
    %imwrite(seg2,fullfile(ot,[currentName,'.tif']),'writemode', 'append')
    imwrite(adapT,fullfile(ot,[currentName,'.tif']),'writemode', 'append')
    imwrite(AFmask,fullfile(ot,[currentName,'.tif']),'writemode', 'append')
    imwrite(filt,fullfile(ot,[currentName,'.tif']),'writemode', 'append')
disp(currentName)
end

