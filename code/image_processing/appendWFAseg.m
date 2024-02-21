ot = '/dcs04/lieber/marmaypag/spatialDLPFC_SCZ_LIBD4100/processed-data/image_processing/samui/';
dt = '/dcs04/lieber/marmaypag/spatialDLPFC_SCZ_LIBD4100/processed-data/VistoSeg/image_processing/WFAseg/';

copyfile('/dcs04/lieber/marmaypag/spatialDLPFC_SCZ_LIBD4100/processed-data/VistoSeg/captureAreas/*1.tif', ot)

myfiles = dir([dt,'*WFAseg.mat']);
for i = 1:numel(myfiles)
    fname=fullfile(myfiles(i).folder, myfiles(i).name);
    load(fname)
    imwrite(BW,[fname1(1:end-4),'.tif'],'writemode', 'append')
