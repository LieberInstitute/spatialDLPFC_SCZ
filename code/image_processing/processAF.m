dt = '/dcs04/lieber/marmaypag/spatialDLPFC_SCZ_LIBD4100/processed-data/VistoSeg/captureAreas/';
dt1 = '/dcs04/lieber/marmaypag/spatialDLPFC_SCZ_LIBD4100/processed-data/image_processing/AFseg/';

ot = '/dcs04/lieber/marmaypag/spatialDLPFC_SCZ_LIBD4100/plots/image_processing/image_histograms';
myfiles = dir([dt,'*1.mat']);
addpath(genpath('/dcs04/lieber/marmaypag/spatialDLPFC_SCZ_LIBD4100/code/image_processing'))
for i = 2:64

load(fullfile(dt,myfiles(i).name), 'AF');
[counts, x]=hist(AF(:),256);
BW = imbinarize(AF);

figure('visible', 'off')
    ax1 = subplot(2,2,1);
    imshow(AF,[])
    ax2 = subplot(2,2,3);
    imshow(BW)
    ax3 = subplot(2,2,2);
    plot(x,counts)
    ax4 = subplot(2,2,4);
    plot(x,counts)
    saveas(gcf,fullfile(ot,[myfiles(i).name(1:end-4),'AF.png']))
    save(fullfile(dt1,[myfiles(i).name(1:end-4), '_AFseg.mat']),'BW')

    close all
    disp(myfiles(i).name)
end
