dt = '/dcs04/lieber/marmaypag/spatialDLPFC_SCZ_LIBD4100/processed-data/VistoSeg/captureAreas/';
ot = '/dcs04/lieber/marmaypag/spatialDLPFC_SCZ_LIBD4100/plots/image_processing/image_histograms';
myfiles = dir([dt,'*1.mat']);

for i = 1:64

load(fullfile(dt,myfiles(i).name), 'WFA');
[counts, x]=hist(WFA(:),256);

lim = 13;
% Plot the smoothed histogram
[thresh, lnP] = triangle_threshold_right_tail(counts(lim:end));
level = x(lim)+x(thresh);
BW=imbinarize(WFA,level);
%imshow(BW)

figure('visible', 'off')
    ax1 = subplot(2,2,1);
    imshow(WFA,[])
    ax2 = subplot(2,2,3);
    imshow(BW)
    ax3 = subplot(2,2,2);
    plot(x,counts)
    ax4 = subplot(2,2,4);
    plot(x(lim:end),counts(lim:end))
    hold(ax4, 'on')
    plot(ax4,x(lnP(1:2)), lnP(3:4), 'Color', 'r')
    xline(level,'g')
    hold(ax4,'off')
    saveas(gcf,fullfile(ot,[myfiles(i).name(1:end-4), 'update.png']))
    save(fullfile(dt,[myfiles(i).name(1:end-4), '_WFAseg.mat']),'BW')
    close all
end
