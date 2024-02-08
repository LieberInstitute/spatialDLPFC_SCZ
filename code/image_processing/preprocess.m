slide = 'V12F14-053';
%slide = 'V12F14-057';

dt = '/dcs04/lieber/marmaypag/spatialDLPFC_SCZ_LIBD4100/processed-data/VistoSeg/captureAreas/';
ot = '/dcs04/lieber/marmaypag/spatialDLPFC_SCZ_LIBD4100/plots/image_processing/image_histograms';
myfiles = dir([dt,slide,'*1.mat']);

for i = 1:4

load(fullfile(dt,myfiles(i).name), 'WFA');
[counts, x]=imhist(WFA);

% Assuming 'counts' and 'x' contain the histogram data
windowSize = 3; % Choose a suitable window size for smoothing

% Apply moving average filter to counts
smoothed_counts = movmean(counts, windowSize);

% Plot the smoothed histogram
level=triangle_th(smoothed_counts(35:end),numel(smoothed_counts(35:end))) + x(34);
BW=imbinarize(WFA,level);
%imshow(BW)

[x0,x1,x11,x22,y0,y1,y11,y22]=extract_coords(smoothed_counts(30:end),x(30:end),level);
figure('visible', 'off')
    ax1 = subplot(2,2,1);
    imshow(WFA,[])
    ax2 = subplot(2,2,3);
    imshow(BW)
    ax3 = subplot(2,2,2);
    plot(x,counts)
    ax4 = subplot(2,2,4);
    plot(x(35:end),smoothed_counts(35:end))
    hold(ax4, 'on')
    plot(ax4,[x11,x22], [y11,y22], 'Color', 'r')
    plot(ax4,[x0,x1], [y0,y1], 'Color', 'r')
    plot(ax4,[x0,x0], [0,y0], 'Color', 'g')
    hold(ax4,'off')
    saveas(gcf,fullfile(ot,[myfiles(i).name(1:end-4), '.png']))
    close all
end

