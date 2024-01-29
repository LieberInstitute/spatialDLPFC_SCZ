% net = denoisingNetwork('DnCNN'); %Load the pretrained denoising convolutional neural network, 'DnCNN'.
% dWFA = denoiseImage(WFA,net); %denoise the image

dt = '/dcs04/lieber/marmaypag/spatialDLPFC_SCZ_LIBD4100/processed-data/VistoSeg/captureAreas/';
ot = '/dcs04/lieber/marmaypag/spatialDLPFC_SCZ_LIBD4100/plots/image_processing/image_histograms';
myfiles = dir([dt,'*1.mat']);

for i = 1:length(myfiles)
    % Read the image
    load(fullfile(dt,myfiles(i).name),'WFA');

[counts, x]=imhist(WFA);
level=triangle_th(counts(5:end),numel(x));
BW=imbinarize(WFA,level);
[x0,x1,x11,x22,y0,y1,y11,y22]=extract_coords(counts,x,level);

    figure('visible', 'off')
    ax1 = subplot(2,2,1);
    imshow(WFA,[])
    ax2 = subplot(2,2,3);
    imshow(BW)
    ax3 = subplot(2,2,2);
    plot(x,counts)
    ax4 = subplot(2,2,4);
    plot(x,counts)
    xlim(ax4, [0.005 1])
    ylim(ax4, [0 10000000]) 
    hold(ax4, 'on')
    plot(ax4,[x11,x22], [y11,y22], 'Color', 'r')
    plot(ax4,[x0,x1], [y0,y1], 'Color', 'r')
    plot(ax4,[x0,x0], [0,y0], 'Color', 'g')
    hold(ax4,'off')
    saveas(gcf,fullfile(ot,[myfiles(i).name(1:end-4), '.png']))
    close all
end