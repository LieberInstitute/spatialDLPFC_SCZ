% net = denoisingNetwork('DnCNN'); %Load the pretrained denoising convolutional neural network, 'DnCNN'.
% dWFA = denoiseImage(WFA,net); %denoise the image

dt = '/dcs04/lieber/marmaypag/spatialDLPFC_SCZ_LIBD4100/processed-data/VistoSeg/captureAreas/';
ot = '/dcs04/lieber/marmaypag/spatialDLPFC_SCZ_LIBD4100/plots/image_processing/image_histograms';
myfiles = dir([dt,'*1.mat']);

for i = 1:length(myfiles)
    % Read the image
    load(fullfile(dt,myfiles(i).name),'WFA');

[lehisto, x]=imhist(WFA);
level=triangle_th(lehisto(5:end),numel(x));
 BW=imbinarize(WFA,level);

    figure('visible', 'off')
    ax1 = subplot(2,2,1);
    imshow(WFA,[])
    ax2 = subplot(2,2,3);
    imshow(BW)
    ax3 = subplot(2,2,2);
    plot(x,lehisto)
    ax4 = subplot(2,2,4);
    plot(x,lehisto)
    xlim(ax4, [0.005 1])
    ylim(ax4, [0 6000000]) 
    hold(ax4, 'on')
    xline(ax4,level, 'Color', 'r')
    hold(ax4,'off')
    saveas(gcf,fullfile(ot,[myfiles(i).name(1:end-4), '.png']))
    close all
end