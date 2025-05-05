fpath = '/dcs04/lieber/marmaypag/spatialDLPFC_SCZ_LIBD4100/processed-data/VistoSeg/captureAreas/V13M06-280_A1.tif';
fname = 'V13M06-280_A1';

ot = '/dcs04/lieber/marmaypag/spatialDLPFC_SCZ_LIBD4100/plots/image_processing/image_histograms';
dt1 = '/dcs04/lieber/marmaypag/spatialDLPFC_SCZ_LIBD4100/processed-data/image_processing/';

addpath(genpath('/dcs04/lieber/marmaypag/spatialDLPFC_SCZ_LIBD4100/code/image_processing/'))
%% preprocess
WFA = im2double(imread(fpath, 5));

[counts, x]=hist(WFA(:),256);
lim = 13;
% Plot the smoothed histogram
[thresh, lnP] = triangle_threshold_right_tail(counts(lim:end));
level = x(lim)+x(thresh)+0.1;
BW=imbinarize(WFA,level);

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
    xline(level+0.1,'g')
    hold(ax4,'off')
    saveas(gcf,fullfile(ot,[fname,'.png']))
    save(fullfile(dt1,'WFAseg',[fname, '_WFAseg1.mat']),'BW')
    close all

%% process_pnnSeg

WFA1 = WFA.*(BW);
noise_idx = imbinarize(nonzeros(WFA1), 'adaptive', 'sensitivity', 0.5);
pnn_idx = find(BW);
pnn = BW;
pnn(pnn_idx(~noise_idx)) =0; 
save(fullfile(dt1,'WFAseg',[fname, '_WFAseg3.mat']),'pnn')

%% processAF
AF = im2double(imread(fpath,1));
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
    saveas(gcf,fullfile(ot,[fname,'AF.png']))
    save(fullfile(dt1,'AFseg',[fname, '_AFseg.mat']),'BW')

%% maskAF
pnnM = pnn;
pnnM(BW) = 0;

stats = struct2table(regionprops(pnnM, WFA, 'MeanIntensity', 'MaxIntensity', 'Area'));
stats.Sample = repmat(fname, height(stats), 1);
if max([stats.MeanIntensity]) < 1
filteredPnnM = bwpropfilt(pnnM, "Area", [0 30000]);
else
filteredPnnM = pnnM;
end
%imshow([pnn, ones(size(pnn,1),50), pnnM, ones(size(pnn,1),50), filteredPnnM])
final = [pnn, ones(size(pnn,1),50), pnnM, ones(size(pnn,1),50), filteredPnnM];
imwrite(imresize(final, 0.5),fullfile(ot,[fname,'AFfilt.png']));
save(fullfile(dt1,'WFAseg',[fname,'_AFseg_mask.mat']),'filteredPnnM')