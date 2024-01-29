dt = '/dcs04/lieber/marmaypag/spatialDLPFC_SCZ_LIBD4100/processed-data/VistoSeg/captureAreas/';
myfiles = dir([dt,'*1.mat']);

% Specify colors for density curves
colors = hsv(length(myfiles)); % Use the HSV color space for better distinguishability
tb = table('Size', [0, 6], 'VariableNames', {'Sample', 'Mean', 'Median', 'Mode', 'Kurtosis', 'Width'}, 'VariableTypes', {'string', 'single', 'single',  'single',  'single',  'single'});

% Plot density curves on the same axes

for i = 1:length(myfiles)
    % Read the image
    load(fullfile(dt,myfiles(i).name));
    
    % Compute the histogram
    [counts, bi] = imhist(WFA);
    
    % Normalize counts to create probability density function (PDF)
    pdf = counts / sum(counts);
    
    % Plot the density curve with the specified color
    % hold on;
    figure('visible', 'off')
    ax1 = subplot(1,2,1);
    imshow(WFA,[])
    ax2 = subplot(1,2,2);
    plot(bi,counts, 'Color', colors(i, :));
    xlim(ax2, [0.005 1])
    ylim(ax2, [0 6000000])     
    saveas(gcf,fullfile(dt,[myfiles(i).name, '.png']))
    close all
N = cellstr(myfiles(i).name(1:end-4));   
m1 = mean(WFA(:));
m2 = median(WFA(:));
m3 = mode(WFA(:));
k = kurtosis(WFA(:));
width = fwhm(bi,counts);
temp = table(N,m1,m2,m3,k,width, 'VariableNames', {'Sample', 'Mean', 'Median', 'Mode', 'Kurtosis', 'Width'});
tb = [tb; temp];
end
writetable(tb, fullfile(dt,'image_metrics.csv'));

% hold off;
% title('Density Curves of Image Histograms');
% xlabel('Pixel Intensity');
% ylabel('Probability Density');
% %grid on;
% xlim([5 300])
% saveas(gcf,fullfile(dt,'raw_histograms.png'))



