dt = '/dcs04/lieber/marmaypag/spatialDLPFC_SCZ_LIBD4100/processed-data/VistoSeg/captureAreas/';
ot = '/dcs04/lieber/marmaypag/spatialDLPFC_SCZ_LIBD4100/plots/image_processing/ALLhistos/';
myfiles = dir([dt,'*1.mat']);

% Specify colors for density curves
colors = hsv(length(myfiles)); % Use the HSV color space for better distinguishability
tb = table('Size', [0, 6], 'VariableNames', {'Sample', 'DAPI', 'NeuN', 'WFA', 'Claudin5'}, 'VariableTypes', {'string', 'single', 'single',  'single',  'single'});

% Plot density curves on the same axes

for i = 1:length(myfiles)
    % Read the image
    load(fullfile(dt,myfiles(i).name));
    
   D =  max(DAPI(:));
   N = max(NeuN(:));
   W = max(WFA(:));
   C = max(Claudin5(:));

   temp = table(myfiles(i).name(1:end-4),D,N,W,C, 'VariableNames', {'Sample', 'DAPI', 'NeuN', 'WFA', 'Claudin5'});
   tb = [tb; temp];

end
    % Compute the histogram
    [Wcounts, Wbi] = hist(WFA(:),256);
    [Ccounts, Cbi] = hist(Claudin5(:),256);
    [Dcounts, Dbi] = hist(DAPI(:),256);
    [Ncounts, Nbi] = hist(NeuN(:),256);
    
    % Plot the density curve with the specified color
    % hold on;
    figure('visible', 'off')
    ax1 = subplot(2,2,1);
    plot(Dbi, Dcounts)
    ax2 = subplot(2,2,2);
    plot(Nbi,Ncounts)
    ax3 = subplot(2,2,3);
    plot(Wbi, Wcounts)
    ax4 = subplot(2,2,4);
    plot(Cbi,Ccounts) 
    saveas(gcf,fullfile(ot,[myfiles(i).name, '.png']))
    close all
end

