dt = '/dcs04/lieber/marmaypag/spatialDLPFC_SCZ_LIBD4100/processed-data/image_processing/WFAseg/';
dt1 = '/dcs04/lieber/marmaypag/spatialDLPFC_SCZ_LIBD4100/processed-data/image_processing/AFseg/';
dt2 = '/dcs04/lieber/marmaypag/spatialDLPFC_SCZ_LIBD4100/processed-data/VistoSeg/captureAreas/';


ot = '/dcs04/lieber/marmaypag/spatialDLPFC_SCZ_LIBD4100/plots/image_processing/AFmasking/';
myfiles = dir([dt,'*3.mat']);
addpath(genpath('/dcs04/lieber/marmaypag/spatialDLPFC_SCZ_LIBD4100/code/image_processing'))


%sum(Area)<300000 and intensity<2.5
% V12D07-334_C1.png
% V13M06-281_D1.png
% V13M06-343_A1.png
% V12F14-057_A1.png
% V13F27-294_B1.png
% V13M06-340_C1.png

%remove max(Area) and intensity<1
%V13M06-340_D1
%V13M06-282_C1
%V13M06-280_D1
%V13M06-279_D1
%V13M06-279_C1
%V13M06-279_B1
%V13F27-336_A1
%V13F27-296_B1
%final_table = table([], [], [], [], 'VariableNames', {'MeanIntensity', 'MaxIntensity', 'Area', 'Sample'});

for i = 1:64
load(fullfile(dt,myfiles(i).name));
load(fullfile(dt1,[myfiles(i).name(1:end-11), 'AFseg.mat']));
load(fullfile(dt2,[myfiles(i).name(1:end-12), '.mat']), 'WFA');

pnnM = pnn;
pnnM(BW) = 0;

stats = struct2table(regionprops(pnnM, WFA, 'MeanIntensity', 'MaxIntensity', 'Area'));
stats.Sample = repmat({myfiles(i).name(1:end-12)}, height(stats), 1);
if max([stats.MeanIntensity]) < 1
filteredPnnM = bwpropfilt(pnnM, "Area", [0 30000]);
else
filteredPnnM = pnnM;
end
%imshow([pnn, ones(size(pnn,1),50), pnnM, ones(size(pnn,1),50), filteredPnnM])
final = [pnn, ones(size(pnn,1),50), pnnM, ones(size(pnn,1),50), filteredPnnM];
imwrite(imresize(final, 0.5),fullfile(ot,[myfiles(i).name(1:end-11),'AFfilt.png']));
save(fullfile(dt,[myfiles(i).name(1:end-11),'AFseg_mask.mat']),'filteredPnnM')

disp(myfiles(i).name)
end

%writetable(final_table, fullfile(dt1,'data.csv'));
