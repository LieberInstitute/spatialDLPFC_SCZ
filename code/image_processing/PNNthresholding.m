% net = denoisingNetwork('DnCNN'); %Load the pretrained denoising convolutional neural network, 'DnCNN'.
% dWFA = denoiseImage(WFA,net); %denoise the image

dt = '/dcs04/lieber/marmaypag/spatialDLPFC_SCZ_LIBD4100/processed-data/VistoSeg/captureAreas/';
myfiles = dir([dt,'*1.mat']);

for i = 1:length(myfiles)
    % Read the image
    load(fullfile(dt,myfiles(i).name),'WFA');

[lehisto, x]=imhist(WFA);
level=triangle_th(lehisto,x);
BW=imbinarize(WFA,level);
end