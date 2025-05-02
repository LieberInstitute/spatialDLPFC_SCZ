% code to convert the stitched mat file (o/p of inform stitch) to tiff file
% the output will have just 1 channel because of memory issues 
% output tiff shape: 73988x18600x1

fname = '/dcs04/lieber/marmaypag/spatialDLPFC_SCZ_LIBD4100/raw-data/images/RealPNN/Inform_0%_overlap/V12D07-334.mat';
img = load(fname);
O = fieldnames(img);

% Save the first channel, change the index (3rd dimension) to 1
channelToSave = 1; % saving just the AF channel

% Get the image data for the specified channel
imageData = img.(O{1})(:,:,channelToSave);

% Rescale the data to the range [0, 255]
maxValue = max(imageData(:));
minValue = min(imageData(:));
scaledData = uint8((imageData - minValue) / (maxValue - minValue) * 255);

% Create a TIFF file and save the image data
outputFilename = '/dcs04/lieber/marmaypag/spatialDLPFC_SCZ_LIBD4100/code/RealPNN/Troubleshoot/stitched_slide_V12D07-334.tif';
imwrite(scaledData, outputFilename);
