fname = 'path and name to your tif file';
DAPI = imread(fname,1);
Claudin5 = imread(fname,2);
NeuN = imread(fname,3);
WFA = imread(fname,4);
AF = imread(fname,5);
save(fullfile('/dcs04/lieber/marmaypag/spatialDLPFC_SCZ_LIBD4100/processed-data/RealPNN/capture_area_segmentations/all_channels_segemented/Test/', 'V12F14-053_A1_thresholded.mat'), 'DAPI', 'Claudin5','NeuN','WFA','AF', '-v7.3')


% Specify the folder where your TIFF images are located
tif_folder = '/dcs04/lieber/marmaypag/spatialDLPFC_SCZ_LIBD4100/processed-data/RealPNN/segmented_channels_stitched/slide16/';

% List all TIFF files in the folder
tif_files = dir(fullfile(tif_folder, '*.tif'));

% Loop through each TIFF file
for i = 1:length(tif_files)
    % Read the current TIFF file
    fname = fullfile(tif_folder, tif_files(i).name);
    % Read individual channels from the TIFF file
    DAPI = imread(fname, 1);
    Claudin5 = imread(fname, 2);
    NeuN = imread(fname, 3);
    WFA = imread(fname, 4);
    AF = imread(fname, 5);
    % Generate the dynamic part of the filename based on the loop index
    dynamic_part = char('A' + i - 1);
    % Create the full filename for saving the MAT file
    mat_filename = fullfile('/dcs04/lieber/marmaypag/spatialDLPFC_SCZ_LIBD4100/processed-data/RealPNN/segmented_channels_stitched/slide16/', ...
        ['V13F27-336_' dynamic_part '1_thresholded.mat']);
    % Save the data to MAT file
    save(mat_filename, 'DAPI', 'Claudin5', 'NeuN', 'WFA', 'AF', '-v7.3');
    % Display a message indicating the processing of the current file
    disp(['Processed: ' tif_files(i).name]);
end
