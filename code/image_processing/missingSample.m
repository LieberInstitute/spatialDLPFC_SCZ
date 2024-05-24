dt = '/dcs04/lieber/marmaypag/spatialDLPFC_SCZ_LIBD4100/processed-data/VistoSeg/captureAreas/';
myfiles = dir([dt,'*1.mat']);
filenames = {myfiles.name}';

split_filenames = arrayfun(@(x) strsplit(filenames{x}, '_'), 1:numel(filenames), 'UniformOutput', false);
% Assuming split_filenames contains the split filenames
slides = cellfun(@(x) x{1}, split_filenames, 'UniformOutput', false); % Extracting the second cell elements
arrays = cellfun(@(x) x{2}, split_filenames, 'UniformOutput', false); % Extracting the second cell elements

% Count occurrences of each number
[gc,grps] = groupcounts(slides');
disp(grps{gc == min(gc)})

[gc1,grps1] = groupcounts(arrays');
disp(grps1{gc1 == min(gc1)})

%% create mat file for missing sample
fpath = '/dcs04/lieber/marmaypag/spatialDLPFC_SCZ_LIBD4100/processed-data/VistoSeg/captureAreas/';
fname = 'V13M06-280_A1';

WFA = im2double(imread(fullfile(fpath,[fname,'.tif']), 5));
AF = im2double(imread(fullfile(fpath,[fname,'.tif']), 1));
DAPI = im2double(imread(fullfile(fpath,[fname,'.tif']), 3));
NeuN = im2double(imread(fullfile(fpath,[fname,'.tif']), 4));
Claudin5 = im2double(imread(fullfile(fpath,[fname,'.tif']), 2));

save(fullfile(fpath,[fname,'.mat']), 'AF', 'Claudin5', 'DAPI', 'NeuN', 'WFA', '-v7.3')

