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

