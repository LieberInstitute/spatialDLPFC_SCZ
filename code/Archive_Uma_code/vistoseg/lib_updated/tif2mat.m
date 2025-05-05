function tif2mat(fname)
%fname = 'path and name to your tif file';
DAPI = imread(fname,1);
Claudin5 = imread(fname,2);
NeuN = imread(fname,3);
WFA = imread(fname,4);
AF = imread(fname,5);
save(fullfile(path you want to save in single quotes, filename ending with .mat in single quotes), 'DAPI', 'Claudin5','NeuN','WFA','AF', '-v7.3')