function tif2mat(fname) %function to convert segmented tif file to mat format
%fname is the path and name of the segmented tif file
img = load(fname)
O = fieldnames(img);
%imwrite(mat2gray(img.(O{1})),[fname(1:end-4),'.tif'])
for i = 1:numel(O)
imwrite(mat2gray(img.(O{i})),[fname(1:end-4),'.tif'],'writemode', 'append')
end