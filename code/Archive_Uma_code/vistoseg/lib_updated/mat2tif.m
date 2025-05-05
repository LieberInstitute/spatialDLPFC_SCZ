function mat2tif(matFileName, tifFileName)

% matFileName = '/dcs04/lieber/marmaypag/Visium_troubleshooting_LIBD001/VSPG/20230501_VSPG_HPC_Real_CNG_S2/InForm_0%_overlap_A1/V12J03-091_A1.mat'
% tiffFileName = '/dcs04/lieber/marmaypag/Visium_troubleshooting_LIBD001/VSPG/20230501_VSPG_HPC_Real_CNG_S2/InForm_0%_overlap_A1/V12J03-091_A1.tif'

fname = 'path to matfile';
img = load(fname);
O = fieldnames(img);

% imwrite(mat2gray(img.(O{1})),[fname(1:end-4),'.tif'])

for i = 2:numel(O)
imwrite(mat2gray(img.(O{i})),[fname(1:end-4),'_A1.tif'],'writemode', 'append')
end