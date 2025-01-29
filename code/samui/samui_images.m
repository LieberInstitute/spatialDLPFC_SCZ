rawImgs = '/dcs04/lieber/marmaypag/spatialDLPFC_SCZ_LIBD4100/processed-data/VistoSeg/captureAreas/';
segImgs = '/dcs04/lieber/marmaypag/spatialDLPFC_SCZ_LIBD4100/processed-data/image_processing/DAPI_NeuN_WFA_Segs/';
outImgs = '/dcs04/lieber/marmaypag/spatialDLPFC_SCZ_LIBD4100/processed-data/samui/';
rawfiles = dir(fullfile(rawImgs,'/*1.mat'));

for i = 1:5
    fname = rawfiles(i).name(1:end-4);
    img = load(fullfile(rawfiles(i).folder, rawfiles(i).name));
    imwrite(mat2gray(img.DAPI), fullfile(outImgs, [fname, '.tif']), 'Compression', 'none')
    imwrite(mat2gray(img.NeuN), fullfile(outImgs, [fname, '.tif']), 'WriteMode','append', 'Compression', 'none')
    imwrite(mat2gray(img.WFA), fullfile(outImgs, [fname, '.tif']), 'WriteMode','append', 'Compression', 'none')
    imwrite(mat2gray(img.Claudin5), fullfile(outImgs, [fname, '.tif']),'WriteMode', 'append', 'Compression', 'none')
      disp('done img')
    segs = load(fullfile(segImgs, [fname, '_segs.mat']));
    imwrite(segs.DAPI, fullfile(outImgs, [fname, '.tif']), 'WriteMode','append', 'Compression', 'none')
    imwrite(segs.NeuN, fullfile(outImgs, [fname, '.tif']), 'WriteMode','append', 'Compression', 'none')
    imwrite(segs.WFA, fullfile(outImgs, [fname, '.tif']), 'WriteMode','append', 'Compression', 'none')
    imwrite(segs.Claudin5, fullfile(outImgs, [fname, '.tif']), 'WriteMode','append', 'Compression', 'none')
        disp('done segs')
end