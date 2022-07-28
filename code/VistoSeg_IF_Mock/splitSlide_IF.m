function splitSlide_IF(fname)

img = load(fname);
O = fieldnames(img);
%N = 4; %number of capture areas
disp(['The IF image has ',num2str(numel(O)),' channels'])

[y,~,~] = size(img.(O{1}));

tic
disp('Splitting whole slide into individual capture areas')

for i = 1:numel(O)
    tic
    I1.(O{i}) = img.(O{i})(1:round(y/2),:);
    I2.(O{i}) = img.(O{i})(round(y/2)+1:round(y/2)*2,:);
    toc
end

save([fname(1:end-4),'_A1.mat'],'-struct','I1', '-v7.3');
save([fname(1:end-4),'_B1.mat'],'-struct','I2', '-v7.3');

imwrite(mat2gray(I1.(O{1})),[fname(1:end-4),'_A1.tif'])
imwrite(mat2gray(I2.(O{1})),[fname(1:end-4),'_B1.tif'])

for i = 2:numel(O)
imwrite(mat2gray(I1.(O{i})),[fname(1:end-4),'_A1.tif'],'writemode', 'append')
imwrite(mat2gray(I2.(O{i})),[fname(1:end-4),'_B1.tif'],'writemode', 'append')

end
