function [count,prop,countC, inten] = countSpots_centroid(BW, img, R, tbl, posPath)
  
count = [];
prop = [];
countC = [];
inten = [];

O = fieldnames(BW);

    nSpots = size(tbl, 1);
    disp([num2str(nSpots),' Visium spots detected'])
    
    disp('Building spot grid')
    crow = round(table2array(tbl(:, 5)));
    ccol = round(table2array(tbl(:, 6)));
    mask = zeros(size(BW.(O{1})));
    mask(sub2ind(size(mask),crow,ccol)) = 1;
    mask = bwdist(mask) <= R;
    mask = bwlabel(mask);
    
for C = 1:numel(O)
    [BW.(O{C}),numROI] = bwlabel(BW.(O{C}));
	disp([num2str(numROI),' ', O{C},' ROIS detected'])
    count.(O{C}) = zeros(nSpots, 1);
    prop.(O{C}) = zeros(nSpots, 1);
    countC.(O{C}) = zeros(nSpots, 1);
    inten.(O{C}) = zeros(nSpots, 1);
end

tic
disp('counting dots per Visium spot')

for C = 1:numel(O)
points = struct2table(regionprops(BW.(O{C}), 'Centroid')).Centroid;
    for i = 1:nSpots 
        idx = mask(crow(i), ccol(i));
        spot = regionprops(mask==idx);
        signal = struct2table(regionprops(mask==idx & BW.(O{C})>0, img.(O{C}), 'Area', 'MeanIntensity'));
        isincircle = sum((points - [ccol(i) crow(i)]).^2,2)<= R^2;
        %check
%         [tempx,tempy] = find(mask == idx);
%         temp = BW.(O{C})(min(tempx):max(tempx),min(tempy):max(tempy));
%         imshow(temp)
%         viscircles(size(temp)/2, repelem(R, 1), 'Color', 'r');
        count.(O{C})(i) = length([signal.Area]);
        prop.(O{C})(i) = sum([signal.Area])/spot.Area;
        inten.(O{C})(i) = mean([signal.MeanIntensity]);
        countC.(O{C})(i) = length(find(isincircle));
        if mod(i,100) == 0
        disp([num2str(i),' spots finished in time ', num2str(toc),'s'])
        end

    end
    disp([num2str(C),'Channel done', num2str(toc),'s'])  
end

for C = 1:numel(O)
 temp = [count.(O{C}) prop.(O{C}) inten.(O{C}) countC.(O{C})];  
 tbl = [tbl array2table(temp, 'VariableNames', {['N',O{C}],['P',O{C}],['I',O{C}],['CN',O{C}]})];        
end

if ~exist(posPath, 'dir')
   mkdir(posPath);
end
    
    disp('writing table')
    writetable(tbl, fullfile(posPath, 'tissue_spot_counts1.csv'), 'Delimiter', ',');

end
