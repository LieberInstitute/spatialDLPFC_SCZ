%% what rows to exclude from tissue_positions_list
load("Stitched_PNN_segmented_A1.mat")
size(WFA)
% ans =
% 
%        17450       16740

%% exclude rows from tissue_positions_list to make it compatible with countNuclei
load("Stitched_PNN_segmented_A1.mat")
posTable = readtable('/dcs04/lieber/marmaypag/spatialDLPFC_SCZ_LIBD4100/processed-data/PNN_segmented_VistoSeg/tissue_positions_list.csv');

 toDel = posTable.Var5<16740;
 posTable(toDel,:) = [];
 toDel = posTable.Var6<16740;
 posTable(toDel,:) = [];     
 max(posTable.Var5)          
writetable(posTable, '/dcs04/lieber/marmaypag/spatialDLPFC_SCZ_LIBD4100/processed-data/tissue_positions_list_modified.csv', 'Delimiter', ',');
