%filename =
%'/dcs04/lieber/marmaypag/spatialDLPFC_SCZ_LIBD4100/raw-data/Mock_tile_annotation_DLAlgorithm/20220712_VIF_MockPNN_Strong_NTC_Scan1/';
%P = 7;
temp = dir(filename);
myfiles = temp(5:end);

disp('Extracting coordinates')
tic
loc = cellfun(@(x) strsplit(x,'_'), {myfiles.name}',  'UniformOutput', false);
temp = cellfun(@(x) strsplit(x{P},','), loc,  'UniformOutput', false);
X = cellfun(@(x) str2double(x{1}(2:end)), temp);
X = sort(unique(X));
Y = cellfun(@(x) str2double(x{2}(1:end-1)), temp);
Y = sort(unique(Y));
