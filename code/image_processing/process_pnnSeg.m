dt = '/dcs04/lieber/marmaypag/spatialDLPFC_SCZ_LIBD4100/processed-data/VistoSeg/captureAreas/';
dt1 = '/dcs04/lieber/marmaypag/spatialDLPFC_SCZ_LIBD4100/processed-data/image_processing/WFAseg/';
WFAfiles = dir([dt,'*.mat']);
%fname = 'V12F14-053_A1';
%fname = 'V13F27-294_C1';
%fname = 'V13F27-293_A1';


for i = 1:numel(WFAfiles)
    fname1=fullfile(WFAfiles(i).folder, WFAfiles(i).name);
    load(fname1, 'WFA')
    
    currentName = WFAfiles(i).name(1:end-4);

    fname=fullfile(dt1, [currentName, '_WFAseg.mat']);
    load(fname)

WFA1 = WFA.*(BW);
noise_idx = imbinarize(nonzeros(WFA1));
noise_idx = imbinarize(nonzeros(WFA1), 'adaptive');
%noise_idx = imbinarize(nonzeros(WFA1), 'adaptive', 'sensitivity', 0.5);
pnn_idx = find(BW);
pnn = BW;
pnn(pnn_idx(~noise_idx)) =0; 

%imshow(pnn)
%imshow(BW) 
save(fullfile(dt1,[currentName, '_WFAseg2.mat']),'pnn')
disp(i)
end
