 %mask = '/dcl02/lieber/ajaffe/SpatialTranscriptomics/LIBD/spatialDLPFC/Images/Liebert_Institute_OTS-20-7748_rush_posterior_2_nuclei.mat';
%jsonname = '/dcl02/lieber/ajaffe/SpatialTranscriptomics/LIBD/spatialDLPFC/outputs/NextSeq/DLPFC_Br3942_post_manual_alignment/outs/spatial/scalefactors_json.json';
%posname = '/dcl02/lieber/ajaffe/SpatialTranscriptomics/LIBD/spatialDLPFC/outputs/NextSeq/DLPFC_Br3942_post_manual_alignment/outs/spatial/tissue_positions_list.csv';

function [count,prop,countC,inten] = countNuclei(mask,img,jsonname,posname,~) 

disp('loading data')
 %img = load(img);
 %BW = load(mask);
 
img = load(img, 'DAPI');
BW = load(mask, 'DAPI');

O = fieldnames(BW);
   
[posPath,~] = fileparts(posname);

w = jsondecode(fileread(jsonname));
R = ceil(w.spot_diameter_fullres/2);
tbl = readtable(posname) ;
count = [];
prop = [];

  % if size(tbl, 2) > 6
  %     b = 7;
  %     for C = 1:numel(O)
  %     count.(O{C}) = table2array(tbl(:, b));
  %     prop.(O{C}) = table2array(tbl(:, b+1));
  %     countC.(O{C}) = table2array(tbl(:, b+2));
  %     inten.(O{C}) = table2array(tbl(:, b+3));
  %     b = b+4;
  %     end
  %	if size(tbl, 2) > 18
  %			disp("file error")
  %  else
       tbl.Properties.VariableNames = {'barcode','tissue','row','col','imagerow','imagecol'};
	   %tbl.Properties.VariableNames = {'barcode','tissue','row','col','imagerow','imagecol','NDAPI','PDAPI','IDAPI','CNDAPI','NNeuN','PNeuN','INeuN','CNNeuN','NWFA','PWFA','IWFA','CNWFA'}
       [count,prop,countC,inten] = countSpots_centroid(BW, img, R, tbl, posPath);
        
  %  end
end

