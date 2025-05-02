ntcR = '/dcs04/lieber/marmaypag/spatialDLPFC_SCZ_LIBD4100/processed-data/VistoSeg/captureAreas/V13M06-342_D1.mat';
sczR = '/dcs04/lieber/marmaypag/spatialDLPFC_SCZ_LIBD4100/processed-data/VistoSeg/captureAreas/V13M06-343_D1.mat';

load(ntcR,'WFA');
ntc = WFA;
load(sczR,'WFA');
scz = WFA;

imwrite(ntc, '/dcs04/lieber/marmaypag/spatialDLPFC_SCZ_LIBD4100/plots/image_processing/sfig_workflow/ntc_raw.png')
imwrite(ntc(9000:10000, 1000:2000), '/dcs04/lieber/marmaypag/spatialDLPFC_SCZ_LIBD4100/plots/image_processing/sfig_workflow/ntc_rawInset.png')
imwrite(scz, '/dcs04/lieber/marmaypag/spatialDLPFC_SCZ_LIBD4100/plots/image_processing/sfig_workflow/scz_raw.png')
imwrite(scz(9000:10000, 8000:9000), '/dcs04/lieber/marmaypag/spatialDLPFC_SCZ_LIBD4100/plots/image_processing/sfig_workflow/scz_rawInset.png')

%% segs
ntcS = '/dcs04/lieber/marmaypag/spatialDLPFC_SCZ_LIBD4100/processed-data/image_processing/DAPI_NeuN_WFA_Segs/V13M06-342_D1_segs.mat';
sczS = '/dcs04/lieber/marmaypag/spatialDLPFC_SCZ_LIBD4100/processed-data/image_processing/DAPI_NeuN_WFA_Segs/V13M06-343_D1_segs.mat';

load(ntcS,'WFA');
ntcS = WFA;
load(sczS,'WFA');
sczS = WFA;

imwrite(ntcS, '/dcs04/lieber/marmaypag/spatialDLPFC_SCZ_LIBD4100/plots/image_processing/sfig_workflow/ntc_seg.png')
imwrite(ntcS(9000:10000, 1000:2000), '/dcs04/lieber/marmaypag/spatialDLPFC_SCZ_LIBD4100/plots/image_processing/sfig_workflow/ntc_segInset.png')
imwrite(sczS, '/dcs04/lieber/marmaypag/spatialDLPFC_SCZ_LIBD4100/plots/image_processing/sfig_workflow/scz_seg.png')
imwrite(sczS(9000:10000, 8000:9000), '/dcs04/lieber/marmaypag/spatialDLPFC_SCZ_LIBD4100/plots/image_processing/sfig_workflow/scz_segInset.png')

%% pre filt
ntcpS = '/dcs04/lieber/marmaypag/spatialDLPFC_SCZ_LIBD4100/processed-data/image_processing/WFAseg/V13M06-342_D1_WFAseg.mat';
sczpS = '/dcs04/lieber/marmaypag/spatialDLPFC_SCZ_LIBD4100/processed-data/image_processing/WFAseg/V13M06-343_D1_WFAseg.mat';

load(ntcpS);
ntcp =BW;
load(sczpS);
sczp =BW;

imwrite(ntcp, '/dcs04/lieber/marmaypag/spatialDLPFC_SCZ_LIBD4100/plots/image_processing/sfig_workflow/ntc_segP.png')
imwrite(ntcp(9000:10000, 1000:2000), '/dcs04/lieber/marmaypag/spatialDLPFC_SCZ_LIBD4100/plots/image_processing/sfig_workflow/ntc_segPInset.png')
imwrite(sczp, '/dcs04/lieber/marmaypag/spatialDLPFC_SCZ_LIBD4100/plots/image_processing/sfig_workflow/scz_segP.png')
imwrite(sczp(9000:10000, 8000:9000), '/dcs04/lieber/marmaypag/spatialDLPFC_SCZ_LIBD4100/plots/image_processing/sfig_workflow/scz_segPInset.png')

%% histograms
addpath(genpath('/dcs04/lieber/marmaypag/spatialDLPFC_SCZ_LIBD4100/code/image_processing/'))
[counts, x]=hist(ntc(:),256);
[thresh, lnP] = triangle_threshold_right_tail(counts(13:end));
level = x(13)+x(thresh);

 plot(x(13:end),counts(13:end))
 %plot(x,counts)
 %set(gca, 'YScale', 'log');
 hold on
% plot(x(lnP(1:2)), lnP(3:4), 'Color', 'r')
 plot([0.1,0.8], lnP(3:4), 'Color', 'r')
 xline(level,'g')
 hold off
saveas(gcf, '/dcs04/lieber/marmaypag/spatialDLPFC_SCZ_LIBD4100/plots/image_processing/sfig_workflow/ntcHist.png');

[counts, x]=hist(scz(:),256);
[thresh, lnP] = triangle_threshold_right_tail(counts(13:end));
level = x(13)+x(thresh);

 plot(x(13:end),counts(13:end))
 %plot(x,counts)
 %set(gca, 'YScale', 'log');
 hold on
% plot(x(lnP(1:2)), lnP(3:4), 'Color', 'r')
 plot([0.08,0.8], lnP(3:4), 'Color', 'r')
 xline(level,'g')
 hold off
saveas(gcf, '/dcs04/lieber/marmaypag/spatialDLPFC_SCZ_LIBD4100/plots/image_processing/sfig_workflow/sczHist.png');
