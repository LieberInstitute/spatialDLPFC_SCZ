setwd('/dcs04/lieber/marmaypag/spatialDLPFC_SCZ_LIBD4100/')
# Load Packages ----
suppressPackageStartupMessages({
  library(tidyverse)
  library(spatialLIBD)
  library(limma)
  library(sessioninfo)
  library(here)
  library(ggrepel)
})

## Load SPG spe object ----
raw_spe <- readRDS(here("processed-data/rds/01_build_spe","fnl_spe_kept_spots_only.rds"))
colData(raw_spe) |> names()

## Subset to neuropil+ only object ----
spg_spe <- raw_spe[, raw_spe$neuropil_pos == TRUE]
# dim: 36601 147965 

spg_spe <- raw_spe[, raw_spe$neuropil_pos &!raw_spe$neun_pos &!raw_spe$pnn_pos &!raw_spe$vasc_pos]
# dim: 36601 124725

# Create pseudobulked data ----
sce_pseudo <-
  registration_pseudobulk(
    spg_spe,
    var_registration = "fnl_spd",
    var_sample_id = "sample_id",
    covars = c("dx", "age", "sex", "lot_num", "slide_id"),
    min_ncells = 10)
	
#2026-09-29 12:42:09.431836 make pseudobulk object
#2026-09-29 12:42:18.720315 dropping 8 pseudo-bulked samples that are below 'min_ncells'.
#2026-09-29 12:42:18.903045 drop lowly expressed genes
#2026-09-29 12:42:19.484341 normalize expression

pb_neuropil <- sce_pseudo[,sce_pseudo$registration_variable %in%sprintf("spd%02d", c(2, 3, 5, 6, 7))]
ncol(pb_neuropil)
#[1] 314

dx_mod <-registration_model(pb_neuropil,covars = c("age", "sex"), var_registration = "dx")
dx_block_cor <- registration_block_cor(pb_neuropil, registration_model = dx_mod, var_sample_id = "sample_id")
dx_res <- registration_stats_enrichment(pb_neuropil,
  block_cor = dx_block_cor, 
  covars = c("age", "sex", "PRECAST_07"),
  var_registration = "dx",
  gene_ensembl = "gene_id",
  gene_name = "gene_name"
)

dx_res |> arrange(fdr_scz) |> head(10)
dx_res |> arrange(t_stat_scz) |> head(10)
dx_res |> arrange(desc(t_stat_scz)) |> head(10)
dx_res |> arrange(desc(logFC_scz)) |> head(10)

write_csv(dx_res |> arrange(fdr_scz),here("processed-data/image_processing/revisions","neuropil-dx_DEG-GM.csv"))

dx_res_old <- read_csv(here("code/analysis/dx_deg_spg_neuropil","neuropil-dx_DEG-GM.csv"), show_col_types = FALSE)

#sum(dx_res$fdr_scz < 0.05)
#[1] 86
#sum(dx_res_old$fdr_scz < 0.05)
#[1] 145
#
#sum(dx_res$fdr_scz < 0.1)
#[1] 199
#sum(dx_res_old$fdr_scz < 0.1)
#[1] 352

# Significant SCZ genes in each result
sig_new <- dx_res |>dplyr::filter(fdr_scz < 0.05) |>dplyr::pull(gene)
sig_new <- dx_res |>dplyr::filter(fdr_scz < 0.1) |>dplyr::pull(gene)
sig_old <- dx_res_old |>dplyr::filter(fdr_scz < 0.05) |>dplyr::pull(gene)

# Overlapping genes
overlap <- intersect(sig_new, sig_old)

# Number overlapping
length(overlap)
#[1] 79 #[1] 126

# Genes unique to dx_res
only_new <- setdiff(sig_new, sig_old)

# Genes unique to dx_res_old
only_old <- setdiff(sig_old, sig_new)

length(only_new)
length(only_old)
# Only in dx_res
dx_res |> dplyr::filter(gene %in% only_new) |> dplyr::select(ensembl, gene, fdr_scz, logFC_scz)
                        ensembl   gene    fdr_scz  logFC_scz
ENSG00000189241 ENSG00000189241 TSPYL1 0.04950127  0.1044001
ENSG00000146386 ENSG00000146386 ABRACL 0.04039126 -0.3084595
ENSG00000136205 ENSG00000136205   TNS3 0.04827705 -0.2252806
ENSG00000108826 ENSG00000108826 MRPL27 0.04352233 -0.1700578
ENSG00000079432 ENSG00000079432    CIC 0.03876601  0.1769005
ENSG00000124440 ENSG00000124440  HIF3A 0.04352233  0.4033422
ENSG00000154723 ENSG00000154723 ATP5PF 0.02204537 -0.2033518

# Only in dx_res_old
as.data.frame(dx_res_old) |> dplyr::filter(gene %in% only_old) |> dplyr::select(ensembl, gene, fdr_scz, logFC_scz)
           ensembl      gene    fdr_scz   logFC_scz
1  ENSG00000110876    SELPLG 0.01574031 -0.60340864
2  ENSG00000122873     CISD1 0.02197387 -0.14821489
3  ENSG00000058668    ATP2B4 0.02422767  0.13522670
4  ENSG00000061938      TNK2 0.02607730  0.17482156
5  ENSG00000166428      PLD4 0.02607730 -0.58534613
6  ENSG00000178209      PLEC 0.02695247  0.15662281
7  ENSG00000261325 LINC02192 0.02695247 -0.50666130
8  ENSG00000184515      BEX5 0.02695247 -0.14093938
9  ENSG00000096060     FKBP5 0.02919621  0.39654282
10 ENSG00000162694     EXTL2 0.03117328 -0.18904735
11 ENSG00000198168      SVIP 0.03117328  0.13086871
12 ENSG00000004779   NDUFAB1 0.03117328 -0.22713346
13 ENSG00000130770   ATP5IF1 0.03194747 -0.20201829
14 ENSG00000111913    RIPOR2 0.03205030  0.21713413
15 ENSG00000119718    EIF2B2 0.03293845 -0.34660509
16 ENSG00000153012      LGI2 0.03312187 -0.33800351
17 ENSG00000256193 LINC00507 0.03312187 -0.37110544
18 ENSG00000112695    COX7A2 0.03314551 -0.23209958
19 ENSG00000149196   HIKESHI 0.03314551 -0.21746668
20 ENSG00000141504      SAT2 0.03314551 -0.15451513
21 ENSG00000154582      ELOC 0.03390761 -0.16928283
22 ENSG00000174886   NDUFA11 0.03631277 -0.16656474
23 ENSG00000167863    ATP5PD 0.03652187 -0.19286133
24 ENSG00000149131  SERPING1 0.03834174  0.35913584
25 ENSG00000106028     SSBP1 0.03962278 -0.17777925
26 ENSG00000104904      OAZ1 0.03962278 -0.13287256
27 ENSG00000169991     IFFO2 0.03991267  0.30919373
28 ENSG00000158828     PINK1 0.03991267 -0.10940919
29 ENSG00000204472      AIF1 0.03991267 -0.53181349
30 ENSG00000179091      CYC1 0.03991267 -0.19591849
31 ENSG00000126432     PRDX5 0.03991267 -0.15405501
32 ENSG00000159199   ATP5MC1 0.03991267 -0.18737962
33 ENSG00000178982     EIF3K 0.03991267 -0.16823664
34 ENSG00000280071    GATD3B 0.03991267 -0.41833865
35 ENSG00000111652    COPS7A 0.04076500 -0.14924565
36 ENSG00000160209      PDXK 0.04141464  0.10532424
37 ENSG00000143420      ENSA 0.04352098 -0.09783222
38 ENSG00000237651   C2orf74 0.04352098 -0.20226138
39 ENSG00000154518   ATP5MC3 0.04352098 -0.20957710
40 ENSG00000162188      GNG3 0.04352098 -0.15183646
41 ENSG00000100804     PSMB5 0.04352098 -0.19178693
42 ENSG00000198408       OGA 0.04418034  0.09865297
43 ENSG00000182175      RGMA 0.04418034  0.21349199
44 ENSG00000184602       SNN 0.04418034 -0.13692908
45 ENSG00000108518      PFN1 0.04418034 -0.17454953
46 ENSG00000134030      CTIF 0.04418034  0.15511317
47 ENSG00000181029   TRAPPC5 0.04418034 -0.13048372
48 ENSG00000072954   TMEM38A 0.04418034 -0.16209097
49 ENSG00000205542    TMSB4X 0.04418034 -0.18298762
50 ENSG00000177606       JUN 0.04463076  0.21005151
51 ENSG00000183648    NDUFB1 0.04470150 -0.20818864
52 ENSG00000187193      MT1X 0.04470150  0.46050145
53 ENSG00000134248   LAMTOR5 0.04586871 -0.18275589
54 ENSG00000196104    SPOCK3 0.04586871 -0.20162633
55 ENSG00000172354      GNB2 0.04635766  0.10743959
56 ENSG00000103035     PSMD7 0.04769753 -0.10501327
57 ENSG00000152234   ATP5F1A 0.04784107 -0.13749101
58 ENSG00000176907      TCIM 0.04786324 -0.41438813
59 ENSG00000177556     ATOX1 0.04803953 -0.18555480
60 ENSG00000111775    COX6A1 0.04836189 -0.19297669
61 ENSG00000101132     PFDN4 0.04836189 -0.15887835
62 ENSG00000164032     H2AFZ 0.04849556 -0.19116748
63 ENSG00000228716      DHFR 0.04849556 -0.28908813
64 ENSG00000119705     SLIRP 0.04864584 -0.21605767
65 ENSG00000196465     MYL6B 0.04955784 -0.14833001
66 ENSG00000184752   NDUFA12 0.04993860 -0.23649474

##### with new_sig<0.1
           ensembl      gene    fdr_scz  logFC_scz
1  ENSG00000110876    SELPLG 0.01574031 -0.6034086
2  ENSG00000166428      PLD4 0.02607730 -0.5853461
3  ENSG00000261325 LINC02192 0.02695247 -0.5066613
4  ENSG00000162694     EXTL2 0.03117328 -0.1890474
5  ENSG00000153012      LGI2 0.03312187 -0.3380035
6  ENSG00000141504      SAT2 0.03314551 -0.1545151
7  ENSG00000106028     SSBP1 0.03962278 -0.1777793
8  ENSG00000126432     PRDX5 0.03991267 -0.1540550
9  ENSG00000178982     EIF3K 0.03991267 -0.1682366
10 ENSG00000111652    COPS7A 0.04076500 -0.1492456
11 ENSG00000182175      RGMA 0.04418034  0.2134920
12 ENSG00000072954   TMEM38A 0.04418034 -0.1620910
13 ENSG00000196104    SPOCK3 0.04586871 -0.2016263
14 ENSG00000172354      GNB2 0.04635766  0.1074396
15 ENSG00000103035     PSMD7 0.04769753 -0.1050133
16 ENSG00000176907      TCIM 0.04786324 -0.4143881
17 ENSG00000177556     ATOX1 0.04803953 -0.1855548
18 ENSG00000164032     H2AFZ 0.04849556 -0.1911675
19 ENSG00000228716      DHFR 0.04849556 -0.2890881


dx_res <- dx_res1 |>
  mutate(
    sig_05 = fdr_scz <= 0.05,
    sig_10 = fdr_scz <= 0.10,
    gene_cat = case_when(
      sig_10 == FALSE ~ "grey",
      sig_10 == TRUE & logFC_scz >= 0 ~ "red",
      sig_10 == TRUE & logFC_scz < 0 ~ "blue"
    )
  ) |>
  arrange(
    desc(p_value_scz)
  )

# -----------------------------
# Settings
# -----------------------------
p_cutoff  <- 0.05
lfc_cutoff <- 0.15

genes_to_label <- c(
  "SST", "CORT", "VGF", "ALDH1A1", "PVALB",
  "CCK", "SNAP25", "NRN1", "GAD2", "BDNF",
  "GAD1", "NPY", "GAP43", "CRHBP", "GFAP",
  "NRXN2", "NTRK2", "CAMK2A", "CACNA1B"
)

# -----------------------------
# Prepare data
# -----------------------------
plot_df <- dx_res %>%
  mutate(
    direction = case_when(
      logFC_scz > 0 ~ "Up",
      logFC_scz < 0 ~ "Down",
      TRUE ~ "Neutral"
    ),

    # intensity based on p-value
    alpha_val = case_when(
      p_value_scz <= p_cutoff ~ 1,
      TRUE ~ 0.25
    )
  )

label_df <- plot_df %>%
  filter(gene %in% genes_to_label)


# -----------------------------
# Volcano plot
# -----------------------------
p <- ggplot(
  plot_df,
  aes(
    x = logFC_scz,
    y = -log10(p_value_scz)
  )
) +

  # p-value threshold
  geom_hline(
    yintercept = -log10(p_cutoff),
    color = "grey50",
    linewidth = 0.6
  ) +

  # fold-change thresholds
  geom_vline(
    xintercept = c(-lfc_cutoff, lfc_cutoff),
    linetype = "dashed",
    color = "grey50",
    linewidth = 0.6
  ) +

  # background / nonsignificant genes
  geom_point(
    data = filter(plot_df, p_value_scz > p_cutoff),
    color = "grey70",
    size = 1.6,
    alpha = 0.45
  ) +

  # downregulated
  geom_point(
    data = filter(
      plot_df,
      p_value_scz <= p_cutoff,
      logFC_scz < 0
    ),
    color = "blue",
    size = 1.7,
    alpha = 0.8
  ) +

  # upregulated
  geom_point(
    data = filter(
      plot_df,
      p_value_scz <= p_cutoff,
      logFC_scz > 0
    ),
    color = "red",
    size = 1.7,
    alpha = 0.8
  ) +

  # gene labels
  geom_label_repel(
    data = label_df,
    aes(
      label = gene,
      color = direction
    ),
    size = 3.2,
    fill = "white",
    box.padding = 0.25,
    point.padding = 0.2,
    segment.color = "grey40",
    segment.size = 0.4,
    min.segment.length = 0,
    max.overlaps = Inf,
    show.legend = FALSE
  ) +

  scale_color_manual(
    values = c(
      "Down" = "blue",
      "Up" = "red",
      "Neutral" = "grey50"
    )
  ) +

  labs(
    x = "log2(Fold Change)",
    y = expression(-log[10](p))
  ) +

  theme_classic(base_size = 14) +

  theme(
    legend.position = "none",

    axis.title = element_text(size = 14),
    axis.text = element_text(size = 12),

    axis.line = element_line(
      color = "black",
      linewidth = 0.7
    ),

    axis.ticks = element_line(
      color = "black"
    ),

    panel.border = element_rect(
      color = "black",
      fill = NA,
      linewidth = 0.7
    )
  )

p