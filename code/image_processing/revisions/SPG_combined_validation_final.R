setwd("/dcs04/lieber/marmaypag/spatialDLPFC_SCZ_LIBD4100/")

suppressPackageStartupMessages({
    library(here)
    library(SpatialExperiment)
    library(SummarizedExperiment)
    library(Matrix)
    library(dplyr)
    library(tidyr)
    library(ggplot2)
	library(readxl)
	library(biomaRt)
	library(babelgene)
})


# ============================================================
# 1. Load data
# ============================================================

domain_levels <- c("spd07","spd06","spd02","spd05", "spd03", "spd01", "spd04")
domain_labels <- c("SpD07-L1/M", "SpD06-L2/3", "SpD02-L3/4", "SpD05-L5","SpD03-L6","SpD01-WMtz", "SpD04-WM")

#domain_levels <- c("spd04", "spd01", "spd03", "spd05","spd02", "spd06", "spd07")
#domain_labels <- c("SpD04-WM", "SpD01-WMtz", "SpD03-L6", "SpD05-L5", "SpD02-L3/4", "SpD06-L2/3", "SpD07-L1/M")
domain_colors <- c(
  "SpD07-L1/M" = "#FFA500",
  "SpD06-L2/3" = "#3388FF",
  "SpD02-L3/4" = "#E5E3E5",
  "SpD05-L5"   = "#40FF40",
  "SpD03-L6"   = "#FF2028",
  "SpD01-WMtz" = "#605860",
  "SpD04-WM"   = "#F000F0"
)

spe_ntc <- readRDS(here( "processed-data", "image_processing", "enrichment", "spe_ntc.rds"))
spe_ntc$slide_id <- sapply(strsplit(spe_ntc$sample_id, "_"), `[`, 1)
spe_ntc$brnum_dx <- paste(spe_ntc$brnum, toupper(spe_ntc$dx), sep = "_")

df <- as.data.frame(colData(spe_ntc))
df$lib_size <- as.numeric(Matrix::colSums(counts(spe_ntc)))

our_gene_symbols <- as.character(SummarizedExperiment::rowData(spe_ntc)$gene_name)

# ==========================================================
# functions
# ==========================================================
## function I

add_spg_gene_expression <- function(spg_name, gene_list, df, spe_ntc, threshold = "q75") {
  gene_list <- unique(toupper(trimws(as.character(gene_list))))
  gene_list <- gene_list[!is.na(gene_list) & gene_list != ""]
  if (length(gene_list) == 0) {stop("gene_list is empty.")}

  # Convert q75 -> 0.75, q85 -> 0.85, etc.
  threshold <- tolower(threshold)
  if (!grepl("^q[0-9]+$", threshold)) {
	  stop("threshold must be formatted as 'q75', 'q85', etc.")
  }
  quantile_prob <- as.numeric(sub("q", "", threshold)) / 100
  if (quantile_prob <= 0 || quantile_prob >= 1) {
    stop("Quantile must be between 0 and 100.")
  }

  prefix <- paste0(spg_name, "_top", length(gene_list))
  raw_col  <- paste0(prefix, "_raw_sum")
  norm_col <- paste0(prefix, "_norm_sum")
  high_col <- paste0(prefix, "_", threshold, "_high")

  # Match genes with spe_ntc
  matched_idx <- which(toupper(our_gene_symbols) %in% gene_list)
  matched_gene_ids <- rownames(spe_ntc)[matched_idx]
  matched_gene_symbols <- unique(toupper(our_gene_symbols[matched_idx]))
  missing_genes <- setdiff(gene_list, matched_gene_symbols)

  cat("\n============================\n")
  cat("SPG:", spg_name, "\n")
  cat("Gene set:", prefix, "\n")
  cat("Requested genes:", length(gene_list), "\n")
  cat("Matched genes:", length(matched_gene_symbols), "\n")
  cat("Missing genes:", length(missing_genes), "\n")

  if (length(missing_genes) > 0) {cat("Missing:", paste(missing_genes, collapse = ", "), "\n")}
  if (length(matched_gene_ids) == 0) {stop("No genes matched spe_ntc.")}
  if (nrow(df) != ncol(spe_ntc)) {stop("df rows do not match spe_ntc spots.")} # Check spot alignment
  if (!identical(rownames(df), colnames(spe_ntc))) {stop("df rows are not aligned with spe_ntc columns.")}

  # Calculate raw counts
  gene_counts <- SummarizedExperiment::assay(spe_ntc, "counts")[matched_gene_ids, , drop = FALSE]
  gene_raw_sum <- Matrix::colSums(gene_counts)

  # Calculate library size if missing
  if (!"lib_size" %in% colnames(df)) {
    df$lib_size <- as.numeric(Matrix::colSums(SummarizedExperiment::assay(spe_ntc, "counts")))
  }

  df[[raw_col]] <- as.numeric(gene_raw_sum)
  df[[norm_col]] <- ifelse(!is.na(df$lib_size) & df$lib_size > 0,df[[raw_col]] / df$lib_size,NA_real_)
  q_value <- as.numeric(quantile(df[[norm_col]],probs = quantile_prob,na.rm = TRUE))
  df[[high_col]] <- !is.na(df[[norm_col]]) &df[[norm_col]] > q_value

  cat("Threshold:", threshold, "\n")
  cat("Global quantile cutoff:", q_value, "\n")
  cat("High spots:", sum(df[[high_col]], na.rm = TRUE), "\n")

  cat("\nClassification:\n")
  print(table(df[[high_col]], useNA = "ifany"))

  cat("============================\n")
  return(df)
}

## function II

make_spg_prop_df <- function(df, spg_name, top_n, threshold = "q75",
                             domain_levels, domain_labels) {

  prefix <- paste0(spg_name, "_top", top_n)
  norm_col <- paste0(prefix, "_norm_sum")
  high_col <- paste0(prefix, "_", tolower(threshold), "_high")

  required_cols <- c("brnum_dx", "sample_id", "PRECAST_07",
                     norm_col, high_col)

  missing_cols <- setdiff(required_cols, colnames(df))

  if (length(missing_cols) > 0) {
    stop("Missing columns in df: ", paste(missing_cols, collapse = ", "))
  }

  plot_df <- df %>%
    dplyr::filter(
      !is.na(brnum_dx),
      !is.na(sample_id),
      PRECAST_07 %in% domain_levels,
      !is.na(.data[[norm_col]]),
      !is.na(.data[[high_col]])
    ) %>%
    dplyr::mutate(PRECAST_07 = as.character(PRECAST_07))

  # Total spots and high spots per sample
  sample_totals <- plot_df %>%
    dplyr::group_by(brnum_dx, sample_id) %>%
    dplyr::summarise(
      sample_spots = dplyr::n(),
      sample_high = sum(.data[[high_col]]),
      sample_high_prop = sample_high / sample_spots,
      .groups = "drop"
    )

  # Domain-level summaries per sample
  prop_df <- plot_df %>%
    dplyr::group_by(brnum_dx, sample_id, PRECAST_07) %>%
    dplyr::summarise(
      n_spots = dplyr::n(),
      n_high = sum(.data[[high_col]]),
      .groups = "drop"
    ) %>%
    dplyr::left_join(
      sample_totals,
      by = c("brnum_dx", "sample_id")
    ) %>%
    dplyr::mutate(
      prop_original = dplyr::if_else(
        sample_high > 0, n_high / sample_high, NA_real_
      ),
      prop_within_domain = n_high / n_spots,
      enrichment = dplyr::if_else(
        sample_high_prop > 0,
        prop_within_domain / sample_high_prop,
        NA_real_
      )
    )

  # Complete missing domains within existing sample combinations
  prop_df <- prop_df %>%
    tidyr::complete(
      tidyr::nesting(brnum_dx, sample_id),
      PRECAST_07 = domain_levels
    ) %>%
    dplyr::mutate(
      n_spots = tidyr::replace_na(n_spots, 0L),
      n_high = tidyr::replace_na(n_high, 0L),
      prop_original = dplyr::if_else(
        n_spots == 0, NA_real_, prop_original
      ),
      prop_within_domain = dplyr::if_else(
        n_spots == 0, NA_real_, prop_within_domain
      ),
      enrichment = dplyr::if_else(
        n_spots == 0, NA_real_, enrichment
      ),
      PRECAST_07 = factor(
        PRECAST_07,
        levels = domain_levels,
        labels = domain_labels
      )
    )

  return(prop_df)
}

## function III
make_heatmap <- function(data, value_col, title, fill_label, limits = NULL) {
  ggplot(data, aes(x = brnum_dx, y = PRECAST_07, fill = .data[[value_col]])) +
  geom_tile(color = "grey60", linewidth = 0.3) +
  scale_fill_gradient(low = "white", high = "black", limits = limits, oob = scales::squish,na.value = "red") +
  labs(title = title, x = "Sample", y = "Spatial domain", fill = fill_label) +
  theme_minimal(base_size = 12) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 7),
      	axis.text.y = element_text(size = 10),
      	panel.grid = element_blank(),
      	plot.title = element_text(face = "bold"))
}

# ============================================================
# PNN
# ============================================================
pnn <- read_excel(here("raw-data/images/SPG_Spot_Valid/PNN/Supple_Table_DataS4_PNN.xlsx"),sheet = "PNN Energy",na = c("", "NA"))
pnn <- as.data.frame(pnn)
pnn$mouse_gene <- trimws(as.character(pnn$gene_acronym))

#  Convert mouse genes to human orthologs
mouse_genes <- unique(na.omit(pnn$mouse_gene))
orthologs <- babelgene::orthologs(genes = mouse_genes,species = "mouse",human = FALSE)

head(orthologs)

# Keep mouse-human mapping
ortholog_map <- orthologs %>%
  transmute(mouse_gene = symbol, human_gene = human_symbol) %>%
  distinct()
pnn_human <- pnn %>% left_join(ortholog_map, by = "mouse_gene")
matched_pnn_genes <- intersect(unique(na.omit(pnn_human$human_gene)),rowData(spe_ntc)$gene_name)

cat("Original mouse genes:", n_distinct(pnn$mouse_gene), "\n")
cat("Mapped human genes:", n_distinct(pnn_human$human_gene, na.rm = TRUE), "\n")
cat("Matched in spe_ntc:", length(matched_pnn_genes), "\n")

#Original mouse genes: 17633 
#Mapped human genes: 13948 
#Matched in spe_ntc: 13601 

ortholog_map_unique <- ortholog_map %>%
    filter(!is.na(human_gene), human_gene != "") %>%
    group_by(mouse_gene) %>%
    filter(n_distinct(human_gene) == 1) %>%
    ungroup() %>%
    distinct(mouse_gene, human_gene)
	
pnn_human <- pnn %>%left_join(ortholog_map_unique, by = "mouse_gene", relationship = "many-to-one")
matched_pnn_genes <- intersect(unique(na.omit(pnn_human$human_gene)),rowData(spe_ntc)$gene_name)

cat("Original mouse genes:", n_distinct(pnn$mouse_gene), "\n")
cat("Mapped human genes:", n_distinct(pnn_human$human_gene, na.rm = TRUE), "\n")
cat("Matched in spe_ntc:", length(matched_pnn_genes), "\n")

#Original mouse genes: 17633 
#Mapped human genes: 13785 
#Matched in spe_ntc: 13477 

pnn_human %>% dplyr::filter(toupper(mouse_gene) == "PVALB" | toupper(human_gene) == "PVALB") %>% dplyr::select(mouse_gene, human_gene, spearman_rho, FDR, p_value)

pnn_ranked <- pnn_human %>% dplyr::filter(FDR < 0.01,p_value < 0.05,spearman_rho > 0,!is.na(human_gene),human_gene != "") %>%
 	 		  arrange(desc(spearman_rho), FDR) %>%
  			  distinct(human_gene, .keep_all = TRUE) %>%
  			  mutate(rank = row_number())
pnn_ranked %>% tibble::as_tibble() %>% print(n = 145, width = Inf)
pnn_ranked %>% dplyr::filter(human_gene == "PVALB") %>% dplyr::select(rank, mouse_gene, human_gene, spearman_rho, FDR, p_value)

pnn_top <- pnn_ranked %>% slice_head(n = 200)
#PNN_genes <- pnn_top$human_gene		  		  

# Generate all three heatmaps
# ==========================================================
df = add_spg_gene_expression("PNN", pnn_top$human_gene, df, spe_ntc, threshold = "q75")
prop_df = make_spg_prop_df(df,"PNN",200,threshold = "q75", domain_levels, domain_labels) 
prop_df_PNN = prop_df

#SPG: PNN 
#Gene set: PNN_top200 
#Requested genes: 200 
#Matched genes: 197 
#Missing genes: 3 
#Missing: ADSS1, LOC728743, LYZL6 

# ============================================================
# vasc
# ============================================================
vasc = read_excel(here("raw-data/images/SPG_Spot_Valid/Vasculature/Supple_Table_2_Vas.xlsx"), sheet="Post Mortem Vascular Subcluster")
vasc = as.data.frame(vasc)
colnames(vasc) <- c("Column1", colnames(vasc)[1:7])
gene_idx <- match(vasc$gene,rowData(spe_ntc)$gene_name)
# Check that all genes were found
missing_genes <- vasc$gene[is.na(gene_idx)]

dim(vasc)
dim(vasc)[1]-length(missing_genes)
#7055 of 7096 reported genes only found in our dataset

vasc_ranked <- vasc %>%
  dplyr::filter(p_val_adj < 0.05, avg_log2FC > 1) %>%
  dplyr::arrange(desc(avg_log2FC)) %>%
  dplyr::distinct(gene, .keep_all = TRUE) %>%
  mutate(rank = row_number())
  
vasc_ranked %>% dplyr::filter(gene == "CLDN5") %>% dplyr::select(rank, gene, avg_log2FC, p_val_adj)
vasc_top <- vasc_ranked %>% slice_head(n = 200)

df = add_spg_gene_expression("Vasc", vasc_top$gene, df, spe_ntc, threshold = "q75")
prop_df = make_spg_prop_df(df,"Vasc",200,threshold = "q75", domain_levels, domain_labels) 
prop_df_Vasc = prop_df

#SPG: Vasc 
#Gene set: Vasc_top200 
#Requested genes: 200 
#Matched genes: 200 
#Missing genes: 0 

# ============================================================
# Neuron
# ============================================================
neun1 = read_excel(here("raw-data/images/SPG_Spot_Valid/NeuN/TableS13_marker_stats_supp_table.xlsx"), sheet="marker_stats_supp_table")
neun1 = as.data.frame(neun1)
neun1$padj = as.numeric(neun1$log.p.value)
celltypes = c("Inhib", "Excit", "Excit_L3/4/5", "Excit_L3", "Excit_L4", "Excit_L6", "Excit_L5/6", "Excit_L5", "Excit_L2/3")
gene_idx <- match(neun1$symbol,rowData(spe_ntc)$gene_name)
missing_genes <- neun1$symbol[is.na(gene_idx)]

neun_ranked <- neun1 %>%
  filter(cellType.target %in% celltypes) %>%
  mutate(rank_marker = as.numeric(rank_marker), std.logFC = as.numeric(std.logFC)) %>%
  filter(!is.na(symbol), symbol != "") %>%
  arrange(rank_marker, desc(std.logFC)) %>%
  distinct(symbol, .keep_all = TRUE)

neun_ranked %>% dplyr::filter(symbol == "PCP4") %>% dplyr::select(rank_marker, symbol, logFC, log.FDR, log.p.value)
neun_top <- neun_ranked %>%slice_head(n = 200)

df = add_spg_gene_expression("Neuron", neun_top$symbol, df, spe_ntc, threshold = "q75")
prop_df = make_spg_prop_df(df,"Neuron",200,threshold = "q75", domain_levels, domain_labels) 
prop_df_Neuron = prop_df

#SPG: Neuron 
#Gene set: Neuron_top200 
#Requested genes: 200 
#Matched genes: 200 
#Missing genes: 0 

# ============================================================
# Neuropil
# ============================================================
neuropil <- read_excel(here("raw-data/images/SPG_Spot_Valid/Neuropil/Supple_Table_10_Neuropil.xlsx"),sheet = "10_human_synase_markers",skip = 1)
neuropil <- neuropil[, 9:15]
colnames(neuropil) <- c("gene", "p_val", "avg_log2FC","pct.1", "pct.2", "p_val_adj", "cluster")
neuropil <- neuropil %>%
 	 		mutate(gene = toupper(trimws(as.character(gene))),
   		 	avg_log2FC = as.numeric(avg_log2FC),
    		p_val_adj = as.numeric(p_val_adj))
gene_idx <- match(neuropil$gene,rowData(spe_ntc)$gene_name)
missing_genes <- neuropil$gene[is.na(gene_idx)]

# Select excitatory and inhibitory synaptic markers
neuropil_ranked <- neuropil %>%
  filter(cluster %in% c("Synapse_ExPFC", "Synapse_In", "LO-synapse", "N-synapse"),
    !is.na(gene), gene != "", p_val_adj < 0.05, avg_log2FC > 0) %>%
  arrange(desc(avg_log2FC), p_val_adj) %>%
  distinct(gene, .keep_all = TRUE) %>%
  mutate(rank = row_number())

neuropil_top <- neuropil_ranked %>% slice_head(n = 200)

df = add_spg_gene_expression("Neuropil", neuropil_top$gene, df, spe_ntc, threshold = "q75")
prop_df = make_spg_prop_df(df,"Neuropil",200,threshold = "q75", domain_levels, domain_labels) 
prop_df_Neuropil = prop_df

#SPG: Neuropil 
#Gene set: Neuropil_top200 
#Requested genes: 200 
#Matched genes: 191 
#Missing genes: 9 
#Missing: MTCO1P12, KIF5C-AS1, MTND4P12, MTCO3P12, MTND2P28, MTCO2P12, RNU4-2, RNY1, AL162581.1 



# ============================================================
# FINAL MULTIPANEL FIGURE: domain strip | image | genes | validation
# ============================================================
suppressPackageStartupMessages({
  library(patchwork)
  library(ggpattern)  # install.packages('ggpattern') if needed
})

spg_levels <- c('Neuron', 'Neuropil', 'PNN', 'Vasc')
#sample_levels <- sort(unique(as.character(df$sample_id[!is.na(df$sample_id)])))
sample_levels <- sort(unique(as.character(df$brnum_dx[!is.na(df$brnum_dx)])))

# Work with original domain IDs until the final plotting step.
# Explicitly complete all sample-domain combinations to show absent tissue.
image_df <- df %>%
  #filter(!is.na(sample_id), PRECAST_07 %in% domain_levels) %>%
  filter(!is.na(brnum_dx), !is.na(sample_id), PRECAST_07 %in% domain_levels) %>%
  dplyr::select(brnum_dx, sample_id, PRECAST_07, neun_pos, neuropil_pos, pnn_pos, vasc_pos) %>%
  #dplyr::select(sample_id, PRECAST_07, neun_pos, neuropil_pos, pnn_pos, vasc_pos) %>%
  mutate(PRECAST_07 = as.character(PRECAST_07)) %>%
  pivot_longer(c(neun_pos, neuropil_pos, pnn_pos, vasc_pos),
               names_to = 'SPG', values_to = 'positive') %>%
  mutate(SPG = recode(SPG, neun_pos = 'Neuron', neuropil_pos = 'Neuropil',
                      pnn_pos = 'PNN', vasc_pos = 'Vasc')) %>%
  #group_by(sample_id, SPG, PRECAST_07) %>%
  group_by(brnum_dx, sample_id, SPG, PRECAST_07) %>%
  summarise(n_spots = n(), n_high = sum(positive %in% TRUE), .groups = 'drop') %>%
  tidyr::complete(SPG = spg_levels, brnum_dx = sample_levels, #sample_id = sample_levels,
                  PRECAST_07 = domain_levels,
                  fill = list(n_spots = 0L, n_high = 0L)) %>%
  group_by(brnum_dx, sample_id, SPG) %>%
  mutate(sample_high = sum(n_high),
         prop_original = if_else(n_spots > 0 & sample_high > 0,
                                 n_high / sample_high, NA_real_),
         prop_within_domain = if_else(n_spots > 0,
                                      n_high / n_spots, NA_real_)) %>%
  ungroup() %>%
  mutate(PRECAST_07 = factor(PRECAST_07, levels = domain_levels,
                             labels = domain_labels))

gene_df <- bind_rows(
  Neuron = prop_df_Neuron,
  Neuropil = prop_df_Neuropil,
  PNN = prop_df_PNN,
  Vasc = prop_df_Vasc,
  .id = 'SPG'
) %>%
  mutate(PRECAST_07 = factor(as.character(PRECAST_07), levels = domain_labels))

# Correlation is one observation per sample x spatial domain.
# Require both measurements to be defined and domain present in both.
cor_df <- gene_df %>%
  transmute(brnum_dx, sample_id, PRECAST_07, SPG,
            gene_prop = prop_original, gene_n_spots = n_spots) %>%
  inner_join(
    image_df %>% transmute(brnum_dx, sample_id, PRECAST_07, SPG,
                            image_prop = prop_original,
                            image_n_spots = n_spots),
    by = c('brnum_dx', 'sample_id', 'PRECAST_07', 'SPG')
  ) %>%
  filter(gene_n_spots > 0, image_n_spots > 0,
         is.finite(gene_prop), is.finite(image_prop))

cor_results <- cor_df %>%
  group_by(SPG) %>%
  summarise(n = n(),
            pearson_r = if (n() >= 3 && sd(gene_prop) > 0 && sd(image_prop) > 0)
              cor(gene_prop, image_prop) else NA_real_,
            pearson_p = if (n() >= 3 && sd(gene_prop) > 0 && sd(image_prop) > 0)
              cor.test(gene_prop, image_prop)$p.value else NA_real_,
            .groups = 'drop')
print(cor_results)

# Plot domain labels top-to-bottom: WM -> L1/M.
# A discrete ggplot y-axis displays factor levels bottom-to-top,
# so reverse the plotting levels for the requested top-to-bottom order.
plot_levels <- rev(domain_labels)
plot_samples <- sample_levels

format_heat_data <- function(x) {
  x %>% mutate(
    #sample_id = factor(as.character(sample_id), levels = plot_samples),
	brnum_dx = factor(as.character(brnum_dx), levels = plot_samples),
    PRECAST_07 = factor(as.character(PRECAST_07), levels = plot_levels),
    missing_domain = n_spots == 0,
    undefined_prop = n_spots > 0 & !is.finite(prop_original)
  )
}
image_plot_df <- format_heat_data(image_df)
gene_plot_df <- format_heat_data(gene_df)

# Missing tissue: white crosshatched squares.
# Present but undefined (sample_high=0): light gray, not crosshatched.
heatmap_panel <- function(dat, title = NULL, show_y = FALSE,
                          show_x = FALSE, show_legend = FALSE) {
  present <- dat %>% filter(!missing_domain)
  absent <- dat %>% filter(missing_domain)

  p <- ggplot() +
    geom_tile(data = present,
              aes(x = brnum_dx, y = PRECAST_07, fill = prop_original),
			  #aes(x = sample_id, y = PRECAST_07, fill = prop_original),
              color = 'grey80', linewidth = 0.12) +
    ggpattern::geom_tile_pattern(
      data = absent,
      aes(x = brnum_dx, y = PRECAST_07),
	  #aes(x = sample_id, y = PRECAST_07),
      fill = 'white', colour = 'grey65', linewidth = 0.12,
      pattern = 'crosshatch', pattern_fill = 'grey65',
      pattern_colour = 'grey65', pattern_angle = 45,
      pattern_density = 0.25, pattern_spacing = 0.035,
      pattern_size = 0.12
    ) +
    scale_fill_gradient(low = 'white', high = 'black', limits = c(0, 1),
                        na.value = 'grey90', name = 'Proportion') +
    scale_x_discrete(drop = FALSE) +
    scale_y_discrete(drop = FALSE) +
    labs(title = title, x = NULL, y = NULL) +
    theme_minimal(base_size = 9) +
    theme(panel.grid = element_blank(),
          plot.title = element_text(face = 'bold', size = 15, hjust = 0.5),
          axis.text.x = if (show_x) element_text(angle = 90, hjust = 1,
                                                  vjust = 0.5, size = 13) else element_blank(),
          axis.ticks.x = element_blank(),
          axis.text.y = if (show_y) element_text(size = 7) else element_blank(),
          axis.ticks.y = element_blank(),
          legend.position = if (show_legend) 'bottom' else 'none',
          legend.key.width = grid::unit(1.2, 'cm'),
          plot.margin = margin(3, 3, 3, 3))
  p
}

# Narrow color annotation at the left of EACH row, aligned with domains.
strip_df <- data.frame(
  PRECAST_07 = factor(domain_labels, levels = plot_levels),
  x = 'Domain',
  domain = domain_labels
)
domain_strip <- function(title = NULL) {
  ggplot(strip_df, aes(x = x, y = PRECAST_07, fill = domain)) +
    geom_tile(color = 'white', linewidth = 0.2) +
    scale_fill_manual(values = domain_colors, guide = 'none') +
    scale_y_discrete(drop = FALSE) +
    labs(title = title, x = NULL, y = NULL) +
    theme_void(base_size = 9) +
    theme(plot.title = element_text(face = 'bold', size = 9, hjust = 0.5),
          plot.margin = margin(3, 0, 3, 0))
}

validation_panel <- function(spg, title = NULL, show_x = FALSE,
                             show_legend = FALSE) {
  sub <- cor_df %>% filter(SPG == spg) %>%
    mutate(PRECAST_07 = factor(as.character(PRECAST_07),
                               levels = domain_labels))
  stats <- cor_results %>% filter(SPG == spg)
  label <- if (nrow(stats) == 1 && is.finite(stats$pearson_r)) {
    paste0('Pearson r = ', sprintf('%.2f', stats$pearson_r),
           '\np = ', format.pval(stats$pearson_p, digits = 2, eps = 0.001),
           '\nn = ', stats$n)
  } else {
    paste0('Pearson r = NA\nn = ', nrow(sub))
  }
  p <- ggplot(sub, aes(x = gene_prop, y = image_prop)) +
    geom_point(aes(color = PRECAST_07), alpha = 0.75, size = 1.2) +
    scale_color_manual(values = domain_colors, breaks = domain_labels,
                       drop = FALSE, name = 'Spatial domain') +
    scale_x_continuous(limits = c(0, 0.83), expand = expansion(mult = 0.03)) +
    scale_y_continuous(limits = c(0, 0.83), expand = expansion(mult = 0.03)) +
    annotate('text', x = 0.04, y = 0.7, label = label,
             hjust = 0, vjust = 1, size = 4) +
    labs(title = title,
        # x = if (show_x) 'Marker genes enrichment ' else NULL,
		x=NULL,
         y = 'Image based classifications') +
    theme_bw(base_size = 9) +
    theme(plot.title = element_text(face = 'bold', size = 15, hjust = 0.5),
          axis.text = element_text(size = 11),
          axis.title = element_text(size = 13),
          legend.position = if (show_legend) 'bottom' else 'none',
          legend.text = element_text(size = 6),
          panel.grid.minor = element_blank(),
          plot.margin = margin(8, 8, 8, 8))
  if (nrow(sub) >= 3 && sd(sub$gene_prop) > 0 && sd(sub$image_prop) > 0) {
    p <- p + geom_smooth(method = 'lm', se = TRUE, color = 'black',
                         linewidth = 0.45, na.rm = TRUE)
  }
  p
}

library(cowplot)

legend_plot <- ggplot(
  data.frame(x = c(0, 1), y = c(1, 1)),
  aes(x = x, y = y, fill = x)
) +
  geom_tile() +
  scale_fill_gradient(
    low = "white",
    high = "black",
    limits = c(0, 1),
    name = "Proportion",
    guide = guide_colorbar(
      direction = "horizontal",
      title.position = "top",
      barwidth = grid::unit(3, "cm"),
      barheight = grid::unit(1, "cm")
    )
  ) +
  theme_void() +
  theme(
    legend.position = "bottom",
	legend.title.align = 0,
	legend.title = element_text(size = 13,hjust = 0,margin = margin(r = 12)),
    legend.text = element_text(size = 9),
    legend.margin = margin(0, 0, 0, 0)
  )

heatmap_legend <- cowplot::get_legend(legend_plot)
# One patchwork row per microenvironment, with equal column widths.

# ============================================================
# LEFT-SIDE ANNOTATIONS
# SPG name | Domain labels | Thin color strip
# ============================================================

spg_display_names <- c(
  "Neuron"   = "Neuronal",
  "Neuropil" = "Neuropil",
  "PNN"      = "PNN",
  "Vasc"     = "Vasculature"
)

strip_df <- data.frame(
  PRECAST_07 = factor(domain_labels, levels = plot_levels),
  domain = domain_labels
)

# 1. Vertical microenvironment name
spg_name_panel <- function(spg) {
  ggplot() +
    annotate(
      "text",
      x = 1, y = 1,
      label = spg_display_names[[spg]],
      angle = 90,
      fontface = "bold",
      size = 5.5
    ) +
    xlim(0.5, 1.5) +
    ylim(0.5, 1.5) +
    theme_void() +
    theme(
      plot.margin = margin(3, 0, 3, 0)
    )
}

# 2. Domain names immediately before color strip
domain_labels_panel <- function() {
  ggplot(
    strip_df,
    aes(x = 1, y = 4 + (as.numeric(PRECAST_07) - 4) * 0.80)
  ) +
    geom_text(
      aes(label = domain),
      hjust = 1,
      size = 3.6
    ) +
    scale_x_continuous(
      limits = c(0, 1),
      expand = c(0, 0)
    ) +
    scale_y_discrete(
      limits = plot_levels,
      drop = FALSE,
      expand = c(0, 0)
    ) +
    coord_cartesian(clip = "off") +
    theme_void() +
    theme(
      plot.margin = margin(3, 0, 3, 0)
    )
}

# 3. Very thin domain color annotation
domain_strip <- function() {
  ggplot(
    strip_df,
    aes(x = 1, y = PRECAST_07, fill = domain)
  ) +
    geom_tile(
      width = 1,
      height = 1,
      color = "white",
      linewidth = 0.15
    ) +
    scale_fill_manual(
      values = domain_colors,
      guide = "none"
    ) +
    scale_x_continuous(
      limits = c(0.5, 1.5),
      expand = c(0, 0)
    ) +
    scale_y_discrete(
      limits = plot_levels,
      drop = FALSE,
      expand = c(0, 0)
    ) +
    theme_void() +
    theme(
      plot.margin = margin(3, 0, 3, 0)
    )
}

# ============================================================
# Ensure heatmap domain rows align with annotations
# ============================================================

# In your existing heatmap_panel(), replace:
# scale_y_discrete(drop = FALSE)
#
# with:
# scale_y_discrete(
#   limits = plot_levels,
#   drop = FALSE,
#   expand = c(0, 0)
# )

# ============================================================
# FINAL MULTIPANEL FIGURE
# ============================================================

row_plots <- lapply(seq_along(spg_levels), function(i) {

  spg <- spg_levels[i]

  top <- i == 1
  bottom <- i == length(spg_levels)

  spg_name_panel(spg) +

    domain_labels_panel() +

    domain_strip() +

    heatmap_panel(
      image_plot_df %>% filter(SPG == spg),
      title = if (top)
        "SPG-associated image-based classifications"
      else NULL,
      show_y = FALSE,
      show_x = bottom,
      show_legend = FALSE
    ) +

    heatmap_panel(
      gene_plot_df %>% filter(SPG == spg),
      title = if (top)
        "SPG-associated marker genes enrichment"
      else NULL,
      show_y = FALSE,
      show_x = bottom,
      show_legend = FALSE
    ) +

    validation_panel(
      spg,
      title = if (top) "Validation" else NULL,
      show_x = bottom,
      show_legend = FALSE
    ) +

    patchwork::plot_layout(
      ncol = 6,
      widths = c(
        0.35,  # Vertical SPG name
        1,  # Spatial domain labels
        0.12,  # Thin domain color strip
        5,     # Image-based heatmap
        5,     # Gene-expression heatmap
        2    # Validation
      )
    )
})

final_plot <- patchwork::wrap_plots(
  row_plots,
  ncol = 1
)

# ============================================================
# SAVE
# ============================================================
final_plot_with_legend <- cowplot::ggdraw() +
  cowplot::draw_plot(final_plot, 0, 0, 1, 1) +
  cowplot::draw_grob(
    heatmap_legend,
    x = -0.001,
    y = 0.0003,
    width = 0.1,
    height = 0.1
  ) +
  cowplot::draw_label(
      "Marker genes enrichment",
      x = 0.91,
      y = 0.050,
      size = 15,
      hjust = 0.5
    )

	out_dir <- here('plots', 'image_processing', 'revisions')
	dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
	
ggsave(
  file.path(
    out_dir,
    "SPG_Combined_Image_Gene_Validation.pdf"
  ),
  plot = final_plot_with_legend ,
  width = 18,
  height = 16,
  limitsize = FALSE,
  device = grDevices::cairo_pdf
)

print(final_plot)

write.csv(cor_results,
          file.path(out_dir, 'SPG_Combined_Pearson_Results.csv'),
          row.names = FALSE)


#################
#dot plots
#################

make_spg_dotplot <- function(df, pos_col, exp_col, property) {
  total_pos <- sum(!is.na(df$PRECAST_07) & df[[pos_col]] == TRUE,na.rm = TRUE)
  props <- df %>%
    filter(!is.na(PRECAST_07), !is.na(.data[[pos_col]])) %>%
    group_by(PRECAST_07) %>%
    summarise(n_spots = n(),
      n_pos = sum(.data[[pos_col]] == TRUE),
      SPGprop = n_pos / total_pos,
      SPDprop = n_pos / n_spots,
      sum_exp = sum(.data[[exp_col]], na.rm = TRUE),
      mean_exp = mean(.data[[exp_col]], na.rm = TRUE),
      .groups = "drop") %>%
    mutate(PRECAST_07 = factor(PRECAST_07, levels = domain_levels, labels = domain_labels))
	print(props %>%dplyr::select(PRECAST_07, mean_exp, sum_exp))

	summary(props$mean_exp)
  p <- ggplot(props, aes(x = PRECAST_07, y = SPGprop, size = SPDprop)) +
    # Horizontal domain color strip
       geom_tile(data = props, aes(x = PRECAST_07, y = -0.035, fill = PRECAST_07),height = 0.025,width = 0.9,inherit.aes = FALSE) +
       scale_fill_manual(values = domain_colors, guide = "none") +
       ggnewscale::new_scale_fill() +
    # Dotplot
       geom_point(aes(fill = mean_exp),alpha = 0.7,shape = 21) +
       scale_fill_gradientn(colours = c("white", "grey98", "grey95", "grey85", "grey60", "grey30", "black"),
      						#values = scales::rescale(c(0, 100, 200, 500, 1000, 2500, 5000), from = c(0, 5000)),
							values = scales::rescale(c(0, 0.01, 0.02, 0.04, 0.06, 0.08, 0.1), from = c(0, 0.1)),
      					  	#limits = c(0, 5000), breaks = c(0, 1000, 2500, 5000),
							limits = c(0, 0.1), breaks = c(0, 0.03, 0.06, 0.1),
      					    labels = scales::label_number(), oob = scales::squish) +
	   scale_size_continuous(limits = c(0, 0.75),  breaks = c(0.10, 0.30, 0.50, 0.70), range = c(5, 18), labels = scales::percent) +
       annotate("text", x = 4, y = Inf, label = property, hjust = 0.5, vjust = 1.5, size = 12, fontface = "bold") + # Title inside top-right corner
       coord_cartesian(ylim = c(-0.06, 0.75),clip = "off") +
       labs(x = "", y = "SPG spot distribution across SpDs", fill = "Mean normalized\nexpression",size = "SPG prevalence\nwithin SpD") +
    theme_minimal() +
    theme(axis.text = element_text(size = 18), axis.text.x = element_text(angle = 90,hjust = 1), axis.title = element_text(size = 18), 
		 legend.title = element_text(size = 16), legend.text = element_text(size = 16),
	  	 panel.grid.major.x = element_line(color = "grey85",linewidth = 0.4),
	     panel.grid.major.y = element_line(color = "grey85",linewidth = 0.4),
	     panel.grid.minor = element_blank())
  return(p)
}

## ------------------------------------------------------------
## Generate plots
## ------------------------------------------------------------

p_wfa <- make_spg_dotplot(df, pos_col = "pnn_pos", exp_col = "PNN_top200_norm_sum", property = "PNN") +
		 guides(size = "none", fill = "none")

p_dapi <- make_spg_dotplot(df, pos_col = "neuropil_pos", exp_col = "Neuropil_top200_norm_sum", property = "Neuropil") +
		  theme(axis.title = element_blank(), axis.text = element_blank(), axis.text.x = element_blank(), axis.ticks = element_blank()) +
		  guides(size = "none")

p_neun <- make_spg_dotplot(df, pos_col = "neun_pos", exp_col = "Neuron_top200_norm_sum",property = "Neuronal") +
		  theme(axis.title.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank()) +
		  guides(size = "none", fill = "none")

p_claudin5 <- make_spg_dotplot(df, pos_col = "vasc_pos", exp_col = "Vasc_top200_norm_sum", property = "Vasculature") +
			  theme(axis.title.y = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank()) + 
			  guides(fill = "none")

#pdf(here("plots","image_processing","revisions","dot_plots_mean1.pdf"),width = 8,height = 8)
library(patchwork)

p_all <- (p_neun | p_dapi) /
         (p_wfa | p_claudin5)

p_all <- wrap_plots(A = p_neun, B = p_dapi, C = p_wfa, D = p_claudin5,
		    design = "A#B
					  ###
					  C#D",
		    widths = c(1, 0.03, 1),heights = c(1, 0.02, 1))

ggsave(here("plots", "image_processing", "revisions", "dot_plots_mean.png"), p_all, width = 16, height = 12, dpi = 300)

## ------------------------------------------------------------
## paired t-test
## ------------------------------------------------------------
library(purrr)
image_df_pnn <- image_df %>% dplyr::filter(SPG == "PNN")
  
paired_SPG_test <- function(data, domain1, domain2) {
  tmp <- image_df_pnn %>%
    filter(PRECAST_07 %in% c(domain1, domain2)) %>%
    dplyr::select(brnum_dx, PRECAST_07, prop_original) %>%
    pivot_wider(names_from = PRECAST_07,values_from = prop_original) %>%
    drop_na()
  test <- t.test(tmp[[domain1]],tmp[[domain2]],paired = TRUE)
  tibble(
    domain1 = domain1,
    domain2 = domain2,
    n_pairs = nrow(tmp),
    mean_domain1 = mean(tmp[[domain1]]),
    mean_domain2 = mean(tmp[[domain2]]),
    mean_difference = mean(tmp[[domain1]] - tmp[[domain2]]),
    t = unname(test$statistic),
    df = unname(test$parameter),
    p_value = test$p.value
  )
}



paired_domain_test <- function(data, domain1, domain2) {
  tmp <- image_df_pnn %>%
    filter(PRECAST_07 %in% c(domain1, domain2)) %>%
    dplyr::select(brnum_dx, PRECAST_07, prop_within_domain) %>%
    pivot_wider(names_from = PRECAST_07,values_from = prop_within_domain) %>%
    drop_na()
  test <- t.test(tmp[[domain1]],tmp[[domain2]],paired = TRUE)
  tibble(
    domain1 = domain1,
    domain2 = domain2,
    n_pairs = nrow(tmp),
    mean_domain1 = mean(tmp[[domain1]]),
    mean_domain2 = mean(tmp[[domain2]]),
    mean_difference = mean(tmp[[domain1]] - tmp[[domain2]]),
    t = unname(test$statistic),
    df = unname(test$parameter),
    p_value = test$p.value
  )
}


domains <- levels(image_df_pnn$PRECAST_07)

spd02_tests <- setdiff(domains, "SpD02-L3/4") %>% 
			   map_dfr(~paired_SPG_test(image_df_pnn,"SpD02-L3/4",.x)) %>%
  			   mutate(p_adj = p.adjust(p_value, method = "BH"))


spd05_tests <- setdiff(domains, "SpD05-L5") %>%
 	 		   map_dfr(~paired_SPG_test(image_df_pnn,"SpD05-L5",.x)) %>%
 			   mutate(p_adj = p.adjust(p_value, method = "BH"))


spd02_tests1 <- setdiff(domains, "SpD02-L3/4") %>% 
			   map_dfr(~paired_domain_test(image_df_pnn,"SpD02-L3/4",.x)) %>%
  			   mutate(p_adj = p.adjust(p_value, method = "BH"))


spd05_tests1 <- setdiff(domains, "SpD05-L5") %>%
 	 		   map_dfr(~paired_domain_test(image_df_pnn,"SpD05-L5",.x)) %>%
 			   mutate(p_adj = p.adjust(p_value, method = "BH"))
			   
as.data.frame(spd02_tests)
as.data.frame(spd05_tests)
as.data.frame(spd02_tests1)
as.data.frame(spd05_tests1)