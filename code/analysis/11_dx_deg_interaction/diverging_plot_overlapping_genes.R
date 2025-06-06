# Load library ----
suppressPackageStartupMessages({
  library(here)
  library(tidyverse)
  library(ggstats)
  library(sessioninfo)
})

# Load Data ----
## Load spd annotation ----
spd_anno_df <- read_csv(
  here(
    "processed-data/man_anno",
    "spd_labels_k7.csv"
  )
) |>
  mutate(
    anno_lab = factor(
      paste0(gsub("spd", "SpD", spd), "-", label),
      levels = c(
        "SpD07-L1/M",
        "SpD06-L2/3",
        "SpD02-L3/4",
        "SpD05-L5",
        "SpD03-L6",
        "SpD01-WMtz",
        "SpD04-WM"
      )
    )
  )

## Load 172 DEGs ----
## 172 genes ----
spd_deg_172 <- read_csv(
  here(
    "code/analysis/10_dx_deg_adjust_spd",
    "172_prelim_fdr010.csv"
  ),
  col_types = cols()
)


## load Layer_specific DEGs ----
spd_files <- list.files(
  "processed-data/rds/11_dx_deg_interaction", ,
  pattern = "layer_specific_logFC_.*\\.csv",
  full.names = TRUE
)

names(spd_files) <- str_extract(
  spd_files,
  "(?<=layer_specific_logFC_).*?(?=\\.csv)"
)

spd_deg_list <- map(
  spd_files,
  ~ read_csv(.x, col_types = cols()) |>
    # mutate(gene = AnnotationDbi::mapIds(
    #   org.Hs.eg.db::org.Hs.eg.db,
    #   keys = gene_id,
    #   column = "SYMBOL",
    #   keytype = "ENSEMBL",
    #   multiVals = "first"
    # )) |>
    filter(
      P.Value < 0.05
    ) |>
    mutate(
      dx_172 = ifelse(
        gene_id %in% spd_deg_172$ensembl,
        TRUE,
        FALSE
      )
    )
)

spd_deg_list |> map(~ .x |> nrow())

# Based on directionality
n_sig_gene_per_spd <- map_dfr(
  spd_deg_list,
  ~ .x |>
    # group_by(gene) |>
    summarise(
      up_novel = sum(logFC > 0 & dx_172 == FALSE),
      up_overlapping = sum(logFC > 0 & dx_172 == TRUE),
      down_novel = sum(logFC < 0 & dx_172 == FALSE),
      down_overlapping = sum(logFC < 0 & dx_172 == TRUE)
    ),
  .id = "spd"
)

n_sig_gene_per_spd$spd <- spd_anno_df$anno_lab[match(n_sig_gene_per_spd$spd, spd_anno_df$spd)]


n_sig_gene_per_spd <- n_sig_gene_per_spd[rev(order(spd_anno_df$label)), ]

long_df <- n_sig_gene_per_spd |>
  pivot_longer(
    cols = c("up_novel", "up_overlapping", "down_novel", "down_overlapping"),
    names_to = "direction",
    values_to = "n_genes"
  ) |>
  mutate(
    direction = factor(
      direction,
      levels = c(
        "up_novel",
        "up_overlapping",
        "down_overlapping",
        "down_novel"
      )
    )
  )

ggplot(long_df) +
  aes(y = spd, fill = direction, weight = n_genes) +
  geom_diverging() +
  # geom_diverging_text() +
  labs(
    title = "Diverging Bar Plot of Significant Genes by Direction",
    x = "Number of Genes",
    y = "SpD Annotation",
    fill = "Direction"
  ) +
  scale_fill_manual(
    labels = c(
      "up_novel" = "Upregulated (Novel)",
      "up_overlapping" = "Upregulated (Overlapping)",
      "down_overlapping" = "Downregulated (Overlapping)",
      "down_novel" = "Downregulated (Novel)"
    ),
    breaks = c(
      "down_novel",
      "down_overlapping",
      "up_overlapping",
      "up_novel"
    ),
    values = c(
      "down_novel" = "#6f6fff",
      "down_overlapping" = "blue",
      "up_overlapping" = "red",
      "up_novel" = "#ef8080"
    )
  ) +
  theme_minimal() +
  theme(
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 14, face = "bold"),
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 10),
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5)
  ) +
  scale_x_continuous(
    labels = label_number_abs(),
    limits = symmetric_limits,
    breaks = seq(-1700, 1700, 300)
  ) +
  scale_y_discrete(limits = rev(levels(long_df$spd)))

ggsave(
  here("plots/11_dx_deg_interaction", "diverging_barplot_genes_nomial_p05.pdf"),
  width = 8,
  height = 6,
  dpi = 300
)


# Add a column to classify the direction of change
# data$Direction <- ifelse(data$LogFC > 0, "Upregulated", "Downregulated")

# Create the diverging bar plot
# ggplot(data, aes(x = reorder(Gene, LogFC), y = LogFC, fill = Direction)) +
#   geom_bar(stat = "identity") +
#   scale_fill_manual(values = c("Upregulated" = "blue", "Downregulated" = "red")) +
#   coord_flip() +
#   labs(
#     title = "Diverging Bar Plot of Genes",
#     x = "Genes",
#     y = "Log Fold Change"
#   ) +
#   theme_minimal()


# ggstats diverging plot example ----
# d <- Titanic |> as.data.frame()

# ggplot(d) +
#   aes(y = Class, fill = Sex, weight = Freq) +
#   geom_diverging() +
#   geom_diverging_text()
