# Load library ----
suppressPackageStartupMessages({
  library(here)
  library(tidyverse)
  library(ggrepel)
  library(UpSetR)
})

# Load data ----
# NOTE: These are layer-restricted DEG results
spd_files <- list.files(
  "processed-data/rds/11_dx_deg_interaction",
  # pattern = "layer_specific_logFC_.*\\.csv", # Previous file names
  pattern = "layer_restricted_logFC_.*\\.csv",
  full.names = TRUE
)

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

spd_deg_list <-
  spd_files |>
  set_names(
    # Retrieve annotated names
    nm = spd_anno_df[
      match(
        str_extract(
          spd_files,
          "(?<=layer_restricted_logFC_).*?(?=\\.csv)"
        ),
        spd_anno_df$spd
      ),
      "anno_lab",
      drop = TRUE
    ]
  ) |>
  map(
    ~ read_csv(.x)
  )


# Descriptive Statistics ----
## Nomial p-value < 0.05 ----
nom_p_vec <- spd_deg_list |>
  map_int(
    ~ sum(.x$P.Value < 0.05, na.rm = TRUE)
  )

# All layer-restricted DEGs. including same genes in multiple SpDs
sum(nom_p_vec)

# Unique genes symbols amopng All layer-restricted DEGs
spd_deg_list |>
  map(
    ~ .x |>
      filter(P.Value < 0.05) |>
      pull(gene_id)
  ) |>
  unlist() |>
  unique() |>
  length()

# Total number of up and down in these two domains. 
spd_deg_list[c("SpD01-WMtz", "SpD04-WM")] |>
  imap_dfr(
    ~ .x |>
      filter(P.Value < 0.05) |>
      mutate(spd = .y)
  ) |>
  mutate(
    direction = case_when(
      logFC > 0 ~ "up",
      logFC < 0 ~ "down"
    )
  ) |>
  group_by(direction) |>
  summarise(
    n = n()
  )



# SpD01-WMtz SpD02-L3/4   SpD03-L6   SpD04-WM   SpD05-L5 SpD06-L2/3
#       1581        284        283       2521        468        782
#   SpD07-L1
#       1059


## FDR < 0.1 ----
fdr_10_vec <- spd_deg_list |>
  map_int(
    ~ sum(.x$adj.P.Val < 0.1, na.rm = TRUE)
  )

# total number of genes
sum(fdr_10_vec)

# SpD01-WMtz SpD02-L3/4   SpD03-L6   SpD04-WM   SpD05-L5 SpD06-L2/3
#         30          0          0        894          1          2
#   SpD07-L1
#         17
# Total number of unique genes
spd_deg_list |>
  map(
    ~ .x |>
      filter(adj.P.Val < 0.1) |>
      pull(gene_id)
  ) |>
  unlist() |>
  unique() |>
  length()



cbind(
  nom_p_vec,
  fdr_10_vec
)

#            nom_p_vec fdr_10_vec
# SpD01-WMtz      1581         30
# SpD02-L3/4       284          0
# SpD03-L6         283          0
# SpD04-WM        2521        894
# SpD05-L5         468          1
# SpD06-L2/3       782          2
# SpD07-L1        1059         17


# Upset plots ----
## Nominal p-value < 0.05 ----
nom_p_geneList <- spd_deg_list |>
  map(
    ~ .x$gene_id[.x$P.Value < 0.05]
  )

pdf(
  here(
    "plots/11_dx_deg_interaction",
    "upset_layer_gene_nominal_p_value.pdf"
  ),
  width = 2.02,
  height = 1.31,
  onefile = FALSE
)
upset(
  fromList(nom_p_geneList),
  nsets = 7,
  nintersects = 7,
  sets = c(
    "SpD07-L1/M", "SpD06-L2/3", "SpD02-L3/4", "SpD05-L5", "SpD03-L6",
    "SpD01-WMtz", "SpD04-WM"
  ) |> rev(),
  keep.order = TRUE,
  order.by = c("degree"),
  decreasing = FALSE,
  sets.x.label = "# of Layer-restricted DEGs",
  mainbar.y.label = "# of Layer-specific DEGs",
  scale.sets = "identity",
  text.scale = c(0.3),
  point.size = 0.5,
  line.size = 0.2,
  mb.ratio = c(0.5, 0.5),
  scale.intersections = "identity"
) #|> ggsave(
#   filename = here(
#     "plots/11_dx_deg_interaction",
#     "upset_layer_gene_nominal_p_value.pdf"
#   ),
#   plot = ggplotify::as.ggplot(_),
#   width = 2.02,
#   height = 1.31,
#   unit = "in"
# )
dev.off()

## FDR p-value < 0.10 ----

fdr_10_geneList <- spd_deg_list |>
  map(
    ~ .x$gene_id[.x$adj.P.Val < 0.10]
  )


pdf(
  here(
    "plots/11_dx_deg_interaction",
    "upset_layer_gene_fdr_10.pdf"
  ),
  width = 10,
  height = 8
)
upset(
  fromList(fdr_10_geneList),
  nsets = 7,
  # sets = sort(
  #   names(enrich_list),
  #   decreasing = TRUE
  # ),
  keep.order = TRUE,
  nintersects = NA,
  order.by = c("freq"),
  # main.bar.color = "blue",
  # sets.bar.color = "red",
  # mainbar.y.label = "Enriched Genes",
  sets.x.label = "Spatial Domains"
)
dev.off()
