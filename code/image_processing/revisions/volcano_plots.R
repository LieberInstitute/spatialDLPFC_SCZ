# ================= Facet-based 1x4 Volcano (Final; no FDR dashed line) =================
# - Unified sans font (Liberation/DejaVu/Arial) to avoid kerning/fallback issues
# - Clean facet strips with title-only by default; switch to "full info" via a single flag
# - Pure red/blue color scheme; nominal = semi-transparent; stroke == fill
# - y-guide: keep solid line at nominal p = 0.05 (no FDR-derived dashed line)
# - Layout polish: thin panel borders, short y-ticks, no major grids, fixed +/- LFC cut lines
# - Exports: PDF/PNG/TIFF (+SVG if svglite installed); writes session info
# - All comments in English; ready for GitHub sharing
# =======================================================================================

setwd('/dcs04/lieber/marmaypag/spatialDLPFC_SCZ_LIBD4100/')
suppressPackageStartupMessages({
  library(tidyverse)
  library(readxl)
  library(glue)
  library(systemfonts)
  library(ragg)
  library(ggrepel)
  library(sessioninfo)
})


stamp   <- format(Sys.time(), "%Y%m%d_%H%M")
out_dir <- here("plots/image_processing/revisions",paste0("Fig6A_volcanos_", stamp, "_Sans_noFDRline"))
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

pdf_out  <- file.path(out_dir, "Fig6A_Volcano_FACET_1x4_6.45x2.0in_noFDRline.pdf")
png_out  <- file.path(out_dir, "Fig6A_Volcano_FACET_1x4_6.45x2.0in_noFDRline.png")


## ---- 2) Fonts & figure sizing -----------------------------------------------------------
# Prefer modern TrueType/OpenType; avoid legacy ArialMT to prevent kerning glitches
BASE_FONT <- "Liberation Sans"
avail <- unique(systemfonts::system_fonts()$family)
if (!BASE_FONT %in% avail) {
  if ("DejaVu Sans" %in% avail)      BASE_FONT <- "DejaVu Sans"
  else if ("Arial" %in% avail)       BASE_FONT <- "Arial"
}
message("Using font family: ", BASE_FONT)

# Ensure outlines (not dingbats) in PDFs
options(ggplot2.useDingbats = FALSE)

FIG_W_IN <- 6.45
FIG_H_IN <- 2.00

## ---- 3) Appearance knobs (colors, widths, spacing) -------------------------------------
AXIS_LINE_WIDTH <- 0.2
AXIS_LINE_COLOR <- "black"
PANEL_BORDER    <- AXIS_LINE_WIDTH

# Pure red / blue and neutrals
COL_UP    <- "#FF0000"
COL_DOWN  <- "#0000FF"
COL_NS    <- "#BFBFBF"
COL_UP_TX <- "#B00000"
COL_DN_TX <- "#0017B0"

THRESH_COL <- "grey45"
LFC_CUT    <- 0.14

BASE_SIZE        <- 6
AXIS_TEXT_PT     <- BASE_SIZE - 1
PANEL_SPACING_MM <- 1.8

# Facet strip styling
STRIP_SIZE_PT <- 6
STRIP_LH     <- 1.05
STRIP_FILL   <- "grey90"

# Points & legend
PT_SIZE <- 0.25; PT_STROKE <- 0.08
LEG_PT  <- 2.8

# Gene label aesthetics
GENE_LABEL_SIZE <- 1.1
LABEL_BORDER    <- 0.08
LABEL_PAD_LINES <- 0.04
LABEL_R_LINES   <- 0.06

set.seed(123)

## ---- 4) Input files ---------------------------------------------------------------------
in_neuropil    <- here('processed-data/image_processing/revisions/neuropil-dx_DEG-GM.csv')
in_neuronal    <- here('processed-data/image_processing/revisions/neun-dx_DEG-GM.csv')
in_pnn         <- here('processed-data/image_processing/revisions/pnn-dx_DEG-GM.csv')
in_vasculature <- here('processed-data/image_processing/revisions/vasc-dx_DEG.csv')

## ---- 5) Loader & preprocessing ----------------------------------------------------------
load_deg <- function(path) {
  readr::read_csv(path, show_col_types = FALSE) %>%
    mutate(
      p_value_scz = pmax(p_value_scz, .Machine$double.xmin),
      neglog10p   = -log10(p_value_scz),
      sig_FDR10   = fdr_scz < 0.10,
      sig_nominal = p_value_scz < 0.05,
      sig_category = case_when(
        sig_FDR10 & logFC_scz >  0 ~ "FDR <0.10 Up",
        sig_FDR10 & logFC_scz <= 0 ~ "FDR <0.10 Down",
        sig_nominal & !sig_FDR10   ~ "Nominal p <0.05",
        TRUE                       ~ "Not significant"
      ),
      sig_category = factor(
        sig_category,
        levels = c("FDR <0.10 Up","FDR <0.10 Down","Nominal p <0.05","Not significant")
      )
    ) %>%
    filter(is.finite(logFC_scz), is.finite(neglog10p))
}

neu  <- load_deg(in_neuropil)    %>% mutate(.panel = "Neuropil")
pnn  <- load_deg(in_pnn)         %>% mutate(.panel = "PNN")
neur <- load_deg(in_neuronal)    %>% mutate(.panel = "Neuronal")
vasc <- load_deg(in_vasculature) %>% mutate(.panel = "Vasculature")

all_df <- bind_rows(neu, pnn, neur, vasc) %>%
  mutate(.panel = factor(.panel, levels = c("Neuropil","PNN","Neuronal","Vasculature")))

## ---- 6) Strip labels (full-info vs title-only; toggle via flag) -------------------------
# Compute per-panel counts for the "full-info" strip (kept for easy switching)
panel_counts <- all_df %>%
  group_by(.panel) %>%
  summarise(
    N_total   = n(),
    N_FDR     = sum(sig_FDR10, na.rm = TRUE),
    N_up      = sum(sig_FDR10 & logFC_scz > 0,  na.rm = TRUE),
    N_down    = sum(sig_FDR10 & logFC_scz <= 0, na.rm = TRUE),
    N_nominal = sum(sig_nominal,   na.rm = TRUE),
    .groups = "drop"
  ) %>%
  # Use "up/down" words instead of arrows to avoid glyph fallback in some fonts
  mutate(strip_lab = glue("{.panel}\nN={N_total}, FDR<0.1: {N_FDR} ({N_up} up, {N_down} down); p<0.05: {N_nominal}"))

# Map for full strip vs simple strip
strip_map_full   <- setNames(panel_counts$strip_lab, panel_counts$.panel)
strip_map_simple <- setNames(panel_counts$.panel,     panel_counts$.panel)

# >>> Toggle here <<<
USE_SIMPLE_STRIP <- TRUE  # TRUE: title-only (recommended by default); FALSE: show counts

## ---- 7) Label sets per panel -------------------------------------------------------------
labels_neuropil <- c(
  "VGF","CORT","SST","BDNF","SNAP25","GAD1","GAD2","PVALB","NPY","CCK",
  "CRH","CRHBP","CACNA1B","CAMK2A","NRN1","GFAP","GAP43","ALDH1A1","NRXN2","NTRK2"
)
labels_neuronal <- c(
  "MAPK3","NOS1","NOS2","GRIN1","CORT","SST","CRH","CRHBP","PVALB","GAD1",
  "GRIK1","GRIK4","GRM7","PTPRS","SLITRK1","GFAP"
)
labels_pnn <- c(
  "NTRK2","GABRA3","GABRA5","GABRA1","GABRG2","CNTN2","KCNA1","NRXN3","GRM5",
  "SYNJ2","SPP1","GFAP"
)
labels_vasc <- c(
  "UBB","FKBP5","MAG","MBP","PLP1","GRIN2B","GRIN1","FGFR1","EPHA4","CEBPD",
  "TF","CLDN11","ABCG2","GFAP","CNTNAP1","C3","PSEN1","NRXN2","A2M"
)

label_tbl <- bind_rows(
  tibble(.panel = "Neuropil",    gene = labels_neuropil),
  tibble(.panel = "PNN",         gene = labels_pnn),
  tibble(.panel = "Neuronal",    gene = labels_neuronal),
  tibble(.panel = "Vasculature", gene = labels_vasc)
)

lab_df <- all_df %>%
  semi_join(label_tbl, by = c(".panel","gene")) %>%
  mutate(
    label_col = case_when(
      logFC_scz >  0 ~ COL_UP_TX,
      logFC_scz <  0 ~ COL_DN_TX,
      TRUE           ~ "#5F5F5F"
    )
  )

## ---- 8) Shared limits & axis breaks ------------------------------------------------------
xr       <- range(all_df$logFC_scz, na.rm = TRUE)
xmax     <- max(1.0, max(abs(xr)))     # ensure ±1.0 visible
xlim     <- c(-xmax, xmax)
x_breaks <- c(-1.0, -0.5, 0, 0.5, 1.0)
y_upper  <- max(all_df$neglog10p, na.rm = TRUE) + 0.25

## ---- 8.5) Significance grouping for aesthetics ------------------------------------------
all_df <- all_df %>%
  mutate(
    sig_group = case_when(
      sig_FDR10 & logFC_scz >  0 ~ "FDR <0.10 Up",
      sig_FDR10 & logFC_scz <= 0 ~ "FDR <0.10 Down",
      !sig_FDR10 & sig_nominal & logFC_scz >  0 ~ "Nominal p <0.05 Up",
      !sig_FDR10 & sig_nominal & logFC_scz <= 0 ~ "Nominal p <0.05 Down",
      TRUE ~ "Not significant"
    ),
    sig_group = factor(
      sig_group,
      levels = c("Not significant",
                 "Nominal p <0.05 Down","Nominal p <0.05 Up",
                 "FDR <0.10 Down","FDR <0.10 Up")
    )
  )

## ---- 9) Plot (no FDR dashed line; keep nominal p=0.05 & +/-LFC lines) -------------------
legend_levels <- c(
  "FDR <0.10 Up","FDR <0.10 Down",
  "Nominal p <0.05 Up","Nominal p <0.05 Down",
  "Not significant"
)
legend_colors <- c(
  "FDR <0.10 Up" = COL_UP,
  "FDR <0.10 Down" = COL_DOWN,
  "Nominal p <0.05 Up" = COL_UP,
  "Nominal p <0.05 Down" = COL_DOWN,
  "Not significant" = COL_NS
)
legend_alpha  <- c(
  "FDR <0.10 Up" = 0.95,
  "FDR <0.10 Down" = 0.95,
  "Nominal p <0.05 Up" = 0.35,
  "Nominal p <0.05 Down" = 0.35,
  "Not significant" = 0.12
)

p <- ggplot(all_df, aes(x = logFC_scz, y = neglog10p)) +
  # Reference lines (keep nominal p=0.05 & +/- LFC; remove any FDR-based horizontal guide)
  geom_hline(yintercept = -log10(0.05), linetype = "solid",  linewidth = 0.28, colour = THRESH_COL) +
  geom_vline(xintercept = c(-LFC_CUT, LFC_CUT), linetype = "dashed", linewidth = 0.24, colour = THRESH_COL) +
  # Points
  geom_point(
    data = all_df %>% arrange(sig_group),
    aes(fill = sig_group, colour = sig_group, alpha = sig_group),
    shape = 21, size = PT_SIZE, stroke = PT_STROKE, show.legend = TRUE
  ) +
  # Labels
  geom_label_repel(
    data = lab_df, aes(label = gene, colour = I(label_col)),
    family = BASE_FONT, size = GENE_LABEL_SIZE, label.size = LABEL_BORDER,
    label.padding = unit(LABEL_PAD_LINES, "lines"), label.r = unit(LABEL_R_LINES, "lines"),
    max.overlaps = Inf, force = 0.8, direction = "both",
    segment.size = 0.20, segment.color = "grey25", segment.alpha = 0.9,
    box.padding = 0.12, point.padding = 0.12, show.legend = FALSE,
    fill = alpha("white", 0.96), seed = 123
  ) +
  # Scales
  scale_fill_manual(values = legend_colors, breaks = legend_levels, name = "Significance") +
  scale_colour_manual(values = legend_colors, breaks = legend_levels, name = "Significance") +
  scale_alpha_manual(values = legend_alpha,  breaks = legend_levels,  name = "Significance") +
  scale_x_continuous(
    limits = xlim, breaks = x_breaks,
    labels = scales::number_format(accuracy = 0.1),
    expand = expansion(mult = c(0.03, 0.03))
  ) +
  scale_y_continuous(limits = c(0, y_upper), expand = expansion(mult = c(0, 0))) +
  coord_cartesian(clip = "off") +
  labs(x = "log2(Fold Change)", y = expression(-log[10](italic(p)))) +
  # Facets with toggleable strip labels
  facet_wrap(
    ~ .panel, nrow = 1,
    labeller = as_labeller(if (USE_SIMPLE_STRIP) strip_map_simple else strip_map_full)
  ) +
  # Theme
  theme_minimal(base_size = BASE_SIZE, base_family = BASE_FONT) +
  theme(
    panel.grid.minor   = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.border       = element_rect(colour = AXIS_LINE_COLOR, fill = NA, linewidth = PANEL_BORDER),
    panel.spacing.x    = unit(PANEL_SPACING_MM, "mm"),
    axis.line.y         = element_blank(),
    axis.ticks.y        = element_line(colour = AXIS_LINE_COLOR, linewidth = AXIS_LINE_WIDTH),
    axis.ticks.length.y = unit(0.6, "mm"),
    axis.ticks.x        = element_line(colour = AXIS_LINE_COLOR, linewidth = AXIS_LINE_WIDTH),
    axis.title.y = element_text(margin = margin(r = 1.0, unit = "mm")),
    axis.text.y  = element_text(size = AXIS_TEXT_PT),
    axis.title.x = element_text(margin = margin(t = 1.0, unit = "mm")),
    axis.text.x  = element_text(size = AXIS_TEXT_PT),
    legend.position    = "bottom",
    legend.key.height  = unit(2.8, "mm"),
    legend.key.width   = unit(7.0, "mm"),
    strip.placement    = "inside",
    strip.background.x = element_rect(fill = STRIP_FILL, color = AXIS_LINE_COLOR, linewidth = AXIS_LINE_WIDTH),
    strip.text.x       = element_text(size = STRIP_SIZE_PT, lineheight = STRIP_LH,
                                      family = BASE_FONT, color = "black")
  ) +
  guides(
    fill   = guide_legend(title.position = "top", nrow = 1, override.aes = list(size = LEG_PT, stroke = 0.3)),
    alpha  = guide_legend(title.position = "top", nrow = 1),
    colour = "none"
  )

## ---- 10) Export --------------------------------------------------------------------------
# PDF (Cairo backend for robust text rendering/embedding)
ggsave(filename = pdf_out, plot = p,
       width = FIG_W_IN, height = FIG_H_IN, units = "in",
       device = cairo_pdf, dpi = 600, limitsize = FALSE, bg = "white")

# PNG/TIFF via ragg (consistent rasterization across platforms)
ggsave(filename = png_out, plot = p,
       width = FIG_W_IN, height = FIG_H_IN, units = "in",
       device = ragg::agg_png, dpi = 600, limitsize = FALSE, bg = "white")


## ---- 11) Session info --------------------------------------------------------------------
writeLines(capture.output(sessioninfo::session_info()), con = sess_out)

message("Done! Exported (no FDR dashed line) to:\n  - ", pdf_out, "\n  - ", png_out, "\n  - ", tif_out,
        if (requireNamespace('svglite', quietly = TRUE)) paste0("\n  - ", svg_out) else "",
        "\nSession info: ", sess_out,
        "\nFolder: ", out_dir)
