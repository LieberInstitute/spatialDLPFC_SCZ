## export saved manuscript Dx statistics; no models are fitted.
## usage: Rscript code/analysis/export_spg_dx_statistics.R [output.xlsx]
spg_dx_tables <- function(root) {
  sources <- c(
    Vascular = "processed-data/spg_pb_de/test_SPD_pseudo_vasc_pos.csv",
    Neuropil = "code/analysis/dx_deg_spg_neuropil/neuropil-dx_DEG-GM.csv",
    Neuronal = "code/analysis/dx_deg_spg_neun/neun-dx_DEG-GM.csv"
  )
  columns <- c(gene_name = "gene", gene_id = "ensembl",
               logFC_DE = "logFC_scz", t_DE = "t_stat_scz",
               P.Value_DE = "p_value_scz", fdr_DE = "fdr_scz")
  tables <- lapply(sources, function(path) {
    x <- read.csv(file.path(root, path), check.names = FALSE,
                  stringsAsFactors = FALSE)
    stopifnot(all(columns %in% names(x)), nrow(x) > 0,
              !anyNA(x$ensembl), !anyDuplicated(x$ensembl))
    x <- x[, unname(columns)]
    names(x) <- names(columns)
    stopifnot(all(vapply(x[3:6], is.numeric, logical(1))),
              all(is.finite(as.matrix(x[3:6]))),
              all(x$P.Value_DE >= 0 & x$P.Value_DE <= 1),
              all(x$fdr_DE >= 0 & x$fdr_DE <= 1))
    x[order(x$gene_name, x$gene_id, na.last = TRUE), ]
  })
  tables$Notes <- data.frame(
    Field = c("Scope", "Lookup", "Missing gene", "gene_name", "gene_id",
              "logFC_DE", "t_DE", "P.Value_DE", "fdr_DE", "Model",
              "Vascular domains", "Neuropil and neuronal domains", names(sources)),
    Description = c(
      "Original manuscript diagnosis statistics. All tested genes retained; no significance filter.",
      "Use the gene_name column filter or Excel Find to look up receptors.",
      "A gene absent from a sheet was not present in that source result table; absence is not a zero effect.",
      "Gene symbol, copied from the source gene column.",
      "Ensembl gene identifier, copied from the source ensembl column.",
      "Log2 fold change, SCZ minus NTC. Positive: higher in SCZ; negative: lower in SCZ.",
      "Moderated t-statistic for SCZ minus NTC.",
      "Unadjusted diagnosis p-value, copied from p_value_scz.",
      "Benjamini-Hochberg FDR across the complete gene family in the source table, copied unchanged from fdr_scz.",
      "Manuscript limma on logcounts, with age, sex, SpD adjustment and donor blocking.",
      "All seven SpDs, including white matter.",
      "Gray matter: SpD02, SpD03, SpD05, SpD06 and SpD07.",
      unname(sources)), stringsAsFactors = FALSE)
  tables
}

export_spg_dx_statistics <- function(root, output) {
  stopifnot(requireNamespace("openxlsx", quietly = TRUE))
  tables <- spg_dx_tables(root)
  wb <- openxlsx::createWorkbook()
  header <- openxlsx::createStyle(fontColour = "#FFFFFF", fgFill = "#24364B",
                                  textDecoration = "bold", halign = "center")
  for (sheet in names(tables)) {
    x <- tables[[sheet]]
    openxlsx::addWorksheet(wb, sheet, gridLines = FALSE)
    openxlsx::writeDataTable(wb, sheet, x, tableStyle = "TableStyleMedium2")
    openxlsx::addStyle(wb, sheet, header, rows = 1, cols = seq_len(ncol(x)), gridExpand = TRUE)
    openxlsx::freezePane(wb, sheet, firstActiveRow = 2, firstActiveCol = 3)
    openxlsx::setColWidths(wb, sheet, cols = seq_len(ncol(x)),
                          widths = if (sheet == "Notes") c(32, 120) else c(25, 25, 18, 18, 18, 18))
    if (sheet != "Notes") {
      for (cols in list(3:4, 5:6)) {
        fmt <- if (cols[1] == 3) "0.0000" else "0.000E+00"
        openxlsx::addStyle(wb, sheet, openxlsx::createStyle(numFmt = fmt),
                          rows = 2:(nrow(x) + 1), cols = cols, gridExpand = TRUE)
      }
    }
  }
  dir.create(dirname(output), recursive = TRUE, showWarnings = FALSE)
  openxlsx::saveWorkbook(wb, output, overwrite = TRUE)
  message("Exported: ", output)
}

if (sys.nframe() == 0L) {
  script <- sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE)[1])
  root <- normalizePath(file.path(dirname(script), "../.."))
  args <- commandArgs(trailingOnly = TRUE)
  stopifnot(length(args) <= 1)
  output <- if (length(args)) args[1] else file.path(root, "processed-data/spg_pb_de/SPG_Dx_statistics_unfiltered.xlsx")
  export_spg_dx_statistics(root, output)
}
