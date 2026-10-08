## vascular PTN SCZ-NTC effect per SpD, with the manuscript SPG model and all seven SpDs.
## usage: Rscript ptn_by_spd.R [project_root] [outdir]
args <- commandArgs(trailingOnly = TRUE)
root <- if (length(args) >= 1) args[1] else "."
out <- if (length(args) >= 2) args[2] else file.path(root, "processed-data", "18_mediation")
suppressPackageStartupMessages(library(SummarizedExperiment))
x <- readRDS(file.path(root, "processed-data/rds/PB_dx_spg/pseudo_vasc_pos_donor_spd.rds"))
cd <- as.data.frame(colData(x)); E <- assay(x, "logcounts")
g <- which(sub("[.].*", "", rowData(x)$gene_id) == "ENSG00000105894")
d <- data.frame(dx = factor(toupper(cd$dx), c("NTC", "SCZ")), age = as.numeric(cd$age),
                sex = factor(cd$sex), spd = factor(cd$registration_variable), donor = cd$brnum)
## donor correlation as in spatialLIBD::registration_block_cor (covariate design without SpD).
rho <- limma::duplicateCorrelation(E, model.matrix(~ dx + age + sex, d), block = d$donor)$consensus.correlation
des <- cbind(model.matrix(~ 0 + spd:dx, d), age = d$age, sexM = as.numeric(d$sex == "M"))
colnames(des) <- make.names(colnames(des))
stopifnot(qr(des)$rank == ncol(des))
fit <- limma::lmFit(E, des, block = d$donor, correlation = rho)
contrast <- function(w) {
  v <- setNames(rep(0, ncol(des)), colnames(des))
  for (s in names(w)) {
    v[paste0("spd", s, ".dxSCZ")] <- w[[s]]; v[paste0("spd", s, ".dxNTC")] <- -w[[s]]
  }
  v
}
gm <- c("spd02", "spd03", "spd05", "spd06", "spd07"); wm <- c("spd01", "spd04")
L <- lapply(setNames(levels(d$spd), levels(d$spd)), function(s) contrast(setNames(list(1), s)))
L$GM_mean <- contrast(setNames(as.list(rep(1 / 5, 5)), gm))
L$WM_mean <- contrast(setNames(as.list(rep(1 / 2, 2)), wm))
L$WM_minus_GM <- L$WM_mean - L$GM_mean
fc <- limma::eBayes(limma::contrasts.fit(fit, do.call(cbind, L)))
n <- table(d$spd, d$dx)
res <- data.frame(contrast = names(L), logFC = fc$coefficients[g, ], p = fc$p.value[g, ],
                  n_NTC = c(n[, "NTC"], NA, NA, NA), n_SCZ = c(n[, "SCZ"], NA, NA, NA))
dir.create(out, recursive = TRUE, showWarnings = FALSE)
write.table(res, file.path(out, "ptn_vasc_by_spd.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
print(res, row.names = FALSE, digits = 3)
