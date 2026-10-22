## export model-consistent plotting values from the unchanged fitted artifacts.
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 1L)
run <- args[1]
options(digits = 17)
hits <- read.csv(file.path(run, "report", "mediation_hits_nominal_p05.csv"), check.names = FALSE)
results <- read.delim(gzfile(file.path(run, "primary", "all_pairs.tsv.gz")), check.names = FALSE)
stopifnot(nrow(hits) == 110L, length(unique(hits$outcome_gene_id)) == 90L)
key <- function(s, g) paste(s, g, sep = ":")
selected <- with(results, historical_p < .05 & c_p < .05 & b_p < .05 & cprime_p >= .05)
stopifnot(setequal(key(results$screen_id[selected], results$gene_id[selected]), key(hits$screen, hits$outcome_gene_id)))
screens <- unique(hits$screen)
donors <- sort(unique(unlist(lapply(screens, function(s) readRDS(file.path(run, "primary", paste0(s, ".rds")))$samples$donor))))
donor_codes <- setNames(sprintf("D%03d", seq_along(donors)), donors)
clean <- function(y, fit) {
    nuisance <- setdiff(colnames(fit$design), c("DxNTC", "DxSCZ"))
    x <- fit$design[, nuisance, drop = FALSE]
    centered <- sweep(x, 2L, colMeans(x))
    y - fit$coefficients[rownames(y), nuisance, drop = FALSE] %*% t(centered)
}
verify_gls <- function(y, fit, donor, rho) {
    covariance <- outer(donor, donor, "==") * rho
    diag(covariance) <- 1
    x <- fit$design
    vinvx <- solve(covariance, x)
    reconstructed <- t(solve(crossprod(x, vinvx), crossprod(vinvx, t(y))))
    expected <- fit$coefficients[rownames(y), , drop = FALSE]
    stopifnot(max(abs(reconstructed - expected)) < 1e-8)
    max(abs(reconstructed - expected))
}
rows <- list(); checks <- list()
for (s in screens) {
    a <- readRDS(file.path(run, "primary", paste0(s, ".rds")))
    h <- hits[hits$screen == s, ]
    rr <- results[match(key(s, h$outcome_gene_id), key(results$screen_id, results$gene_id)), ]
    stopifnot(identical(a$signature, unique(rr$run_signature)), !anyDuplicated(a$samples$key))
    stopifnot(identical(rownames(a$design), a$samples$key), identical(colnames(a$tdge), a$samples$key),
              identical(colnames(a$sdge), a$samples$key))
    b <- unclass(a$baseline); j <- unclass(a$joint); m <- unclass(a$source_fit)
    stopifnot(identical(b$design, a$design), identical(m$design, a$design),
              identical(rownames(j$design), a$samples$key))
    y <- a$tdge[h$outcome_gene_id, , drop = FALSE]
    med <- a$sdge[a$screen$mediator_id, , drop = FALSE]
    stopifnot(max(abs(as.numeric(med) - a$samples$M_logCPM)) < 1e-12)
    stopifnot(max(abs(j$design[, "M"] - a$samples$M)) < 1e-12)
    stopifnot(max(abs(as.numeric(scale(as.numeric(med))) - a$samples$M)) < 1e-12)
    errors <- c(verify_gls(y, b, a$samples$donor, a$rho[["target"]]),
                verify_gls(y, j, a$samples$donor, a$rho[["target"]]),
                verify_gls(med, m, a$samples$donor, a$rho[["source"]]))
    before <- clean(y, b); after <- clean(y, j); mediator <- as.numeric(clean(med, m))
    stopifnot(max(abs(rowMeans(before) - rowMeans(y))) < 1e-12,
              max(abs(rowMeans(after) - rowMeans(y))) < 1e-12,
              abs(mean(mediator) - mean(med)) < 1e-12)
    contrast <- function(fit, ids) fit$coefficients[ids, "DxSCZ"] - fit$coefficients[ids, "DxNTC"]
    stopifnot(max(abs(contrast(b, h$outcome_gene_id) - rr$c_beta)) < 1e-12,
              max(abs(contrast(j, h$outcome_gene_id) - rr$cprime_beta)) < 1e-12,
              max(abs(j$coefficients[h$outcome_gene_id, "M"] - rr$b_beta)) < 1e-12,
              max(abs(contrast(m, a$screen$mediator_id) - rr$a_beta)) < 1e-12)
    stopifnot(ncol(y) == unique(h$n_observations), length(unique(a$samples$donor)) == unique(h$n_donors))
    for (i in seq_len(nrow(h))) {
        rows[[length(rows) + 1L]] <- data.frame(screen = s, gene_id = h$outcome_gene_id[i],
            donor = unname(donor_codes[a$samples$donor]), Dx = a$samples$Dx, SpD = a$samples$SpD,
            outcome_baseline = before[i, ], outcome_adjusted = after[i, ], mediator_adjusted = mediator)
    }
    checks[[s]] <- data.frame(screen = s, pairs = nrow(h), observations = ncol(y),
        donors = length(unique(a$samples$donor)), source_rho = a$rho[["source"]], target_rho = a$rho[["target"]],
        max_GLS_coefficient_error = max(errors), artifact_md5 = unname(tools::md5sum(file.path(run, "primary", paste0(s, ".rds")))))
}
out <- file.path(run, "candidate_boxplots")
dir.create(out, showWarnings = FALSE)
con <- gzfile(file.path(out, "plot_values.tsv.gz"), "wt")
write.table(do.call(rbind, rows), con, sep = "\t", quote = FALSE, row.names = FALSE)
close(con)
write.table(do.call(rbind, checks), file.path(out, "model_verification.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
cat("Exported", nrow(hits), "pairs; saved GLS coefficients independently reconstructed for all five screens.\n")
