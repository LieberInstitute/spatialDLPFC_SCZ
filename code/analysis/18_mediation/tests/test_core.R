## deterministic contract and statistical-behavior checks; run under conda_R/4.5.
arg <- grep("^--file=", commandArgs(), value = TRUE)[1]
base <- dirname(dirname(normalizePath(sub("^--file=", "", arg))))
source(file.path(base, "R", "core.R")); source(file.path(base, "config.R"))
expect_error <- function(f) assert(inherits(try(f(), silent = TRUE), "try-error"), "Expected error was not raised")
config$min_donors <- 4L; config$min_per_dx <- 2L
n <- 40L; d <- expand.grid(SpD = c("spd02", "spd03"), donor = paste0("d", seq_len(n)), stringsAsFactors = FALSE)
d$Dx <- ifelse(as.integer(sub("d", "", d$donor)) <= n / 2, "NTC", "SCZ")
d$age <- rep(seq(30, 69), each = 2); d$sex <- rep(rep(c("M", "F"), 20), each = 2)
set.seed(2718); d$age <- rep(runif(n, 30, 70), each = 2)
d$section <- paste0("section_", d$donor); d$slide_id <- rep(rep(c("A", "B", "C", "D"), 10), each = 2)
d$rin <- rep(runif(n, 6, 9), each = 2); d$ncells <- 30L; d$key <- paste(d$donor, d$SpD, sep = "|")
o <- list(samples = d)
z <- match_samples(o, o, config)
shuffled <- o; shuffled$samples <- d[sample(nrow(d)), ]
assert(identical(z$samples$key, match_samples(shuffled, o, config)$samples$key), "Matching depends on row order")
bad <- o; bad$samples$Dx[1] <- "SCZ"
expect_error(function() match_samples(bad, o, config))
expect_error(function() unique_keys(c("a", "a"), "test"))
bad_d <- d; bad_d$age[1] <- NA_real_; expect_error(function() make_design(bad_d, config))
bad_d <- d; bad_d$M <- 1; expect_error(function() make_design(bad_d, config, "M"))
base <- data.frame(historical_p = rep(.01, 4), c_p = rep(.01, 4), a_p = c(.01, .01, .2, .01),
 b_q = rep(.01, 4), cprime_p = c(.1, .01, .1, .1), c_beta = rep(1, 4), cprime_beta = c(1.2, .5, .5, -.5), b_q_global = rep(.01, 4))
strict <- config; strict$require_mediator_gate <- TRUE
z <- classify_pairs(base, strict)
assert(classify_pairs(base, config)$screen_hit[3], "Optional mediator gate still enforced")
assert(z$screen_hit[1] && !z$coefficient_shrinkage[1] && !z$higher_priority[1], "Significance loss incorrectly implies shrinkage")
assert(z$attenuated_still_significant[2] && !z$screen_hit[2], "Retained significance classification wrong")
assert(!z$screen_hit[3] && !z$mediator_gate[3], "Failed mediator gate ignored")
assert(z$screen_hit[4] && z$direction_reversal[4] && !z$higher_priority[4], "Sign reversal incorrectly prioritized")
p <- c(.0001, .02, rep(.9, 98)); assert(p.adjust(p, "BH")[2] > .05, "Full-family FDR test failed")
assert(!identical(hash(d), hash(d[-1, ])), "Sample change did not change signature")
assert(!identical(hash(make_design(d, config)), hash(make_design(d, config, "rin"))), "Formula change did not change signature")
## synthetic mediator signal: fit all genes and verify known coefficient direction and shrinkage.
d <- z <- o$samples
set.seed(1234)
d$M <- rep(rnorm(n), each = 2) + 1.2 * (d$Dx == "SCZ") + rnorm(nrow(d), 0, .15)
ng <- 400L
mu <- matrix(exp(6 + rnorm(ng, 0, .3)), ng, nrow(d))
mu[1:8, ] <- mu[1:8, ] * matrix(exp(.65 * d$M), 8, nrow(d), byrow = TRUE)
counts <- matrix(rnbinom(length(mu), mu = mu, size = 60), ng)
rownames(counts) <- paste0("gene", seq_len(ng)); colnames(counts) <- d$key
obj <- list(counts = counts)
des <- make_design(d, config); dge <- prepare_dge(obj, d$key, des)
f0 <- fit_voom(dge, des, d); f1 <- fit_voom(dge, make_design(d, config, "M"), d)
t0 <- coef_table(f0); t1 <- coef_table(f1); tb <- coef_table(f1, "M")
ix <- match(paste0("gene", 1:8), t0$gene_id)
assert(all(t0$beta[ix] > 0) && all(tb$beta[ix] > 0), "SCZ or mediator coefficient direction wrong")
assert(mean(abs(t1$beta[ix])) < mean(abs(t0$beta[ix])), "Known mediator signal did not attenuate diagnosis coefficients")
assert(nrow(tb) == nrow(dge) && isTRUE(all.equal(tb$q, p.adjust(tb$p, "BH"))), "M-Y FDR family mismatch")
assert(identical(rownames(f0$coefficients), rownames(f1$coefficients)), "Gene universe changed in nested models")
## an independent negative-control mediator must not systematically attenuate the signal.
d$M_null <- rnorm(nrow(d))
fnull <- fit_voom(dge, make_design(d, config, "M_null"), d)
tnull <- coef_table(fnull)
assert(mean(abs(tnull$beta[ix] - t0$beta[ix])) < mean(abs(t1$beta[ix] - t0$beta[ix])), "Null mediator attenuated as much as true signal")
cat("All contract, classification, multiplicity, and synthetic-model checks passed.\n")
