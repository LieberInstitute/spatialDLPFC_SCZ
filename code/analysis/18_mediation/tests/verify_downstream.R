## independently reconstruct GO universes, term memberships, hypergeometric and BH results.
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) %in% c(1L, 2L))
run <- args[1]
mode <- if (length(args) == 2L) args[2] else "eligible"
stopifnot(mode %in% c("eligible", "all_tested"))
out <- file.path(run, if (mode == "eligible") "downstream_nominal" else "downstream_nominal_all_tested")
r <- read.delim(gzfile(file.path(out, "GO_results_all.tsv.gz")), check.names = FALSE)
source <- read.delim(gzfile(file.path(run, "primary", "all_pairs.tsv.gz")))
hits <- read.csv(file.path(run, "report", "mediation_hits_nominal_p05.csv"))
key <- function(s, g) paste(s, g, sep = ":")
pass <- with(source, historical_p < .05 & c_p < .05 & b_p < .05 & cprime_p >= .05)
stopifnot(setequal(key(source$screen_id[pass], source$gene_id[pass]), key(hits$screen, hits$outcome_gene_id)))
stopifnot(nrow(hits) == 110L, length(unique(hits$outcome_gene_id)) == 90L)
expected <- with(r, phyper(Count - 1, term_background,
                           background_annotated - term_background,
                           query_annotated, lower.tail = FALSE))
stopifnot(isTRUE(all.equal(r$pvalue, expected, tolerance = 1e-12)))
families <- split(seq_len(nrow(r)), paste(r$screen, r$ontology))
for (ix in families) {
    stopifnot(isTRUE(all.equal(r$p_adjust[ix], p.adjust(r$pvalue[ix], "BH"), tolerance = 1e-12)))
}
stopifnot(isTRUE(all.equal(r$p_adjust_all_screens_ontologies,
                          p.adjust(r$pvalue, "BH"), tolerance = 1e-12)))
membership <- read.delim(file.path(out, "gene_overlap_membership.tsv"), check.names = FALSE)
screens <- setdiff(names(membership), c("gene_id", "gene_symbol"))
stopifnot(!anyDuplicated(membership$gene_id), sum(as.matrix(membership[screens])) == nrow(hits))
for (s in screens) {
    stopifnot(setequal(membership$gene_id[membership[[s]] == 1], hits$outcome_gene_id[hits$screen == s]))
}
summary <- read.delim(file.path(out, "GO_group_summary.tsv"))
universes <- read.delim(file.path(out, "GO_eligible_universes.tsv"))
mapping <- read.delim(file.path(out, "GO_gene_mapping.tsv"), colClasses = "character")
annotation <- read.delim(gzfile(file.path(out, "GO_annotation_snapshot.tsv.gz")), colClasses = "character")
term_sets <- lapply(split(annotation, annotation$ontology), function(a) split(a$entrez_id, a$go_id))
for (i in seq_len(nrow(summary))) {
    screen <- summary$screen[i]
    ont <- summary$ontology[i]
    ix <- r$screen == screen & r$ontology == ont
    rr <- r[ix, ]
    eligible <- source$screen_id == screen
    if (mode == "eligible") eligible <- eligible & source$historical_p < .05 & source$c_p < .05
    genes <- unique(source$gene_id[eligible])
    stopifnot(setequal(genes, universes$gene_id[universes$screen == screen]), length(genes) == summary$background_ensembl[i])
    universe <- unique(mapping$entrez_id[mapping$gene_id %in% genes])
    selected <- hits$outcome_gene_id[hits$screen == screen]
    candidates <- unique(mapping$entrez_id[mapping$gene_id %in% selected])
    sets <- term_sets[[ont]]
    background <- intersect(universe, unique(unlist(sets, use.names = FALSE)))
    query <- intersect(candidates, background)
    sizes <- vapply(sets, function(x) length(intersect(x, background)), integer(1))
    eligible_terms <- names(sizes)[sizes >= 10 & sizes <= 500]
    overlaps <- lapply(sets[eligible_terms], intersect, y = query)
    overlaps <- overlaps[lengths(overlaps) > 0]
    stopifnot(setequal(names(overlaps), rr$go_id), length(eligible_terms) == summary$all_eligible_terms[i])
    stopifnot(all(rr$background_annotated == length(background)), all(rr$query_annotated == length(query)))
    stopifnot(all(rr$term_background == sizes[rr$go_id]), all(rr$Count == lengths(overlaps[rr$go_id])))
    for (j in seq_len(nrow(rr))) stopifnot(setequal(strsplit(rr$entrez_ids[j], ";", fixed = TRUE)[[1]], overlaps[[rr$go_id[j]]]))
    stopifnot(isTRUE(all.equal(rr$p_adjust_all_eligible_terms, p.adjust(rr$pvalue, "BH", n = length(eligible_terms)), tolerance = 1e-12)))
    stopifnot(sum(ix) == summary$tested_terms[i], sum(rr$p_adjust < .05) == summary$BH05_terms[i])
    stopifnot(length(background) == summary$background_annotated[i], length(query) == summary$query_annotated[i])
}
cat("Verified", mode, ":", nrow(r), "hypergeometric tests,", length(families),
    "BH families, pooled and expanded corrections, input gates, overlap, GO memberships and denominators.\n")
