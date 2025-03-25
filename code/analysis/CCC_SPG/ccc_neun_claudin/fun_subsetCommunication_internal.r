subsetCommunication_internal <- function (net, LR, cells.level, slot.name = "net", sources.use = NULL, 
    targets.use = NULL, signaling = NULL, pairLR.use = NULL, 
    thresh = 0.05, datasets = NULL, ligand.pvalues = NULL, ligand.logFC = NULL, 
    ligand.pct.1 = NULL, ligand.pct.2 = NULL, receptor.pvalues = NULL, 
    receptor.logFC = NULL, receptor.pct.1 = NULL, receptor.pct.2 = NULL) 
{
    if (!is.data.frame(net)) {
        prob <- net$prob
        pval <- net$pval
        prob[pval >= thresh] <- 0
        net <- reshape2::melt(prob, value.name = "prob")
        colnames(net)[1:3] <- c("source", "target", "interaction_name")
        net.pval <- reshape2::melt(pval, value.name = "pval")
        net$pval <- net.pval$pval
        net <- subset(net, prob > 0)
    }
    if (!("ligand" %in% colnames(net))) {
        col.use <- intersect(c("interaction_name_2", "pathway_name", 
            "ligand", "receptor", "annotation", "evidence"), 
            colnames(LR))
        pairLR <- dplyr::select(LR, col.use)
        idx <- match(net$interaction_name, rownames(pairLR))
        net <- cbind(net, pairLR[idx, ])
    }
    if (!is.null(signaling)) {
        pairLR.use <- data.frame()
        for (i in 1:length(signaling)) {
            pairLR.use.i <- searchPair(signaling = signaling[i], 
                pairLR.use = LR, key = "pathway_name", matching.exact = T, 
                pair.only = T)
            pairLR.use <- rbind(pairLR.use, pairLR.use.i)
        }
    }
    if (!is.null(pairLR.use)) {
        net <- tryCatch({
            subset(net, interaction_name %in% pairLR.use$interaction_name)
        }, error = function(e) {
            subset(net, pathway_name %in% pairLR.use$pathway_name)
        })
    }
    if (!is.null(datasets)) {
        if (!("datasets" %in% colnames(net))) {
            stop("Please run `identifyOverExpressedGenes` and `netMappingDEG` before selecting 'datasets'")
        }
        net <- net[net$datasets %in% datasets, , drop = FALSE]
    }
    if (!is.null(ligand.pvalues)) {
        if (!("ligand.pvalues" %in% colnames(net))) {
            stop("Please run `identifyOverExpressedGenes` and `netMappingDEG` before using the threshold 'ligand.pvalues'")
        }
        net <- net[net$ligand.pvalues <= ligand.pvalues, , drop = FALSE]
    }
    if (!is.null(ligand.logFC)) {
        if (!("ligand.logFC" %in% colnames(net))) {
            stop("Please run `identifyOverExpressedGenes` and `netMappingDEG` before using the threshold 'ligand.logFC'")
        }
        if (ligand.logFC >= 0) {
            net <- net[net$ligand.logFC >= ligand.logFC, , drop = FALSE]
        }
        else {
            net <- net[net$ligand.logFC <= ligand.logFC, , drop = FALSE]
        }
    }
    if (!is.null(ligand.pct.1)) {
        if (!("ligand.pct.1" %in% colnames(net))) {
            stop("Please run `identifyOverExpressedGenes` and `netMappingDEG` before using the threshold 'ligand.pct.1'")
        }
        net <- net[net$ligand.pct.1 >= ligand.pct.1, , drop = FALSE]
    }
    if (!is.null(ligand.pct.2)) {
        if (!("ligand.pct.2" %in% colnames(net))) {
            stop("Please run `identifyOverExpressedGenes` and `netMappingDEG` before using the threshold 'ligand.pct.2'")
        }
        net <- net[net$ligand.pct.2 >= ligand.pct.2, , drop = FALSE]
    }
    if (!is.null(receptor.pvalues)) {
        if (!("receptor.pvalues" %in% colnames(net))) {
            stop("Please run `identifyOverExpressedGenes` and `netMappingDEG` before using the threshold 'receptor.pvalues'")
        }
        net <- net[net$receptor.pvalues <= receptor.pvalues, 
            , drop = FALSE]
    }
    if (!is.null(receptor.logFC)) {
        if (!("receptor.logFC" %in% colnames(net))) {
            stop("Please run `identifyOverExpressedGenes` and `netMappingDEG` before using the threshold 'receptor.logFC'")
        }
        if (receptor.logFC >= 0) {
            net <- net[net$receptor.logFC >= receptor.logFC, 
                , drop = FALSE]
        }
        else {
            net <- net[net$receptor.logFC <= receptor.logFC, 
                , drop = FALSE]
        }
    }
    if (!is.null(receptor.pct.1)) {
        if (!("receptor.pct.1" %in% colnames(net))) {
            stop("Please run `identifyOverExpressedGenes` and `netMappingDEG` before using the threshold 'receptor.pct.1'")
        }
        net <- net[net$receptor.pct.1 >= receptor.pct.1, , drop = FALSE]
    }
    if (!is.null(receptor.pct.2)) {
        if (!("receptor.pct.2" %in% colnames(net))) {
            stop("Please run `identifyOverExpressedGenes` and `netMappingDEG` before using the threshold 'receptor.pct.2'")
        }
        net <- net[net$receptor.pct.2 >= receptor.pct.2, , drop = FALSE]
    }
    net <- net[rowSums(is.na(net)) != ncol(net), , drop = FALSE]
    if (nrow(net) == 0) {
        # browser()
        warning("No significant signaling interactions are inferred based on the input!")
    }
    if (slot.name == "netP") {
        col.use <- intersect(c("source", "target", "pathway_name", 
            "prob", "pval", "annotation"), colnames(net))
        net <- dplyr::select(net, col.use)
        net$source_target <- paste(net$source, net$target, sep = "sourceTotarget")
        net.pval <- net %>% group_by(source_target, pathway_name) %>% 
            summarize(pval = mean(pval), .groups = "drop")
        net <- net %>% group_by(source_target, pathway_name) %>% 
            summarize(prob = sum(prob), .groups = "drop")
        a <- stringr::str_split(net$source_target, "sourceTotarget", 
            simplify = T)
        net$source <- as.character(a[, 1])
        net$target <- as.character(a[, 2])
        net <- dplyr::select(net, -source_target)
        net$pval <- net.pval$pval
    }
    if (!is.null(sources.use)) {
        if (is.numeric(sources.use)) {
            sources.use <- cells.level[sources.use]
        }
        net <- subset(net, source %in% sources.use)
    }
    if (!is.null(targets.use)) {
        if (is.numeric(targets.use)) {
            targets.use <- cells.level[targets.use]
        }
        net <- subset(net, target %in% targets.use)
    }
    net <- BiocGenerics::as.data.frame(net, stringsAsFactors = FALSE)
    if (nrow(net) == 0) {
        warning("No significant signaling interactions are inferred!")
    }
    else {
        rownames(net) <- 1:nrow(net)
    }
    if (slot.name == "net") {
        if (("ligand.logFC" %in% colnames(net)) & ("datasets" %in% 
            colnames(net))) {
            col.use <- intersect(c("source", "target", "ligand", 
                "receptor", "prob", "pval", "interaction_name", 
                "interaction_name_2", "pathway_name", "annotation", 
                "evidence", "datasets", "ligand.logFC", "ligand.pct.1", 
                "ligand.pct.2", "ligand.pvalues", "receptor.logFC", 
                "receptor.pct.1", "receptor.pct.2", "receptor.pvalues"), 
                colnames(net))
            net <- net[, col.use]
        }
        else if ("ligand.logFC" %in% colnames(net)) {
            col.use <- intersect(c("source", "target", "ligand", 
                "receptor", "prob", "pval", "interaction_name", 
                "interaction_name_2", "pathway_name", "annotation", 
                "evidence", "ligand.logFC", "ligand.pct.1", "ligand.pct.2", 
                "ligand.pvalues", "receptor.logFC", "receptor.pct.1", 
                "receptor.pct.2", "receptor.pvalues"), colnames(net))
            net <- net[, col.use]
        }
        else {
            col.use <- intersect(c("source", "target", "ligand", 
                "receptor", "prob", "pval", "interaction_name", 
                "interaction_name_2", "pathway_name", "annotation", 
                "evidence"), colnames(net))
            net <- net[, col.use]
        }
    }
    else if (slot.name == "netP") {
        col.use <- intersect(c("source", "target", "pathway_name", 
            "prob", "pval"), colnames(net))
        net <- net[, col.use]
    }
    return(net)
}