IntegrateSpaData_debug <-function (PRECASTObj, species = "Human", custom_housekeep = NULL, 
    covariates_use = NULL, seuList = NULL, subsample_rate = 1, 
    sample_seed = 1) 
{
    verbose <- PRECASTObj@parameterList$verbose
    if (!inherits(PRECASTObj, "PRECASTObj")) 
        stop("IntegrateSpaData: Check the argument: PRECASTObj!  PRECASTObj must be a PRECASTObj object.")
    if (is.null(PRECASTObj@seulist)) 
        stop("IntegrateSpaData: Check the argument: PRECASTObj! The slot seulist in PRECASTObj is NULL!")
    if (length(subsample_rate) > 1 | subsample_rate < 0 | subsample_rate > 
        1) 
        stop("subsample_rate must be a real between 0 and 1")
    if (!tolower(species) %in% c("human", "mouse", "unknown")) 
        stop("IntegrateSpaData: Check the argument: species! it must be one of 'Human', 'Mouse' and 'Unknown'!")
    defAssay_vec <- sapply(PRECASTObj@seulist, DefaultAssay)
    if (any(defAssay_vec != defAssay_vec[1])) 
        warning("IntegrateSpaData: there are different default assays in PRECASTObj@seulist that will be used to integrating!")
    n_r <- length(defAssay_vec)
    if (is.null(seuList) && (!is.null(PRECASTObj@seuList))) {
        message("Use PRECASTObj@seuList as input to remove unwanted variation since it is preserved in PRECASTObj.")
        seuList <- PRECASTObj@seuList
    }
    if (!is.null(seuList)) {
        message("seuList is not NULL. Filter the spots in seuList but not in PRECASTObj!")
        if (!is.null(names(seuList)) && !is.null(names(PRECASTObj@seulist))) {
            if (!all(names(seuList) == names(PRECASTObj@seulist))) {
                stop("The names of seuList must be the same as the names of PRECASTObj@seulist.")
            }
        }
        seuList <- lapply(seq_along(PRECASTObj@seulist), function(j) {
            seu1 <- PRECASTObj@seulist[[j]]
            seu <- seuList[[j]]
            if (length(setdiff(colnames(seu1), colnames(seu)))) {
                stop("The spot's name in PRECASTObj@seulist[[", 
                  j, "]] must be a subset of the spot's name in that of seuList!")
            }
            seu <- seu[, colnames(seu1)]
            return(seu)
        })
        seuList <- lapply(seuList, NormalizeData, verbose = FALSE)
        XList <- lapply(1:n_r, function(r) Matrix::t(GetAssayData(seuList[[r]], 
            assay = defAssay_vec[r], slot = "data")))
    }
    else {
        XList <- lapply(1:n_r, function(r) Matrix::t(GetAssayData(PRECASTObj@seulist[[r]], 
            assay = defAssay_vec[r], slot = "data")))
    }
    browser()
    if (!is.null(covariates_use)) {
        covariateList <- lapply(PRECASTObj@seulist, function(x) x@meta.data[covariates_use])
    }
    else {
        covariateList <- NULL
    }
    if (tolower(species) == "mouse") {
        for (r in 1:n_r) {
            colnames(XList[[r]]) <- firstup(colnames(XList[[r]]))
        }
        if (!is.null(custom_housekeep)) {
            custom_housekeep <- firstup(custom_housekeep)
        }
    }
    if (tolower(species) == "human") {
        for (r in 1:n_r) {
            colnames(XList[[r]]) <- toupper(colnames(XList[[r]]))
        }
        if (!is.null(custom_housekeep)) {
            custom_housekeep <- toupper(custom_housekeep)
        }
    }
    barcodes_all <- lapply(XList, row.names)
    if (any(duplicated(unlist(barcodes_all)))) {
        for (r in 1:n_r) {
            row.names(XList[[r]]) <- paste0(row.names(XList[[r]]), 
                r)
        }
    }
    genelist <- colnames(XList[[1]])
    lower_species <- tolower(species)
    houseKeep <- switch(lower_species, human = {
        intersect(toupper(genelist), PRECAST::Human_HK_genes$Gene)
    }, mouse = {
        intersect(firstup(genelist), PRECAST::Mouse_HK_genes$Gene)
    }, unknown = {
        character()
    })
    houseKeep <- c(houseKeep, custom_housekeep)
    houseKeep <- intersect(houseKeep, colnames(XList[[1]]))
    nvec <- sapply(XList, nrow)
    if (sum(nvec) > 80000) {
        subsample_rate <- 50000/sum(nvec)
        message("The total number of spots exceeds 8e4, thus the subsampling schema will be used to speed up computation.")
    }
    if (subsample_rate < 1 && subsample_rate > 0) 
        message("IntegrateSRTData: the subsampling schema will be used to speed up computation since subsample_rate is smaller than 1.")
    tstart <- Sys.time()
    if (length(houseKeep) < 5) {
        if (verbose) {
            message("Using only PRECAST results to obtain the batch corrected gene expressions since species is unknown or the genelist in PRECASTObj has less than 5 overlapp with the housekeeping genes of given species.")
            message("Start integration...")
        }
        hX <- PRECAST:::get_correct_mean_exp(XList, PRECASTObj@resList$hV, 
            covariateList = covariateList, subsample_rate = subsample_rate, 
            sample_seed = sample_seed)
    }
    else {
        if (verbose) {
            message("Using bouth housekeeping gene and PRECAST results to obtain the batch corrected gene expressions.")
            message("Start integration...")
        }
        hX <- get_correct_exp(XList, PRECASTObj@resList$Rf, houseKeep = houseKeep, 
            q_unwanted = min(10, length(houseKeep)), covariateList = covariateList, 
            subsample_rate = subsample_rate, sample_seed = sample_seed)
    }
    .logDiffTime(sprintf(paste0("%s Data integration finished!"), 
        "*****"), t1 = tstart, verbose = verbose)
    if (verbose) 
        message("Put the data into a new Seurat object...")
    tstart <- Sys.time()
    meta_data <- data.frame(batch = factor(get_sampleID(XList)), 
        cluster = factor(unlist(PRECASTObj@resList$cluster)))
    row.names(meta_data) <- row.names(hX)
    rm(XList)
    count <- sparseMatrix(i = 1, j = 1, x = 0, dims = dim(t(hX)))
    row.names(count) <- colnames(hX)
    colnames(count) <- row.names(hX)
    seuInt <- CreateSeuratObject(counts = count, assay = "PRE_CAST", 
        meta.data = meta_data)
    if (inherits(seuInt[["PRE_CAST"]], "Assay5")) {
        seuInt <- SetAssayData(object = seuInt, slot = "data", 
            assay = "PRE_CAST", new.data = t(hX))
    }
    else {
        seuInt[["PRE_CAST"]]@data <- t(hX)
    }
    rm(hX)
    seuInt <- Add_embed(matlist2mat(PRECASTObj@resList$hZ), seuInt, 
        embed_name = "PRECAST", assay = "PRE_CAST")
    posList <- lapply(PRECASTObj@seulist, function(x) cbind(x$row, 
        x$col))
    seuInt <- Add_embed(matlist2mat(posList), seuInt, embed_name = "position", 
        assay = "PRE_CAST")
    Idents(seuInt) <- factor(meta_data$cluster)
    .logDiffTime(sprintf(paste0("%s New Seurat object is generated!"), 
        "*****"), t1 = tstart, verbose = verbose)
    return(seuInt)
}