subsetCommunication <- function(
    object = NULL, net = NULL, slot.name = "net", sources.use = NULL,
    targets.use = NULL, signaling = NULL, pairLR.use = NULL,
    thresh = 0.05, datasets = NULL, ligand.pvalues = NULL, ligand.logFC = NULL,
    ligand.pct.1 = NULL, ligand.pct.2 = NULL, receptor.pvalues = NULL,
    receptor.logFC = NULL, receptor.pct.1 = NULL, receptor.pct.2 = NULL) {
  if (!is.null(pairLR.use)) {
    if (!is.data.frame(pairLR.use)) {
      stop("pairLR.use should be a data frame with a signle column named either 'interaction_name' or 'pathway_name' ")
    } else if ("pathway_name" %in% colnames(pairLR.use)) {
      message("slot.name is set to be 'netP' when pairLR.use contains signaling pathways")
      slot.name <- "netP"
    }
  }
  if (!is.null(pairLR.use) & !is.null(signaling)) {
    stop("Please do not assign values to 'signaling' when using 'pairLR.use'")
  }
  if (object@options$mode == "single") {
    if (is.null(net)) {
      net <- slot(object, "net")
    }
    LR <- object@LR$LRsig
    cells.level <- levels(object@idents)
    df.net <- subsetCommunication_internal(net, LR, cells.level,
      slot.name = slot.name, sources.use = sources.use,
      targets.use = targets.use, signaling = signaling,
      pairLR.use = pairLR.use, thresh = thresh, datasets = datasets,
      ligand.pvalues = ligand.pvalues, ligand.logFC = ligand.logFC,
      ligand.pct.1 = ligand.pct.1, ligand.pct.2 = ligand.pct.2,
      receptor.pvalues = receptor.pvalues, receptor.logFC = receptor.logFC,
      receptor.pct.1 = receptor.pct.1, receptor.pct.2 = receptor.pct.2
    )
  } else if (object@options$mode == "merged") {
    if (is.null(net)) {
      net0 <- slot(object, "net")
      df.net <- vector("list", length(net0))
      names(df.net) <- names(net0)
      for (i in 1:length(net0)) {
        net <- net0[[i]]
        LR <- object@LR[[i]]$LRsig
        cells.level <- levels(object@idents[[i]])
        # print("here")
        df.net[[i]] <- subsetCommunication_internal(net,
          LR, cells.level,
          slot.name = slot.name, sources.use = sources.use,
          targets.use = targets.use, signaling = signaling,
          pairLR.use = pairLR.use, thresh = thresh, datasets = datasets,
          ligand.pvalues = ligand.pvalues, ligand.logFC = ligand.logFC,
          ligand.pct.1 = ligand.pct.1, ligand.pct.2 = ligand.pct.2,
          receptor.pvalues = receptor.pvalues, receptor.logFC = receptor.logFC,
          receptor.pct.1 = receptor.pct.1, receptor.pct.2 = receptor.pct.2
        )
      }
    } else {
      LR <- data.frame()
      for (i in 1:length(object@LR)) {
        LR <- rbind(LR, object@LR[[i]]$LRsig)
      }
      LR <- unique(LR)
      cells.level <- levels(object@idents$joint)
      df.net <- subsetCommunication_internal(net, LR, cells.level,
        slot.name = slot.name, sources.use = sources.use,
        targets.use = targets.use, signaling = signaling,
        pairLR.use = pairLR.use, thresh = thresh, datasets = datasets,
        ligand.pvalues = ligand.pvalues, ligand.logFC = ligand.logFC,
        ligand.pct.1 = ligand.pct.1, ligand.pct.2 = ligand.pct.2,
        receptor.pvalues = receptor.pvalues, receptor.logFC = receptor.logFC,
        receptor.pct.1 = receptor.pct.1, receptor.pct.2 = receptor.pct.2
      )
    }
  }
  return(df.net)
}
