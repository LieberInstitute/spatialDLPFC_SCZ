extract_diff_pathway <- function(
    object.CT, object.DIS,
    measure = "weight", # The lot name of cell chat object
    signaling, # Name of a pathway
    thresh = 0.05) {
  measure <- match.arg(measure)


  # Identify the significant pathways
  pairLR.DIS <- searchPair(
    signaling = signaling, pairLR.use = object.DIS@LR$LRsig,
    key = "pathway_name", matching.exact = T, pair.only = T
  )
  pairLR.CT <- searchPair(
    signaling = signaling, pairLR.use = object.CT@LR$LRsig,
    key = "pathway_name", matching.exact = T, pair.only = T
  )


  pairLR <- merge(pairLR.DIS, pairLR.CT)
  rownames(pairLR) <- pairLR$interaction_name

  # if (is.null(signaling.name)) {
  signaling.name <- signaling
  # }

  net.DIS <- object.DIS@net
  net.CT <- object.CT@net

  pairLR.use.name.DIS <- dimnames(net.DIS$prob)[[3]]
  pairLR.use.name.CT <- dimnames(net.CT$prob)[[3]]

  pairLR.name <- intersect(
    intersect(rownames(pairLR), pairLR.use.name.DIS), pairLR.use.name.CT
  )

  pairLR <- pairLR[pairLR.name, ]
  # browser()
  if (measure == "weight") {
    prob.DIS <- net.DIS$prob
    pval.DIS <- net.DIS$pval
    prob.CT <- net.CT$prob
    pval.CT <- net.CT$pval

    prob.DIS[pval.DIS > thresh] <- 0
    prob.CT[pval.CT > thresh] <- 0


    if (length(pairLR.name) > 1) {
      pairLR.name.use.DIS <- pairLR.name[apply(prob.DIS[, , pairLR.name], 3, sum) != 0]
      pairLR.name.use.CT <- pairLR.name[apply(prob.CT[, , pairLR.name], 3, sum) != 0]
    } else {
      pairLR.name.use.DIS <- pairLR.name[sum(prob.DIS[, , pairLR.name]) != 0]
      pairLR.name.use.CT <- pairLR.name[sum(prob.CT[, , pairLR.name]) != 0]
    }
    # pairLR.name.use <- intersect(pairLR.name.use.DIS, pairLR.name.use.CT)
    pairLR.name.use <- union(pairLR.name.use.DIS, pairLR.name.use.CT)

    if (length(pairLR.name.use) == 0) {
      stop(paste0("There is no significant communication of ", signaling.name))
    } else {
      pairLR <- pairLR[pairLR.name.use, ]
    }

    # browser()

    prob.DIS <- prob.DIS[, , pairLR.name.use]
    pval.DIS <- pval.DIS[, , pairLR.name.use]
    prob.CT <- prob.CT[, , pairLR.name.use]
    pval.CT <- pval.CT[, , pairLR.name.use]

    if (length(dim(prob.DIS)) == 2) {
      prob.DIS <- replicate(1, prob.DIS, simplify = "array")
      pval.DIS <- replicate(1, pval.DIS, simplify = "array")
      prob.CT <- replicate(1, prob.CT, simplify = "array")
      pval.CT <- replicate(1, pval.CT, simplify = "array")
    }
  }


  # BOYI: NOTE
  # Calculate the difference.
  # prob <-(prob-min(prob))/(max(prob)-min(prob))
  prob.sum.DIS <- apply(prob.DIS, c(1, 2), sum)
  prob.sum.CT <- apply(prob.CT, c(1, 2), sum)
  prob.sum <- prob.sum.DIS - prob.sum.CT
  # if (!is.null(cell.order)) {
  #   prob.sum <- prob.sum[cell.order, cell.order]
  # }

  return(prob.sum)
}
