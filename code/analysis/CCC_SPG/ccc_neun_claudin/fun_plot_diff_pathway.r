my_netVisual_aggregate <- function(
    object.DIS, object.CT, signaling,
    edge.transparency = FALSE, signaling.name = NULL, # signaling.name doesn't seem to have particular use.
     color.use = NULL,
    thresh = 0.05, vertex.receiver = NULL, sources.use = NULL, targets.use = NULL,
    idents.use = NULL, top = 1, remove.isolate = FALSE,
    vertex.weight = 1, vertex.weight.max = NULL, vertex.size.max = NULL,
    measure = c("weight", "count"),
    layout = c("circle", "chord"),
    weight.scale = TRUE, edge.weight.max = NULL, edge.width.max = 8,
    pt.title = 12, title.space = 6, vertex.label.cex = 0.8, title.cex = 1.1,
    alpha.image = 0.15, point.size = 1.5,
    group = NULL, cell.order = NULL, small.gap = 1, big.gap = 10, scale = FALSE, reduce = -1, show.legend = FALSE, legend.pos.x = 20, legend.pos.y = 20,
    ...) {
  layout <- match.arg(layout)
  measure <- match.arg(measure)
  if (is.null(vertex.size.max)) {
    if (length(unique(vertex.weight)) == 1) {
      vertex.size.max <- 5
    } else {
      vertex.size.max <- 15
    }
  }
  browser()
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

  if (is.null(signaling.name)) {
    signaling.name <- signaling
  }

  net.DIS <- object.DIS@net
  net.CT <- object.CT@net

  pairLR.use.name.DIS <- dimnames(net.DIS$prob)[[3]]
  pairLR.use.name.CT <- dimnames(net.CT$prob)[[3]]

  pairLR.name <- intersect(
    intersect(rownames(pairLR), pairLR.use.name.DIS), pairLR.use.name.CT)

  pairLR <- pairLR[pairLR.name, ]
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

  # prob <-(prob-min(prob))/(max(prob)-min(prob))
  prob.sum.DIS <- apply(prob.DIS, c(1, 2), sum)
  prob.sum.CT <- apply(prob.CT, c(1, 2), sum)
  prob.sum <- prob.sum.DIS - prob.sum.CT
  if (!is.null(cell.order)) {
    prob.sum <- prob.sum[cell.order, cell.order]
  }
  if (layout == "circle") {
    # gg <- my_netVisual_circle(prob.sum, sources.use = sources.use, targets.use = targets.use, idents.use = idents.use,
    #                        remove.isolate = remove.isolate, top = top, color.use = color.use, vertex.weight = vertex.weight,
    #                        vertex.weight.max = vertex.weight.max, vertex.size.max = vertex.size.max, weight.scale = weight.scale,
    #                        edge.weight.max = edge.weight.max, edge.width.max=edge.width.max, title.cex=title.cex,
    #                        title.name = paste0(signaling.name, " signaling pathway network"), vertex.label.cex = vertex.label.cex,...)
    stop("Not implemented.")
    gg <- my_netVisual_circle(prob.sum,
      sources.use = sources.use, targets.use = targets.use, idents.use = idents.use,
      remove.isolate = remove.isolate, top = top, color.use = color.use, vertex.weight = vertex.weight,
      vertex.weight.max = vertex.weight.max, vertex.size.max = vertex.size.max, weight.scale = weight.scale,
      edge.weight.max = edge.weight.max, edge.width.max = edge.width.max, title.cex = title.cex,
      title.name = paste0(signaling, " signaling pathway network"), vertex.label.cex = vertex.label.cex, ...
    )
  } else if (layout == "chord") {
    # gg <- my_netVisual_chord_cell_internal(prob.sum, color.use = color.use, sources.use = sources.use, targets.use = targets.use, remove.isolate = remove.isolate,
    #                                     group = group, cell.order = cell.order,
    #                                     lab.cex = vertex.label.cex,small.gap = small.gap, big.gap = big.gap,
    #                                     scale = scale, reduce = reduce, title.cex = title.cex,
    #                                     title.name = paste0(signaling.name, " signaling pathway network"), show.legend = show.legend, legend.pos.x = legend.pos.x, legend.pos.y= legend.pos.y)
    gg <- my_netVisual_chord_cell_internal(prob.sum,
      edge.transparency = edge.transparency, color.use = color.use, sources.use = sources.use, targets.use = targets.use, remove.isolate = remove.isolate,
      group = group, cell.order = cell.order,
      lab.cex = vertex.label.cex, small.gap = small.gap, big.gap = big.gap,
      scale = scale, reduce = reduce, title.cex = title.cex,
      title.name = paste0(signaling, " signaling pathway network"), show.legend = show.legend, legend.pos.x = legend.pos.x, legend.pos.y = legend.pos.y
    )
  }

  return(gg)
}
