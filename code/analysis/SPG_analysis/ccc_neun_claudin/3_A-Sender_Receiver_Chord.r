# NOTE: Code adapted from https://github.com/dssikdar/C2Cv0

library(circlize)

# my_netVisual_aggregate2 <- function(
#   net.diff, edge.transparency = FALSE, signaling.name = NULL, color.use = NULL, thresh = 0.05, vertex.receiver = NULL, sources.use = NULL, targets.use = NULL, idents.use = NULL, top = 1, remove.isolate = FALSE,
#    vertex.weight = 1, vertex.weight.max = NULL, vertex.size.max = NULL,
#                                     measure = c("weight", "count"),
#                                     layout = c("circle", "chord"),
#                                     weight.scale = TRUE, edge.weight.max = NULL, edge.width.max = 8,
#                                     pt.title = 10, title.space = 6, vertex.label.cex = 0.8, title.cex = 1.1,
#                                     alpha.image = 0.15, point.size = 1.5,
#                                     group = NULL, cell.order = NULL, small.gap = 1, big.gap = 10, scale = FALSE, reduce = -1, show.legend = FALSE, legend.pos.x = 20, legend.pos.y = 20,
#                                     ...) {
#   gg <- my_netVisual_chord_cell_internal(net.diff,
#     edge.transparency = edge.transparency, color.use = color.use, sources.use = sources.use, targets.use = targets.use, remove.isolate = remove.isolate,
#     group = group, cell.order = cell.order,
#     lab.cex = vertex.label.cex, small.gap = small.gap, big.gap = big.gap,
#     scale = scale, reduce = reduce, title.cex = title.cex,
#     ti
#   tle.name = "", show.legend = show.legend, legend.pos.x = legend.pos.x, legend.pos.y = legend.pos.y
#   )
#   return(gg)
# }


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
  # browser()
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




my_netVisual_chord_cell_internal <- function(net, edge.transparency = FALSE, color.use = NULL, group = NULL, cell.order = NULL,
                                             sources.use = NULL, targets.use = NULL,
                                             lab.cex = 0.8, small.gap = 1, big.gap = 10, annotationTrackHeight = c(0.06),
                                             remove.isolate = FALSE, link.visible = TRUE, scale = FALSE, directional = 1, link.target.prop = TRUE, reduce = -1,
                                             transparency = 0.4, link.border = NA, title.cex = 1,
                                             title.name = NULL, show.legend = FALSE, legend.pos.x = 20, legend.pos.y = 20, ...) {
  if (inherits(x = net, what = c("matrix", "Matrix"))) {
    cell.levels <- union(rownames(net), colnames(net))
    net <- reshape2::melt(net, value.name = "prob")
    colnames(net)[1:2] <- c("source", "target")
  } else if (is.data.frame(net)) {
    if (all(c("source", "target", "prob") %in% colnames(net)) == FALSE) {
      stop("The input data frame must contain three columns named as source, target, prob")
    }
    cell.levels <- as.character(union(net$source, net$target))
  }
  if (!is.null(cell.order)) {
    cell.levels <- cell.order
  }
  net$source <- as.character(net$source)
  net$target <- as.character(net$target)

  # keep the interactions associated with sources and targets of interest
  if (!is.null(sources.use)) {
    if (is.numeric(sources.use)) {
      sources.use <- cell.levels[sources.use]
    }
    net <- subset(net, source %in% sources.use)
  }
  if (!is.null(targets.use)) {
    if (is.numeric(targets.use)) {
      targets.use <- cell.levels[targets.use]
    }
    net <- subset(net, target %in% targets.use)
  }
  # remove the interactions with zero values
  net <- subset(net, prob != 0)
  if (dim(net)[1] <= 0) {
    message("No interaction between those cells")
  }
  # create a fake data if keeping the cell types (i.e., sectors) without any interactions
  if (!remove.isolate) {
    cells.removed <- setdiff(cell.levels, as.character(union(net$source, net$target)))
    if (length(cells.removed) > 0) {
      net.fake <- data.frame(cells.removed, cells.removed, 1e-10 * sample(length(cells.removed), length(cells.removed)))
      colnames(net.fake) <- colnames(net)
      net <- rbind(net, net.fake)
      link.visible <- net[, 1:2]
      link.visible$plot <- FALSE
      if (nrow(net) > nrow(net.fake)) {
        link.visible$plot[1:(nrow(net) - nrow(net.fake))] <- TRUE
      }
      # directional <- net[, 1:2]
      # directional$plot <- 0
      # directional$plot[1:(nrow(net) - nrow(net.fake))] <- 1
      # link.arr.type = "big.arrow"
      # message("Set scale = TRUE when remove.isolate = FALSE")
      scale <- TRUE
    }
  }

  df <- net
  cells.use <- union(df$source, df$target)

  # define grid order
  order.sector <- cell.levels[cell.levels %in% cells.use]

  # define grid color
  if (is.null(color.use)) {
    color.use <- scPalette(length(cell.levels))
    names(color.use) <- cell.levels
  } else if (is.null(names(color.use))) {
    names(color.use) <- cell.levels
  }

  grid.col <- color.use[order.sector]
  names(grid.col) <- order.sector

  # set grouping information
  if (!is.null(group)) {
    group <- group[names(group) %in% order.sector]
  }
  # browser()
  # define edge color
  edge.color <- ifelse(df$prob > 0, "#2166ac", "#b2182b")

  ###############################################
  # Create a custom color palette from red to blue
  # range_value <- max(abs(net$prob))
  # # Create a custom color palette from red to blue
  # color_ramp <- colorRamp(c("blue", "white", "red"))
  # NOTE: this normalization could be problematic
  # edge.color <- rgb(color_ramp((net$prob + range_value) / (2 * range_value)), maxColorValue = 255)
  ###############################################

  rgb_matrix <- col2rgb(edge.color) / 255
  if (edge.transparency) {
    max_prob <- max(abs(df$prob))
    for (i in 1:dim(df)[1]) {
      edge.color[i] <- rgb(rgb_matrix[1, i], rgb_matrix[2, i], rgb_matrix[3, i], abs(df$prob[i] / max_prob))
    }
  }
  # edge.color <- color.use[as.character(df$source)]
  # return (rgb_matrix)
  if (directional == 0 | directional == 2) {
    link.arr.type <- "triangle"
  } else {
    link.arr.type <- "big.arrow"
  }

  circos.clear()
  chordDiagram(df,
    order = order.sector,
    col = edge.color,
    grid.col = grid.col,
    transparency = transparency,
    link.border = link.border,
    directional = directional,
    direction.type = c("diffHeight", "arrows"),
    link.arr.type = link.arr.type, # link.border = "white",
    annotationTrack = "grid",
    annotationTrackHeight = annotationTrackHeight,
    preAllocateTracks = list(track.height = max(strwidth(order.sector))),
    small.gap = small.gap,
    big.gap = big.gap,
    link.visible = link.visible,
    scale = TRUE,
    group = group,
    link.target.prop = link.target.prop,
    diffHeight = mm_h(4), target.prop.height = mm_h(2),
    reduce = reduce,
    ...
  )
  circos.track(track.index = 1, panel.fun = function(x, y) {
    xlim <- get.cell.meta.data("xlim")
    xplot <- get.cell.meta.data("xplot")
    ylim <- get.cell.meta.data("ylim")
    sector.name <- get.cell.meta.data("sector.index")
    circos.text(mean(xlim), ylim[1], sector.name, facing = "inside", niceFacing = TRUE, adj = c(0.5, 0), cex = lab.cex)
  }, bg.border = NA)


  # https://jokergoo.github.io/circlize_book/book/legends.html
  if (show.legend) {
    lgd <- ComplexHeatmap::Legend(at = names(grid.col), type = "grid", legend_gp = grid::gpar(fill = grid.col), title = "Cell State")
    ComplexHeatmap::draw(lgd, x = unit(1, "npc") - unit(legend.pos.x, "mm"), y = unit(legend.pos.y, "mm"), just = c("right", "bottom"))
  }

  if (!is.null(title.name)) {
    title(title.name, cex = 1)
    # text(0, 1.3, title.name, cex=title.cex)
  }
  circos.clear()
  gg <- recordPlot()
  return(gg)
}


# Deprecated funciton ----
# myVisual_diffInteraction <- function(object, celltype_name, comparison = c(1, 2), measure = c("count", "weight", "count.merged", "weight.merged"), node.color = c("input", "output"), color.use = NULL, color.edge = c("#b2182b", "#2166ac"), title.name = NULL, title.cex = 1.1, sources.use = NULL, targets.use = NULL, remove.isolate = FALSE, top = 1,
#                                      weight.scale = FALSE, vertex.weight = 20, vertex.weight.max = NULL, vertex.size.max = 15, vertex.label.cex = 1, vertex.label.color = "black",
#                                      edge.weight.max = NULL, edge.width.max = 8, alpha.edge = 0.6, label.edge = FALSE, edge.label.color = "black", edge.label.cex = 0.8,
#                                      edge.curved = 0.2, shape = "circle", layout = in_circle(), margin = 0.2,
#                                      arrow.width = 1, arrow.size = 0.2, title = TRUE, sender = NULL, receiver = NULL) {
#   options(warn = -1)
#   measure <- match.arg(measure)
#   node.color <- match.arg(node.color)

#   # get the matrix
#   obj1 <- object@net[[comparison[1]]][[measure]]
#   obj2 <- object@net[[comparison[2]]][[measure]]
#   melt_1 <- reshape2::melt(obj1, value.name = "count")
#   melt_2 <- reshape2::melt(obj2, value.name = "count")
#   sum1 <- sum(melt_1$count)
#   sum2 <- sum(melt_2$count)

#   # balance
#   net.diff <- obj2 * (sum1 / sum2) - obj1
#   # net.diff <- obj2 - obj1

#   # reorder the celltypes
#   net.diff <- net.diff[celltype_name, celltype_name]

#   if (!is.null(receiver)) {
#     net.diff[, setdiff(celltype_name, receiver)] <- 0 # Keep only edges to specified receivers
#   }

#   # title name
#   if (measure %in% c("count", "count.merged")) {
#     if (is.null(title.name)) {
#       title.name <- "Differential number of interactions"
#     }
#   } else if (measure %in% c("weight", "weight.merged")) {
#     if (is.null(title.name)) {
#       title.name <- "Differential interaction strength"
#     }
#   }
#   if (!title) {
#     title.name <- NULL
#   }
#   net <- net.diff

#   # from source code, only show those celltypes we are interested in
#   if ((!is.null(sources.use)) | (!is.null(targets.use))) {
#     df.net <- reshape2::melt(net, value.name = "value")
#     colnames(df.net)[1:2] <- c("source", "target")
#     # keep the interactions associated with sources and targets of interest
#     if (!is.null(sources.use)) {
#       if (is.numeric(sources.use)) {
#         sources.use <- rownames(net.diff)[sources.use]
#       }
#       df.net <- subset(df.net, source %in% sources.use)
#     }
#     if (!is.null(targets.use)) {
#       if (is.numeric(targets.use)) {
#         targets.use <- rownames(net.diff)[targets.use]
#       }
#       df.net <- subset(df.net, target %in% targets.use)
#     }
#     cells.level <- rownames(net.diff)
#     df.net$source <- factor(df.net$source, levels = cells.level)
#     df.net$target <- factor(df.net$target, levels = cells.level)
#     df.net$value[is.na(df.net$value)] <- 0
#     net <- tapply(df.net[["value"]], list(df.net[["source"]], df.net[["target"]]), sum)
#     net[is.na(net)] <- 0
#   }

#   # from source code, remove the celltypes which the sum of input or output is zero
#   if (remove.isolate) {
#     idx1 <- which(Matrix::rowSums(net) == 0)
#     idx2 <- which(Matrix::colSums(net) == 0)
#     idx <- intersect(idx1, idx2)
#     net <- net[-idx, ]
#     net <- net[, -idx]
#   }

#   net[abs(net) < stats::quantile(abs(net), probs = 1 - top)] <- 0

#   # start generate figures
#   g <- graph_from_adjacency_matrix(net, mode = "directed", weighted = T)
#   edge.start <- igraph::ends(g, es = igraph::E(g), names = FALSE)
#   coords <- layout_(g, layout)
#   if (nrow(coords) != 1) {
#     coords_scale <- scale(coords)
#   } else {
#     coords_scale <- coords
#   }

#   # color the node by input or output
#   if (is.null(color.use)) {
#     if (node.color == "input") {
#       s <- colSums(net)
#     } else if (node.color == "output") {
#       s <- rowSums(net)
#     }
#     s <- s / max(abs(s) + 1e-10)
#     color.use <- ifelse(s[igraph::V(g)] > 0, rgb(0.698, 0.094, 0.168, abs(s[igraph::V(g)])), rgb(0.129, 0.4, 0.674, abs(s[igraph::V(g)])))
#   }

#   ## Original coloring
#   # if (is.null(color.use)) {
#   #   color.use = scPalette(length(igraph::V(g)))
#   # }

#   # the rest are all from source code
#   if (is.null(vertex.weight.max)) {
#     vertex.weight.max <- max(vertex.weight)
#   }
#   vertex.weight <- vertex.weight / vertex.weight.max * vertex.size.max + 5

#   loop.angle <- ifelse(coords_scale[igraph::V(g), 1] > 0, -atan(coords_scale[igraph::V(g), 2] / coords_scale[igraph::V(g), 1]), pi - atan(coords_scale[igraph::V(g), 2] / coords_scale[igraph::V(g), 1]))
#   igraph::V(g)$size <- vertex.weight
#   igraph::V(g)$color <- color.use[igraph::V(g)]
#   igraph::V(g)$frame.color <- color.use[igraph::V(g)]
#   igraph::V(g)$label.color <- vertex.label.color
#   igraph::V(g)$label.cex <- vertex.label.cex
#   if (label.edge) {
#     igraph::E(g)$label <- igraph::E(g)$weight
#     igraph::E(g)$label <- round(igraph::E(g)$label, digits = 1)
#   }
#   igraph::E(g)$arrow.width <- arrow.width
#   igraph::E(g)$arrow.size <- arrow.size
#   igraph::E(g)$label.color <- edge.label.color
#   igraph::E(g)$label.cex <- edge.label.cex
#   # igraph::E(g)$color<- grDevices::adjustcolor(igraph::V(g)$color[edge.start[,1]],alpha.edge)
#   igraph::E(g)$color <- ifelse(igraph::E(g)$weight > 0, color.edge[1], color.edge[2])
#   igraph::E(g)$color <- grDevices::adjustcolor(igraph::E(g)$color, alpha.edge)

#   igraph::E(g)$weight <- abs(igraph::E(g)$weight)

#   if (is.null(edge.weight.max)) {
#     edge.weight.max <- max(igraph::E(g)$weight)
#   }
#   if (weight.scale == TRUE) {
#     # E(g)$width<-0.3+edge.width.max/(max(E(g)$weight)-min(E(g)$weight))*(E(g)$weight-min(E(g)$weight))
#     igraph::E(g)$width <- 0.3 + igraph::E(g)$weight / edge.weight.max * edge.width.max
#   } else {
#     igraph::E(g)$width <- 0.3 + edge.width.max * igraph::E(g)$weight
#   }
#   if (sum(edge.start[, 2] == edge.start[, 1]) != 0) {
#     igraph::E(g)$loop.angle[which(edge.start[, 2] == edge.start[, 1])] <- loop.angle[edge.start[which(edge.start[, 2] == edge.start[, 1]), 1]]
#   }
#   radian.rescale <- function(x, start = 0, direction = 1) {
#     c.rotate <- function(x) (x + start) %% (2 * pi) * direction
#     c.rotate(scales::rescale(x, c(0, 2 * pi), range(x)))
#   }
#   label.locs <- radian.rescale(x = 1:length(igraph::V(g)), direction = -1, start = 0)
#   label.dist <- vertex.weight / max(vertex.weight) + 2
#   plot(g,
#     edge.curved = edge.curved, vertex.shape = shape, vertex.frame.color = "black", vertex.frame.width = 0.2, margin = margin, layout = coords_scale, vertex.label.dist = label.dist,
#     vertex.label.degree = label.locs, vertex.label.family = "Helvetica", edge.label.family = "Helvetica"
#   ) # "sans"
#   if (!is.null(title.name)) {
#     text(0, 1.45, title.name, cex = title.cex)
#   }
#   # https://www.andrewheiss.com/blog/2016/12/08/save-base-graphics-as-pseudo-objects-in-r/
#   # grid.echo()
#   # gg <-  grid.grab()
#   gg <- recordPlot()
#   return(gg)
# }
