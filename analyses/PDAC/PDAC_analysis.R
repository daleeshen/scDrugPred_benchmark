library(pROC)
library(caret)
library(PRROC)
library(dplyr)
library(clusterProfiler)
library(org.Hs.eg.db)

#### modify ClusterGVis main function ####
my_visCluster<-function (object = NULL, htColList = list(col_range = c(-2, 
                                                                       0, 2), col_color = c("#08519C", "white", "#A50F15")), border = TRUE, 
                         plotType = c("line", "heatmap", "both"), msCol = c("#0099CC", 
                                                                            "grey90", "#CC3333"), lineSize = 0.1, lineCol = "grey90", 
                         addMline = TRUE, mlineSize = 2, mlineCol = "#CC3333", ncol = 4, 
                         ctAnnoCol = NULL, setMd = "median", textboxPos = c(0.5, 
                                                                            0.8), textboxSize = 8, panelArg = c(2, 0.25, 4, "grey90", 
                                                                                                                NA), ggplotPanelArg = c(2, 0.25, 4, "grey90", NA), annoTermData = NULL, 
                         annoTermMside = "right", termAnnoArg = c("grey95", "grey50"), 
                         addBar = FALSE, barWidth = 8, textbarPos = c(0.8, 0.8), 
                         goCol = NULL, goSize = NULL, byGo = "anno_link", annoKeggData = NULL, 
                         annoKeggMside = "right", keggAnnoArg = c("grey95", "grey50"), 
                         addKeggBar = FALSE, keggCol = NULL, keggSize = NULL, byKegg = "anno_link", 
                         wordWrap = TRUE, addNewLine = TRUE, addBox = FALSE, boxCol = NULL, 
                         boxArg = c(0.1, "grey50"), addPoint = FALSE, pointArg = c(19, 
                                                                                   "orange", "orange", 1), addLine = TRUE, lineSide = "right", 
                         markGenes = NULL, markGenesSide = "right", genesGp = c("italic", 
                                                                                10, NA), termTextLimit = c(10, 18), mulGroup = NULL, 
                         lgdLabel = NULL, showRowNames = FALSE, subgroupAnno = NULL, 
                         annnoblockText = TRUE, annnoblockGp = c("white", 8), addSampleAnno = TRUE, 
                         sampleGroup = NULL, sampleCol = NULL, sampleOrder = NULL, 
                         clusterOrder = NULL, sampleCellOrder = NULL, heatmapAnnotation = NULL, 
                         columnSplit = NULL, clusterColumns = FALSE, pseudotimeCol = NULL, 
                         gglist = NULL, rowAnnotationObj = NULL, ...) 
{
  if (!requireNamespace("ComplexHeatmap", quietly = TRUE)) {
    stop("Package 'ComplexHeatmap' is required. Please install it.")
  }
  ComplexHeatmap::ht_opt(message = FALSE)
  if (is.null(htColList[["col_range"]])) {
    col_range <- c(-2, 0, 2)
  }
  else {
    col_range <- htColList[["col_range"]]
  }
  if (is.null(htColList[["col_color"]])) {
    col_color <- c("#08519C", "white", "#A50F15")
  }
  else {
    col_color <- htColList[["col_color"]]
  }
  if (!requireNamespace("circlize", quietly = TRUE)) {
    stop("Package 'circlize' is required. Please install it.")
  }
  col_fun <- circlize::colorRamp2(col_range, col_color)
  plotType <- match.arg(plotType)
  if (plotType == "line") {
    if (object$type %in% c("scRNAdata", "monocle", "wgcna")) {
      data <- dplyr::arrange(data.frame(object$long.res), 
                             .data[["cluster"]])
    }
    else {
      if (object$type %in% c("mfuzz", "TCseq")) {
        data <- dplyr::arrange(data.frame(object$long.res), 
                               .data[["cluster"]], .data[["membership"]])
      }
      else {
        data <- dplyr::arrange(data.frame(object$long.res), 
                               .data[["cluster"]])
      }
    }
    data$gene <- factor(data$gene, levels = unique(data$gene))
    if (!is.null(sampleOrder)) {
      data$cell_type <- factor(data$cell_type, levels = sampleOrder)
    }
    line <- ggplot2::ggplot(data, ggplot2::aes(x = cell_type, 
                                               y = norm_value))
    if (object$type %in% c("mfuzz", "TCseq")) {
      line <- line + ggplot2::geom_line(ggplot2::aes(color = membership, 
                                                     group = gene), linewidth = lineSize) + ggplot2::scale_color_gradient2(low = msCol[1], 
                                                                                                                           mid = msCol[2], high = msCol[3], midpoint = 0.5)
    }
    else {
      line <- line + ggplot2::geom_line(ggplot2::aes(group = gene), 
                                        color = lineCol, linewidth = lineSize)
    }
    if (addMline == TRUE) {
      if (object$type == "wgcna") {
        linec <- unique(data$modulecol)
        names(linec) <- linec
        line <- line + ggplot2::geom_line(stat = "summary", 
                                          fun = "median", linewidth = mlineSize, ggplot2::aes(group = 1, 
                                                                                              color = modulecol)) + ggplot2::scale_color_manual(values = linec)
      }
      else {
        line <- line + ggplot2::geom_line(stat = "summary", 
                                          fun = "median", colour = mlineCol, linewidth = mlineSize, 
                                          ggplot2::aes(group = 1))
      }
    }
    else {
      line <- line
    }
    line1 <- line + ggplot2::theme_classic(base_size = 14) + 
      ggplot2::ylab("Normalized expression") + ggplot2::xlab("") + 
      ggplot2::theme(axis.ticks.length = ggplot2::unit(0.1, 
                                                       "cm"), axis.text.x = ggplot2::element_text(angle = 45, 
                                                                                                  hjust = 1, color = "black"), strip.background = ggplot2::element_blank()) + 
      ggplot2::facet_wrap(~cluster_name, ncol = ncol, 
                          scales = "free")
    return(line1)
  }
  else {
    data <- dplyr::arrange(data.frame(object$wide.res, check.names = FALSE), 
                           as.numeric(as.character(cluster)))
    if (object$type %in% c("mfuzz", "TCseq")) {
      mat <- dplyr::select(dplyr::arrange(data, as.numeric(as.character(cluster))), 
                           -gene, -cluster, -membership)
    }
    else if (object$type == "wgcna") {
      mat <- dplyr::select(dplyr::arrange(data, as.numeric(as.character(cluster))), 
                           -gene, -cluster, -modulecol)
    }
    else if (object$type == "scRNAdata") {
      mat <- dplyr::select(dplyr::arrange(data, as.numeric(as.character(cluster))), 
                           -gene, -cluster)
    }
    else if (object$type == "monocle") {
      mat <- dplyr::select(dplyr::arrange(data, as.numeric(as.character(cluster))), 
                           -gene, -cluster)
    }
    else {
      mat <- dplyr::select(dplyr::arrange(data, as.numeric(as.character(cluster))), 
                           -gene, -cluster)
    }
    rownames(mat) <- data$gene
    if (object$geneMode == "all" | ncol(mat) > 20) {
      use_raster <- TRUE
    }
    else {
      use_raster <- FALSE
    }
    if (!is.null(sampleOrder)) {
      mat <- mat[, sampleOrder]
    }
    cl.info <- dplyr::arrange(dplyr::mutate(data.frame(table(data$cluster)), 
                                            Var1 = as.numeric(as.character(Var1))), Var1)
    cluster.num <- nrow(cl.info)
    subgroup <- unlist(lapply(seq_len(nrow(cl.info)), function(x) {
      nm <- rep(as.character(cl.info$Var1[x]), cl.info$Freq[x])
      paste("C", nm, sep = "")
    }))
    if (!is.null(clusterOrder)) {
      subgroup <- factor(subgroup, levels = paste("C", 
                                                  clusterOrder, sep = ""))
      cluster_row_slices <- FALSE
    }
    else {
      cluster_row_slices <- TRUE
    }
    if (object$geneMode == "all" & object$type == "scRNAdata") {
      celltype <- vapply(strsplit(colnames(mat), split = "\\|"), 
                         function(x) {
                           x[2]
                         }, character(1))
      cell.num.info <- table(celltype)[unique(celltype)]
      if (is.null(sampleCellOrder)) {
        column_split <- factor(rep(names(cell.num.info), 
                                   cell.num.info), levels = unique(celltype))
      }
      else {
        column_split <- factor(rep(names(cell.num.info), 
                                   cell.num.info), levels = sampleCellOrder)
      }
      if (is.null(sampleCol)) {
        block.col <- seq_len(length(cell.num.info))
      }
      else {
        block.col <- sampleCol
      }
    }
    else {
      if (is.null(sampleGroup)) {
        sample.info <- colnames(mat)
        if (is.null(heatmapAnnotation)) {
          if (object$geneType == "branched") {
            if (ncol(mat) == 200) {
              column_split <- rep(c("branch1", "branch2"), 
                                  each = 100)
            }
            else {
              column_split <- rep(levels(object$pseudotime), 
                                  rep(100, ncol(mat)/100))
            }
          }
          else {
            column_split <- NULL
          }
        }
        else {
          column_split <- columnSplit
        }
      }
      else {
        sample.info <- sampleGroup
        column_split <- factor(sampleGroup, levels = unique(sampleGroup))
      }
      if (object$type != "monocle") {
        sample.info <- factor(sample.info, levels = unique(sample.info))
        if (is.null(sampleCol)) {
          if (!requireNamespace("circlize", quietly = TRUE)) {
            stop("Package 'circlize' is required. Please install it.")
          }
          scol <- circlize::rand_color(n = length(sample.info))
          names(scol) <- sample.info
        }
        else {
          scol <- sampleCol
          names(scol) <- sample.info
        }
      }
      else {
        sample.info <- factor(object$pseudotime, levels = unique(object$pseudotime))
        if (is.null(pseudotimeCol)) {
          if (object$geneType == "branched") {
            if (length(unique(object$pseudotime)) == 
                3) {
              pseudotimeCol <- c("red", "grey80", "blue")
            }
            else {
              if (!requireNamespace("circlize", quietly = TRUE)) {
                stop("Package 'circlize' is required. Please install it.")
              }
              pseudotimeCol <- circlize::rand_color(n = length(unique(object$pseudotime)))
            }
          }
          else {
            pseudotimeCol <- c("blue", "red")
          }
        }
        else {
          pseudotimeCol <- pseudotimeCol
        }
        if (is.null(sampleCol)) {
          if (object$type != "monocle") {
            if (!requireNamespace("circlize", quietly = TRUE)) {
              stop("Package 'circlize' is required. Please install it.")
            }
            scol <- circlize::rand_color(n = length(sample.info))
            names(scol) <- sample.info
          }
          else {
            if (object$geneType == "branched") {
              if (length(unique(object$pseudotime)) == 
                  3) {
                scol <- rep(pseudotimeCol, table(object$pseudotime)[unique(object$pseudotime)])
                names(scol) <- sample.info
              }
              else {
                scol <- rep(pseudotimeCol, table(object$pseudotime)[unique(object$pseudotime)])
                names(scol) <- sample.info
              }
            }
            else {
              scol <- (grDevices::colorRampPalette(pseudotimeCol))(100)
              names(scol) <- sample.info
            }
          }
        }
        else {
          scol <- sampleCol
          names(scol) <- sample.info
        }
      }
    }
    if (addSampleAnno == TRUE) {
      if (object$geneMode == "all" & object$type == "scRNAdata") {
        topanno <- ComplexHeatmap::HeatmapAnnotation(cluster = ComplexHeatmap::anno_block(gp = grid::gpar(fill = block.col), 
                                                                                          labels = NULL), show_annotation_name = FALSE)
      }
      else {
        if (is.null(heatmapAnnotation)) {
          topanno <- ComplexHeatmap::HeatmapAnnotation(sample = sample.info, 
                                                       col = list(sample = scol), gp = grid::gpar(col = ifelse(object$type == 
                                                                                                                 "monocle", NA, "white")), show_legend = ifelse(object$type == 
                                                                                                                                                                  "monocle", FALSE, TRUE), show_annotation_name = FALSE)
        }
        else {
          topanno <- heatmapAnnotation
        }
      }
    }
    else {
      topanno <- NULL
    }
    if (is.null(ctAnnoCol)) {
      if (!requireNamespace("circlize", quietly = TRUE)) {
        stop("Package 'circlize' is required. Please install it.")
      }
      colanno <- circlize::rand_color(n = cluster.num)
    }
    else {
      colanno <- ctAnnoCol
    }
    names(colanno) <- seq_len(cluster.num)
    align_to <- split(seq_len(nrow(mat)), subgroup)
    anno.block <- ComplexHeatmap::anno_block(align_to = align_to, 
                                             panel_fun = function(index, nm) {
                                               npos <- as.numeric(unlist(strsplit(nm, split = "C"))[2])
                                               grid::grid.rect(gp = grid::gpar(fill = colanno[npos], 
                                                                               col = NA))
                                               if (annnoblockText == TRUE) {
                                                 grid::grid.text(label = paste("Num: ", length(index), 
                                                                               sep = ""), rot = 90, gp = grid::gpar(col = annnoblockGp[1], 
                                                                                                                    fontsize = as.numeric(annnoblockGp[2])))
                                               }
                                             }, which = "row")
    if (!is.null(markGenes)) {
      rowGene <- rownames(mat)
      annoGene <- markGenes
      gene.col <- dplyr::filter(dplyr::select(data, gene, 
                                              cluster), gene %in% annoGene)
      gene.col <- purrr::map_df(seq_len(cluster.num), 
                                function(x) {
                                  tmp <- dplyr::mutate(dplyr::filter(gene.col, 
                                                                     as.numeric(cluster) == x), col = colanno[x])
                                })
      gene.col <- gene.col[match(annoGene, gene.col$gene), 
      ]
      if (is.na(genesGp[3])) {
        gcol <- gene.col$col
      }
      else {
        gcol <- genesGp[3]
      }
      index <- match(annoGene, rowGene)
      geneMark <- ComplexHeatmap::anno_mark(at = index, 
                                            labels = annoGene, which = "row", side = markGenesSide, 
                                            labels_gp = grid::gpar(fontface = genesGp[1], 
                                                                   fontsize = as.numeric(genesGp[2]), col = gcol))
    }
    else {
      geneMark <- NULL
    }
    right_annotation <- ComplexHeatmap::rowAnnotation(gene = geneMark, 
                                                      cluster = anno.block)
    if (object$type == "monocle" | object$geneMode == "all" | 
        ncol(mat) > 20) {
      show_column_names <- FALSE
    }
    else {
      show_column_names <- TRUE
    }
    if (object$geneType == "non-branched") {
      rg <- range(as.numeric(as.character(sample.info)))
      if (!requireNamespace("circlize", quietly = TRUE)) {
        stop("Package 'circlize' is required. Please install it.")
      }
      col_fun2 <- circlize::colorRamp2(c(rg[1], rg[2]), 
                                       pseudotimeCol)
      lgd <- ComplexHeatmap::Legend(col_fun = col_fun2, 
                                    title = "pseudotime")
      lgd_list <- list(lgd)
    }
    else if (object$geneType == "branched") {
      if (length(levels(sample.info)) == 3) {
        lgd <- ComplexHeatmap::Legend(labels = levels(sample.info), 
                                      legend_gp = grid::gpar(fill = pseudotimeCol), 
                                      title = "branch")
      }
      else {
        lgd <- ComplexHeatmap::Legend(labels = levels(sample.info), 
                                      legend_gp = grid::gpar(fill = pseudotimeCol), 
                                      title = "branch")
      }
      lgd_list <- list(lgd)
    }
    else {
      lgd_list <- NULL
    }
    if (plotType == "heatmap") {
      if (!is.null(rowAnnotationObj)) {
        left_annotation_ht <- rowAnnotationObj
      }
      else {
        left_annotation_ht <- NULL
      }
      htf <- ComplexHeatmap::Heatmap(as.matrix(mat), name = "Z-score", 
                                     cluster_columns = clusterColumns, show_row_names = showRowNames, 
                                     border = border, column_split = column_split, 
                                     row_split = subgroup, cluster_row_slices = cluster_row_slices, 
                                     column_names_side = "top", show_column_names = show_column_names, 
                                     top_annotation = topanno, left_annotation = left_annotation_ht, 
                                     right_annotation = right_annotation, col = col_fun, 
                                     use_raster = use_raster, ...)
      ComplexHeatmap::draw(htf, merge_legend = TRUE, annotation_legend_list = lgd_list)
    }
    else {
      rg <- range(mat)
      if (!is.null(gglist)) {
        anno_ggplot2 <- ComplexHeatmap::anno_zoom(align_to = align_to, 
                                                  which = "row", panel_fun = function(index, 
                                                                                      nm) {
                                                    g <- gglist[[nm]]
                                                    g <- grid::grid.grabExpr(grid::grid.draw(g))
                                                    grid::pushViewport(grid::viewport())
                                                    grid::grid.rect()
                                                    grid::grid.draw(g)
                                                    grid::popViewport()
                                                  }, size = grid::unit(as.numeric(ggplotPanelArg[1]), 
                                                                       "cm"), gap = grid::unit(as.numeric(ggplotPanelArg[2]), 
                                                                                               "cm"), width = grid::unit(as.numeric(ggplotPanelArg[3]), 
                                                                                                                         "cm"), side = "right", link_gp = grid::gpar(fill = ggplotPanelArg[4], 
                                                                                                                                                                     col = ggplotPanelArg[5]))
      }
      else {
        anno_ggplot2 <- NULL
      }
      panel_fun <- function(index, nm) {
        if (addBox == TRUE & addLine != TRUE) {
          xscale <- c(-0.1, 1.1)
        }
        else {
          xscale <- c(-0.1, 1.1)
          panel_scale <- c(0.1, 0.9)
        }
        grid::pushViewport(grid::viewport(xscale = xscale, 
                                          yscale = c(0, 1)))
        grid::grid.rect()
        if (object$geneMode == "all" & object$type == 
            "scRNAdata") {
          mulGroup <- cell.num.info
          grid::grid.lines(x = c(0, 1), y = rep(0.5, 
                                                2), gp = grid::gpar(col = "black", lty = "dashed"))
          cu <- cumsum(mulGroup)
          seqn <- data.frame(st = c(1, cu[seq(1, (length(cu) - 
                                                    1), 1)] + 1), sp = c(cu[1], cu[seq(2, length(cu), 
                                                                                       1)]))
        }
        else {
          if (is.null(mulGroup)) {
            mulGroup <- ncol(mat)
            seqn <- data.frame(st = 1, sp = ncol(mat))
          }
          else {
            mulGroup <- mulGroup
            grid::grid.lines(x = c(0, 1), y = rep(0.5, 
                                                  2), gp = grid::gpar(col = "black", lty = "dashed"))
            cu <- cumsum(mulGroup)
            seqn <- data.frame(st = c(1, cu[seq(1, length(cu) - 
                                                  1, 1)] + 1), sp = c(cu[1], cu[seq(2, length(cu), 
                                                                                    1)]))
          }
        }
        if (object$geneMode == "all" && object$type == 
            "scRNAdata") {
          cell.ave <- purrr::map_dfr(seq_len(nrow(seqn)), 
                                     function(x) {
                                       tmp <- seqn[x, ]
                                       tmpmat <- mat[index, seq(tmp$st, tmp$sp, 
                                                                1)]
                                       rg <- base::range(mat[index, ])
                                       if (setMd == "mean") {
                                         mdia <- base::mean(base::rowMeans(tmpmat))
                                       }
                                       else if (setMd == "median") {
                                         mdia <- stats::median(base::apply(tmpmat, 
                                                                           1, stats::median))
                                       }
                                       else {
                                         message("supply mean/median !")
                                       }
                                       res <- data.frame(x = x, val = mdia)
                                       return(res)
                                     })
          if (addLine == TRUE) {
            grid::grid.lines(x = scales::rescale(cell.ave$x, 
                                                 to = c(0, 1)), y = scales::rescale(cell.ave$val, 
                                                                                    to = c(0.1, 0.9)), gp = grid::gpar(lwd = 3, 
                                                                                                                       col = mlineCol))
          }
        }
        else {
          base::lapply(seq_len(nrow(seqn)), function(x) {
            tmp <- seqn[x, ]
            tmpmat <- mat[index, seq(tmp$st, tmp$sp, 
                                     1)]
            if (setMd == "mean") {
              mdia <- base::colMeans(tmpmat)
            }
            else if (setMd == "median") {
              mdia <- base::apply(tmpmat, 2, stats::median)
            }
            else {
              message("supply mean/median !")
            }
            pos <- scales::rescale(seq_len(ncol(tmpmat)), 
                                   to = c(0, 1))
            if (is.null(boxCol)) {
              boxCol <- rep("grey90", ncol(tmpmat))
            }
            else {
              boxCol <- boxCol
            }
            if (addBox == TRUE) {
              lapply(seq_len(ncol(tmpmat)), function(x) {
                ComplexHeatmap::grid.boxplot(scales::rescale(tmpmat[, 
                                                                    x], to = c(0, 1), from = c(rg[1] - 
                                                                                                 0.5, rg[2] + 0.5)), pos = pos[x], 
                                             direction = "vertical", box_width = as.numeric(boxArg[1]), 
                                             outline = FALSE, gp = grid::gpar(col = boxArg[2], 
                                                                              fill = boxCol[x]))
              })
            }
            if (addPoint == TRUE) {
              grid::grid.points(x = scales::rescale(seq_len(ncol(tmpmat)), 
                                                    to = c(0.1, 0.9)), y = scales::rescale(mdia, 
                                                                                           to = c(0, 1), from = c(rg[1] - 0.5, 
                                                                                                                  rg[2] + 0.5)), pch = as.numeric(pointArg[1]), 
                                gp = grid::gpar(fill = pointArg[2], 
                                                col = pointArg[3]), size = grid::unit(as.numeric(pointArg[4]), 
                                                                                      "char"))
            }
            if (addLine == TRUE) {
              grid::grid.lines(x = scales::rescale(seq_len(ncol(tmpmat)), 
                                                   to = c(0.1, 0.9)), y = scales::rescale(mdia, 
                                                                                          to = c(0, 1), from = c(rg[1] - 0.5, 
                                                                                                                 rg[2] + 0.5)), gp = grid::gpar(lwd = 3, 
                                                                                                                                                col = mlineCol[x]))
            }
          })
        }
        text <- paste("Gene size:", nrow(mat[index, ]), sep = " ")
        grid::grid.text(text, x = textboxPos[1], y = textboxPos[2], 
                        gp = grid::gpar(fontsize = textboxSize, fontface = "italic"))
        grid::popViewport()
      }
      if (!is.null(subgroupAnno)) {
        align_to <- split(seq_len(nrow(mat)), subgroup)
        align_to <- align_to[subgroupAnno]
      }
      else {
        align_to <- subgroup
      }
      anno <- ComplexHeatmap::anno_link(align_to = align_to, 
                                        which = "row", panel_fun = panel_fun, size = grid::unit(as.numeric(panelArg[1]), 
                                                                                                "cm"), gap = grid::unit(as.numeric(panelArg[2]), 
                                                                                                                        "cm"), width = grid::unit(as.numeric(panelArg[3]), 
                                                                                                                                                  "cm"), side = lineSide, link_gp = grid::gpar(fill = panelArg[4], 
                                                                                                                                                                                               col = panelArg[5]))
      if (!is.null(annoTermData)) {
        termanno <- annoTermData
        if (ncol(termanno) == 2) {
          colnames(termanno) <- c("id", "term")
        }
        else if (ncol(termanno) == 3) {
          colnames(termanno) <- c("id", "term", "pval")
        }
        else if (ncol(termanno) == 4) {
          colnames(termanno) <- c("id", "term", "pval", 
                                  "ratio")
        }
        else {
          message("No more than 4 columns!")
        }
        if (is.null(goCol)) {
          if (!requireNamespace("circlize", quietly = TRUE)) {
            stop("Package 'circlize' is required. Please install it.")
          }
          gocol <- circlize::rand_color(n = nrow(termanno))
        }
        else {
          gocol <- goCol
        }
        if (is.null(goSize)) {
          gosize <- rep(12, nrow(termanno))
        }
        else {
          if (goSize == "pval") {
            termanno.tmp <- purrr::map_df(unique(termanno$id), 
                                          function(x) {
                                            tmp <- dplyr::mutate(dplyr::filter(termanno, 
                                                                               id == x), size = scales::rescale(-log10(pval), 
                                                                                                                to = termTextLimit))
                                          })
            gosize <- termanno.tmp$size
          }
          else {
            gosize <- goSize
          }
        }
        termanno <- dplyr::mutate(dplyr::ungroup(termanno), 
                                  col = gocol, fontsize = gosize)
        term.list <- lapply(seq_len(length(unique(termanno$id))), 
                            function(x) {
                              tmp <- termanno[which(termanno$id == unique(termanno$id)[x]), 
                              ]
                              df <- data.frame(text = tmp$term, col = tmp$col, 
                                               fontsize = tmp$fontsize)
                              return(df)
                            })
        names(term.list) <- unique(termanno$id)
        if (!is.null(subgroupAnno)) {
          align_to2 <- split(seq_along(subgroup), subgroup)
          align_to2 <- align_to2[subgroupAnno]
          term.list <- term.list[subgroupAnno]
        }
        else {
          align_to2 <- subgroup
          term.list <- term.list
        }
        textbox <- ComplexHeatmap::anno_textbox(align_to2, 
                                                term.list, word_wrap = wordWrap, add_new_line = addNewLine, 
                                                side = annoTermMside, background_gp = grid::gpar(fill = termAnnoArg[1], 
                                                                                                 col = termAnnoArg[2]), by = byGo)
        if (ncol(termanno) - 2 > 2) {
          anno_gobar <- function(data = NULL, barWidth = 0.1, 
                                 align_to = NULL, panelArg = panelArg, ...) {
            if (ncol(data) - 2 == 3) {
              data <- dplyr::mutate(data, bary = -log10(pval))
            }
            else {
              data <- dplyr::mutate(data, bary = ratio)
            }
            ComplexHeatmap::anno_zoom(align_to = align_to, 
                                      which = "row", panel_fun = function(index, 
                                                                          nm) {
                                        grid::pushViewport(grid::viewport(xscale = c(0, 
                                                                                     1), yscale = c(0, 1)))
                                        grid::grid.rect()
                                        tmp <- dplyr::filter(data, id == nm)
                                        grid::grid.segments(x0 = rep(0, nrow(tmp)), 
                                                            x1 = scales::rescale(rev(tmp$bary), 
                                                                                 to = c(0.1, 0.9)), y0 = scales::rescale(seq_len(nrow(tmp)), 
                                                                                                                         to = c(0.1, 0.9)), y1 = scales::rescale(seq_len(nrow(tmp)), 
                                                                                                                                                                 to = c(0.1, 0.9)), gp = grid::gpar(lwd = barWidth, 
                                                                                                                                                                                                    col = rev(tmp$col), lineend = "butt"))
                                        grid.textbox <- utils::getFromNamespace("grid.textbox", 
                                                                                "ComplexHeatmap")
                                        text <- nm
                                        grid.textbox(text, x = textbarPos[1], 
                                                     y = textbarPos[2], gp = grid::gpar(fontsize = textboxSize, 
                                                                                        fontface = "italic", col = unique(tmp$col), 
                                                                                        ...))
                                        grid::popViewport()
                                      }, size = grid::unit(as.numeric(panelArg[1]), 
                                                           "cm"), gap = grid::unit(as.numeric(panelArg[2]), 
                                                                                   "cm"), width = grid::unit(as.numeric(panelArg[3]), 
                                                                                                             "cm"), side = "right", link_gp = grid::gpar(fill = termAnnoArg[1], 
                                                                                                                                                         col = termAnnoArg[2]), ...)
          }
          baranno <- anno_gobar(data = termanno, align_to = align_to2, 
                                panelArg = panelArg, barWidth = barWidth)
        }
        if (addBar == TRUE) {
          baranno
        }
        else {
          baranno <- NULL
        }
      }
      else {
        textbox <- NULL
        baranno <- NULL
      }
      if (!is.null(annoKeggData)) {
        termanno <- annoKeggData
        if (ncol(termanno) == 2) {
          colnames(termanno) <- c("id", "term")
        }
        else if (ncol(termanno) == 3) {
          colnames(termanno) <- c("id", "term", "pval")
        }
        else if (ncol(termanno) == 4) {
          colnames(termanno) <- c("id", "term", "pval", 
                                  "ratio")
        }
        else {
          message("No more than 4 columns!")
        }
        if (is.null(keggCol)) {
          if (!requireNamespace("circlize", quietly = TRUE)) {
            stop("Package 'circlize' is required. Please install it.")
          }
          gocol <- circlize::rand_color(n = nrow(termanno))
        }
        else {
          gocol <- keggCol
        }
        if (is.null(keggSize)) {
          gosize <- rep(12, nrow(termanno))
        }
        else {
          if (keggSize == "pval") {
            termanno.tmp <- purrr::map_df(unique(termanno$id), 
                                          function(x) {
                                            tmp <- dplyr::mutate(dplyr::filter(termanno, 
                                                                               id == x), size = scales::rescale(-log10(pval), 
                                                                                                                to = termTextLimit))
                                          })
            gosize <- termanno.tmp$size
          }
          else {
            gosize <- keggSize
          }
        }
        termanno <- dplyr::mutate(dplyr::ungroup(termanno), 
                                  col = gocol, fontsize = gosize)
        term.list <- lapply(seq_len(length(unique(termanno$id))), 
                            function(x) {
                              tmp <- termanno[which(termanno$id == unique(termanno$id)[x]), 
                              ]
                              df <- data.frame(text = tmp$term, col = tmp$col, 
                                               fontsize = tmp$fontsize)
                              return(df)
                            })
        names(term.list) <- unique(termanno$id)
        if (!is.null(subgroupAnno)) {
          align_to2 <- split(seq_along(subgroup), subgroup)
          align_to2 <- align_to2[subgroupAnno]
          term.list <- term.list[subgroupAnno]
        }
        else {
          align_to2 <- subgroup
          term.list <- term.list
        }
        textbox.kegg <- ComplexHeatmap::anno_textbox(align_to2, 
                                                     term.list, word_wrap = wordWrap, add_new_line = addNewLine, 
                                                     side = annoKeggMside, background_gp = grid::gpar(fill = keggAnnoArg[1], 
                                                                                                      col = keggAnnoArg[2]), by = byKegg)
        if (ncol(termanno) - 2 > 2) {
          anno_keggbar <- function(data = NULL, barWidth = 0.1, 
                                   align_to = NULL, panelArg = panelArg, ...) {
            if (ncol(data) - 2 == 3) {
              data <- dplyr::mutate(data, bary = -log10(pval))
            }
            else {
              data <- dplyr::mutate(data, bary = ratio)
            }
            ComplexHeatmap::anno_zoom(align_to = align_to, 
                                      which = "row", panel_fun = function(index, 
                                                                          nm) {
                                        grid::pushViewport(grid::viewport(xscale = c(0, 
                                                                                     1), yscale = c(0, 1)))
                                        grid::grid.rect()
                                        tmp <- dplyr::filter(data, id == nm)
                                        grid::grid.segments(x0 = rep(0, nrow(tmp)), 
                                                            x1 = scales::rescale(rev(tmp$bary), 
                                                                                 to = c(0.1, 0.9)), y0 = scales::rescale(seq_len(nrow(tmp)), 
                                                                                                                         to = c(0.1, 0.9)), y1 = scales::rescale(seq_len(nrow(tmp)), 
                                                                                                                                                                 to = c(0.1, 0.9)), gp = grid::gpar(lwd = barWidth, 
                                                                                                                                                                                                    col = rev(tmp$col), lineend = "butt"))
                                        grid.textbox <- utils::getFromNamespace("grid.textbox", 
                                                                                "ComplexHeatmap")
                                        text <- nm
                                        grid.textbox(text, x = textbarPos[1], 
                                                     y = textbarPos[2], gp = grid::gpar(fontsize = textboxSize, 
                                                                                        fontface = "italic", col = unique(tmp$col), 
                                                                                        ...))
                                        grid::popViewport()
                                      }, size = grid::unit(as.numeric(panelArg[1]), 
                                                           "cm"), gap = grid::unit(as.numeric(panelArg[2]), 
                                                                                   "cm"), width = grid::unit(as.numeric(panelArg[3]), 
                                                                                                             "cm"), side = "right", link_gp = grid::gpar(fill = keggAnnoArg[1], 
                                                                                                                                                         col = keggAnnoArg[2]), ...)
          }
          baranno.kegg <- anno_keggbar(data = termanno, 
                                       align_to = align_to2, panelArg = panelArg, 
                                       barWidth = barWidth)
        }
        if (addKeggBar == TRUE) {
          baranno.kegg
        }
        else {
          baranno.kegg <- NULL
        }
      }
      else {
        textbox.kegg <- NULL
        baranno.kegg <- NULL
      }
      if (lineSide == "right") {
        if (markGenesSide == "right") {
          right_annotation2 <- ComplexHeatmap::rowAnnotation(gene = geneMark, 
                                                             cluster = anno.block, line = anno, anno_ggplot2 = anno_ggplot2, 
                                                             textbox = textbox, bar = baranno, textbox.kegg = textbox.kegg, 
                                                             baranno.kegg = baranno.kegg)
          if (!is.null(rowAnnotationObj)) {
            left_annotation <- rowAnnotationObj
          }
          else {
            left_annotation <- NULL
          }
        }
        else {
          right_annotation2 <- ComplexHeatmap::rowAnnotation(cluster = anno.block, 
                                                             line = anno, anno_ggplot2 = anno_ggplot2, 
                                                             textbox = textbox, bar = baranno, textbox.kegg = textbox.kegg, 
                                                             baranno.kegg = baranno.kegg)
          if (!is.null(rowAnnotationObj)) {
            left_annotation <- do.call(ComplexHeatmap::rowAnnotation, 
                                       modifyList(list(gene = geneMark), rowAnnotationObj))
          }
          else {
            left_annotation <- ComplexHeatmap::rowAnnotation(gene = geneMark)
          }
        }
      }
      else {
        if (markGenesSide == "right") {
          right_annotation2 <- ComplexHeatmap::rowAnnotation(gene = geneMark, 
                                                             cluster = anno.block, anno_ggplot2 = anno_ggplot2, 
                                                             textbox = textbox, bar = baranno, textbox.kegg = textbox.kegg, 
                                                             baranno.kegg = baranno.kegg)
          if (!is.null(rowAnnotationObj)) {
            left_annotation <- do.call(ComplexHeatmap::rowAnnotation, 
                                       modifyList(list(line = anno), rowAnnotationObj))
          }
          else {
            left_annotation <- ComplexHeatmap::rowAnnotation(line = anno)
          }
        }
        else {
          right_annotation2 <- ComplexHeatmap::rowAnnotation(cluster = anno.block, 
                                                             anno_ggplot2 = anno_ggplot2, textbox = textbox, 
                                                             bar = baranno, textbox.kegg = textbox.kegg, 
                                                             baranno.kegg = baranno.kegg)
          if (!is.null(rowAnnotationObj)) {
            left_annotation <- do.call(ComplexHeatmap::rowAnnotation, 
                                       modifyList(list(gene = geneMark, line = anno), 
                                                  rowAnnotationObj))
          }
          else {
            left_annotation <- ComplexHeatmap::rowAnnotation(line = anno, 
                                                             gene = geneMark)
          }
        }
      }
      if (object$type == "monocle" | object$geneMode == 
          "all" | ncol(mat) > 20) {
        show_column_names <- FALSE
      }
      else {
        show_column_names <- TRUE
      }
      htf <- ComplexHeatmap::Heatmap(as.matrix(mat), name = "Z-score", 
                                     cluster_columns = clusterColumns, show_row_names = showRowNames, 
                                     border = border, column_split = column_split, 
                                     top_annotation = topanno, right_annotation = right_annotation2, 
                                     left_annotation = left_annotation, column_names_side = "top", 
                                     show_column_names = show_column_names, row_split = subgroup, 
                                     cluster_row_slices = cluster_row_slices, col = col_fun, 
                                     use_raster = use_raster, ...)
      if (is.null(mulGroup)) {
        ComplexHeatmap::draw(htf, merge_legend = TRUE, 
                             annotation_legend_list = lgd_list)
      }
      else {
        if (is.null(lgdLabel)) {
          lgdLabel <- paste("group", seq_len(length(mulGroup)), 
                            sep = "")
        }
        else {
          lgdLabel <- lgdLabel
        }
        lgd_list2 <- ComplexHeatmap::Legend(labels = lgdLabel, 
                                            type = "lines", legend_gp = grid::gpar(col = mlineCol, 
                                                                                   lty = 1))
        if (!is.null(lgd_list)) {
          lgd_list_com <- ComplexHeatmap::packLegend(lgd_list, 
                                                     lgd_list2)
        }
        else {
          lgd_list_com <- lgd_list2
        }
        ComplexHeatmap::draw(htf, annotation_legend_list = lgd_list_com, 
                             merge_legend = TRUE)
      }
    }
  }
}

############## END modification ############## 

pal_all_eval_df <- data.frame()
tram_all_eval_df <- data.frame()
pal_meta <- as.data.frame(pal_dmso@meta.data)
tram_meta <- as.data.frame(tram_dmso@meta.data)

evaluate_drug_tool_flexible <- function(meta_df, res_matrix, drug_name, tool_name, method = "median", fixed_threshold = NULL, positive_class = "resistant") {
  common_ids <- intersect(rownames(meta_df), res_matrix$ID)
  if(length(common_ids) == 0) stop("未找到匹配的 ID")
  y_true <- as.factor(meta_df[common_ids,"condition"])
  matched_res <- res_matrix[match(common_ids,res_matrix$ID),]
  y_pred_score <- as.numeric(matched_res[[drug_name]])
  roc_obj <- pROC::roc(y_true,y_pred_score,quiet = TRUE)
  auc_val <- as.numeric(pROC::auc(roc_obj))
  pos_scores <- y_pred_score[y_true == positive_class]
  neg_scores <- y_pred_score[y_true != positive_class]
  pr_obj <- PRROC::pr.curve(scores.class0 = pos_scores,scores.class1 = neg_scores,curve = FALSE)
  ap_val <- pr_obj$auc.integral
  threshold <- 0.5
  if(method == "youden") threshold <- pROC::coords(roc_obj,"best",ret = "threshold",transpose = TRUE)[1]
  if(method == "median") threshold <- median(y_pred_score,na.rm = TRUE)
  if(method == "mean") threshold <- mean(y_pred_score,na.rm = TRUE)
  if(method == "fixed" && !is.null(fixed_threshold)) threshold <- fixed_threshold
  other_class <- setdiff(levels(y_true),positive_class)
  y_pred_class <- ifelse(y_pred_score >= threshold,positive_class,other_class[1])
  y_pred_class <- factor(y_pred_class,levels = levels(y_true))
  cm <- caret::confusionMatrix(data = y_pred_class,reference = y_true,positive = positive_class,mode = "everything")
  data.frame(drug = drug_name,AUC = auc_val,AP = ap_val,ACC = as.numeric(cm$overall["Accuracy"]),F1 = as.numeric(cm$byClass["F1"]),Recall = as.numeric(cm$byClass["Recall"]),Precision = as.numeric(cm$byClass["Precision"]),Tool = tool_name,stringsAsFactors = FALSE)
}

pal_common_ids_scdr <- intersect(rownames(pal_meta),colnames(pal_scdr_res))
pal_scdr_pred_df <- data.frame(ID = pal_common_ids_scdr,condition = pal_meta[pal_common_ids_scdr,"condition"],paclitaxel = as.numeric(pal_scdr_res["paclitaxel",pal_common_ids_scdr]),stringsAsFactors = FALSE)
rownames(pal_scdr_pred_df) <- pal_scdr_pred_df$ID
pal_scDr_eval_df <- evaluate_drug_tool_flexible(pal_meta,pal_scdr_pred_df,"paclitaxel","scDr",method = "median",positive_class = "resistant")
pal_all_eval_df <- rbind(pal_all_eval_df,pal_scDr_eval_df)

tram_common_ids_scdr <- intersect(rownames(tram_meta),colnames(tram_scdr_res))
tram_scdr_pred_df <- data.frame(ID = tram_common_ids_scdr,condition = tram_meta[tram_common_ids_scdr,"condition"],trametinib = as.numeric(tram_scdr_res["trametinib",tram_common_ids_scdr]),stringsAsFactors = FALSE)
rownames(tram_scdr_pred_df) <- tram_scdr_pred_df$ID
tram_scDr_eval_df <- evaluate_drug_tool_flexible(tram_meta,tram_scdr_pred_df,"trametinib","scDr",method = "median",positive_class = "resistant")
tram_all_eval_df <- rbind(tram_all_eval_df,tram_scDr_eval_df)

pal_beyondcell_all_res <- read.table("./pred_res/paclitaxel/beyondcell/paclitaxel.beyondcell_res.tsv",header = TRUE,check.names = FALSE)
rownames(pal_beyondcell_all_res) <- pal_beyondcell_all_res$cell_id
pal_common_ids_beyondcell <- intersect(rownames(pal_meta),pal_beyondcell_all_res$cell_id)
pal_beyondcell_res.1 <- data.frame(ID = pal_common_ids_beyondcell,condition = pal_meta[pal_common_ids_scdr,"condition"],paclitaxel = as.numeric(pal_beyondcell_all_res[pal_common_ids_beyondcell,5]))
pal_beyondcell_res.2 <- data.frame(ID = pal_common_ids_beyondcell,condition = pal_meta[pal_common_ids_scdr,"condition"],paclitaxel = as.numeric(pal_beyondcell_all_res[pal_common_ids_beyondcell,6]))
pal_beyondcell_res.3 <- data.frame(ID = pal_common_ids_beyondcell,condition = pal_meta[pal_common_ids_scdr,"condition"],paclitaxel = as.numeric(pal_beyondcell_all_res[pal_common_ids_beyondcell,7]))
pal_beyondcell_res.4 <- data.frame(ID = pal_common_ids_beyondcell,condition = pal_meta[pal_common_ids_scdr,"condition"],paclitaxel = as.numeric(pal_beyondcell_all_res[pal_common_ids_beyondcell,8]))
pal_beyondcell_eval.1 <- evaluate_drug_tool_flexible(pal_meta,pal_beyondcell_res.1,"paclitaxel","beyondcell",method = "median",positive_class = "sensitive")
pal_beyondcell_eval.2 <- evaluate_drug_tool_flexible(pal_meta,pal_beyondcell_res.2,"paclitaxel","beyondcell",method = "median",positive_class = "sensitive")
pal_beyondcell_eval.3 <- evaluate_drug_tool_flexible(pal_meta,pal_beyondcell_res.3,"paclitaxel","beyondcell",method = "median",positive_class = "sensitive")
pal_beyondcell_eval.4 <- evaluate_drug_tool_flexible(pal_meta,pal_beyondcell_res.4,"paclitaxel","beyondcell",method = "median",positive_class = "sensitive")
pal_beyondcell_eval_df <- pal_beyondcell_eval.1
pal_all_eval_df <- rbind(pal_all_eval_df,pal_beyondcell_eval_df)

tram_beyondcell_all_res <- read.table("./pred_res/trametinib/beyondcell/trametinib.beyondcell_res.tsv",header = TRUE,check.names = FALSE)
rownames(tram_beyondcell_all_res) <- tram_beyondcell_all_res$cell_id
tram_common_ids_beyondcell <- intersect(rownames(tram_meta),tram_beyondcell_all_res$cell_id)
tram_beyondcell_res <- data.frame(ID = tram_common_ids_beyondcell,condition = tram_meta[tram_common_ids_beyondcell,"condition"],trametinib = as.numeric(tram_beyondcell_all_res[tram_common_ids_beyondcell,5]))
tram_beyondcell_eval_df <- evaluate_drug_tool_flexible(tram_meta,tram_beyondcell_res,"trametinib","beyondcell",method = "median",positive_class = "sensitive")
tram_all_eval_df <- rbind(tram_all_eval_df,tram_beyondcell_eval_df)

pal_sciduc_res_ccle <- read.csv("./pred_res/paclitaxel/scIDUC/scIDUC_res_CCLE.csv",header = TRUE,check.names = FALSE,row.names = 1)
pal_sciduc_res_gdsc <- read.csv("./pred_res/paclitaxel/scIDUC/scIDUC_res_GDSC.csv",header = TRUE,check.names = FALSE,row.names = 1)
names(pal_sciduc_res_ccle) <- c("ID","paclitaxel")
names(pal_sciduc_res_gdsc) <- c("ID","paclitaxel")
pal_sciduc_eval_ccle <- evaluate_drug_tool_flexible(pal_meta,pal_sciduc_res_ccle,"paclitaxel","scIDUC",method = "median",positive_class = "resistant")
pal_sciduc_eval_gdsc <- evaluate_drug_tool_flexible(pal_meta,pal_sciduc_res_gdsc,"paclitaxel","scIDUC",method = "median",positive_class = "resistant")
pal_sciduc_eval_df <- pal_sciduc_eval_ccle
pal_all_eval_df <- rbind(pal_all_eval_df,pal_sciduc_eval_df)

tram_sciduc_res_ccle <- read.csv("./pred_res/trametinib/scIDUC/scIDUC_res_CCLE.csv",header = TRUE,check.names = FALSE,row.names = 1)
tram_sciduc_res_gdsc <- read.csv("./pred_res/trametinib/scIDUC/scIDUC_res_GDSC.csv",header = TRUE,check.names = FALSE,row.names = 1)
names(tram_sciduc_res_ccle) <- c("ID","trametinib")
names(tram_sciduc_res_gdsc) <- c("ID","trametinib")
tram_sciduc_eval_ccle <- evaluate_drug_tool_flexible(tram_meta,tram_sciduc_res_ccle,"trametinib","scIDUC",method = "median",positive_class = "resistant")
tram_sciduc_eval_gdsc <- evaluate_drug_tool_flexible(tram_meta,tram_sciduc_res_gdsc,"trametinib","scIDUC",method = "median",positive_class = "resistant")
tram_sciduc_eval_df <- tram_sciduc_eval_ccle
tram_all_eval_df <- rbind(tram_all_eval_df,tram_sciduc_eval_df)

rownames(pal_dreep_res) <- pal_dreep_res$cells
pal_common_ids_dreep <- intersect(rownames(pal_meta),pal_dreep_res$cells)
pal_dreep_res.1 <- data.frame(ID = pal_common_ids_dreep,condition = pal_meta[pal_common_ids_dreep,"condition"],paclitaxel = as.numeric(pal_dreep_res[pal_common_ids_dreep,3]))
pal_dreep_res.2 <- data.frame(ID = pal_common_ids_dreep,condition = pal_meta[pal_common_ids_dreep,"condition"],paclitaxel = as.numeric(pal_dreep_res[pal_common_ids_dreep,4]))
pal_dreep_res.3 <- data.frame(ID = pal_common_ids_dreep,condition = pal_meta[pal_common_ids_dreep,"condition"],paclitaxel = as.numeric(pal_dreep_res[pal_common_ids_dreep,5]))
pal_dreep_eval.1 <- evaluate_drug_tool_flexible(pal_meta,pal_dreep_res.1,"paclitaxel","DREEP",method = "median",positive_class = "resistant")
pal_dreep_eval.2 <- evaluate_drug_tool_flexible(pal_meta,pal_dreep_res.2,"paclitaxel","DREEP",method = "median",positive_class = "resistant")
pal_dreep_eval.3 <- evaluate_drug_tool_flexible(pal_meta,pal_dreep_res.3,"paclitaxel","DREEP",method = "median",positive_class = "resistant")
pal_dreep_eval_df <- pal_dreep_eval.3
pal_all_eval_df <- rbind(pal_all_eval_df,pal_dreep_eval_df)

rownames(tram_dreep_res) <- tram_dreep_res$cells
tram_common_ids_dreep <- intersect(rownames(tram_meta),tram_dreep_res$cells)
tram_dreep_res.1 <- data.frame(ID = tram_common_ids_dreep,condition = tram_meta[tram_common_ids_dreep,"condition"],trametinib = as.numeric(tram_dreep_res[tram_common_ids_dreep,3]))
tram_dreep_res.2 <- data.frame(ID = tram_common_ids_dreep,condition = tram_meta[tram_common_ids_dreep,"condition"],trametinib = as.numeric(tram_dreep_res[tram_common_ids_dreep,4]))
tram_dreep_res.3 <- data.frame(ID = tram_common_ids_dreep,condition = tram_meta[tram_common_ids_dreep,"condition"],trametinib = as.numeric(tram_dreep_res[tram_common_ids_dreep,5]))
tram_dreep_eval.1 <- evaluate_drug_tool_flexible(tram_meta,tram_dreep_res.1,"trametinib","DREEP",method = "median",positive_class = "resistant")
tram_dreep_eval.2 <- evaluate_drug_tool_flexible(tram_meta,tram_dreep_res.2,"trametinib","DREEP",method = "median",positive_class = "resistant")
tram_dreep_eval.3 <- evaluate_drug_tool_flexible(tram_meta,tram_dreep_res.3,"trametinib","DREEP",method = "median",positive_class = "resistant")
tram_dreep_eval_df <- tram_dreep_eval.3
tram_all_eval_df <- rbind(tram_all_eval_df,tram_dreep_eval_df)

pal_cadrres_all_res.1 <- read.csv("./pred_res/paclitaxel/cadrres-sc/pal_cadrres-wo-sample-bias_test_pred.csv",header = TRUE,check.names = FALSE)
pal_cadrres_all_res.2 <- read.csv("./pred_res/paclitaxel/cadrres-sc/pal_cadrres-wo-sample-bias-weight_test_pred.csv",header = TRUE,check.names = FALSE)
names(pal_cadrres_all_res.1)[1] <- "cell_id"
names(pal_cadrres_all_res.2)[1] <- "cell_id"
rownames(pal_cadrres_all_res.1) <- pal_cadrres_all_res.1$cell_id
rownames(pal_cadrres_all_res.2) <- pal_cadrres_all_res.2$cell_id
pal_common_ids_cadrres <- intersect(rownames(pal_meta),pal_cadrres_all_res.1$cell_id)
pal_cadrres_res.1 <- data.frame(ID = pal_common_ids_cadrres,condition = pal_meta[pal_common_ids_cadrres,"condition"],paclitaxel = as.numeric(pal_cadrres_all_res.1[pal_common_ids_cadrres,"11"]))
pal_cadrres_res.2 <- data.frame(ID = pal_common_ids_cadrres,condition = pal_meta[pal_common_ids_cadrres,"condition"],paclitaxel = as.numeric(pal_cadrres_all_res.2[pal_common_ids_cadrres,"11"]))
pal_cadrres_eval.1 <- evaluate_drug_tool_flexible(pal_meta,pal_cadrres_res.1,"paclitaxel","CaDRReS-sc",method = "median",positive_class = "sensitive")
pal_cadrres_eval.2 <- evaluate_drug_tool_flexible(pal_meta,pal_cadrres_res.2,"paclitaxel","CaDRReS-sc",method = "median",positive_class = "sensitive")
pal_cadrres_eval_df <- pal_cadrres_eval.1
pal_all_eval_df <- rbind(pal_all_eval_df,pal_cadrres_eval_df)

pal_precily_all_res <- read.csv("./pred_res/paclitaxel/precily/paclitaxel_precily_results.csv",header = TRUE,check.names = FALSE)
rownames(pal_precily_all_res) <- pal_precily_all_res$DRUGS
pal_common_ids_precily <- intersect(rownames(pal_meta),colnames(pal_precily_all_res))
pal_precily_res <- data.frame(ID = pal_common_ids_precily,condition = pal_meta[pal_common_ids_precily,"condition"],paclitaxel = as.numeric(pal_precily_all_res["Paclitaxel",pal_common_ids_precily]))
pal_precily_eval_df <- evaluate_drug_tool_flexible(pal_meta,pal_precily_res,"paclitaxel","precily",method = "median",positive_class = "resistant")
pal_all_eval_df <- rbind(pal_all_eval_df,pal_precily_eval_df)

tram_precily_all_res <- read.csv("./pred_res/trametinib/precily/trametinib_precily_results.csv",header = TRUE,check.names = FALSE)
rownames(tram_precily_all_res) <- tram_precily_all_res$DRUGS
tram_common_ids_precily <- intersect(rownames(tram_meta),colnames(tram_precily_all_res))
tram_precily_res <- data.frame(ID = tram_common_ids_precily,condition = tram_meta[tram_common_ids_precily,"condition"],trametinib = as.numeric(tram_precily_all_res["Trametinib",tram_common_ids_precily]))
tram_precily_eval_df <- evaluate_drug_tool_flexible(tram_meta,tram_precily_res,"trametinib","precily",method = "median",positive_class = "resistant")
tram_all_eval_df <- rbind(tram_all_eval_df,tram_precily_eval_df)

pal_drugformer_files <- list.files("./pred_res/paclitaxel/DrugFormer/paclitaxel_pred_res/",pattern = "epoch_4_fold_.*_res\\.csv",full.names = TRUE)
pal_drugformer_all_res <- do.call(rbind,lapply(pal_drugformer_files,read.csv))
rownames(pal_drugformer_all_res) <- pal_drugformer_all_res$sample_id
pal_common_ids_drugformer <- intersect(rownames(pal_meta),pal_drugformer_all_res$sample_id)
pal_drugformer_res <- data.frame(ID = pal_common_ids_drugformer,condition = pal_meta[pal_common_ids_drugformer,"condition"],paclitaxel = as.numeric(pal_drugformer_all_res$pred_score))
pal_drugformer_eval_df <- evaluate_drug_tool_flexible(pal_meta,pal_drugformer_res,"paclitaxel","DrugFormer",method = "median",positive_class = "sensitive")
pal_all_eval_df <- rbind(pal_all_eval_df,pal_drugformer_eval_df)

tram_drugformer_files <- list.files("./pred_res/trametinib/DrugFormer/trametinib_pred_res/",pattern = "epoch_4_fold_.*_res\\.csv",full.names = TRUE)
tram_drugformer_all_res <- do.call(rbind,lapply(tram_drugformer_files,read.csv))
rownames(tram_drugformer_all_res) <- tram_drugformer_all_res$sample_id
tram_common_ids_drugformer <- intersect(rownames(tram_meta),tram_drugformer_all_res$sample_id)
tram_drugformer_res <- data.frame(ID = tram_common_ids_drugformer,condition = tram_meta[tram_common_ids_drugformer,"condition"],trametinib = as.numeric(tram_drugformer_all_res[tram_common_ids_drugformer,"pred_score"]),stringsAsFactors = FALSE)
tram_drugformer_eval_df <- evaluate_drug_tool_flexible(tram_meta,tram_drugformer_res,"trametinib","DrugFormer",method = "median",positive_class = "sensitive")
tram_all_eval_df <- rbind(tram_all_eval_df,tram_drugformer_eval_df)

pal_scad_all_res.1 <- read.table("./pred_res/paclitaxel/SCAD/best_param_results_norm/markus_pal_Paclitaxel/bimodal/best_models/FX/test/mbS32_mbT32/hdim512_zdim256_lr0.001_epoch30_lambda13_mbS32_mbT32_seed42_predicted_y.tsv",header = FALSE)
pal_scad_all_res.2 <- read.table("./pred_res/paclitaxel/SCAD/best_param_results_ppi/markus_pal_Paclitaxel/bimodal/best_models/FX/test/mbS32_mbT32/hdim512_zdim128_lr0.001_epoch10_lambda15_mbS32_mbT32_seed42_predicted_y.tsv",header = FALSE)
pal_scad_all_res.3 <- read.table("./pred_res/paclitaxel/SCAD/best_param_results_tp4k/markus_pal_Paclitaxel/bimodal/best_models/FX/test/mbS32_mbT32/hdim1024_zdim128_lr0.001_epoch20_lambda11.0_mbS32_mbT32_seed42_predicted_y.tsv",header = FALSE)
names(pal_scad_all_res.1) <- names(pal_scad_all_res.2) <- names(pal_scad_all_res.3) <- c("splits","cell_id","condition","pred_score")
rownames(pal_scad_all_res.1) <- pal_scad_all_res.1$cell_id
rownames(pal_scad_all_res.2) <- pal_scad_all_res.2$cell_id
rownames(pal_scad_all_res.3) <- pal_scad_all_res.3$cell_id
pal_common_ids_scad <- intersect(rownames(pal_meta),rownames(pal_scad_all_res.1))
pal_scad_res.1 <- data.frame(ID = pal_common_ids_scad,condition = pal_meta[pal_common_ids_scad,"condition"],paclitaxel = as.numeric(pal_scad_all_res.1[pal_common_ids_scad,"pred_score"]))
pal_scad_res.2 <- data.frame(ID = pal_common_ids_scad,condition = pal_meta[pal_common_ids_scad,"condition"],paclitaxel = as.numeric(pal_scad_all_res.2[pal_common_ids_scad,"pred_score"]))
pal_scad_res.3 <- data.frame(ID = pal_common_ids_scad,condition = pal_meta[pal_common_ids_scad,"condition"],paclitaxel = as.numeric(pal_scad_all_res.3[pal_common_ids_scad,"pred_score"]))
pal_scad_eval.1 <- evaluate_drug_tool_flexible(pal_meta,pal_scad_res.1,"paclitaxel","SCAD",method = "median",positive_class = "sensitive")
pal_scad_eval.2 <- evaluate_drug_tool_flexible(pal_meta,pal_scad_res.2,"paclitaxel","SCAD",method = "median",positive_class = "sensitive")
pal_scad_eval.3 <- evaluate_drug_tool_flexible(pal_meta,pal_scad_res.3,"paclitaxel","SCAD",method = "median",positive_class = "sensitive")
pal_scad_eval_df <- pal_scad_eval.3
pal_all_eval_df <- rbind(pal_all_eval_df,pal_scad_eval_df)

tram_scad_all_res.1 <- read.table("./pred_res/trametinib/SCAD/best_param_results_norm/markus_tram_Trametinib/bimodal/best_models/FX/test/mbS32_mbT32/hdim512_zdim128_lr0.001_epoch100_lambda12_mbS32_mbT32_seed42_predicted_y.tsv",header = FALSE)
tram_scad_all_res.2 <- read.table("./pred_res/trametinib/SCAD/best_param_results_ppi/markus_tram_Trametinib/bimodal/best_models/FX/test/mbS32_mbT32/hdim512_zdim128_lr0.001_epoch10_lambda12_mbS32_mbT32_seed42_predicted_y.tsv",header = FALSE)
tram_scad_all_res.3 <- read.table("./pred_res/trametinib/SCAD/best_param_results_tp4k/markus_tram_Trametinib/bimodal/best_models/FX/test/mbS8_mbT8/hdim512_zdim256_lr0.001_epoch10_lambda13_mbS8_mbT8_seed42_predicted_y.tsv",header = FALSE)
names(tram_scad_all_res.1) <- names(tram_scad_all_res.2) <- names(tram_scad_all_res.3) <- c("splits","cell_id","condition","pred_score")
rownames(tram_scad_all_res.1) <- tram_scad_all_res.1$cell_id
rownames(tram_scad_all_res.2) <- tram_scad_all_res.2$cell_id
rownames(tram_scad_all_res.3) <- tram_scad_all_res.3$cell_id
tram_common_ids_scad <- intersect(rownames(tram_meta),rownames(tram_scad_all_res.1))
tram_scad_res.1 <- data.frame(ID = tram_common_ids_scad,condition = tram_meta[tram_common_ids_scad,"condition"],trametinib = as.numeric(tram_scad_all_res.1[tram_common_ids_scad,"pred_score"]))
tram_scad_res.2 <- data.frame(ID = tram_common_ids_scad,condition = tram_meta[tram_common_ids_scad,"condition"],trametinib = as.numeric(tram_scad_all_res.2[tram_common_ids_scad,"pred_score"]))
tram_scad_res.3 <- data.frame(ID = tram_common_ids_scad,condition = tram_meta[tram_common_ids_scad,"condition"],trametinib = as.numeric(tram_scad_all_res.3[tram_common_ids_scad,"pred_score"]))
tram_scad_eval.1 <- evaluate_drug_tool_flexible(tram_meta,tram_scad_res.1,"trametinib","SCAD",method = "median",positive_class = "sensitive")
tram_scad_eval.2 <- evaluate_drug_tool_flexible(tram_meta,tram_scad_res.2,"trametinib","SCAD",method = "median",positive_class = "sensitive")
tram_scad_eval.3 <- evaluate_drug_tool_flexible(tram_meta,tram_scad_res.3,"trametinib","SCAD",method = "median",positive_class = "sensitive")
tram_scad_eval_df <- tram_scad_eval.1
tram_all_eval_df <- rbind(tram_all_eval_df,tram_scad_eval_df)

pal_scdeal_files <- list.files("./pred_res/paclitaxel/scDEAL/",pattern = "\\.csv$",full.names = TRUE)
pal_scdeal_all_evals <- data.frame()
for(pal_scdeal_file in pal_scdeal_files) {
  pal_scdeal_res_file <- read.csv(pal_scdeal_file)
  rownames(pal_scdeal_res_file) <- pal_scdeal_res_file[,1]
  pal_common_ids_scdeal <- intersect(rownames(pal_meta),pal_scdeal_res_file[,1])
  pal_scdeal_res.x <- data.frame(ID = pal_common_ids_scdeal,condition = pal_meta[pal_common_ids_scdeal,"condition"],paclitaxel = as.numeric(pal_scdeal_res_file[pal_common_ids_scdeal,"sens_preds_pret"]))
  pal_scdeal_eval.x <- evaluate_drug_tool_flexible(pal_meta,pal_scdeal_res.x,"paclitaxel","scDEAL",method = "median",positive_class = "sensitive")
  pal_scdeal_eval.x$file_source <- basename(pal_scdeal_file)
  pal_scdeal_all_evals <- rbind(pal_scdeal_all_evals,pal_scdeal_eval.x)
}
pal_scdeal_eval_df <- pal_scdeal_all_evals[which.max(pal_scdeal_all_evals$AUC),]
pal_all_eval_df <- rbind(pal_all_eval_df,pal_scdeal_eval_df[,1:8])

tram_scdeal_files <- list.files("./pred_res/trametinib/scDEAL/",pattern = "\\.csv$",full.names = TRUE)
tram_scdeal_all_evals <- data.frame()
for(tram_scdeal_file in tram_scdeal_files) {
  tram_scdeal_res_file <- read.csv(tram_scdeal_file)
  rownames(tram_scdeal_res_file) <- tram_scdeal_res_file[,1]
  tram_common_ids_scdeal <- intersect(rownames(tram_meta),tram_scdeal_res_file[,1])
  tram_scdeal_res.x <- data.frame(ID = tram_common_ids_scdeal,condition = tram_meta[tram_common_ids_scdeal,"condition"],trametinib = as.numeric(tram_scdeal_res_file[tram_common_ids_scdeal,"sens_preds_pret"]),stringsAsFactors = FALSE)
  tram_scdeal_eval.x <- evaluate_drug_tool_flexible(tram_meta,tram_scdeal_res.x,"trametinib","scDEAL",method = "median",positive_class = "sensitive")
  tram_scdeal_eval.x$file_source <- basename(tram_scdeal_file)
  tram_scdeal_all_evals <- rbind(tram_scdeal_all_evals,tram_scdeal_eval.x)
}
tram_scdeal_eval_df <- tram_scdeal_all_evals[which.max(tram_scdeal_all_evals$AUC),]
tram_all_eval_df <- rbind(tram_all_eval_df,tram_scdeal_eval_df[,1:8])

markus_all_eval_df <- rbind(pal_all_eval_df,tram_all_eval_df)
markus_all_eval_df$Tool[markus_all_eval_df$Tool == "beyondcell"] <- "Beyondcell"
markus_all_eval_df$Tool[markus_all_eval_df$Tool == "precily"] <- "Precily"

pal_scdeal_best_pred_df <- read.csv(file.path("./pred_res/paclitaxel/scDEAL",pal_scdeal_eval_df$file_source),header = TRUE,check.names = FALSE,row.names = 1)
pal_dmso$scDEAL_pred <- pal_scdeal_best_pred_df[match(rownames(pal_dmso@meta.data),rownames(pal_scdeal_best_pred_df)),"sens_preds_pret"]
tram_scdeal_best_pred_df <- read.csv(file.path("./pred_res/trametinib/scDEAL",tram_scdeal_eval_df$file_source),header = TRUE,check.names = FALSE,row.names = 1)
tram_dmso$scDEAL_pred <- tram_scdeal_best_pred_df[match(rownames(tram_dmso@meta.data),rownames(tram_scdeal_best_pred_df)),"sens_preds_pret"]

run_cluster_enrichment <- function(seurat_obj, ident_col, ident_1, ident_2, organism = "hsa", species_db = org.Hs.eg.db) {
  Idents(seurat_obj) <- ident_col
  degs <- FindMarkers(seurat_obj,ident.1 = ident_1,ident.2 = ident_2,logfc.threshold = 0.25)
  degs$gene <- rownames(degs)
  degs$direction <- ifelse(degs$avg_log2FC > 0,"Up-regulated","Down-regulated")
  degs_sig <- degs[degs$p_val < 0.05,]
  if(nrow(degs_sig) == 0) stop("No significant DEGs found (p_val < 0.05).")
  gene_df <- bitr(degs_sig$gene,fromType = "SYMBOL",toType = "ENTREZID",OrgDb = species_db)
  gene_with_dir <- merge(gene_df,degs_sig[,c("gene","direction")],by.x = "SYMBOL",by.y = "gene")
  bp_res <- compareCluster(ENTREZID ~ direction,data = gene_with_dir,fun = "enrichGO",OrgDb = species_db,ont = "BP",pvalueCutoff = 0.05)
  kegg_res <- compareCluster(ENTREZID ~ direction,data = gene_with_dir,fun = "enrichKEGG",organism = organism,pvalueCutoff = 0.05)
  list(degs = degs_sig,go_res = bp_res,kegg_res = kegg_res)
}

pal_dmso@meta.data$scDEAL_pred_cond <- ifelse(pal_dmso@meta.data$scDEAL_pred > median(pal_dmso@meta.data$scDEAL_pred,na.rm = TRUE),"sensitive","resistant")
tram_dmso@meta.data$scDEAL_pred_cond <- ifelse(tram_dmso@meta.data$scDEAL_pred > median(tram_dmso@meta.data$scDEAL_pred,na.rm = TRUE),"sensitive","resistant")
Idents(pal_dmso) <- "scDEAL_pred_cond"
Idents(tram_dmso) <- "scDEAL_pred_cond"

cc.genes <- Seurat::cc.genes
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
pal_dmso <- CellCycleScoring(pal_dmso,s.features = s.genes,g2m.features = g2m.genes,set.ident = FALSE)
tram_dmso <- CellCycleScoring(tram_dmso,s.features = s.genes,g2m.features = g2m.genes,set.ident = FALSE)
pal_dmso <- ScaleData(pal_dmso,vars.to.regress = c("S.Score","G2M.Score"))
tram_dmso <- ScaleData(tram_dmso,vars.to.regress = c("S.Score","G2M.Score"))

pal_enrich_results <- run_cluster_enrichment(pal_dmso,"scDEAL_pred_cond","resistant","sensitive")
tram_enrich_results <- run_cluster_enrichment(tram_dmso,"scDEAL_pred_cond","resistant","sensitive")

library(clusterProfiler)
library(org.Hs.eg.db)
library(pheatmap)
library(ggpubr)
library(msigdbr)
library(enrichplot)
library(GSVA)
library(UCell)
library(monocle)
library(ClusterGVis)

run_cluster_enrichment <- function(seurat_obj, ident_col, ident_1, ident_2, organism = "hsa", species_db = org.Hs.eg.db) {
  Idents(seurat_obj) <- ident_col
  degs <- FindMarkers(seurat_obj, ident.1 = ident_1, ident.2 = ident_2, logfc.threshold = 0.25)
  degs$gene <- rownames(degs)
  degs$direction <- ifelse(degs$avg_log2FC > 0, "Up-regulated", "Down-regulated")
  degs_sig <- degs[degs$p_val < 0.05,]
  if(nrow(degs_sig) == 0) stop("No significant DEGs found (p_val < 0.05).")
  gene_df <- bitr(degs_sig$gene, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = species_db)
  gene_with_dir <- merge(gene_df, degs_sig[,c("gene","direction")], by.x = "SYMBOL", by.y = "gene")
  bp_res <- compareCluster(ENTREZID ~ direction, data = gene_with_dir, fun = "enrichGO", OrgDb = species_db, ont = "BP", pvalueCutoff = 0.05)
  kegg_res <- compareCluster(ENTREZID ~ direction, data = gene_with_dir, fun = "enrichKEGG", organism = organism, pvalueCutoff = 0.05)
  list(degs = degs_sig, go_res = bp_res, kegg_res = kegg_res)
}

pal_dmso@meta.data$scDEAL_pred_cond <- ifelse(pal_dmso@meta.data$scDEAL_pred > median(pal_dmso@meta.data$scDEAL_pred, na.rm = TRUE), "sensitive", "resistant")
tram_dmso@meta.data$scDEAL_pred_cond <- ifelse(tram_dmso@meta.data$scDEAL_pred > median(tram_dmso@meta.data$scDEAL_pred, na.rm = TRUE), "sensitive", "resistant")
Idents(pal_dmso) <- "scDEAL_pred_cond"
Idents(tram_dmso) <- "scDEAL_pred_cond"

cc.genes <- Seurat::cc.genes
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
pal_dmso <- CellCycleScoring(pal_dmso, s.features = s.genes, g2m.features = g2m.genes, set.ident = FALSE)
tram_dmso <- CellCycleScoring(tram_dmso, s.features = s.genes, g2m.features = g2m.genes, set.ident = FALSE)

pal_df_phase <- pal_dmso@meta.data %>% count(scDEAL_pred_cond, Phase) %>% group_by(scDEAL_pred_cond) %>% mutate(prop = n/sum(n))
tram_df_phase <- tram_dmso@meta.data %>% count(scDEAL_pred_cond, Phase) %>% group_by(scDEAL_pred_cond) %>% mutate(prop = n/sum(n))

pal_dmso <- ScaleData(pal_dmso, vars.to.regress = c("S.Score","G2M.Score"))
tram_dmso <- ScaleData(tram_dmso, vars.to.regress = c("S.Score","G2M.Score"))
pal_enrich_results <- run_cluster_enrichment(pal_dmso, "scDEAL_pred_cond", "resistant", "sensitive")
tram_enrich_results <- run_cluster_enrichment(tram_dmso, "scDEAL_pred_cond", "resistant", "sensitive")

png("./figs/pal_go_bp_dotplot.png", units = "in", res = 300, width = 8, height = 8)
print(dotplot(pal_enrich_results$go_res) + theme(axis.text.y = element_text(size = 16), axis.text.x = element_text(size = 16)))
dev.off()

png("./figs/tram_go_bp_dotplot.png", units = "in", res = 300, width = 8, height = 8)
print(dotplot(tram_enrich_results$go_res) + theme(axis.text.y = element_text(size = 16), axis.text.x = element_text(size = 16)))
dev.off()

m_t2g <- msigdbr(species = "Homo sapiens", category = "H") %>% dplyr::select(gs_name, gene_symbol)
c2_set <- msigdbr(species = "Homo sapiens", category = "C2") %>% dplyr::select(gs_name, gene_symbol)

Idents(pal_dmso) <- "scDEAL_pred_cond"
pal_deg_all <- FindMarkers(pal_dmso, ident.1 = "resistant", ident.2 = "sensitive", logfc.threshold = 0, min.pct = 0)
pal_deg_gene_list <- pal_deg_all$avg_log2FC
names(pal_deg_gene_list) <- rownames(pal_deg_all)
pal_deg_gene_list <- sort(pal_deg_gene_list, decreasing = TRUE)
pal_gsea_res_h <- GSEA(pal_deg_gene_list, TERM2GENE = m_t2g, pvalueCutoff = 0.99)
pal_gsea_res_c2 <- GSEA(pal_deg_gene_list, TERM2GENE = c2_set, pvalueCutoff = 0.99)

Idents(tram_dmso) <- "scDEAL_pred_cond"
tram_deg_all <- FindMarkers(tram_dmso, ident.1 = "resistant", ident.2 = "sensitive", logfc.threshold = 0, min.pct = 0)
tram_deg_gene_list <- tram_deg_all$avg_log2FC
names(tram_deg_gene_list) <- rownames(tram_deg_all)
tram_deg_gene_list <- sort(tram_deg_gene_list, decreasing = TRUE)
tram_gsea_res_h <- GSEA(tram_deg_gene_list, TERM2GENE = m_t2g, pvalueCutoff = 0.99)
tram_gsea_res_c2 <- GSEA(tram_deg_gene_list, TERM2GENE = c2_set, pvalueCutoff = 0.99)

plot_gsea_with_stats <- function(gsea_results, gene_set_id, base_size = 20, rel_heights = c(1.5,0.5,1)) {
  res_df <- as.data.frame(gsea_results)
  if(!(gene_set_id %in% res_df$ID)) stop("Gene Set ID not found")
  nes <- round(res_df[res_df$ID == gene_set_id,"NES"],2)
  p_raw <- res_df[res_df$ID == gene_set_id,"pvalue"]
  p_exp <- floor(log10(p_raw))
  p_base <- round(p_raw/10^p_exp,2)
  plot_title <- bquote(atop(bold(.(gene_set_id)), (NES == .(nes) ~ "," ~ p == .(p_base) %*% 10^{.(p_exp)})))
  p <- enrichplot::gseaplot2(gsea_results, geneSetID = gene_set_id, title = NULL, base_size = base_size, rel_heights = rel_heights)
  p[[1]] <- p[[1]] + ggtitle(plot_title) + theme(plot.title = element_text(hjust = 0.5))
  p
}

png("./figs/pal_GSEA_g2m_checkpoint.png", units = "in", res = 300, width = 9, height = 8)
p <- plot_gsea_with_stats(pal_gsea_res_h, "HALLMARK_G2M_CHECKPOINT")
aplot::plot_list(p, ncol = 1, heights = c(10.5,0,0))
dev.off()

png("./figs/pal_GSEA_mitotic_spindle.png", units = "in", res = 300, width = 9, height = 8)
p <- plot_gsea_with_stats(pal_gsea_res_h, "HALLMARK_MITOTIC_SPINDLE")
aplot::plot_list(p, ncol = 1, heights = c(10.5,0,0))
dev.off()

png("./figs/tram_GSEA_mitotic_spindle.png", units = "in", res = 300, width = 9, height = 8)
p <- plot_gsea_with_stats(tram_gsea_res_h, "HALLMARK_MITOTIC_SPINDLE")
aplot::plot_list(p, ncol = 1, heights = c(10.5,0,0))
dev.off()

png("./figs/tram_GSEA_g2m_checkpoint.png", units = "in", res = 300, width = 9, height = 8)
p <- plot_gsea_with_stats(tram_gsea_res_h, "HALLMARK_G2M_CHECKPOINT")
aplot::plot_list(p, ncol = 1, heights = c(10.5,0,0))
dev.off()

png("./figs/tram_GSEA_KEGG_RAS_PI3K.png", units = "in", res = 300, width = 12, height = 8)
p <- plot_gsea_with_stats(tram_gsea_res_c2, "KEGG_MEDICUS_REFERENCE_GF_RTK_RAS_PI3K_SIGNALING_PATHWAY", base_size = 17)
aplot::plot_list(p, ncol = 1, heights = c(10.5,0,0))
dev.off()

png("./figs/tram_GSEA_KEGG_MAPK.png", units = "in", res = 300, width = 9, height = 8)
p <- plot_gsea_with_stats(tram_gsea_res_c2, "KEGG_MAPK_SIGNALING_PATHWAY")
aplot::plot_list(p, ncol = 1, heights = c(10.5,0,0))
dev.off()

png("./figs/tram_GSEA_MYC_TARGETS.png", units = "in", res = 300, width = 9, height = 8)
p <- plot_gsea_with_stats(tram_gsea_res_h, "HALLMARK_MYC_TARGETS_V1")
aplot::plot_list(p, ncol = 1, heights = c(10.5,0,0))
dev.off()

png("./figs/tram_GSEA_hallmark_KRAS.png", units = "in", res = 300, width = 9, height = 8)
p <- plot_gsea_with_stats(tram_gsea_res_h, "HALLMARK_KRAS_SIGNALING_UP")
aplot::plot_list(p, ncol = 1, heights = c(10.5,0,0))
dev.off()

plot_pretty_vln <- function(obj, features, group_by, colors = c("#66c2a5","#fc8d62"), stat_method = "wilcox.test") {
  plot_data <- FetchData(obj, vars = c(features,group_by))
  plot_data_long <- plot_data %>% pivot_longer(cols = all_of(features), names_to = "Gene", values_to = "Expression")
  ggplot(plot_data_long, aes(x = .data[[group_by]], y = Expression, fill = .data[[group_by]])) +
    geom_violin(trim = TRUE, color = "black", linewidth = 0.5) +
    geom_boxplot(width = 0.1, linewidth = 0.5, color = "black", outlier.shape = NA, fill = "white") +
    facet_wrap(~Gene, scales = "free_y", ncol = 3) + theme_bw() +
    theme(panel.grid = element_blank(), strip.background = element_blank(), text = element_text(size = 24), axis.title.x = element_blank(), legend.position = "none", legend.title = element_blank()) +
    scale_fill_manual(values = colors) + labs(y = "Expression") +
    stat_compare_means(method = stat_method, comparisons = list(c("sensitive","resistant")), size = 15, vjust = 0.5, label = "p.signif", label.x = 1.5, symnum.args = list(cutpoints = c(0,0.0001,0.001,0.01,0.05,1), symbols = c("****","***","**","*","ns"))) +
    scale_y_continuous(expand = expansion(mult = c(0.05,0.1)))
}

for(gene in c("CCNB1","TOP2A","BUB1","BAX","E2F1","PCNA")) {
  png(paste0("./figs/pal_",gene,".png"), units = "in", res = 300, width = 4, height = 6)
  print(plot_pretty_vln(pal_dmso, features = gene, group_by = "scDEAL_pred_cond", colors = c(my_colors[3],my_colors[1])))
  dev.off()
}

Idents(tram_dmso) <- factor(Idents(tram_dmso), levels = c("sensitive","resistant"))
png("./figs/tram_tar_MAP2K1.png", units = "in", res = 300, height = 6, width = 5)
print(VlnPlot(tram_dmso, features = "MAP2K1") + scale_fill_manual(values = c("skyblue","tomato")) + annotate("segment",x = 1,xend = 2,y = 2.3,yend = 2.3) + annotate("text",x = 1.5,y = 2.4,label = "p == 4.4 %*% 10^-6",parse = TRUE,size = 5) + ylim(0,2.5) + ggtitle("MAP2K1") + theme(text = element_text(size = 18),axis.text = element_text(size = 18),axis.title.x = element_blank()) + NoLegend())
dev.off()

png("./figs/tram_tar_MAP2K2.png", units = "in", res = 300, height = 6, width = 5)
print(VlnPlot(tram_dmso, features = "MAP2K2") + scale_fill_manual(values = c("skyblue","tomato")) + annotate("segment",x = 1,xend = 2,y = 2.8,yend = 2.8) + annotate("text",x = 1.5,y = 2.9,label = "NS",size = 5) + ylim(0,3) + ggtitle("MAP2K2") + theme(text = element_text(size = 18),axis.text = element_text(size = 18),axis.title.x = element_blank()) + NoLegend())
dev.off()

Idents(pal_dmso) <- factor(Idents(pal_dmso), levels = c("sensitive","resistant"))
for(gene in c("TUBB1","BCL2","NR1I2")) {
  ymax <- if(gene == "NR1I2") 1.25 else 1.5
  line_y <- if(gene == "NR1I2") 1.1 else 1.3
  text_y <- if(gene == "NR1I2") 1.2 else 1.4
  png(paste0("./figs/pal_tar_",gene,".png"), units = "in", res = 300, height = 6, width = 5)
  print(VlnPlot(pal_dmso, features = gene) + scale_fill_manual(values = c("skyblue","tomato")) + annotate("segment",x = 1,xend = 2,y = line_y,yend = line_y) + annotate("text",x = 1.5,y = text_y,label = "NS",size = 5) + ylim(0,ymax) + ggtitle(gene) + theme(text = element_text(size = 18),axis.text = element_text(size = 18),axis.title.x = element_blank()) + NoLegend())
  dev.off()
}

m_df <- msigdbr(species = "Homo sapiens", category = "H")
pal_gene_sets <- split(m_df$gene_symbol,m_df$gs_name)
pal_cnts <- as.matrix(GetAssayData(pal_dmso, slot = "data"))
pal_ssgsea_res <- gsva(ssgseaParam(expr = pal_cnts, geneSets = pal_gene_sets, minSize = 10, maxSize = 500))
pal_pathways <- list(Paclitaxel_Resistance = c("ABCB1","ABCC1","TUBB1","MAPT","BCL2","NR1I2"), Microtubule_Stability = c("MAP4","STMN1","MAP1A"))
pal_u_scores <- ScoreSignatures_UCell(pal_cnts, features = pal_pathways)
pal_dmso <- AddMetaData(pal_dmso, metadata = as.data.frame(pal_u_scores))
pal_ssgsea_df <- as.data.frame(t(pal_ssgsea_res))
pal_pathways_to_keep <- intersect(c("HALLMARK_XENOBIOTIC_METABOLISM","HALLMARK_APOPTOSIS","HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION","HALLMARK_G2M_CHECKPOINT"),colnames(pal_ssgsea_df))
pal_dmso <- AddMetaData(pal_dmso,metadata = pal_ssgsea_df[,pal_pathways_to_keep,drop = FALSE])
pal_plot_df <- pal_dmso@meta.data

png("./figs/pal_ssGSEA_apoptosis.png", units = "in", res = 300, height = 6, width = 5)
print(ggplot(pal_plot_df,aes(x = scDEAL_pred_cond,y = HALLMARK_APOPTOSIS,fill = scDEAL_pred_cond)) + geom_violin(trim = FALSE,alpha = 0.7) + geom_boxplot(width = 0.1,fill = "white",outlier.shape = NA) + scale_fill_manual(values = c("tomato","skyblue")) + theme_classic() + labs(title = "APOPTOSIS (ssGSEA)",y = "Enrichment Score",x = NULL) + theme(plot.title = element_text(hjust = 0.5,face = "bold"),axis.text = element_text(size = 18),text = element_text(size = 18),legend.position = "none") + annotate("segment",x = 1,xend = 2,y = 0.4,yend = 0.4) + annotate("text",x = 1.5,y = 0.42,label = "p = 0.041",size = 5) + ylim(NA,0.45))
dev.off()

png("./figs/pal_ssGSEA_g2m.png", units = "in", res = 300, height = 6, width = 5)
print(ggplot(pal_plot_df,aes(x = scDEAL_pred_cond,y = HALLMARK_G2M_CHECKPOINT,fill = scDEAL_pred_cond)) + geom_violin(trim = FALSE,alpha = 0.7) + geom_boxplot(width = 0.1,fill = "white",outlier.shape = NA) + scale_fill_manual(values = c("tomato","skyblue")) + theme_classic() + labs(title = "G2M CHECKPOINT (ssGSEA)",y = "Enrichment Score",x = NULL) + theme(plot.title = element_text(hjust = 0.5,face = "bold"),axis.text = element_text(size = 18),text = element_text(size = 18),legend.position = "none") + annotate("segment",x = 1,xend = 2,y = 0.55,yend = 0.55) + annotate("text",x = 1.5,y = 0.58,label = "p == 2.2 %*% 10^-16",parse = TRUE,size = 5) + ylim(NA,0.6))
dev.off()

png("./figs/pal_ssGSEA_Microtubule_Stability.png", units = "in", res = 300, height = 6, width = 5)
print(ggplot(pal_plot_df,aes(x = scDEAL_pred_cond,y = Microtubule_Stability_UCell,fill = scDEAL_pred_cond)) + geom_violin(trim = FALSE,alpha = 0.7) + geom_boxplot(width = 0.1,fill = "white",outlier.shape = NA) + scale_fill_manual(values = c("tomato","skyblue")) + theme_classic() + labs(title = "Microtubule Stability",y = "Enrichment Score",x = NULL) + theme(plot.title = element_text(hjust = 0.5,face = "bold"),text = element_text(size = 18),axis.text = element_text(size = 18),legend.position = "none") + annotate("segment",x = 1,xend = 2,y = 0.65,yend = 0.65) + annotate("text",x = 1.5,y = 0.68,label = "p == 2.2 %*% 10^-8",parse = TRUE,size = 5) + ylim(NA,0.7))
dev.off()

png("./figs/pal_ssGSEA_Paclitaxel_Resistance_UCell.png", units = "in", res = 300, height = 6, width = 5)
print(ggplot(subset(pal_plot_df,Paclitaxel_Resistance_UCell <= 0.05),aes(x = scDEAL_pred_cond,y = Paclitaxel_Resistance_UCell,fill = scDEAL_pred_cond)) + geom_violin(trim = FALSE,alpha = 0.7) + geom_boxplot(width = 0.1,fill = "white",outlier.shape = NA) + scale_fill_manual(values = c("tomato","skyblue")) + theme_classic() + labs(title = "Paclitaxel Resistance",y = "Enrichment Score",x = NULL) + theme(plot.title = element_text(hjust = 0.5,face = "bold"),text = element_text(size = 18),axis.text = element_text(size = 18),legend.position = "none") + annotate("segment",x = 1,xend = 2,y = 0.055,yend = 0.055) + annotate("text",x = 1.5,y = 0.058,label = "p = 0.017",size = 5) + ylim(NA,0.06))
dev.off()

tram_gene_sets <- split(m_df$gene_symbol,m_df$gs_name)
tram_cnts <- as.matrix(GetAssayData(tram_dmso, slot = "data"))
tram_ssgsea_res <- gsva(ssgseaParam(expr = tram_cnts,geneSets = tram_gene_sets,minSize = 10,maxSize = 500))
tram_pathways <- list(MAPK_Feedback_Targets = c("DUSP6","DUSP4","SPRY2","SPRY4","ETV4","ETV5"), MEK_Resistance_Bypass = c("MET","AXL","ERBB3","FGFR1","RAF1","BRAF"), MAPK_Kinases = c("MAP2K1","MAP2K2","MAPK1","MAPK3"))
tram_u_scores <- ScoreSignatures_UCell(tram_cnts,features = tram_pathways)
tram_dmso <- AddMetaData(tram_dmso,metadata = as.data.frame(tram_u_scores))
tram_ssgsea_df <- as.data.frame(t(tram_ssgsea_res))
tram_pathways_to_keep <- intersect(c("HALLMARK_KRAS_SIGNALING_UP","HALLMARK_KRAS_SIGNALING_DN","HALLMARK_MYC_TARGETS_V1","HALLMARK_PI3K_AKT_MTOR_SIGNALING","HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION"),colnames(tram_ssgsea_df))
tram_dmso <- AddMetaData(tram_dmso,metadata = tram_ssgsea_df[,tram_pathways_to_keep,drop = FALSE])
tram_plot_df <- tram_dmso@meta.data

png("./figs/tram_ssGSEA_KRAS_UP.png", units = "in", res = 300, height = 6, width = 5)
print(ggplot(tram_plot_df,aes(x = scDEAL_pred_cond,y = HALLMARK_KRAS_SIGNALING_UP,fill = scDEAL_pred_cond)) + geom_violin(trim = FALSE,alpha = 0.7) + geom_boxplot(width = 0.1,fill = "white",outlier.shape = NA) + scale_fill_manual(values = c("tomato","skyblue")) + theme_classic() + labs(title = "KRAS SIGNALING (ssGSEA)",y = "Enrichment Score",x = NULL) + theme(plot.title = element_text(hjust = 0.5,face = "bold"),axis.text = element_text(size = 18),text = element_text(size = 18),legend.position = "none") + annotate("segment",x = 1,xend = 2,y = 0.23,yend = 0.23) + annotate("text",x = 1.5,y = 0.25,label = "p == 2.2 %*% 10^-16",parse = TRUE,size = 5) + ylim(NA,0.26))
dev.off()

png("./figs/tram_ssGSEA_EMT.png", units = "in", res = 300, height = 6, width = 5)
print(ggplot(tram_plot_df,aes(x = scDEAL_pred_cond,y = HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION,fill = scDEAL_pred_cond)) + geom_violin(trim = FALSE,alpha = 0.7) + geom_boxplot(width = 0.1,fill = "white",outlier.shape = NA) + scale_fill_manual(values = c("tomato","skyblue")) + theme_classic() + labs(title = "EMT (ssGSEA)",y = "Enrichment Score",x = NULL) + theme(plot.title = element_text(hjust = 0.5,face = "bold"),axis.text = element_text(size = 18),text = element_text(size = 18),legend.position = "none") + annotate("segment",x = 1,xend = 2,y = 0.27,yend = 0.27) + annotate("text",x = 1.5,y = 0.29,label = "p == 2.2 %*% 10^-16",parse = TRUE,size = 5) + ylim(NA,0.3))
dev.off()

png("./figs/tram_ssGSEA_MAPK_Feedback_Targets_UCell.png", units = "in", res = 300, height = 6, width = 5)
print(ggplot(tram_plot_df,aes(x = scDEAL_pred_cond,y = MAPK_Feedback_Targets_UCell,fill = scDEAL_pred_cond)) + geom_violin(trim = FALSE,alpha = 0.7) + geom_boxplot(width = 0.1,fill = "white",outlier.shape = NA) + scale_fill_manual(values = c("tomato","skyblue")) + theme_classic() + labs(title = "MAPK Feedback Targets",y = "Enrichment Score",x = NULL) + theme(plot.title = element_text(hjust = 0.5,face = "bold"),axis.text = element_text(size = 18),text = element_text(size = 18),legend.position = "none") + annotate("segment",x = 1,xend = 2,y = 0.37,yend = 0.37) + annotate("text",x = 1.5,y = 0.39,label = "p == 8.6 %*% 10^-7",parse = TRUE,size = 5) + ylim(NA,0.4))
dev.off()

pal_count_df <- as(as.matrix(pal_dmso@assays$RNA@counts),"sparseMatrix")
pal_pd <- new("AnnotatedDataFrame",data = pal_dmso@meta.data)
pal_fData <- data.frame(gene_short_name = rownames(pal_count_df),row.names = rownames(pal_count_df))
pal_fd <- new("AnnotatedDataFrame",data = pal_fData)
pal_cds <- newCellDataSet(pal_count_df,phenoData = pal_pd,featureData = pal_fd,expressionFamily = negbinomial.size())
pal_cds <- estimateSizeFactors(pal_cds)
pal_cds <- estimateDispersions(pal_cds)
pal_ordering_genes <- pal_enrich_results$degs %>% filter(p_val < 0.05) %>% slice_max(order_by = abs(avg_log2FC),n = 300) %>% rownames()
pal_cds <- setOrderingFilter(pal_cds,pal_ordering_genes)
pal_cds <- reduceDimension(pal_cds,max_components = 2,method = "DDRTree",ncenter = 40)
pal_cds <- orderCells(pal_cds,root_state = 4)

png("./figs/pal_monocle_state.png", units = "in", res = 300, width = 6, height = 6)
print(plot_cell_trajectory(pal_cds,color_by = "State",cell_size = 3) + scale_color_manual(values = my_colors[1:5]) + theme(text = element_text(size = 22),legend.text = element_text(size = 28)))
dev.off()

png("./figs/pal_monocle_pseudotime.png", units = "in", res = 300, width = 6, height = 6)
print(plot_cell_trajectory(pal_cds,color_by = "Pseudotime",cell_size = 3) + scale_color_gradient(high = my_colors[3],low = my_colors[2],breaks = c(0,5,10)) + theme(text = element_text(size = 22)))
dev.off()

png("./figs/pal_monocle_scdeal_prob.png", units = "in", res = 300, width = 6, height = 6)
print(plot_cell_trajectory(pal_cds,color_by = "scDEAL_pred",cell_size = 3) + scale_color_gradient(high = "tomato",low = "skyblue",breaks = c(0,0.5,1),limits = c(0,1),name = "Sensitivity prob") + theme(text = element_text(size = 22)))
dev.off()

pal_state_genes <- differentialGeneTest(pal_cds,fullModelFormulaStr = "~State",cores = 4)
pal_sig_state_genes <- pal_state_genes %>% filter(qval < 1e-10) %>% arrange(qval)
pal_sig_genes_htmap <- pal_state_genes %>% filter(qval < 1e-10) %>% slice_min(order_by = qval,n = 20) %>% rownames()

tram_count_df <- as(as.matrix(tram_dmso@assays$RNA@counts),"sparseMatrix")
tram_pd <- new("AnnotatedDataFrame",data = tram_dmso@meta.data)
tram_fData <- data.frame(gene_short_name = rownames(tram_count_df),row.names = rownames(tram_count_df))
tram_fd <- new("AnnotatedDataFrame",data = tram_fData)
tram_cds <- newCellDataSet(tram_count_df,phenoData = tram_pd,featureData = tram_fd,expressionFamily = negbinomial.size())
tram_cds <- estimateSizeFactors(tram_cds)
tram_cds <- estimateDispersions(tram_cds)
tram_ordering_genes <- tram_enrich_results$degs %>% filter(p_val < 0.05) %>% slice_max(order_by = abs(avg_log2FC),n = 1000) %>% rownames()
tram_cds <- setOrderingFilter(tram_cds,tram_ordering_genes)
tram_cds <- reduceDimension(tram_cds,max_components = 2,method = "DDRTree")
tram_cds <- orderCells(tram_cds,root_state = 2)

tram_valid_cells <- rownames(subset(pData(tram_cds),State != "9"))
tram_cds_clean <- tram_cds[,tram_valid_cells]
tram_cds_clean <- setOrderingFilter(tram_cds_clean,tram_ordering_genes)
tram_cds_clean <- reduceDimension(tram_cds_clean,max_components = 2,method = "DDRTree")
tram_cds_clean <- orderCells(tram_cds_clean)
pData(tram_cds_clean)$State <- droplevels(pData(tram_cds_clean)$State)
levels(pData(tram_cds_clean)$State) <- seq_along(levels(pData(tram_cds_clean)$State))
keep_logic <- (pData(tram_cds_clean)$State != "3") | (pData(tram_cds_clean)$State == "3" & pData(tram_cds_clean)$scDEAL_pred < 0.3)
tram_cds_subset <- tram_cds_clean[,keep_logic]
tram_cds_subset <- tram_cds_subset[,rownames(pData(tram_cds_subset))]
tram_cds_subset@auxOrderingData[["DDRTree"]]$cell_pair_dist <- NULL
tram_cds_subset <- reduceDimension(tram_cds_subset,max_components = 2,method = "DDRTree")
tram_cds_subset <- orderCells(tram_cds_subset)

png("./figs/tram_monocle_state.png", units = "in", res = 300, width = 6, height = 6)
print(plot_cell_trajectory(tram_cds_subset,color_by = "State",cell_size = 3) + scale_color_manual(values = my_colors[1:6]) + theme(text = element_text(size = 22),legend.text = element_text(size = 28)))
dev.off()

png("./figs/tram_monocle_pseudotime.png", units = "in", res = 300, width = 6, height = 6)
print(plot_cell_trajectory(tram_cds_subset,color_by = "Pseudotime",cell_size = 3) + scale_color_gradient(high = my_colors[3],low = my_colors[2],breaks = c(0,5,10)) + theme(text = element_text(size = 22)))
dev.off()

png("./figs/tram_monocle_scdeal_prob.png", units = "in", res = 300, width = 6, height = 6)
print(plot_cell_trajectory(tram_cds_subset,color_by = "scDEAL_pred",cell_size = 3) + scale_color_gradient(high = "tomato",low = "skyblue",breaks = c(0,0.5,1),limits = c(0,1),name = "Sensitivity prob") + theme(text = element_text(size = 22)))
dev.off()

tram_state_genes <- differentialGeneTest(tram_cds_subset,fullModelFormulaStr = "~State",cores = 8)
tram_sig_state_genes <- tram_state_genes %>% filter(qval < 1e-10) %>% arrange(qval)
tram_sig_genes_htmap <- tram_state_genes %>% filter(qval < 1e-10) %>% slice_min(order_by = qval,n = 20) %>% rownames()

pal_heatmap_res <- plot_pseudotime_heatmap2(pal_cds[pal_sig_genes_htmap,],num_clusters = 3,cores = 8)
tram_heatmap_res <- plot_pseudotime_heatmap2(tram_cds_subset[tram_sig_genes_htmap,],num_clusters = 3,cores = 8)

set.seed(1)
pal_clusterGVis_gene <- sample(pal_heatmap_res$wide.res$gene,20,replace = FALSE)
tram_clusterGVis_gene <- sample(tram_heatmap_res$wide.res$gene,20,replace = FALSE)

pal_clusterGVis_enrich <- enrichCluster(object = pal_heatmap_res,OrgDb = org.Hs.eg.db,type = "BP",organism = "hsa",pvalueCutoff = 0.05,topn = 10)
pal_c3_filtered <- pal_clusterGVis_enrich[pal_clusterGVis_enrich$group == "C3" & !grepl("digestion|digestive",pal_clusterGVis_enrich$Description,ignore.case = TRUE),]
if(nrow(pal_c3_filtered) >= 2) pal_c3_filtered <- pal_c3_filtered[-2,]
pal_c1_df <- pal_clusterGVis_enrich[pal_clusterGVis_enrich$group == "C1",]
pal_c2_df <- pal_clusterGVis_enrich[pal_clusterGVis_enrich$group == "C2",]
pal_clusterGVis_enrich <- bind_rows(head(pal_c1_df,5),head(pal_c2_df,5),head(pal_c3_filtered,5))

tram_clusterGVis_enrich <- enrichCluster(object = tram_heatmap_res,OrgDb = org.Hs.eg.db,type = "BP",organism = "hsa",pvalueCutoff = 0.05,topn = 5)

text_colors <- c("blue4","red4","darkgreen","purple4","darkorange3","grey30")

png("./figs/pal_monocle_clusterGVis_compact_text.png", units = "in", res = 300, width = 12, height = 8)
my_visCluster(object = pal_heatmap_res,plotType = "both",column_names_rot = 45,show_row_dend = FALSE,markGenes = pal_clusterGVis_gene,markGenesSide = "left",pseudotimeCol = c(my_colors[2],my_colors[3]),annoTermData = pal_clusterGVis_enrich,ctAnnoCol = text_colors[1:3],goCol = rep(text_colors[1:3],each = 5),genesGp = c("italic",15,NA),goSize = 16,addBar = FALSE,textboxSize = 12,textboxPos = c(0.5,0.88),lineSide = "left")
dev.off()

png("./figs/tram_monocle_clusterGVis_compact_text.png", units = "in", res = 300, width = 12, height = 8)
my_visCluster(object = tram_heatmap_res,plotType = "both",column_names_rot = 45,show_row_dend = FALSE,markGenes = tram_clusterGVis_gene,markGenesSide = "left",pseudotimeCol = c(my_colors[2],my_colors[3]),annoTermData = tram_clusterGVis_enrich,ctAnnoCol = text_colors[1:3],goCol = rep(text_colors[1:3],each = 5),genesGp = c("italic",15,NA),goSize = 16,textboxSize = 12,textboxPos = c(0.5,0.88),addBar = FALSE,lineSide = "left")
dev.off()