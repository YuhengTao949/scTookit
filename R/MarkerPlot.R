#' Visualize average expression and percentage of expressed genes by cluster
#'
#' This function generates a dot plot visualizing marker gene expression patterns
#' across cell clusters. Dot size represents the percentage of cells expressing a gene,
#' while color represents relative expression (scaled per gene across clusters).
#'
#' @param markerlist A named list of gene vectors. Each element name represents a functional
#'   group (e.g., pathways), and the vector contains associated marker genes.
#' @param seu A Seurat object containing single-cell data.
#' @param assay Name of assay to use (default: "RNA").
#' @param slot Data slot to use (default: "counts").
#' @param facet A logical value that controls whether faceted display is enabled.
#'   - If `TRUE` (default): enables faceting using the \code{names(markerlist)}.
#'   - If `FALSE`: generates a unified plot without faceting.
#'
#' @return A ggplot object showing marker expression patterns across clusters.
#'
#' @importFrom Seurat Idents WhichCells GetAssayData
#' @importFrom Matrix rowMeans rowSums
#' @importFrom ggplot2 ggplot aes geom_point scale_color_gradientn guide_colorbar
#'   scale_size_continuous facet_grid theme_minimal theme element_text element_rect
#'   element_blank xlab ylab
#' @importFrom RColorBrewer brewer.pal
#' @importFrom stats setNames
#'
#' @examples
#' # Example using Seurat object 'pbmc' and marker list
#' Idents(pbmc) <- pbmc$seurat_clusters
#' Idents(pbmc) <- factor(Idents(pbmc), levels = c("0", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12"))
#' marker_list <- list(
#'   T_cell = c("CD3D", "CD8A"),
#'   B_cell = c("CD79A", "MS4A1")
#' )
#' MarkDotplot(markerlist = marker_list, seu = pbmc)
#'
#' @export

MarkDotplot <- function(markerlist = NULL, seu = NULL, assay = "RNA", slot = "counts", facet = FALSE) {

  if (is.null(seu) || is.null(markerlist)) stop("Both seu and markerlist must be provided")

  genes_all <- unlist(markerlist)
  genes_valid <- intersect(genes_all, rownames(seu))
  clusters <- levels(Idents(seu))

  cluster_cells <- lapply(setNames(nm = clusters), function(cl) {
    WhichCells(seu, ident = cl)
  })

  counts_data <- GetAssayData(seu, assay = assay, slot = slot)

  result_data <- do.call(rbind, lapply(names(markerlist), function(func_name) {
    genes_func <- intersect(markerlist[[func_name]], genes_valid)
    if (length(genes_func) == 0) return(NULL)

    do.call(rbind, lapply(clusters, function(cl) {
      cells <- cluster_cells[[cl]]
      mat <- counts_data[genes_func, cells, drop = FALSE]

      avg_expr <- Matrix::rowMeans(mat)
      pct_expr <- Matrix::rowSums(mat > 0) / length(cells) * 100

      data.frame(
        Function = func_name,
        Gene = rep(genes_func, each = 1),
        Group = cl,
        AvgExpr = avg_expr,
        PctExpr = pct_expr,
        stringsAsFactors = FALSE
      )
    }))
  }))

  result_split <- split(result_data$AvgExpr, result_data$Gene)
  max_vals <- vapply(result_split, max, numeric(1))
  max_vals[max_vals <= 0] <- 1
  result_data$AvgExpr <- result_data$AvgExpr / max_vals[result_data$Gene]

  result_data$Function <- factor(result_data$Function, levels = names(markerlist))
  unique_genes <- rev(genes_all[!duplicated(genes_all)])
  result_data$Gene <- factor(result_data$Gene, levels = unique_genes)
  result_data$Group <- factor(result_data$Group, levels = rev(levels(Idents(seu))))

  p <- ggplot(result_data,
              aes(y = Group, x = Gene,
                  color = AvgExpr, size = PctExpr)) +
    geom_point() +
    scale_color_gradientn(
      colours = rev(brewer.pal(n = 10, name = "RdBu")),
      guide = guide_colorbar(ticks.colour = "black", frame.colour = "black"),
      name = "Relative\nExpression"
    ) +
    scale_size_continuous(
      name = "Percentage\nExpressed",
      range = c(0, 6)) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
      panel.border = element_rect(fill = NA, color = "black", size = 1.5),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      strip.background = element_rect(fill = "white", color = "black", size = 1.5),
      strip.text = element_text(face = "bold"),
      panel.background = element_rect(fill = "white", color = "black", size = 1.5),
      plot.background = element_rect(fill = "white")
    ) +
    xlab("Genes") +
    ylab("Cell Clusters")

  if (facet) {
    p <- p + facet_grid(
      . ~ Function,
      scales = "free_x",
      space = "free_x"
    )
  }

  return(p)
}

