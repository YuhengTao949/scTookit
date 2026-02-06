#' Sample Quality Control
#'
#' Calculates various quality control metrics for each sample including
#' mitochondrial, ribosomal, hemoglobin, heat shock protein, and dissociation
#' gene percentages. Adds these metrics as new columns to the meta.data slot
#' of the input Seurat object.
#'
#' @param sc A Seurat object containing single-cell RNA sequencing data
#'
#' @return Modified Seurat object with added QC metrics in meta.data:
#'   - percent_mito_by_sample (if MT genes present)
#'   - percent_ribo_by_sample (if RP genes present)
#'   - percent_hb_by_sample (if HB genes present)
#'   - percent_hsp_by_sample (if HSP genes present)
#'   - percent_dissociation_by_sample (if dissociation genes present)
#'
#' @importFrom Seurat PercentageFeatureSet
#' @export
#'
#' @examples
#' seurat_obj <- calcQCMetrics(seurat_obj)
calcQCMetrics <- function(sc){
  mito_genes <- rownames(sc@assays$RNA)[grep("^MT-", rownames(sc@assays$RNA), ignore.case = TRUE)]
  if (length(mito_genes) > 0) {
    cat("\nCalculating percentage of mitochondrial gene expression \n")
    sc <- PercentageFeatureSet(sc, features = mito_genes, col.name = "percent_mito_by_sample")
  } else {
    stop("\nError: No mitochondrial gene was found in this Seurat object.\n")
  }

  ribo_genes <- rownames(sc@assays$RNA)[grep("^RP[SL]", rownames(sc@assays$RNA), ignore.case = TRUE)]
  if (length(ribo_genes) > 0) {
    cat("Calculating percentage of ribosomal gene expression\n")
    sc <- PercentageFeatureSet(sc, features = ribo_genes, col.name = "percent_ribo_by_sample")
  } else {
    stop("Error: No ribosomal genes were found in this Seurat object.\n")
  }

  hb_genes <- rownames(sc@assays$RNA)[grep("^HB", rownames(sc@assays$RNA), ignore.case = TRUE)]
  if (length(hb_genes) > 0) {
    cat("Calculating percentage of hemoglobin gene expression\n")
    sc <- PercentageFeatureSet(sc, features = hb_genes, col.name = "percent_hb_by_sample")
  } else {
    stop("Error: Red blood cell genes were found in this Seurat object.\n")
  }

  hsp_genes <- rownames(sc@assays$RNA)[grep("^HSP", rownames(sc@assays$RNA), ignore.case = TRUE)]
  if (length(hsp_genes) > 0) {
    cat("Calculating percentage of heat shock protein gene expression\n")
    sc <- PercentageFeatureSet(sc, features = hsp_genes, col.name = "percent_hsp_by_sample")
  } else {
    stop("Error: Heat shock protein genes were found in this Seurat object.\n")
  }

  ds_genes <- c(
    "ATF3", "BTG2", "CEBPB", "CEBPB-AS1", "CEBPD", "CXCL1", "EGR1",
    "FOS", "FOSB", "FOSL1", "FOSL1P1", "FOSL2", "ID3", "IER2",
    "JUN", "JUNB", "JUND", "MT1A", "MT1B", "MT1E", "MT1F", "MT1G",
    "MT1H", "MT1L", "MT1M", "MT1X", "MT2A", "NFKBIA", "NR4A1",
    "PPP1R15A", "SOCS3", "UBC", "ZFP36"
  )
  ds_genes <- ds_genes[ds_genes %in% rownames(sc@assays$RNA)]
  if (length(ds_genes) > 0) {
    cat("Calculating percentage of dissociation gene expression\n")
    sc <- PercentageFeatureSet(sc, features = ds_genes, col.name = "percent_dissociation_by_sample")
  } else {
    stop("Error: Dissociation genes were found in this Seurat object.\n")
  }

  return(sc)
}


#' Get upper cuttoff
#' @export
get_upper_cutoff <- function(x) {
  if(length(x) < 10 || all(x == 0)) return(max(x, na.rm = T))
  return(boxplot.stats(x)$stats[5])
}

#' Generate QC Overview Boxplots for Multiple Samples
#'
#' This function creates boxplots to visualize quality control metrics across all samples
#' in a list of Seurat objects. It provides a convenient way to compare the distribution
#' of QC metrics between different samples.
#'
#' @param seulist A named list of Seurat objects where each object represents a single sample.
#'                The list names are used as sample identifiers.
#' @param metric Character string specifying the QC metric to visualize (e.g., "nCount_RNA",
#'               "percent_mito_by_sample"). Must be a column in the meta.data slot of
#'               all Seurat objects.
#' @param y_limits Optional numeric vector of length 2 specifying the y-axis limits
#'                 (e.g., c(0, 100)). If NULL (default), limits are determined automatically.
#'
#' @return A ggplot object containing the boxplot visualization of the specified QC metric
#'         across all samples.
#'
#' @export
#'
#' @importFrom ggplot2 ggplot aes geom_boxplot labs theme_classic theme coord_cartesian
#' @importFrom ggplot2 element_text
#' @importFrom stats setNames
#'
#' @examples
#' # Create list of Seurat objects
#' seulist <- list(sample1 = seurat_obj1, sample2 = seurat_obj2)
#' names(seulist) <- c("sample1", "sample2")
#'
#' # Generate plot for mitochondrial percentage
#' PlotMultiQC(seulist, "percent_mito_by_sample")
#'
#' # Generate with y-axis limits
#' PlotMultiQC(seulist, "nCount_RNA", y_limits = c(0, 30000))
#'
PlotMultiQC <- function(seulist, metric, y_limits = NULL) {
  sample_names <- names(seulist)
  for (sample_name in sample_names) {
    if (!metric %in% colnames(seulist[[sample_name]]@meta.data)) {
      stop(paste0("Metric '", metric, "' not found in sample: ", sample_name))
    }
  }

  plot_data <- data.frame()
  for (sample_name in sample_names) {
    sample_df <- data.frame(
      Value = seulist[[sample_name]]@meta.data[[metric]],
      Sample = sample_name
    )
    plot_data <- rbind(plot_data, sample_df)
  }

  p <- ggplot(plot_data, aes(x = Sample, y = Value)) +
    geom_boxplot(
      color = "black",
      fill = NA,
      width = 0.7,
      alpha = 0.8,
      outlier.shape = NA
    ) +
    labs(
      title = metric,
      x = "Sample",
      y = metric
    ) +
    theme_classic() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 12, face = "bold"),
      axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
      axis.text.y = element_text(size = 8),
      axis.title = element_text(size = 10)
    )

  if (!is.null(y_limits)) {
    p <- p + coord_cartesian(ylim = y_limits)
  }

  return(p)
}

#' Diagnostic QC Violin Plot
#'
#' Creates a customized violin plot with boxplot and jittered points for visualizing
#' a single quality control metric. Includes automatic threshold line based on
#' distribution characteristics.
#'
#' @param data Data frame containing the QC metric to visualize
#' @param colname Column name (as string) of the QC metric to plot on y-axis
#' @param Breaks Sequence of breaks for y-axis (default: seq(0, 100, by = 5))
#' @param violin_color Color of violin outline (default: "steelblue")
#' @param violin_fill Fill color of violin (default: "lightblue")
#' @param box_fill Fill color of boxplot (default: "white")
#' @param point_color Color of jittered points (default: "darkblue")
#' @param violin_side Side to draw violin ("r" = right, "l" = left, default: "r")
#' @param nudge_amount Horizontal nudge amount for violin position (default: 0.2)
#'
#' @return A ggplot object displaying the QC metric distribution with:
#'   - Half violin plot showing density distribution
#'   - Boxplot showing quartiles
#'   - Jittered points representing individual cells
#'   - Dashed red line indicating automatic threshold cutoff
#'
#' @importFrom ggplot2 ggplot aes geom_boxplot geom_jitter geom_hline labs
#'   theme_minimal theme scale_y_continuous element_blank element_text
#'   position_nudge
#' @importFrom gghalves geom_half_violin
#' @export
PlotQC <- function(data, colname, Breaks = seq(0, 100, by = 5),
                         violin_color = "steelblue", violin_fill = "lightblue",
                         box_fill = "white", point_color = "darkblue",
                         violin_side = "r", nudge_amount = 0.2) {

  upper_cutoff <- get_upper_cutoff(data[, colname])
  cat(colname, "cut off:", upper_cutoff, "\n")

  p <- ggplot(data, aes(x = "", y = .data[[colname]])) +
    geom_half_violin(
      aes(fill = factor(1)),
      side = violin_side,
      position = position_nudge(x = nudge_amount, y = 0),
      color = violin_color,
      fill = violin_fill,
      alpha = 0.7,
      linewidth = 0.5
    ) +
    geom_boxplot(
      fill = box_fill,
      width = 0.15,
      position = position_nudge(x = nudge_amount, y = 0),
      color = "black",
      outlier.shape = NA,
      linewidth = 0.4
    ) +
    geom_jitter(
      width = 0.1,
      height = 0,
      size = 0.2,
      alpha = 0.6,
      color = point_color
    ) +
    geom_hline(
      yintercept = upper_cutoff,
      color = "red4",
      linewidth = 0.5,
      linetype = "dashed"
    ) +
    labs(
      title = colname,
      x = "",
      y = ""
    ) +
    theme_minimal() +
    theme(
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank(),
      legend.position = "none",  # 隐藏图例
      plot.title = element_text(hjust = 0.5, face = "bold", size = 12)
    ) +
    scale_y_continuous(breaks = Breaks)

  return(p)
}


#' Sample-Level Quality Control and Filtering
#'
#' Performs comprehensive quality control by generating diagnostic violin plots for key QC metrics,
#' applies automated threshold-based filtering, and returns a filtered Seurat object.
#'
#' @param seu A Seurat object containing single-cell data with metadata
#' @param plotdir Directory path where the QC plot should be saved
#' @param filename Base filename for output plot (without file extension)
#'
#' @return A filtered Seurat object after applying QC thresholds. Original and remaining cell counts
#'         are printed to the console.
#'
#' @importFrom ggplot2 ggsave
#' @importFrom patchwork wrap_plots
#' @export
#'
#' @examples
#' filtered_seu <- sample_qc(seu = my_seurat,
#'                          plotdir = "qc_plots",
#'                          filename = "sample_qc_report")
sample_qc <- function(seu, plotdir, filename){
  metadata <- seu@meta.data
  data <- metadata[, c("nCount_RNA", "nFeature_RNA", "contamination_score",
                       "doublet_score", "percent_mito_by_sample", "percent_ribo_by_sample",
                       "percent_hb_by_sample", "percent_hsp_by_sample",
                       "percent_dissociation_by_sample")]
  p1 <- PlotQC(data, colname = "nCount_RNA", Breaks = seq(0,100000,by = 10000))
  p2 <- PlotQC(data, colname = "nFeature_RNA", Breaks = seq(0,10000,by = 1000))
  p3 <- PlotQC(data, colname = "contamination_score", Breaks = seq(0, 1, by = 0.05))
  p4 <- PlotQC(data, colname = "doublet_score", Breaks = seq(0, 1, by = 0.05))
  p5 <- PlotQC(data, colname = "percent_mito_by_sample", Breaks = seq(0, 100, by = 5))
  p6 <- PlotQC(data, colname = "percent_ribo_by_sample", Breaks = seq(0, 100, by = 5))
  p7 <- PlotQC(data, colname = "percent_hb_by_sample", Breaks = seq(0, 1, by = 0.05))
  p8 <- PlotQC(data, colname = "percent_hsp_by_sample", Breaks = seq(0, 100, by = 5))
  p9 <- PlotQC(data, colname = "percent_dissociation_by_sample", Breaks = seq(0, 100, by = 5))
  plot <- patchwork::wrap_plots(list(p1, p2, p3, p4, p5, p6, p7, p8, p9), nrow = 3)
  # print(plot)
  ggsave(filename = file.path(plotdir, paste0(filename, ".pdf")), plot, height = 9, width = 9)

  seu_new <- subset(seu, subset =
                      nCount_RNA < get_upper_cutoff(data[,"nCount_RNA"]) &
                      nFeature_RNA < get_upper_cutoff(data[,"nFeature_RNA"]) &
                      contamination_score < get_upper_cutoff(data[,"contamination_score"]) &
                      # doublet_score < get_upper_cutoff(data["doublet_score"])&
                      percent_mito_by_sample < get_upper_cutoff(data[,"percent_mito_by_sample"]) &
                      percent_hb_by_sample < get_upper_cutoff(data[,"percent_hb_by_sample"]) &
                      percent_hsp_by_sample < get_upper_cutoff(data[,"percent_hsp_by_sample"]) &
                      percent_dissociation_by_sample < get_upper_cutoff(data[,"percent_dissociation_by_sample"])

  )
  cat("\nRow cell number:", ncol(seu))
  cat("\nRemaining cell number:", ncol(seu_new))
  return(seu_new)
}
