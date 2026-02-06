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
#' seurat_obj <- sample_qc(seurat_obj)
sample_qc <- function(sc){
  mito_genes <- rownames(sc@assays$RNA)[grep("^MT-", rownames(sc@assays$RNA), ignore.case = TRUE)]
  if (length(mito_genes) > 0) {
    print("Mitochondrial genes:")
    print(mito_genes)
    sc <- PercentageFeatureSet(sc, features = mito_genes, col.name = "percent_mito_by_sample")
  } else {
    print("Warning: No mitochondrial gene was found in this Seurat object.")
  }

  ribo_genes <- rownames(sc@assays$RNA)[grep("^RP[SL]", rownames(sc@assays$RNA), ignore.case = TRUE)]
  if (length(ribo_genes) > 0) {
    print("Ribosomal genes:")
    print(ribo_genes)
    sc <- PercentageFeatureSet(sc, features = ribo_genes, col.name = "percent_ribo_by_sample")
  } else {
    print("Warning: No ribosomal genes were found in this Seurat object.")
  }

  hb_genes <- rownames(sc@assays$RNA)[grep("^HB", rownames(sc@assays$RNA), ignore.case = TRUE)]
  if (length(hb_genes) > 0) {
    print("Red blood cell genes:")
    print(hb_genes)
    sc <- PercentageFeatureSet(sc, features = hb_genes, col.name = "percent_hb_by_sample")
  } else {
    print("Warning: Red blood cell genes were found in this Seurat object.")
  }

  hsp_genes <- rownames(sc@assays$RNA)[grep("^HSP", rownames(sc@assays$RNA), ignore.case = TRUE)]
  if (length(hsp_genes) > 0) {
    print("Heat shock protein genes:")
    print(hsp_genes)
    sc <- PercentageFeatureSet(sc, features = hsp_genes, col.name = "percent_hsp_by_sample")
  } else {
    print("Warning: Heat shock protein genes were found in this Seurat object.")
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
    print("Dissociation genes genes:")
    print(ds_genes)
    sc <- PercentageFeatureSet(sc, features = ds_genes, col.name = "percent_dissociation_by_sample")
  } else {
    print("Warning: Dissociation genes were found in this Seurat object.")
  }

  return(sc)
}


#' Get upper cuttoff
#' @export
get_upper_cutoff <- function(x) {
  if(length(x) < 10 || all(x == 0)) return(max(x, na.rm = T))
  return(boxplot.stats(x)$stats[5])
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
qc_Plot <- function(data, colname, Breaks = seq(0, 100, by = 5),
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
