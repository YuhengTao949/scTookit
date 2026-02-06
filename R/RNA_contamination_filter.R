#' RNA Contamination Filtering
#'
#' Estimates and filters out contamination from ambient RNA using DecontX algorithm.
#' Adds contamination scores and multiple binary classification thresholds to metadata.
#'
#' @param seu A Seurat object containing single-cell RNA sequencing data
#'
#' @return Modified Seurat object with added columns in meta.data:
#'   - contamination_score: Continuous contamination probability (0-1)
#'   - contamination_th_0.1 to contamination_th_0.5: Binary classifications at
#'     different thresholds labeling cells as "pure_cell" or "contaminated_cell"
#'
#' @importFrom Seurat GetAssayData
#' @importFrom SingleCellExperiment colData SingleCellExperiment
#' @importFrom decontX decontX
#' @export
#'
#' @examples
#' seu <- RNA_contamination_filter(seu)
RNA_contamination_filter <- function(seu){

  counts <- GetAssayData(object = seu, slot = "counts")
  sce <- SingleCellExperiment(list(counts = counts))
  sce <- decontX::decontX(sce)

  seu$contamination_score <- SingleCellExperiment::colData(sce)$decontX_contamination

  thresholds <- c(0.1, 0.2, 0.3, 0.4, 0.5)
  for (th in thresholds) {
    col_name <- paste0("contamination_th_", th)
    seu[[col_name]] <- ifelse(seu$contamination_score < th,
                              "pure_cell",
                              "contaminated_cell")
  }
  return(seu)
}
