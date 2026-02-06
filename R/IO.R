#' Read and Filter 10X Genomics Data with Gene Symbol Unification
#'
#' Processes 10X Genomics scRNA-seq data by reading count matrices and feature annotations,
#' ensuring unique gene symbol representation. Filters out low-quality cells and genes,
#' while resolving duplicate gene symbols by retaining the highest expressed transcript.
#'
#' @param filename Character string specifying the sample name. Used as a prefix
#'                for cell barcodes and as the project name in the Seurat object.
#' @param data_path Character string specifying the directory containing 10X data files
#'                 (matrix.mtx, barcodes.tsv, and features file).
#' @param features_filename Character string specifying the name of the features file
#'                         (e.g., "features.tsv.gz" or "genes.tsv") within \code{data_path}.
#'
#' @return A \code{\link[Seurat]{SeuratObject}} containing filtered and gene-unified count data.
#'         Cell barcodes are prefixed with \code{filename}, and gene symbols are used as row names.
#'
#' @importFrom Seurat Read10X CreateSeuratObject
#' @importFrom Matrix rowMeans
#' @importFrom dplyr filter arrange desc distinct
#' @importFrom utils read.table
#' @export
#'
#' @examples
#' \dontrun{
#' seurat_obj <- read_10x_genefilter(
#'   filename = "sample1",
#'   data_path = "/path/to/data",
#'   features_filename = "features.tsv.gz"
#' )
#' }
read_10x_genefilter <- function(filename, data_path, features_filename) {

  scdata <- Seurat::Read10X(data.dir = data_path, gene.column = 1)
  colnames(scdata) <- paste(filename, colnames(scdata), sep = "_")

  feat_path <- file.path(data_path, features_filename)
  if (!file.exists(feat_path)) stop("File not found: ", feat_path)

  genes <- read.table(
    gzfile(feat_path),
    header = FALSE,
    sep = "\t",
    stringsAsFactors = FALSE
  )

  gene_info <- data.frame(
    ensembl_id = rownames(scdata),
    symbol = genes$V2[match(rownames(scdata), genes$V1)],
    mean_expression = Matrix::rowMeans(scdata),
    stringsAsFactors = FALSE
  )

  unique_genes <- gene_info %>%
    dplyr::filter(!is.na(symbol) & symbol != "") %>%
    dplyr::arrange(dplyr::desc(mean_expression)) %>%
    dplyr::distinct(symbol, .keep_all = TRUE)

  scdata_filtered <- scdata[unique_genes$ensembl_id, , drop = FALSE]
  rownames(scdata_filtered) <- unique_genes$symbol

  seu_obj <- Seurat::CreateSeuratObject(
    counts = scdata_filtered,
    project = filename,
    min.cells = 3,
    min.features = 200
  )

  return(seu_obj)
}
