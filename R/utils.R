#' Create nested subdirectories for project organization
#'
#' This function creates standardized subdirectory structures within specified parent directories.
#'
#' @param parent_path Character string specifying the root path where directories will be created (default = "analysis")
#' @param dirs Character vector naming parent directories to create (default = c("data", "figure"))
#' @param subdir Character string naming the target subdirectory to create within each parent directory
#'
#' @return Invisibly returns NULL. Function is called for side effects (directory creation).
#' @details Recursively creates paths of form: `<parent_path>/<dir>/<subdir>`
#' @export
#' @examples
#' create_dirs(parent_path = "projects/analysis1", subdir = "experimentA")
#' create_dirs(parent_path = "reports", dirs = c("output", "docs"), subdir = "2023Q4")
create_dirs <- function(subdir, parent_path = "analysis", dirs = c("data", "figure")) {
  for (dir in dirs) {
    target_dir <- file.path(parent_path, dir, subdir)

    if (dir.exists(target_dir)) {
      warning("Directory already exists, skipping creation: ", target_dir)
    } else {
      dir.create(target_dir, recursive = TRUE)
      message("Successfully created directory: ", target_dir)
    }
  }
}
