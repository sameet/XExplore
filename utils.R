get_dir_name <- function(s) {
  if (dir.exists(dirname(s))) {
    dirname(s)
  } else {
    NULL
  }
}

make_board <- function(dir_name) {
  pins::board_folder(dir_name)
}

make_colors <- function(sce, col_name) {
  sce@meta.data |>
    dplyr::distinct(get(col_name)) |>
    dplyr::pull() |>
    sort() -> use_cells

  my_cols <- colorRampPalette(rev(RColorBrewer::brewer.pal(
    12,
    "Paired"
  )))(length(use_cells))
  names(my_cols) <- use_cells
  my_cols
}
