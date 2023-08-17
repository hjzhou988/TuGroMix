#' Transposes gene expression data from thin (each row is a gene) to fat (each row is a sample)
#' @param data.wide gene expression data frame or matrix in the "thin" format, in which each row is a gene.
#'  The data should be all numeric, and exclude the gene id column.
#' @param gene.id a vector of gene ID strings that must be unique and matches the sequence of the genes.
#'
#'
#' @return returns a data frame with the first column as Model and the rest columns as gene expression data
#' @export
#'
transpose = function(data.thin, gene.id){
  data.fat = base::t(data.thin)
  data.fat = base::as.data.frame(data.fat)
  base::colnames(data.fat) = gene.id
  dplyr::bind_cols(Model = colnames(data.thin), data.fat)
}
