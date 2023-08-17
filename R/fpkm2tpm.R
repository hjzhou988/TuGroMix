#' Converts FPKM data to TPM
#'
#' Converts FPKM (or log2FPKM) data to TPM (or log2TPM) data. We recommend use log2TPM as the input data for biomarker screening,
#' as it has better performance than log2FPKM in our case studies.
#'
#' @param fpkm The input FPKM data frame or matrix. Columns are samples, and genes are rows. All columns should be numeric and gene name/symbol column need to be excluded. It can be in the original scale or log2 scale. If in the log2 scale, need to set parameter fpkm.log = TRUE
#' @param log.fpkm Logical variable indicating whether the input FPKM data in the original scale or log2 scale.  If in the log scale, need to set parameter fpkm.log = TRUE; otherwise set fpkm.log = FALSE
#' @param log.tpm Logical variable indicating whether the output TPM should be in the original scale or log2 scale. If TRUE, the output is in the log2 scale; otherwise, the output is in the original scale.
#' @param offset A value that was added onto the original fpkm to prevent log2(fpkm) from being negative infinity. The default is 0.
#' @return Returns a tpm data frame with the first column as gene ID and the rest columns are mouse models.
#'
#' @export
#'
#'
fpkm2tpm = function(fpkm, gene.id, log.fpkm = T, log.tpm = T, offset = 0){
  # fpkm = cdx.exp[,-1]; gene.id = cdx.exp[[1]]; log.fpkm = T; log.tpm = T

  if (log.fpkm==T){
    fpkm = 2^(fpkm)
    fpkm = fpkm- offset
  }

  fpkm.sum = base::colSums(fpkm)

  tpm = fpkm/fpkm.sum*10^6

  if (log.tpm == T){
    tpm = base::log2(tpm+1)
  }

  tpm  = base::as.data.frame(tpm)

  tpm = base::cbind(Gene.ID = gene.id,tpm)

  tpm
}
