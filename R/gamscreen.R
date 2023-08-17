#' Uses generalized additive mixed model to screen potential biomarker genes
#'
#' Uses GAMM to screen potential biomarker genes one by one. It quotes gamm4() function from gamm4 package to fit growth curves.
#'
#' @importFrom magrittr %>%
#' @importFrom rlang :=
#' @import foreach
#'
#' @param tv.data tv data frame in the "long" format, must contain a "Model"
#' column that matches the sample names of gene expression data. Must have mouse ID variable that is unique at least within the experimental group.
#' @param ref.group The reference group, usually a control/vehicle group for the drug group to compare with. Required.
#' @param model.var The name of the mouse model column. Do not add quotation marks.
#' @param group.var The name of the group column. Do not add quotation marks.
#' @param mouse.id.var The name of the mouse ID column. Do not add quotation marks.
#' @param gene.expr.data gene expression data in the "fat" format, in which every column is a gene, and every row is a mouse model.
#' The first column must be the sample names.
#'
#'
#' @export
#'
gamscreen = function(tv.data, ref.group, model.var, group.var, mouse.id.var,gene.expr.data){
  # tv.data = read.csv("../../Irinotecan_for_validation/TV_data_combined.csv",row.names = 1)
  # expre = readRDS("../../Irinotecan_for_validation/Irinotecan_RSEM_log2TPM_selected_genes.RDS")
  # expre = expre[,c(ncol(expre),1:(ncol(expre)-1))]
  # gene.expr.data = expre
  # gene.expr.data$Model = sub("P[0-9]+$","", rownames(gene.expr.data))

  tv.data = tv.data %>% dplyr::mutate({{group.var}}:=base::as.factor({{group.var}}))
  tv.data = tv.data %>% dplyr::mutate({{group.var}}:=stats::relevel({{group.var}},ref = ref.group))
  tv.data = tv.data %>% dplyr::arrange({{ model.var }},{{ group.var }},{{ mouse.id.var }},Day)
  tv.data = tv.data %>%
    dplyr::group_by({{ model.var }},{{ group.var }},{{ mouse.id.var }}) %>% dplyr::mutate(TV0 = dplyr::first(TV))
  tv.data = tv.data %>% dplyr::mutate(Trt = {{group.var}} )
  # tv.data <- transform(tv.data, fTrt = ordered(Trt,levels = c("Vehicle", "Centuximab")))
  tv.data = tv.data %>% tidyr::unite(col="Mouse.ID.2",{{model.var}},{{group.var}},{{mouse.id.var}},remove = FALSE)
  tv.data = tv.data %>% dplyr::mutate({{model.var}}:= base::gsub("-",".",{{model.var}}))
  # tv.data = tv.data %>% rename(Model = model.var)
  tv.data = tv.data %>% dplyr::rename( Model:={{model.var}}) # renamed to "Model"

  # print(head(tv.data))

  base::colnames(gene.expr.data)[1]="Model"
  genes = base::colnames(gene.expr.data)
  base::colnames(gene.expr.data) = base::gsub(pattern = "-",replacement = ".",genes)
  genes = base::colnames(gene.expr.data)[-1]  # get rid of the first column
  # print(head(genes))
  # tv.data

  GAM_fit = function(gene,logCPM,tv.long,REML){
    # i = 1
  # gene= genes[1];logCPM=gene.expr.data;tv.long=tv.data;REML=FALSE
    dat1 = tv.data %>% dplyr::left_join(logCPM[,c(gene,"Model")],by = "Model")
    # print(head(dat1))
    f.string.2 = base::paste0("log(TV+1)~ Trt + s(",gene,",Day,by=Trt)")
    fit3.ML = NULL
    base::tryCatch(fit3.ML <- gamm4::gamm4(stats::as.formula(f.string.2),random= ~(1+Day|Model/Mouse.ID.2),data = dat1,REML = REML),
             error = function(e){print(paste("Error:",e))})
    if(!is.null(fit3.ML)){
      fit3.ML.sum = base::summary(fit3.ML[["gam"]])
      p_table = fit3.ML.sum$p.table
      # fit3.ML.sum.mer = summary(fit3.ML[["mer"]])
      out = base::list(Gene = gene,
           Coef= p_table[nrow(p_table),1],
           t.val = p_table[nrow(p_table),3],
           p.val = p_table[nrow(p_table),4]
           # AIC = fit3.ML.sum.mer$AICtab[[1]],
           # BIC = fit3.ML.sum.mer$AICtab[[2]],
           # logLik = fit3.ML.sum.mer$AICtab[[3]]
           )
      base::print(out)
    }
  }
  #
  numCores <- parallel::detectCores()
  doParallel::registerDoParallel(numCores-1)

  fit_res = foreach::foreach(i=1:length(genes)) %dopar%{ #
    GAM_fit(genes[i],gene.expr.data,tv.data,FALSE)
  }

  fit_res.df = dplyr::bind_rows(fit_res)

  fit_res.df = fit_res.df %>% dplyr::arrange(p.val)

  fit_res.df
}
# tv.data = read.csv("../../Irinotecan_for_validation/TV_data_combined.csv",row.names = 1)
# tv.data = tv.data %>% arrange(Model)
# tv.data = tv.data[1:535,]
# tv.data.screen = tv.data
# expre = readRDS("../../Irinotecan_for_validation/Irinotecan_RSEM_log2TPM_selected_genes.RDS")
# expre = expre[,c(ncol(expre),1:(ncol(expre)-1))]
# expre = expre[,1:20]
# expre = expre %>% filter(Model %in% tv.data$Model)
# expre.screen = expre
# save(tv.data.screen,expre.screen,file="screen.RData")
# test = gamscreen(tv.data.screen,"Vehicle",Model,Trt,Mouse.ID,expre.screen[,1:5])

