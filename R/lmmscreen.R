#' Uses linear mixed-effects model to screen potential biomarker genes
#'
#' @import foreach
#' @importFrom magrittr %>%
#' @importFrom rlang :=
#'
#' @param tv.data tv data frame in the "long" format, must contain a "Model"
#' column that matches the sample names of gene expression data. Must have mouse ID
#' variable that is unique at least within the experimental group.
#' @param ref.group The reference group, usually a control/vehicle group for the drug group to compare with. Required.
#' @param model.var The name of the mouse model column. Do not add quotation marks.
#' @param group.var The name of the group column. Do not add quotation marks.
#' @param mouse.id.var The name of the mouse ID column. Do not add quotation marks.
#' @param gene.expr.data gene expression data in the "fat" format, in which every column is a gene, and every row is a mouse model.
#' The first column must be the sample names.
#'
#' @details
#' The statistical model for gene screening is log(TV+1)~ log(TV0+1) + Day + gene:Day + Trt:Day + gene:Trt:Day + (1+Day|Model/Mouse.ID.2)
#'
#' @export
#'
lmmscreen = function(tv.data, ref.group, model.var, group.var, mouse.id.var,gene.expr.data){
  tv.data = tv.data %>% dplyr::mutate({{group.var}}:=base::as.factor({{group.var}}))
  tv.data = tv.data %>% dplyr::mutate({{group.var}}:=stats::relevel({{group.var}},ref = ref.group))
  # tv.data = read.csv("../../Irinotecan_for_validation/TV_data_combined.csv",row.names = 1); model.var = "Model";group.var="Trt";mouse.id.var="Mouse.ID";gene.expr.data = readRDS("../../Irinotecan_for_validation/Irinotecan_RSEM_log2TPM_selected_genes.RDS")
  # gene.expr.data$Model = sub("P[0-9]+$","", rownames(gene.expr.data))


  tv.data = tv.data %>% dplyr::arrange({{ model.var }},{{ group.var }},{{ mouse.id.var }},Day)
  tv.data = tv.data %>%
    dplyr::group_by({{ model.var }},{{ group.var }},{{ mouse.id.var }}) %>% dplyr::mutate(TV0 = dplyr::first(TV))
  tv.data = tv.data %>% dplyr::mutate(Trt = {{group.var}} )

  # tv.data <- transform(tv.data, fTrt = ordered(Trt,levels = c("Vehicle", "Centuximab")))
  tv.data = tv.data %>% tidyr::unite(col="Mouse.ID.2",{{model.var}},{{group.var}},{{mouse.id.var}},remove = FALSE)
  tv.data = tv.data %>% dplyr::mutate({{model.var}}:=gsub("-",".",{{model.var}}))
  # tv.data = tv.data %>% rename(Model = model.var)
  tv.data = tv.data %>% dplyr::rename( Model:={{model.var}}) # renamed to "Model"
  # print(head(tv.data))
  base::colnames(gene.expr.data)[1]="Model"

  genes = base::colnames(gene.expr.data)
  base::colnames(gene.expr.data) = base::gsub(pattern = "-",replacement = ".",genes)
  genes = base::colnames(gene.expr.data)[-1]  # get rid of the first column
  # print(head(genes))
  # tv.data
  LMM_fit = function(gene,logCPM,tv.data,REML){
    # gene=genes[766];logCPM=logcpm.t;tv.data = tv.data;REML=FALSE
    # gene=genes[2]; logCPM=logcpm.t;tv.data = tv.data;REML=FALSE
    # print("not_joined")
    dat1 = tv.data %>% dplyr::left_join(logCPM[,c(gene,"Model")],by = "Model")
    # print(head(dat1))
    exp.str.tv0 = base::paste0("log(TV+1)~ log(TV0+1)+Day + ",gene,":Day + Trt:Day + ",gene,":Trt:Day + (1+Day|Model/Mouse.ID.2)")
    fit3.ML.tv0= NULL

    tryCatch(fit3.ML.tv0 <- lmerTest::lmer(stats::as.formula(exp.str.tv0),data = dat1,REML = REML),error=function(e){print(paste("Error:",e))} )
    if(!is.null(fit3.ML.tv0)){
      ms=base::summary(fit3.ML.tv0)
      coefs=ms$coefficients
      out=base::list(Gene = gene,
               Term = rownames(coefs)[nrow(coefs)],
               Coef=coefs[nrow(coefs),1],
               p.val=coefs[nrow(coefs),ncol(coefs)],
               AIC=ms$AICtab[[1]],
               BIC=ms$AICtab[[2]],
               logLik=ms$AICtab[[3]])
      out
    }
  }
  #
  numCores <- parallel::detectCores()
  doParallel::registerDoParallel(numCores-1)

  fit_res = foreach::foreach(i=1:length(genes)) %dopar%{ #
    # i = 2
    LMM_fit(genes[i],gene.expr.data,tv.data,FALSE)
  }

  fit_res.df = dplyr::bind_rows(fit_res)

  fit_res.df = fit_res.df %>% dplyr::arrange(p.val)

  fit_res.df
}

#test = lmmscreen(tv.data.screen,"Vehicle",Model,Trt,Mouse.ID,expre.screen[,1:5],REML=F)
# test.no.ri = lmmscreen(tv.data.screen,"Vehicle",Model,Trt,Mouse.ID,expre.screen[,1:5],REML=F)


