#' Uses eGR to screen potential biomarker genes
#'
#' @import foreach
#' @importFrom magrittr %>%
#'
#' @param tv.data tv data frame in the "long" format, must contain a "Model"
#' column that matches the sample names of gene expression data. Must have mouse ID
#' variable that is unique at least within the experimental group.
#' @param ref.group The reference group, usually a control/vehicle group for the drug group to compare with. Required.
#' @param model.var The name of the mouse model column. Required.
#' @param group.var The name of the group column. Required.
#' @param mouse.id.var The name of the mouse ID column. Required.
#' @param time.var The name of the time column. Required.
#' @param tv.var The name of the TV column. Required.
#' @param gene.expr.data gene expression data in the "fat" format, in which every column is a gene, and every row is a mouse model.
#' The first column must be the sample names.Required.
#'
#' @return a list of eGR difference or ratio result for each model, and gene list ordered by p-value
#'
#' @export
#'
eGRscreen = function(tv.data, ref.group, model.var, group.var, mouse.id.var,time.var,tv.var,gene.expr.data, type = c("difference","ratio")){
  # tv.data = tv.data.screen; ref.group = "Vehicle"; model.var = "Model"; group.var="Trt";mouse.id.var="Mouse.ID"; gene.expr.data = transpose(expres.screen, gene.id = rownames(expres.screen))

  type = match.arg(type)

  tv.data = dplyr::rename(tv.data, Group = {{group.var}},
                          Mouse = {{mouse.id.var}},
                          Model = {{model.var}},
                          Day = {{time.var}},
                          TV = {{tv.var}})

  tv.data = tv.data %>% dplyr::mutate(Group = base::as.factor(Group))
  tv.data = tv.data %>% dplyr::mutate(Group = stats::relevel(Group,ref = ref.group))

  tv.data = tv.data %>% dplyr::arrange(Model,Group,Mouse,Day)

  proj_lst = split(tv.data,tv.data$Model)

   # Calculate eGR diff or ratio
  auc.res = lapply(proj_lst, function(x) get_Model_eGR(x,
                                                       ref.group="Vehicle",
                                                       group.id.var = "Group",
                                                       mouse.id.var = "Mouse",
                                                       time.var = "Day",
                                                       tv.var = "TV",
                                                       ci = F,
                                                       type = type))
  auc.res.df = dplyr::bind_rows(auc.res, .id = "Model")

  # spearman correlation
  base::colnames(gene.expr.data)[1]="Model"
  auc.res.med.expr = auc.res.df %>% dplyr::left_join(gene.expr.data, by = "Model")
  genes = base::names(auc.res.med.expr)[-(1:6)]

  numCores <- parallel::detectCores()
  doParallel::registerDoParallel(numCores-1)
  cor.res = foreach(i = 1:length(genes)) %dopar% {
    # i = 1
    r = stats::cor.test(x = auc.res.med.expr[[3]], # either difference or ratio
                 y = auc.res.med.expr[[genes[i]]],
                 method = "spearman",exact = F)
    base::list(Gene = genes[i],rho = r$estimate, p.val = r$p.value)
  }
  cor.res = dplyr::bind_rows(cor.res)
  cor.res = cor.res %>% dplyr::arrange(p.val)

  list(cor = cor.res,
       eGR = auc.res.df)

}

#test = lmmscreen(tv.data.screen,"Vehicle",Model,Trt,Mouse.ID,expre.screen[,1:5],REML=F)
# test.no.ri = lmmscreen(tv.data.screen,"Vehicle",Model,Trt,Mouse.ID,expre.screen[,1:5],REML=F)


