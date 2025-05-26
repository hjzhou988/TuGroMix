#' get eGR ratio or difference within a trial
#'
#' @param tv.data tumor volume data in long format, containing columns "Group","Mouse","Day","TV"
#' @param ref.group the reference group in eGR ratio or difference
#' @param group.id.var The name of the "group" column of your input data. Required.
#' @param mouse.id.var The name of the "mouse.id" column of your input data. Required.
#' @param time.var The name of the "day" column of your input data. Required.
#' @param tv.var  The name of the "tumor volume" column of your input data. Required.
#' @param ci a logical variable to do bootstrapping for confidence interval
#' @param type calculate eGR "ratio" or "difference"
#'
#' @return a data frame containing eGR ratio/difference for each treatment group
#'
#' @importFrom magrittr %>%
#' @importFrom dplyr bind_rows
#' @export
#'
get_Model_eGR=function(tv.data,ref.group,group.id.var,mouse.id.var,time.var,tv.var,ci = F, type = c("difference","ratio")){
  # tv.data = tv_long %>% filter(Model=="GA0006"); ref.group = "Vehicle"; group.id.var = "Trt"; mouse.id.var = "Mouse.ID"; time.var = "Day"; tv.var = "TV";
  type = match.arg(type)
  tv.data = dplyr::rename(tv.data, Group = {{group.id.var}}, Mouse = {{mouse.id.var}}, Day = {{time.var}}, TV = {{tv.var}})

  tv.data = tv.data[,c("Group","Mouse","Day","TV")]
  tv.data = tv.data[stats::complete.cases(tv.data),]
  # tv.data = tv.data %>% dplyr::mutate(Mouse = base::paste(Group,Mouse,sep="_")) # so that mouse ID are unique in one experiment.

  # tv.data$Mouse = base::paste(tv.data[["Group"]],tv.data[["Mouse"]],sep="_")
  Grp.dat = base::split(tv.data,tv.data$Group)
  all.eGR = base::lapply(Grp.dat, function(g){
    mouses=base::split(g,g$Mouse)
    auc_mouse=base::do.call(base::rbind,base::lapply(mouses,get_eGR))
    # get_group_eGR(tv.data = g,mouse.id.var ="Mouse",time.var = "Day",tv.var = "TV",ci = F)
  })
  all.eGR.df = dplyr::bind_rows(all.eGR,.id = "Group")


  grps=base::unique(tv.data$Group)
  trt_grps=base::sort(grps[!(grps %in% ref.group)])
  ref_grp = base::setdiff(grps,trt_grps)
  eGR.res = base::lapply(trt_grps,function(g) geteGRRatio(all.eGR.df,g,ref_grp,type = type ,ci=ci))
  names(eGR.res) = trt_grps
  eGR.res.df =dplyr::bind_rows(eGR.res, .id = "Group")
  if (length(eGR.res.df)==0) return(NULL) else eGR.res.df
}
