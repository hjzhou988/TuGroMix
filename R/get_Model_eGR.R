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
#' @export
#'
get_Model_eGR=function(tv.data,ref.group,group.id.var,mouse.id.var,time.var,tv.var,ci = F, type = c("difference","ratio")){
  #  tv.data = tv_wide.long;ref.group = "G1";ci = F; type = "ratio";group.id.var="Group";mouse.id.var="Animal.No.";time.var="Day";tv.var="TV"
  # tv.data = proj_lst[[2]];ref.group = 1
  # tv.data = proj_lst[[7]];ref.group = c("Grp-1A","Grp-2A","Grp-3A","Grp-4A","Grp-5A","Grp-6A")
  # tv.data = tv.data.split[[1]];ref.group="G1"



  type = match.arg(type)
  tv.data = tv.data[stats::complete.cases(tv.data),] # important to get rid of NA here.
  tv.data = dplyr::rename(tv.data, Group = {{group.id.var}}, Mouse = {{mouse.id.var}}, Day = {{time.var}}, TV = {{tv.var}})

  # tv.data = tv.data %>% dplyr::mutate(Mouse = base::paste(Group,Mouse,sep="_")) # so that mouse ID are unique in one experiment.
  tv.data$Mouse = base::paste(tv.data[["Group"]],tv.data[["Mouse"]],sep="_")
  mouses=base::split(tv.data,tv.data$Mouse) # it assumes that mouse ID is unique over all groups in one experiment. Shall I check whether there are same IDs in different groups in one experiment?
  auc_mouse=base::do.call(base::rbind,base::lapply(mouses,get_eGR))
  grps=base::unique(tv.data$Group)
  trt_grps=base::sort(grps[!(grps %in% ref.group)])
  ref_grp = base::setdiff(grps,trt_grps)

  eGR.res =dplyr::bind_rows(base::lapply(trt_grps,function(g) geteGRRatio(auc_mouse,g,ref_grp,type = type ,ci=ci)))
  if (length(eGR.res)==0) return(NULL) else eGR.res

}
