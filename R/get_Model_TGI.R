#' get tumor growth inhibition TGI for treatment groups
#'
#' @param tv.data tumor volume data in long format, containing columns "Group","Mouse","Day","TV"
#' @param ref.group the reference group in eGR ratio or difference
#' @param def definition of TGI. 1: TV_trt/TV_ctl, 2: DeltaTV_trt/DeltaTV_ctl, 3: 1-DeltaTV_trt/DeltaTV_ctl, 4:RTV_trt/RTV_ctl, 5: 1-TV_trt/TV_ctl, 6: 1-RTV_trt/RTV_ctl
#' @param type the type of the mean of TGI, choose from "arithmetic" and "geometric".
#' @param group.id.var The name of the "group" column of your input data. Required.
#' @param mouse.id.var The name of the "mouse.id" column of your input data. Required.
#' @param time.var The name of the "day" column of your input data. Required.
#' @param tv.var  The name of the "tumor volume" column of your input data. Required.
#' @param ci a logical variable to do bootstrapping for confidence interval
#' @param bound a number for the confidence interval boundary. Suggest 90 or 95.
#' @return a tibble containing TGI for each treatment group
#'
#' @importFrom magrittr %>%
#' @export
#'
get_Model_TGI = function(tv.data,ref.group,type = c("geometric","arithmetic"),def=6, group.id.var, mouse.id.var,time.var,tv.var,ci = F,bound = 90){
  type = match.arg(type)
  # proj_lst = split(tv_long, tv_long$Model)
  # proj_lst = split(tv.data.screen, tv.data.screen$Model)
  # tv.data = proj_lst2[[1]];ref.group="Vehicle";group.id.var="Trt";mouse.id.var="Mouse.ID";time.var="Day";tv.var="TV"
  # tv.data = readRDS("~/Documents/VCG/2024-8-1_example_data_for_tugromix_tgi_test_syn2.rds");ref.group = "Grp-1";def = 6;group.id.var = "Grp.Name";mouse.id.var = "session_key";time.var = "Day";tv.var = "TV";ci = T
  tv.data = dplyr::rename(tv.data, Group = {{group.id.var}}, Mouse = {{mouse.id.var}}, Day = {{time.var}}, TV = {{tv.var}})

  tv.data = tv.data[,c("Group","Mouse","Day","TV")]
  tv.data = tv.data[stats::complete.cases(tv.data),]

  tv.data$TV=as.numeric(tv.data$TV)

  # tv.data$Mouse = base::paste(tv.data[["Group"]],tv.data[["Mouse"]],sep="_")
  tv.data = tv.data %>% dplyr::arrange(Group, Mouse,Day)
  tv.data = tv.data %>% dplyr::group_by(Group, Mouse) %>% dplyr::mutate(TV0 = dplyr::first(TV))

  tv_var='TV'

  tv.data$RTV=ifelse(tv.data$TV0==0,tv.data$TV,tv.data$TV/tv.data$TV0)
  tv.data$DeltaTV=tv.data$TV-tv.data$TV0
  grps=base::unique(tv.data$Group)
  grps=base::sort(grps[grps!=ref.group])

  # tgi_df=base::as.data.frame(base::do.call(base::rbind,base::lapply(grps,function(g) tgi_group(tv.data,ref.group,g,def=def))))
  # tgi_df$Model=proj
  tgi_df=dplyr::bind_rows(base::lapply(grps,function(g) {
    tgi_group(tv.data,ref.group,g,type = type,def=def,ci = ci,bound = bound)
    })
  )

  if (length(tgi_df)==0) return(NULL) else tgi_df

  #   auc_mouse=base::do.call(base::rbind,base::lapply(mouses,get_eGR))
  # grps=base::unique(tv.data$Group)
  # trt_grps=base::sort(grps[!(grps %in% ref.group)])
  # ref_grp = base::setdiff(grps,trt_grps)
  # eGR.res = base::lapply(trt_grps,function(g) geteGRRatio(auc_mouse,g,ref_grp,type = type ,ci=ci))
  # names(eGR.res) = trt_grps
  # eGR.res.df =dplyr::bind_rows(eGR.res, .id = "Group")
  # if (length(eGR.res.df)==0) return(NULL) else eGR.res.df

}
