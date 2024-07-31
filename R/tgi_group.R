#' Calculate tumor growth inhibition for treatment groups
#'
#' @param tv.data a data frame containing tumor volume data in long format
#' @param ref.group the reference group
#' @param def definition of TGI
#' @param group.id.var The name of the "group" column of your input data. Required.
#' @param mouse.id.var The name of the "mouse.id" column of your input data. Required.
#' @param time.var The name of the "day" column of your input data. Required.
#' @param tv.var  The name of the "tumor volume" column of your input data. Required.
#' @details
#' Definitions of TGI:
#' 1. TV_t/TV_c
#' 2. DeltaTV_t/DeltaTV_c
#' 3. 1 - DeltaTV_t/DeltaTV_c
#' 4. RTV_t/RTV_c
#' 5. 1 - TV_t/TV_c
#' 6. 1 - RTV_t/RTV_c
#'
#'
#'
#' @return a list of TGI information
#'
#'

tgi_group = function(tv.data,ref.group,def=3,group.id.var,mouse.id.var,time.var,tv.var){

  # tv.data = tv.data; ref.group = ref.group; grp = "Dosing";def = 3
  tv.data = dplyr::rename(tv.data, Group = {{group.id.var}}, Mouse = {{mouse.id.var}}, Day = {{time.var}}, TV = {{tv.var}})

  tv.data = tv.data[,c("Group","Mouse","Day","TV")]
  tv.data = tv.data[stats::complete.cases(tv.data),]
  tv.data = tv.data %>% arrange(Group,Mouse,Day)
  tv.data = tv.data %>% group_by(Group,Mouse) %>% mutate(TV0 = first(TV))
  tv.data$TV0 = ifelse(tv.data$TV0==0, 1, tv.data$TV0)
  tv.data$RTV = tv.data$TV/tv.data$TV0
  tv.data$DeltaTV = tv.data$TV - tv.data$TV0


  tv_var='TV'
  if(def %in% c(2,3)){
    tv_var='DeltaTV'
  }else if(def %in% c(4,6)){ # added 6
    tv_var='RTV'
  }

  ds = tv.data # %>% dplyr::filter(Group %in% c(ref.group,grp))
  day_n=base::table(ds$Group,ds$Day)
  day_pct=day_n/day_n[,1]
  xx=base::apply(day_pct,2,function(v) all(v>=0.6))
  sel_day=base::max(base::as.numeric(base::names(xx[xx])))
  ds=dplyr::filter(ds,Day==sel_day)
  ms=base::tapply(ds[,tv_var,drop = T],ds$Group,mean)
  ms1=ms[-which(names(ms)==ref.group)]/ms[ref.group]
  if(def %in% c(3,5,6)) ms1 = 1-ms1 # added 6
  list('Group'=grp,'TGI'=as.numeric(ms1),'Day'=sel_day)
}