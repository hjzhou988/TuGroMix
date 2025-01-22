#' Calculate tumor growth inhibition for treatment groups
#'
#' @param tv.data a data frame containing tumor volume data in long format
#' @param ref.group the reference group
#' @param grp the group to calculate TGI for
#' @param type the type of the mean of TGI, choose from "arithmetic" and "geometric".
#' @param def definition of TGI
#' @param ci whether to do 1000 times bootstrapping to determine the 95 confidence interval for the calculated TGI.
#' @param bound a number for the confidence interval boundary. Suggest 90 or 95.
# #' @param group.id.var The name of the "group" column of your input data. Required.
# #' @param mouse.id.var The name of the "mouse.id" column of your input data. Required.
# #' @param time.var The name of the "day" column of your input data. Required.
# #' @param tv.var  The name of the "tumor volume" column of your input data. Required.
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

tgi_group = function(tv.data,ref.group,grp,type = c("geometric","arithmetic"),def=6, ci = F, bound = 90){ #,group.id.var,mouse.id.var,time.var,tv.var

  type = match.arg(type)
  # tv.data = tv.data; ref.group = ref.group; grp = "Dosing";def = 3
  # tv.data = dplyr::rename(tv.data, Group = {{group.id.var}}, Mouse = {{mouse.id.var}}, Day = {{time.var}}, TV = {{tv.var}})
  #
  # tv.data = tv.data[,c("Group","Mouse","Day","TV")]
  # tv.data = tv.data[stats::complete.cases(tv.data),]
  # tv.data = tv.data %>% arrange(Group,Mouse,Day)
  # tv.data = tv.data %>% group_by(Group,Mouse) %>% mutate(TV0 = first(TV))
  # tv.data$TV0 = ifelse(tv.data$TV0==0, 1, tv.data$TV0)
  # tv.data$RTV = tv.data$TV/tv.data$TV0
  # tv.data$DeltaTV = tv.data$TV - tv.data$TV0
  # tv.data = tv.data; ref.group = ref.group; grp = "Grp-2";def = 6

  tv_var='TV'
  if(def %in% c(2,3)){
    tv_var='DeltaTV'
  }else if(def %in% c(4,6)){ # added 6
    tv_var='RTV'
  }

  # load("data/tv_data.RData")
  # tv_wide_long = wide2long(tv_wide,id.vars = 1:2)
  # tv_wide_long = tv_wide_long %>% group_by(Group, Animal.No.) %>% mutate(TV0 = first(TV))
  # tv_wide_long = tv_wide_long  %>% mutate(RTV = TV/TV0,DeltaTV = TV-TV0)
  # tv.data = tv_wide_long

  # ds = tv.data %>% dplyr::filter(Group %in% c(ref.group,grp))
  ds = subset(tv.data,Group==ref.group | Group == grp)
  day_n=base::table(ds$Group,ds$Day)
  day_pct=day_n/day_n[,1]
  xx=base::apply(day_pct,2,function(v) all(v>=0.6))
  sel_day=base::max(base::as.numeric(base::names(xx[xx])),na.rm = T) # add na.rm = T, just in case it returns NA for sel_day
  ds=dplyr::filter(ds,Day==sel_day)

  if (type=="arithmetic"){
    ms=tapply(ds[[tv_var]],ds$Group,mean,na.rm = T)
  }else if(type=="geometric"){
    ms=tapply(ds[[tv_var]],ds$Group, function(x) exp(mean(log(x+1),na.rm =T)))
  }

  ms1=ms[-which(names(ms)==ref.group)]/ms[ref.group]

  if(def %in% c(3,5,6)) ms1 = 1-ms1 # added 6
  if(ci == F) list('Group'=grp,'TGI'=as.numeric(ms1),'Day'=sel_day) else {
    # tic()
    res = base::lapply(1:1000,function(i){
      # ds.s = ds %>% dplyr::group_by(Group) %>% dplyr::slice_sample(prop = 1, replace = T)
      sampled = tapply(ds[[tv_var]], ds$Group, sample,replace = TRUE)

      # ms=base::tapply(ds.s[,tv_var,drop = T],ds.s$Group,mean)

      if (type=="arithmetic"){
        ms = lapply(sampled, mean, na.rm = T)
      }else if(type=="geometric"){
        ms = lapply(sampled, function(x) exp(mean(log(x+1),na.rm =T)))
      }

      # ms.b=ms[-which(names(ms)==ref.group)]/ms[ref.group]
      ms.b=ms[[-which(names(ms)==ref.group)]]/ms[[ref.group]]
      if(def %in% c(3,5,6)) ms.b = 1-ms.b # added 6
      ms.b
    })
    # toc()
    qts = stats::quantile(base::unlist(res),c((100-bound)/200, (bound + (100-bound)/2)/100))
    out = base::list('Group'=grp,'TGI'=as.numeric(ms1),'3' = qts[[1]],'4'= qts[[2]],'Day'=sel_day)
    names(out)[3]=paste0(bound,"_CI_lower_bound")
    names(out)[4]=paste0(bound,"_CI_upper_bound")
    out
  }
}
