#' get eGR for a single group of tumors
#'
#' @param tv.data tumor volume data in long format, containing columns "Group","Mouse","Day","TV"
#' @param mouse.id.var The name of the "mouse.id" column of your input data. Required.
#' @param time.var The name of the "day" column of your input data. Required.
#' @param tv.var  The name of the "tumor volume" column of your input data. Required.
#' @param ci a logical variable to do bootstrapping for 90% confidence interval
#' @param nrep the number of bootstrapping for confidence interval
#' @return a data frame containing eGR for this group
#'
#' @export
#'
get_group_eGR=function(tv.data,mouse.id.var,time.var,tv.var,ci = F,nrep=1000){
  # tv.data = ind_study [[1]];  mouse.id.var = "Mouse";time.var = "Day";tv.var = "TV";ci=T;nrep = 1000
  tv.data = dplyr::rename(tv.data,  Mouse = {{mouse.id.var}}, Day = {{time.var}}, TV = {{tv.var}})
  tv.data = tv.data[,c("Mouse","Day","TV")]
  tv.data = tv.data[stats::complete.cases(tv.data),]
  mouses=base::split(tv.data,tv.data$Mouse)
  auc_mouse=base::do.call(base::rbind,base::lapply(mouses,get_eGR))

  med.eGR = median(auc_mouse[,"eGR"],na.rm = T)

  N = nrow(auc_mouse)

  cis=c(NA,NA)
  if(ci){
    if (sum(!is.na(auc_mouse[,"eGR"]))<=1)
      cis = c(NA,NA) # at least two mice in a group
    else{
      boot.res = sapply(1:nrep, function(i){
        auc_boot = base::sample(auc_mouse[, "eGR"], size = N,
                                replace = T)
        auc_boot_med = median(auc_boot,na.rm = T)
      })
      cis = stats::quantile(boot.res, c(0.05, 0.95),na.rm = T)
    }
  }
  eGR.res.df = data.frame(Median.eGR = med.eGR, "Lower.Bound.90CI"=cis[1],"Upper.Bound.90CI"=cis[2],row.names = NULL)
  if (sum(!is.na(auc_mouse[,"eGR"]))==0) return(NULL) else eGR.res.df
}
