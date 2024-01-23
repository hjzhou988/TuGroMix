#' get eGR ratio for treatment groups
#'
#' @param df a data frame containing calculated eGR
#' @param grp a specific treatment group
#' @param ref_grp the reference group
#' @param type calculate eGR ratio or difference
#' @param ci a logical value to bootstrap confidence interval or not
#' @param nrep the number of bootstrapping for confidence interval
#'
#' @return a list of eGR information
#'
#'
geteGRRatio=function(df,grp,ref_grp,type = c("difference","ratio"),ci=F,nrep=1000){

  # ref_grp=levels(df$Group)[1]
  # df = auc_mouse; grp = trt_grps; ref_grp= "Vehicle"
  # df = auc_mouse; grp = "Group 02"; ref_grp= "Group 01";ci=T;nrep=1000
  type = match.arg(type)

  auc_v=df[df$Group%in%ref_grp,'eGR']
  auc_t=df[df$Group==grp,'eGR']
  aucs=expand.grid('eGR_tr'=auc_t,'eGR_v'=auc_v)
  if(type=="ratio"){
    medianeGRRatio=base::with(aucs,stats::median(eGR_tr/eGR_v,na.rm =TRUE))
  }else if(type=="difference"){
    medianeGRRatio=base::with(aucs,stats::median(eGR_tr - eGR_v,na.rm =TRUE))
  }

  #medianeGRRatio
  cis=c(NA,NA)
  if(ci){
    if (length(auc_v)==1 | length(auc_t)==1)  cis = c(NA,NA) # at least two mice in each group
    else{
      auc_boot = base::sapply(1:nrep,function(i){
        auc_vi = base::sample(auc_v,replace = T) # if auc_v is a single value, and less than 1, then sampling will only take its own auc_v value.
        auc_ti = base::sample(auc_t,replace = T)
        aucsi = base::expand.grid('eGR_tr'=auc_ti,'eGR_v'=auc_vi)
        if(type=="ratio"){
          base::with(aucsi,stats::median(eGR_tr/eGR_v,na.rm =TRUE))
        }else if(type=="difference"){
          base::with(aucsi,stats::median(eGR_tr - eGR_v,na.rm =TRUE))
        }

      })
      cis=stats::quantile(auc_boot,c(0.05,0.95),na.rm = TRUE)
    }
  }
  base::names(cis)=NULL
  if(type=="ratio"){
    base::list('Median.eGR.Ratio'=medianeGRRatio,"Lower.Bound.90CI"=cis[1],
               "Upper.Bound.90CI"=cis[2],'Reference.Group'=ref_grp)
  }else if(type=="difference"){
    base::list('Median.eGR.Difference'=medianeGRRatio,"Lower.Bound.90CI"=cis[1],
               "Upper.Bound.90CI"=cis[2],'Reference.Group'=ref_grp)
  }

}
