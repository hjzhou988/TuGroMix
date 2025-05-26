#' calculate eGR for a single tumor
#'
#' @param df a single tumor volume data frame in long format, which contains columns  "Group","Mouse","Day", "TV"
#' @return a data frame of eGR information
#' @export
#'
#'
get_eGR=function(df){
  # df = test_mouse_sample
  # df = mouses[[1]]
  df=base::as.data.frame(df)
  df=df[base::order(df$Day),]
  n=base::nrow(df)
  day_v=df$Day
  df$TV = base::as.numeric(df$TV)
  TV_v=base::log(df$TV+1)
  # TV_v = df$TV
  n1=base::length(day_v) # why n and n1? what's the difference?
  AUC=0
  if(n1>=4){
    for(i in 2:n1){
      AUC = AUC+(day_v[i]-day_v[i-1])*(TV_v[i]+TV_v[i-1])/2
    }
    day_span=day_v[n1]-day_v[1]
    AUC = AUC-TV_v[1]*day_span
    eGR = 2*AUC/(day_span*day_span) #tumor growth rate under exponential growth model.
  }else{
    eGR=NA
  }
  auc_df=base::data.frame('Mouse'=base::unique(df$Mouse),'eGR'=eGR,'Day'=day_v[n1]) #,'Group'=base::unique(base::as.character(df$Group))
  auc_df
}
