#' Fit tumor growth curves using generalized additive model (GAM)
#'
#' This function uses generalized additive model to test if there are significant differences in growth rates between treatment groups.
#'
#' @param tv.data Tumor volume data in "long format", that is, tumor volumes are in separate lines, with a Day column. Required.
#' @param ref.group The reference group, usually a control/vehicle group for the drug group to compare with. Required.
#' @param group.id.var The name of the "group" column of your input data. Required.
#' @param mouse.id.var The name of the "mouse.id" column of your input data. Required.
#'
#' @details
#' Generalized additive model fits tumor growth curves with smooth functions (with regularization to avoid over-fitting).
#' p-value is given to test whether there is significant difference between the intercepts of groups (which can be a surrogate for tumor growth rates)
#'
#' @return An S3 object of class "gamfit" which stores coefficient tables and a "gam" object of the fitted GAM model.
#'
#' @export
#'


gamFit <- function(tv.data, ref.group,group.id.var,mouse.id.var){
  # tv.data = tv_long;ref.group="Vehicle";group.id.var="Trt";mouse.id.var = "Mouse.ID.2"
  tv.data$Group = tv.data[[group.id.var]]
  tv.data[["Group"]] = base::as.factor(tv.data[["Group"]])
  tv.data[["Group"]]  = stats::relevel( tv.data[["Group"]] , ref = ref.group)
  f.string.2 = base::paste0("log(TV+1) ~ Group + s(Day,by=Group)")
  f.string.r = base::paste0("~(Day-1|",mouse.id.var,")")
  fit = NULL
  base::tryCatch(fit <- gamm4::gamm4(stats::as.formula(f.string.2),random= as.formula(f.string.r),data = tv.data,REML = T),
           error = function(e){print(paste("Error:",e))})
  gam.sum = base::summary(fit[["gam"]])
  Coef= gam.sum$p.table[-1,,drop=F]
  Coef = base::as.data.frame(Coef)
  colnames(Coef)[4]="p.value"
  Coef$Holm.adjusted.p.value = stats::p.adjust(Coef$p.value,method = "holm")
  Coef$Hommel.adjusted.p.value = stats::p.adjust(Coef$p.value,method = "hommel")
  obj = base::list(coefs = Coef, gam = fit$gam, mer = fit$mer)
  base::class(obj)<-"gamfit"
  obj
}
