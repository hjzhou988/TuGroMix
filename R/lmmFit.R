#' Fit by linear mixed-effects model (LMM)
#'
#' This function uses linear mixed model to test if there are significant differences in growth rates between treatment groups.
#' Need to set reference group. Results will be
#' in the simplest case, random slope for each mouse.
#'
#'
#'
#' @param tv.data Tumor volume data in "long format", that is, tumor volumes are in separate lines, with a Day column. Required.
#' @param ref.group The reference group, usually a control/vehicle group for the drug group to compare with. Required.
#' @param group.id.var The name of the "group" column of your input data. Required.
#' @param mouse.id.var The name of the "mouse.id" column of your input data. Required.
#' @details
#' Assuming tumor grows exponentially, log-transformed tumor growth curve will become a straight line versus time. LMM
#' is comparing the slopes of the lines, and will generate a p-value indicating whether there are significant differences in the slopes.
#'
#' @export
#' @return An S3 object of class "lmmfit" which stores coefficient tables and an S4 object of class "lmerModLmerTest".


lmmFit <- function(tv.data, ref.group, group.id.var, mouse.id.var){
  tv.data = base::as.data.frame(tv.data)
  tv.data$Group = tv.data[[group.id.var]]
  tv.data[["Group"]] = base::as.factor(tv.data[["Group"]])
  tv.data[["Group"]]  = stats::relevel( tv.data[["Group"]] , ref = ref.group)
  f.string = paste0("log(TV+1)~ Day + Day:Group + (Day-1|",mouse.id.var,")")
  lmm1<-lmerTest::lmer(stats::as.formula(f.string), data=tv.data)
  lmm1.sum = base::summary(lmm1)
  lmm_coefs = stats::coef(lmm1.sum)[-1,]
  lmm_coefs = base::as.data.frame(lmm_coefs)
  lmm_coefs$Term = base::gsub(':Group',':',rownames(lmm_coefs))
  # browser()
  lmm_coefs=lmm_coefs[,c(6,1:5)]
  colnames(lmm_coefs)[6]="p.value"
  lmm_coefs$Holm.adjusted.p.value = NA
  lmm_coefs$Holm.adjusted.p.value[-1] = stats::p.adjust(lmm_coefs$p.value[-1],method  = "holm")
  lmm_coefs$Holm.adjusted.p.value[1] = lmm_coefs$p.value[1]
  lmm_coefs$Hommel.adjusted.p.value = NA
  lmm_coefs$Hommel.adjusted.p.value[-1] = stats::p.adjust(lmm_coefs$p.value[-1],method  = "hommel")
  lmm_coefs$Hommel.adjusted.p.value[1] = lmm_coefs$p.value[1]
  # print(lmm_coefs)
  # list(coefs = lmm_coefs, mer = lmm1)
  obj = base::list(coefs = lmm_coefs, mer = lmm1)
  base::class(obj) <- "lmmfit"
  obj
}

# lmmFit()
