#' Plot GAM fitted tumor growth curves
#'
#' This function visualized the GAM fitted tumor growth curves, by using the function is plot_smooths() in package tidymv
#' @param object an object returned by gamFit()
#' @param show.points a logical value to decide whether to show the original data points or not.
#'
#' @export

plot.gamfit<- function(object,show.points = T){ # x is an object of from gamFit
  # dat = res$gam[["model"]]
  if (!inherits(object, "gamfit"))
    stop("plot.gamfit can only be used to plot gamfit objects")

  dat = object$gam[["model"]]
  # x$gam %>% predict_gam() %>% plot(series = "Day")
  # res$gam %>% predict_gam() %>% plot(series = "Day",comparison = "Group")
  if(show.points==T){
    tidymv::plot_smooths(model = object$gam, series = Day,comparison=Group)+
      ggplot2::geom_point(ggplot2::aes(x=Day, y= `log(TV + 1)`, group = Group, color = Group), data = dat,alpha = 0.3)
  }else{
    tidymv::plot_smooths(model = object$gam, series = Day,comparison=Group)
  }
  # par(mfrow=c(2,1))
  # ggplot(dat, aes(x=Day, y= `log(TV + 1)`, group = Group))+geom_point(aes(color = Group))
  # ggplot(dat, aes(x=Day, y= TV.p, group = Group))+geom_line(aes(color = Group))+ylab("log(TV+1)")
  #
  # print(p1)
}
# pdf("test.pdf")
# plot.lmmfit(res)
# dev.off()
# require(mgcv)
# plot.gam(res$gam)
# require(tidymv)

