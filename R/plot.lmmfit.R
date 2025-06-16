#' Plot LMM fitted tumor growth curves
#'
#' This function visualized the LMM fitted tumor growth curves, by using the ggpredict() function from ggeffects package.
#'
#' @param object an object returned by lmmFit()
#' @param show.points a logical value to decide whether to show the original data points or not.
#'
#' @return a ggplot2 object of plot
#' @export

plot.lmmfit<- function(object,show.points = T){ # x is an object of from lmmFit
  if (!inherits(object, "lmmfit"))
    stop("plot.lmmfit can only be used to plot lmmfit objects")
  dat = object$mer@frame
  mydf = ggeffects::ggpredict(object$mer,terms = c("Day","Group"),back_transform = F)
  base::names(mydf)[ncol(mydf)]="Group"
  # dat$predicted = predict(res.lmmfit$mer,re.form = NA)
  # ggplot(dat)+geom_line(aes(x=Day,y=predicted,color = Group))+geom_point(aes(x=x, y = predicted,color = group),data = mydf)
  if (show.points==T){
    p=ggplot2::ggplot(dat)+
      ggplot2::geom_point(ggplot2::aes(x=Day, y= `log(TV + 1)`, group = Group, color = Group),alpha = 0.3)+
      ggplot2::geom_line(ggplot2::aes(x=x,y=predicted,color = Group),data = mydf) +
      ggplot2::geom_ribbon(ggplot2::aes(x=x,y = predicted,ymin = conf.low, ymax = conf.high,fill = Group),data=mydf, alpha = 0.2)
  }else{
    p=ggplot2::ggplot(dat)+
      ggplot2::geom_line(ggplot2::aes(x=x,y=predicted,color = Group),data = mydf) +
      ggplot2::geom_ribbon(ggplot2::aes(x=x,y = predicted,ymin = conf.low, ymax = conf.high,fill = Group),data=mydf, alpha = 0.2)+
      ggplot2::xlab("Day")+ggplot2::ylab("log(TV + 1)")
  }
  base::print(p)
}
# pdf("test.pdf")
# plot.lmmfit(res)
# dev.off()
