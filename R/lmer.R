function (formula, data = NULL, REML = TRUE, control = lmerControl(),
          start = NULL, verbose = 0L, subset, weights, na.action, offset,
          contrasts = NULL, devFunOnly = FALSE)
{
  mc <- mcout <- match.call()
  missCtrl <- missing(control)
  if (!missCtrl && !inherits(control, "lmerControl")) {
    if (!is.list(control))
      stop("'control' is not a list; use lmerControl()")
    warning("passing control as list is deprecated: please use lmerControl() instead",
            immediate. = TRUE)
    control <- do.call(lmerControl, control)
  }
  mc$control <- control
  mc[[1]] <- quote(lme4::lFormula)
  lmod <- eval(mc, parent.frame(1L))
  mcout$formula <- lmod$formula
  lmod$formula <- NULL
  devfun <- do.call(mkLmerDevfun, c(lmod, list(start = start,
                                               verbose = verbose, control = control)))
  if (devFunOnly)
    return(devfun)
  if (identical(control$optimizer, "none"))
    stop("deprecated use of optimizer=='none'; use NULL instead")
  opt <- if (length(control$optimizer) == 0) {
    s <- getStart(start, environment(devfun)$pp)
    list(par = s, fval = devfun(s), conv = 1000, message = "no optimization")
  }
  else {
    optimizeLmer(devfun, optimizer = control$optimizer, restart_edge = control$restart_edge,
                 boundary.tol = control$boundary.tol, control = control$optCtrl,
                 verbose = verbose, start = start, calc.derivs = control$calc.derivs,
                 use.last.params = control$use.last.params)
  }
  cc <- checkConv(attr(opt, "derivs"), opt$par, ctrl = control$checkConv,
                  lbound = environment(devfun)$lower)
  mkMerMod(environment(devfun), opt, lmod$reTrms, fr = lmod$fr,
           mc = mcout, lme4conv = cc)
}
