
function (rho, opt, reTrms, fr, mc, lme4conv = NULL)
{
  if (missing(mc))
    mc <- match.call()
  stopifnot(is.environment(rho), is(pp <- rho$pp, "merPredD"),
            is(resp <- rho$resp, "lmResp"), is.list(opt), "par" %in%
              names(opt), c("conv", "fval") %in% substr(names(opt),
                                                        1, 4), is.list(reTrms), c("flist", "cnms", "Gp",
                                                                                  "lower") %in% names(reTrms), length(rcl <- class(resp)) ==
              1)
  n <- nrow(pp$V)
  p <- ncol(pp$V)
  isGLMM <- (rcl == "glmResp")
  dims <- c(N = nrow(pp$X), n = n, p = p, nmp = n - p, q = nrow(pp$Zt),
            nth = length(pp$theta), nAGQ = rho$nAGQ, compDev = rho$compDev,
            useSc = !(isGLMM && hasNoScale(resp$family)), reTrms = length(reTrms$cnms),
            spFe = 0L, REML = if (rcl == "lmerResp") resp$REML else 0L,
            GLMM = isGLMM, NLMM = (rcl == "nlsResp"))
  storage.mode(dims) <- "integer"
  fac <- as.numeric(rcl != "nlsResp")
  if (trivial.y <- (length(resp$y) == 0)) {
    sqrLenU <- wrss <- pwrss <- NA
  }
  else {
    sqrLenU <- pp$sqrL(fac)
    wrss <- resp$wrss()
    pwrss <- wrss + sqrLenU
  }
  beta <- pp$beta(fac)
  if (!is.null(sc <- attr(pp$X, "scaled:scale"))) {
    warning("auto(un)scaling not yet finished/tested")
    beta2 <- beta
    beta2[names(sc)] <- sc * beta2[names(sc)]
    beta <- beta2
  }
  if (!is.null(attr(pp$X, "scaled:center"))) {
    warning("auto(un)centering not yet implemented")
  }
  sigmaML <- pwrss/n
  if (rcl != "lmerResp") {
    pars <- opt$par
    if (length(pars) > length(pp$theta))
      beta <- pars[-(seq_along(pp$theta))]
  }
  cmp <- c(ldL2 = pp$ldL2(), ldRX2 = pp$ldRX2(), wrss = wrss,
           ussq = sqrLenU, pwrss = pwrss, drsum = if (rcl == "glmResp" &&
                                                      !trivial.y) resp$resDev() else NA, REML = if (rcl ==
                                                                                                    "lmerResp" && resp$REML != 0L && !trivial.y) opt$fval else NA,
           dev = if (rcl == "lmerResp" && resp$REML != 0L || trivial.y) NA else opt$fval,
           sigmaML = sqrt(unname(if (!dims["useSc"] || trivial.y) NA else sigmaML)),
           sigmaREML = sqrt(unname(if (rcl != "lmerResp" || trivial.y) NA else sigmaML *
                                     (dims["n"]/dims["nmp"]))), tolPwrss = rho$tolPwrss)
  if (missing(fr))
    fr <- data.frame(resp$y)
  new(switch(rcl, lmerResp = "lmerMod", glmResp = "glmerMod",
             nlsResp = "nlmerMod"), call = mc, frame = fr, flist = reTrms$flist,
      cnms = reTrms$cnms, Gp = reTrms$Gp, theta = pp$theta,
      beta = beta, u = if (trivial.y)
        rep(NA_real_, nrow(pp$Zt))
      else pp$u(fac), lower = reTrms$lower, devcomp = list(cmp = cmp,
                                                           dims = dims), pp = pp, resp = resp, optinfo = .optinfo(opt,
                                                                                                                  lme4conv))
}
