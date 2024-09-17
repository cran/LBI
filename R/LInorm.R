LInorm = function(x, k, conf.level=0.95, PLOT="", LOCATE=FALSE, Resol=201)
{
  x = x[!is.na(x)]
  n0 = length(x)
  if (!is.numeric(x) | sum(is.infinite(x) > 0) | sum(is.nan(x)) > 0 | n0 < 3 | length(unique(x)) == 1) stop("Check the input!")
  m0 = mean(x)
  v0 = var(x)*(n0 - 1)/n0
  s0 = sqrt(v0)
  maxLL = sum(dnorm(x, mean=m0, sd=s0, log=TRUE))

  if (!missing(k)) {
    logk = log(k)
  } else {
#    logk = n0/2*log(1 + 2*qf(conf.level, 2, n0 - 2)/(n0 - 2))
    logk = n0/2*log(1 + 1*qf(conf.level, 1, n0 - 2)/(n0 - 2)) # two parameters with one nuisance
    logk = min(logk, log(2/(1 - conf.level)))
  }

  O1 = function(th) maxLL - sum(dnorm(x, mean=th, sd=s0, log=TRUE)) - logk
  O2 = function(th) maxLL - sum(dnorm(x, mean=m0, sd=th, log=TRUE)) - logk
  meanLL = uniroot(O1, c(m0 - 10*s0, m0))$root
  meanUL = uniroot(O1, c(m0, m0 + 10*s0))$root
  sdLL = uniroot(O2, c(1e-7, s0))$root
  sdUL = uniroot(O2, c(s0, 100*s0))$root
  varLL = sdLL^2
  varUL = sdUL^2
  Res = cbind(PE = c(m0, s0, v0), LL=c(meanLL, sdLL, varLL), UL=c(meanUL, sdUL, varUL))
  rownames(Res) = c("mean", "sd", "var")
  attr(Res, "n") = n0
  attr(Res, "k") = exp(logk)
  attr(Res, "log(k)") = logk
  attr(Res, "maxLogLik") = maxLL

  if (toupper(trimws(PLOT)) %in% c("1D", "PROFILE")) {
    x1 = seq(m0 - 2*(m0 - meanLL), m0 + 2*(meanUL - m0), length.out=Resol)
    y1 = rep(NA, Resol)
    for (i in 1:Resol) y1[i] = sum(dnorm(x, mean=x1[i], sd=s0, log=TRUE))

    x2 = seq(max(1e-4, s0 - 2*(s0 - sdLL)), max(0, s0 + 2*(sdUL - s0)), length.out=Resol)
    y2 = rep(NA, Resol)
    for (i in 1:Resol) y2[i] = sum(dnorm(x, mean=m0, sd=x2[i], log=TRUE))

    oPar = par(mfrow=c(2, 2))

    plot(x1, y1, type="l", xlab="Mean", ylab="log Likelihood", ylim=c(maxLL - 3*logk, maxLL))
    abline(h=maxLL - logk, col="red")

    plot(x2, y2, type="l", xlab="Standard Deviation", ylab="log Likelihood", ylim=c(maxLL - 3*logk, maxLL))
    abline(h=maxLL - logk, col="red")

    plot(x1, exp(y1), type="l", xlab="Mean", ylab="Likelihood")
    abline(h=exp(maxLL - logk), col="red")

    plot(x2, exp(y2), type="l", xlab="Standard Deviation", ylab="Likelihood")
    abline(h=exp(maxLL - logk), col="red")

    par(oPar)
  } else if (toupper(trimws(PLOT)) %in% c("2D", "CONTOUR")) {
    Xs = seq(m0 - 2*(m0 - meanLL), m0 + 2*(meanUL - m0), length.out=Resol)
    Ys = seq(max(1e-4, s0 - 2*(s0 - sdLL)), s0 + 2*(sdUL - s0), length.out=Resol)
    mLik = matrix(NA, nrow=Resol, ncol=Resol)
    for (i in 1:Resol) for (j in 1:Resol) mLik[i, j] = sum(dnorm(x, mean=Xs[i], sd=Ys[j], log=TRUE))
    mLik[mLik < maxLL - 3*logk] = maxLL - 3*logk
    contour(Xs, Ys, mLik, xlab="Mean", ylab="Standard Deviation")
    contour(Xs, Ys, mLik, nlevels=1, levels=maxLL - logk, col="red", add=TRUE)
    abline(h=s0, v=m0, lty=3)
    points(m0, s0, pch="+", col="red")

    if (LOCATE) {
      print(Res)
      options(locatorBell = FALSE) # locator bell is annoying
      while (TRUE) {
        ans = locator(n=1)
        if (length(ans) < 1) break
        logL = sum(dnorm(x, mean=ans$x, sd=ans$y, log=TRUE))
        print(c(Mean=ans$x, SD=ans$y, logLik=logL))
        flush.console()
      }
      invisible(Res)
    }
  }

  return(Res)
}
