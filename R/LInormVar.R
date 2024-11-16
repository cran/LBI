LInormVar = function(x, k, conf.level=0.95)
{
  x = x[!is.na(x)]
  n0 = length(x)
  if (!is.numeric(x) | sum(is.infinite(x) > 0) | sum(is.nan(x)) > 0 | n0 < 3 | length(unique(x)) == 1) stop("Check the input!")
  n0v0 = sum((x - mean(x))^2)
  v0 = n0v0/n0
  maxLL = -n0*(log(2*pi*v0) + 1)/2

  if (!missing(k)) {
    logk = log(k)
  } else {
    logk = n0/2*log(1 + qf(conf.level, 1, n0 - 2)/(n0 - 2)) # two parameters with one nuisance (mean)
    logk = min(logk, log(2/(1 - conf.level)))
  }

  O2 = function(th) maxLL + (n0*log(2*pi*th) + n0v0/th)/2 - logk

  options(warn=-1)
  varLL = uniroot(O2, c(1e-8, v0))$root
  varUL = uniroot(O2, c(v0, 1e6*v0))$root
  options(warn=0)

  Res = c(PE=v0, LL=varLL, UL=varUL)
  Res = rbind(Res, sqrt(Res))
  rownames(Res) = c("var", "sd")
  attr(Res, "n") = n0
  attr(Res, "k") = exp(logk)
  attr(Res, "log(k)") = logk
  attr(Res, "maxLogLik") = maxLL
  return(Res)
}
