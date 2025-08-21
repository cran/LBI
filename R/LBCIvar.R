LBCIvar = function(x, conf.level=0.95)
{
  x = x[!is.na(x)] ; n0 = length(x)
  if (!is.numeric(x) | sum(is.infinite(x) > 0) | sum(is.nan(x)) > 0 | n0 < 3 | length(unique(x)) == 1) stop("Check the input!")
  n0v0 = sum(x^2) - sum(x)^2/n0
  v0   = n0v0/n0
  maxLL = -n0*(log(2*pi*v0) + 1)/2

  logk0 = n0/2*log(1 + qf(conf.level, 1, n0 - 2)/(n0 - 2)) # initial guess for logk
  logk0 = min(logk0, log(2/(1 - conf.level)))

  O1 = function(th, logk) maxLL + (n0*log(2*pi*th) + n0v0/th)/2 - logk
  O2 = function(logk) {
    varLL = uniroot(O1, interval=c(1e-8, v0), logk=logk)$root
    varUL = uniroot(O1, interval=c(v0, 1e6*v0), logk=logk)$root
    pchisq(n0v0/varLL, n0 - 1) - pchisq(n0v0/varUL, n0 - 1) - conf.level
  }

  options(warn=-1)
  logk = uniroot(O2, interval=c(0.5*logk0, 2*logk0))$root
  varLL = uniroot(O1, interval=c(1e-8, v0), logk=logk)$root   # One more after fixing logk
  varUL = uniroot(O1, interval=c(v0, 1e6*v0), logk=logk)$root # One more after fixing logk
  options(warn=0)

  Res = cbind(PE=v0, LL=varLL, UL=varUL)
  Res = rbind(Res, sqrt(Res))
  rownames(Res) = c("var", "sd")
  attr(Res, "n") = n0
  attr(Res, "k") = exp(logk)
  attr(Res, "log(k)") = logk
  attr(Res, "maxLogLik") = maxLL
  return(Res)
}
