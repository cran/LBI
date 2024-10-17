LBCIvar = function(x, conf.level=0.95)
{
  x = x[!is.na(x)]
  n0 = length(x) ; sxx = sum(x^2) ; sx = sum(x)
  if (!is.numeric(x) | sum(is.infinite(x) > 0) | sum(is.nan(x)) > 0 | n0 < 3 | length(unique(x)) == 1) stop("Check the input!")
  m0 = sx/n0
  v0 = sxx/n0 - m0^2
  s0 = sqrt(v0)
  maxLL = -n0*(log(2*pi*v0) + 1)/2

  logk0 = n0/2*log(1 + qf(conf.level, 1, n0 - 2)/(n0 - 2)) # two parameters with one nuisance (mean)
  logk0 = min(logk0, log(2/(1 - conf.level)))

  O1 = function(th, logk) maxLL + (n0*log(2*pi*th) + (sxx - sx^2/n0)/th)/2 - logk

  O2 = function(logk) {
    varLL = uniroot(O1, interval=c(1e-8, v0), logk=logk)$root
    varUL = uniroot(O1, interval=c(v0, 100*v0), logk=logk)$root
    ltProb = pchisq((n0 - 1)*varLL/v0, n0 - 1)
    rtProb = pchisq((n0 - 1)*varUL/v0, n0 - 1)
    rtProb - ltProb - conf.level
  }

  logk = uniroot(O2, interval=c(logk0, 2*logk0))$root
  varLL = uniroot(O1, interval=c(1e-8, v0), logk=logk)$root   # One more after fixing logk
  varUL = uniroot(O1, interval=c(v0, 100*v0), logk=logk)$root # One more after fixing logk
  sdLL = sqrt(varLL)
  sdUL = sqrt(varUL)
  Res = cbind(PE = c(s0, v0), LL=c(sdLL, varLL), UL=c(sdUL, varUL))
  rownames(Res) = c("sd", "var")
  attr(Res, "n") = n0
  attr(Res, "k") = exp(logk)
  attr(Res, "log(k)") = logk
  attr(Res, "maxLogLik") = maxLL
  return(Res)
}
