LIvar = function(s1, n1, k, conf.level=0.95)
{
  n1v1 = (n1 - 1)*s1^2
  v1 = n1v1/n1
  maxLL = -n1*(log(2*pi*v1) + 1)/2

  if (!missing(k)) {
    logk = log(k)
  } else {
    logk = n1/2*log(1 + qf(conf.level, 1, n1 - 1.2)/(n1 - 1.2)) # 1.2 is empirical
    logk = min(logk, log(2/(1 - conf.level)))
  }

  O2 = function(th) maxLL + (n1*log(2*pi*th) + n1v1/th)/2 - logk

  options(warn=-1)
  varLL = uniroot(O2, c(1e-8, v1))$root
  varUL = uniroot(O2, c(v1, 1e6*v1))$root
  options(warn=0)

  Res = c(PE=v1, LL=varLL, UL=varUL)
  Res = rbind(Res, sqrt(Res))
  rownames(Res) = c("var", "sd")
  attr(Res, "n") = n1
  attr(Res, "k") = exp(logk)
  attr(Res, "log(k)") = logk
  attr(Res, "maxLogLik") = maxLL
  return(Res)
}
