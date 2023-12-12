LInorm = function(x, k, conf.level=0.95) {
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
    logk = n0/2*log(1 + qf(conf.level, 1, n0 - 2)/(n0 - 2))
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
  attr(Res, "2*log(k)") = 2*logk
  return(Res)
}
