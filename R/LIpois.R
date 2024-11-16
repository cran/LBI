LIpois = function(x, k, n=1, conf.level=0.95, eps=1e-8)
{
  if (!is.finite(x)) stop("Check the input!")
  n0 = length(x)
  if (n0 == 1) { # x is mean x
    sx = x*n
  } else {       # x is raw data vector
    sx = sum(x)
    n  = n0
  }

  Lam0 = sx/n    # PE for mean
  maxLL = sx*log(Lam0) - n*Lam0 # excluding factorial part

  if (!missing(k)) {
    logk = log(k)
  } else if (n == 1) {
    logk = log(2/(1 - conf.level))
  } else {
    logk = n/2*log(1 + qf(conf.level, 1, n - 1)/(n - 1))
    logk = min(logk, log(2/(1 - conf.level))) # Pawitan p240 k = 20 -> p < 0.05
  }

  Obj = function(lam, ylevel) {
    ll = ifelse(lam > 0, sx*log(lam) - n*lam, 0)
    maxLL - ll - ylevel
  }
  if (Lam0 > 0) {
    LL = uniroot(Obj, c(eps, Lam0), ylevel=logk)$root
  } else {
    LL = 0
  }
  if (!is.finite(maxLL)) {
    UL = Inf
  } else {
    UL = uniroot(Obj, c(Lam0, 1e9), ylevel=logk)$root
  }
  Res = c(PE=Lam0, LL=LL, UL=UL)
  attr(Res, "n") = n
  attr(Res, "k") = exp(logk)
  attr(Res, "logk") = logk
  attr(Res, "maxLL without factorial") = maxLL
  return(Res)
}
