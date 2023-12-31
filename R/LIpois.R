LIpois = function(x, k, conf.level=0.95, eps=1e-8)
{
  Lam0 = mean(x)
  maxLL = sum(dpois(x, Lam0, log=T))

  n = length(x)
  if (!missing(k)) {
    logk = log(k)
  } else if (n == 1) {
    logk = log(2/(1 - conf.level))    
  } else {
    logk = n/2*log(1 + qf(conf.level, 1, n - 1)/(n - 1))
    logk = min(logk, log(2/(1 - conf.level))) # Pawitan p240 k = 20 -> p < 0.05
  }

  Obj = function(lam, ylevel) {
    ll = ifelse(lam > 0, sum(dpois(x, lam, log=T)), 0)
    return(maxLL - ll - ylevel)
  }
  if (Lam0 > 0) {
    LL = uniroot(Obj, c(eps, Lam0), ylevel=logk)$root
  } else {
    LL = 0
  }
  UL = uniroot(Obj, c(Lam0, 1e9), ylevel=logk)$root
  return(c(PE=Lam0, LL=LL, UL=UL, k=exp(logk)))
}
