LIvtest = function(s1, n1, s2, n2, k, conf.level=0.95)
{
  n1v1 = (n1 - 1)*s1^2 ; v1 = n1v1/n1
  n2v2 = (n2 - 1)*s2^2 ; v2 = n2v2/n2
  R0 = v1/v2
  n = n1 + n2
  maxLL = -(n1*(log(2*pi*v1) + 1) + n2*(log(2*pi*v2) + 1))/2 

  if (!missing(k)) {
    logk = log(k)
    conf.level = pf((n - 1)*(k^(2/n) - 1), 1, n - 1)
  } else {
    logk = n/2*log(1 + qf(conf.level, 1, n - 2.4)/(n - 2.4)) # empirical for two parameters
    logk = min(logk, log(2/(1 - conf.level)))
  }

  Obj = function(r) {
    mLL = function(v2t) (n1*log(2*pi*r*v2t) + n1v1/(r*v2t) + n2*log(2*pi*v2t) + n2v2/v2t)/2
    v2t = nlminb(v2, mLL, lower=1e-8, upper=1e6*v2)
    return(maxLL + v2t$objective - logk)
  }

  options(warn=-1)
  LL = uniroot(Obj, c(0, R0))$root
  UL = uniroot(Obj, c(R0, 1e9))$root
  options(warn=0)

  Res = c(group1=v1, group2=v2, PE=R0, LL=LL, UL=UL)
  Res = rbind(Res, sqrt(Res))
  rownames(Res) = c("var", "sd")
  attr(Res, "logk") = logk
  attr(Res, "maxLL") = maxLL
  attr(Res, "conf.level") = conf.level
  return(Res)
}
