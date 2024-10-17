LIvtest = function(m1, s1, n1, m2, s2, n2, k, conf.level=0.95)
{
  v1 = (n1 - 1)*s1^2/n1 ; sx = n1*m1 ; sxx = n1*v1 + 2*m1*sx - n1*m1^2
  v2 = (n2 - 1)*s2^2/n2 ; sy = n2*m2 ; syy = n2*v2 + 2*m2*sy - n2*m2^2
  R0 = v1/v2 ; n = n1 + n2
  maxLL = -(n1*(log(2*pi*v1) + 1) + n2*(log(2*pi*v2) + 1))/2 

  if (!missing(k)) {
    logk = log(k)
    conf.level = pf((n - 1)*(k^(2/n) - 1), 1, n - 1)
  } else {
    logk = n/2*log(1 + qf(conf.level, 1, n - 1)/(n - 1))
    logk = min(logk, log(2/(1 - conf.level)))
  }

  Obj = function(r) {
    mLL = function(v2t) (n1*log(2*pi*r*v2t) + (sxx - sx^2/n1)/(r*v2t) + n2*log(2*pi*v2t) + (syy - sy^2/n2)/v2t)/2
    v2t = nlminb(v2, mLL, lower=1e-8, upper=100*v2)
    return(maxLL + v2t$objective - logk)
  }

  options(warn=-1)
  LL = uniroot(Obj, c(0, R0))$root
  UL = uniroot(Obj, c(R0, 1e4))$root
  options(warn=0)

  Res = c(g1=v1, g2=v2, PE=R0, LL=LL, UL=UL)
  Res = rbind(Res, sqrt(Res))
  rownames(Res) = c("var", "sd")
  attr(Res, "logk") = logk
  attr(Res, "maxLL") = maxLL
  attr(Res, "conf.level") = conf.level
  return(Res)
}


LIvtest(5.4, 10.5, 3529, 5.1, 8.9, 5190)
LIvtest(5.4, 10.5, 3529, 5.1, 8.9, 5190, exp(1.92106))
LIvtest(65, 2, 10, 62, 3, 10)
LIvtest(65, 2, 10, 62, 3, 10, exp(2.07474))

LIvtest(65, 2, 10, 62, 3, 10, k=15)
LIvtest(65, 2, 10, 62, 3, 10, conf.level=0.974868)

