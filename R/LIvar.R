LIvar = function(x, y, k, conf.level=0.95)
{
  x = x[!is.na(x)] ; n1 = length(x) ; sxx = sum(x^2) ; sx = sum(x)
  y = y[!is.na(y)] ; n2 = length(y) ; syy = sum(y^2) ; sy = sum(y)
  if (!is.numeric(x) | sum(is.infinite(x) > 0) | sum(is.nan(x)) > 0 | n1 < 3 | length(unique(x)) == 1) stop("Check the input!")
  if (!is.numeric(y) | sum(is.infinite(y) > 0) | sum(is.nan(y)) > 0 | n2 < 3 | length(unique(y)) == 1) stop("Check the input!")
  
  m1 = sx/n1 ; v1 = sxx/n1 - m1^2 ; s1 = sqrt(v1)
  m2 = sy/n2 ; v2 = syy/n2 - m2^2 ; s2 = sqrt(v2)
  R0 = v1/v2 ; n = n1 + n2
  maxLL = -(n1*(log(2*pi*v1) + 1) + n2*(log(2*pi*v2) + 1))/2 

  if (!missing(k)) {
    logk = log(k)
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

  return(c(v1=v1, v2=v2, PE=R0, LL=LL, UL=UL, logk=logk, maxLL=maxLL))
}
