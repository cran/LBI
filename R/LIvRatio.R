LIvRatio = function(x, y, k, conf.level=0.95)
{
  x = x[!is.na(x)] ; n1 = length(x) ; n1v1 = sum((x - mean(x))^2) ; v1 = n1v1/n1
  y = y[!is.na(y)] ; n2 = length(y) ; n2v2 = sum((y - mean(y))^2) ; v2 = n2v2/n2
  if (!is.numeric(x) | sum(is.infinite(x) > 0) | sum(is.nan(x)) > 0 | n1 < 3 | length(unique(x)) == 1) stop("Check the input!")
  if (!is.numeric(y) | sum(is.infinite(y) > 0) | sum(is.nan(y)) > 0 | n2 < 3 | length(unique(y)) == 1) stop("Check the input!")

  R0 = v1/v2
  n = n1 + n2
  maxLL = -(n1*(log(2*pi*v1) + 1) + n2*(log(2*pi*v2) + 1))/2

  if (!missing(k)) {
    logk = log(k)
  } else {
    logk = n/2*log(1 + qf(conf.level, 1, n - 2)/(n - 2)) # two parameters
    logk = min(logk, log(2/(1 - conf.level)))
  }

  Obj = function(r) {
    mLL = function(v2t) (n1*log(2*pi*r*v2t) + n1v1/(r*v2t) + n2*log(2*pi*v2t) + n2v2/v2t)/2
    v2t = nlminb(v2, mLL, lower=1e-8, upper=1e6*v2)
    maxLL + v2t$objective - logk
  }

  options(warn=-1)
  LL = uniroot(Obj, c(0, R0))$root
  UL = uniroot(Obj, c(R0, 1e4))$root
  options(warn=0)

  Res = c(group1=v1, group2=v2, PE=R0, LL=LL, UL=UL)
  Res = rbind(Res, sqrt(Res))
  rownames(Res) = c("var", "sd")
  attr(Res, "logk") = logk
  attr(Res, "maxLL") = maxLL
  attr(Res, "conf.level") = conf.level
  return(Res)
}
