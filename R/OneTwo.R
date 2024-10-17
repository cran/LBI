OneTwo = function(x, alpha=0.05)
{
  x = x[!is.na(x)] ; n0 = length(x) ; sxx = sum(x^2) ; sx = sum(x)
  if (!is.numeric(x) | sum(is.infinite(x) > 0) | sum(is.nan(x)) > 0 | n0 < 3 | length(unique(x)) == 1) stop("Check the input!")
  m0 = sx/n0
  v0 = sxx/n0 - m0^2
  s0 = sqrt(v0)
  logL0 = -n0*(log(2*pi*v0) + 1)/2
  r1 = list(n=n0, mean=m0, var=v0, sd=s0, logLik=logL0)
  p1 = 1
  mLL = function(TH, d0) {
    el1 = dnorm(d0, mean=TH[1], sd=TH[2])
    el2 = dnorm(d0, mean=TH[3], sd=TH[4])

    O2 = function(p0) -sum(log(p0*el1 + (1 - p0)*el2))
    rt = optimize(O2, c(0, 1))
    p1 <<- rt$minimum
    return(rt$objective)
  }

  m0LL = min(x, na.rm=T)
  m0UL = max(x, na.rm=T)
  r2 = nlminb(c(0.8*m0, s0, 1.2*m0, s0), mLL, lower=c(m0LL, 1e-4, m0, 1e-4),
              upper=c(m0, s0, m0UL, s0), d0=x)

  el1 = p1*dnorm(x, mean=r2$par[1], sd=r2$par[2])
  el2 = (1 - p1)*dnorm(x, mean=r2$par[3], sd=r2$par[4])
  n1 = sum(el1 > el2)
  n2 = n0 - n1

  r3 = LRT(n0, pFull=4, pReduced=2, -r2$objective, logL0, alpha=alpha)
  r4 = LRT(n0, pFull=4, pReduced=2, -r2$objective, logL0, alpha=alpha, Wilks=T)

  r0 = NULL
  r0$Estimate = data.frame(Model=c("One group", "Two group", ""),
                           GroupNo = c(1, 1, 2),
                           n = c(n0, n1, n2),
                           Mean = c(m0, r2$par[1], r2$par[3]),
                           SD = c(s0, r2$par[2], r2$par[4]),
                           Prior = c(1, p1, 1 - p1))
  r0$Delta = data.frame(Model = c("One group", "Two group", "Difference"),
                        nPara = c(2, 4, 2),
                        logLik = c(logL0, -r2$objective, max(-r2$objective - logL0, 0)))
  r0$Statistic = data.frame(Distribution = c("F", "Chisq"),
                            CutoffLogLik = c(r3["cutoff"], r4["cutoff"]), # cut-off value in log-likelihood
                            StatCrit = c(qf(1 - alpha, 2, n0 - 2), qchisq(1 - alpha, 2)), # statistics critical value
                            Statistic = c(r3["Fval"], r4["Chisq"]),
                            pval = c(r3["pval"], r4["pval"]),
                            Favor = c(ifelse(r3["pval"] < alpha, "Two group", "One group"), ifelse(r4["pval"] < alpha, "Two group", "One group")))
  attr(r0$Statistic, "alpha") = alpha
  return(r0)
}
