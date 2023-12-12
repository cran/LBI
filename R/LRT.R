LRT = function(n, pFull, pReduced, logLikFull, logLikReduced, alpha=0.05, Wilks=FALSE)
{
  if (n - pFull < 1 | pFull < pReduced | alpha <= 0 | alpha >= 0.5) stop("Check the input!")

  dP = pFull - pReduced
  dLL = max(logLikFull - logLikReduced, 0)

  if (Wilks) {
    cut0 = qchisq(1 - alpha, dP)/2
    Chisq = 2*dLL
    pval = 1 - pchisq(Chisq, dP)
  } else {
    cut0 = n/2*log(1 + dP*qf(1 - alpha, dP, n - pFull)/(n - pFull))
    Fval = (exp(dLL)^(2/n) - 1)*(n - pFull)/dP
    pval = 1 - pf(Fval, dP, n - pFull)
  }

  if (dLL > cut0) { Verdict = "Full model is better."
  } else if (dLL < cut0) { Verdict = "Reduced model is better."
  } else { Verdict = "Both models are equivalent." }

  if (Wilks) {
    Res = list(n = n, paraFull = pFull, paraReduced = pReduced, deltaPara = dP,
               cutoff = cut0, deltaLogLik = dLL, Chisq=Chisq, pval=pval,
               Verdict=Verdict)
    if (n < 100)
    attr(Res, "Warning") = "Sample size (n) is too small to use Wilks' theorem!"
  } else {
    Res = list(n = n, paraFull = pFull, paraReduced = pReduced, deltaPara = dP,
               cutoff = cut0, deltaLogLik = dLL, Fval=Fval, pval=pval,
               Verdict=Verdict)
  }
  return(Res)
}
