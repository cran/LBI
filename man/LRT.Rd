\name{LRT}
\alias{LRT}
\title{Likelihood Ratio Test}
\description{
  Likelihood ratio test with given fitting results, sample size, number of parameters, log-likelihoods, and alpha
}
\usage{
  LRT(n, pFull, pReduced, logLikFull, logLikReduced, alpha=0.05, Wilks=FALSE)
}
\arguments{
  \item{n}{number of observations}
  \item{pFull}{number of parameters of full model}
  \item{pReduced}{number of parameters of reduced model}
  \item{logLikFull}{log likelihood of full model}
  \item{logLikReduced}{log likelihood of reduced model}
  \item{alpha}{alpha value for type I error, significance level}
  \item{Wilks}{if TRUE, Wilks theorem (chi-square distribution) will be used, otherwise F distribution will be used.}
}
\details{
  It performs likelihood ratio test with given fitting results. The default test is using F distribution. For small n (i.e. less than 100), you need to use F distribution.
}
\value{
  \item{n}{number of observations}
  \item{paraFull}{number of parameters of full model}
  \item{paraReduced}{number of parameters of reduced model}
  \item{deltaPara}{difference of parameter counts}
  \item{cutoff}{cutoff, threshold, critical value of log-likelihood for the test}
  \item{deltaLogLik}{difference of log likelihood, if negative 0 is used.}
  \item{Chisq or Fval}{statistics according to the used distribution Chi-square of F}
  \item{pval}{p-value of null hypothesis. i.e. the reduced model is better.}
  \item{Verdict}{the model preferred.}
}
\examples{
  LRT(20, 4, 2, -58.085, -60.087)
  LRT(20, 4, 2, -58.085, -60.087, Wilks=TRUE)
  LRT(20, 4, 2, -57.315, -66.159)
  LRT(20, 4, 2, -57.315, -66.159, Wilks=TRUE)
}
\author{Kyun-Seop Bae k@acr.kr}
