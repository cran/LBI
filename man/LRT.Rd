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
  If the residuals are normally distributed, the delta -2 log likelihood (the difference between -2LL, the objective function value of each model) follows exactly an F-distribution, independent of sample size.
  When the distribution of the residuals is not normal (no matter what the distribution of the residuals is), it approaches a chi-square distribution as sample size increases (Wilks' theorem).
  The extreme distribution of the F-distribution (when the degrees of freedom in the denominator go to infinity) is chi-square distribution.
  The p-value from the F-distribution is slightly larger than the p-value from the chi-square distribution, meaning the F-distribution is more conservative. 
  The difference decreases as sample size increases.
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

  r1 = lm(mpg ~ disp + drat + wt, mtcars)
  r2 = lm(mpg ~ disp + drat, mtcars)
  anova(r2, r1)
  LRT(nrow(mtcars), r1$rank, r2$rank, logLik(r1), logLik(r2))  
}
\references{
  \enumerate{
    \item Ruppert D, Cressie N, Carroll RJ. A Transformation/Weighting Model For Estimating Michaelis-Menten Parameters. School of Operations Research and Industrial Engineering, College of Engineering, Cornell University. Technical Report No. 796. May 1988.
    \item Scheff\eacute{} H. The Analysis of Variance. Wiley. 1959.
    \item Wilks SS. The Large-Sample Distribution of the Likelihood Ratio for Testing Composite Hypotheses. \emph{Annals Math. Statist.} 1938;9:60-62
  }
}
\author{Kyun-Seop Bae k@acr.kr}
