\name{LIvRatio}
\alias{LIvRatio}
\title{Likelihood Interval of the ratio of two variances from two groups}
\description{
  Likelihood interval of  the ratio of two variances from two groups assuming normal distribution. Likelihood interval usually gives a narrower interval when the likelihood function is asymmetric.
}
\usage{
  LIvRatio(x, y, k, conf.level=0.95)
}
\arguments{
  \item{x}{observations from the first group, the test group, used for the numerator}
  \item{y}{observations from the second group, the control group, used for the denominator}
  \item{k}{1/k likelihood interval will be provided}  
  \item{conf.level}{approximate confidence level}
}
\details{
  It calculates likelihood interval of the ratio of two variances from two groups. Likelihood interval usually gives a narrower interval when the likelihood function is asymmetric.
}
\value{
  \item{PE}{maximum likelihood estimate on the ratio}
  \item{LL}{lower limit of likelihood interval}
  \item{UL}{upper limit of likelihood interval}
  \item{logk}{log(k) value used for LI}
  \item{maxLL}{maximum log likelihood}
  \item{conf.level}{approximate confidence level}  
}
\examples{
  LIvRatio(mtcars$drat, mtcars$wt)
  var.test(mtcars$drat, mtcars$wt)
  LIvRatio(mtcars$qsec, mtcars$wt)
  var.test(mtcars$qsec, mtcars$wt)
  LIvtest(sd(mtcars$qsec), nrow(mtcars), sd(mtcars$wt), nrow(mtcars))
}
\author{Kyun-Seop Bae k@acr.kr}
