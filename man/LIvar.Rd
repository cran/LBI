\name{LIvar}
\alias{LIvar}
\title{Likelihood Interval of the ratio of two variances from two groups}
\description{
  Likelihood interval of  the ratio of two variances from two groups assuming normal distribution. Likelihood interval usually gives a narrower interval when the likelihood function is asymmetric.
}
\usage{
  LIvar(x, y, k, conf.level=0.95)
}
\arguments{
  \item{x}{observations from the first group, used for the numerator}
  \item{y}{observations from the second group, used for the denominator}
  \item{k}{1/k likelihood interval will be provided}  
  \item{conf.level}{approximate confidence level}
}
\details{
  It calculates likelihood interval of the ratio of two variances from two groups. Likelihood interval usually gives a narrower interval when the likelihood function is asymmetric.
}
\value{
  \item{v1}{variance of the first group}
  \item{v2}{variance of the second group}
  \item{PE}{point estimation for the proportion}
  \item{LL}{lower limit of likelihood interval}
  \item{UL}{upper limit of likelihood interval}
  \item{logk}{log(k) value}
  \item{maxLL}{maximum log likelihood}
}
\examples{
  LIvar(mtcars$drat, mtcars$wt)
  var.test(mtcars$drat, mtcars$wt)
  LIvar(mtcars$qsec, mtcars$wt)
  var.test(mtcars$qsec, mtcars$wt)
}
\author{Kyun-Seop Bae k@acr.kr}
