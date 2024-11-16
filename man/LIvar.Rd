\name{LIvar}
\alias{LIvar}
\title{Likelihood Interval of variance and sd assuming Norml Distribution with sample mean and sample size}
\description{
  Likelihood interval of sd and variance assuming normal distribution. This is estimated likelihood interval, not profile likelihood interval.
}
\usage{
  LIvar(s1, n1, k, conf.level=0.95)
}
\arguments{
  \item{s1}{standard deviation of the sample}
  \item{n1}{sample size}
  \item{k}{1/k likelihood interval will be calculated}
  \item{conf.level}{approximately corresponding confidence level. If k is specified, this is ignored.}
}
\details{
  It calculates likelihood interval of sd and variance assuming normal distribution in one group. The likelihood interval is asymmetric and there is no standard error in the output.
}
\value{
  \item{PE}{maximum likelihood estimate on the population variance}
  \item{LL}{lower limit of likelihood interval}
  \item{UL}{upper limit of likelihood interval}
}
\examples{
  x = c(-5.3, -4.5, -1.0, -0.7, 3.7, 3.9, 4.2, 5.5, 6.8, 7.4, 9.3)
  LInormVar(x)
  LIvar(sd(x), length(x))
}
\author{Kyun-Seop Bae k@acr.kr}
