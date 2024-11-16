\name{LIpois}
\alias{LIpois}
\title{Likelihood Interval of the Mean assuming Poisson Distribution}
\description{
  Likelihood interval of lambda assuming Poisson distribution.
}
\usage{
  LIpois(x, k, n = 1, conf.level = 0.95, eps = 1e-8)
}
\arguments{
  \item{x}{raw data vector or a mean value.  If the length of x is 1, x is considered as a mean.}
  \item{k}{1/k likelihood interval will be calculated.}
  \item{n}{number of observations. If the length of x is 1, x is considered as the mean.}
  \item{conf.level}{approximately corresponding confidence level. If k is specified, this is ignored.}
  \item{eps}{estimated values less than this eps are considered as 0.}
}
\details{
  It calculates likelihood interval of mean(lambda) assuming Poisson distribution. The likelihood interval is asymmetric and there is no standard error in the output.
}
\value{
  \item{PE}{maximum likelihood estimate on the lambda}
  \item{LL}{lower limit of likelihood interval}
  \item{UL}{upper limit of likelihood interval}
  \item{n}{number of observations}
  \item{k}{1/k likelihood interval provided}
  \item{logk}{log(k) of k value}
  \item{maxLL without factorial}{maximum log likelihood without factorial part}
}
\examples{
  LIpois(4, k=15)     # Fisher
  LIpois(4, k=exp(2)) # Edwards
  LIpois(4, k=1/0.15) # Pawitan
  LIpois(4, k=8)      # Rhode
  LIpois(4, n=4)      # Bae
  LIpois(4)           # Bae
#  poisson.test(4)
  LIpois(4, k=32)     # 0.7454614 11.7893612
  LIpois(2.1, n=60)   # 1.750222   2.493533
}
\author{Kyun-Seop Bae k@acr.kr}
