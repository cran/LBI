\name{LInormVar}
\alias{LInormVar}
\title{Likelihood Interval of sd and variance assuming Norml Distribution}
\description{
  Likelihood interval of sd and variance assuming normal distribution. This is estimated likelihood interval, not profile likelihood interval.
}
\usage{
  LInormVar(x, k, conf.level=0.95)
}
\arguments{
  \item{x}{a vector of observation}
  \item{k}{1/k likelihood interval will be calculated}
  \item{conf.level}{approximately corresponding confidence level. If k is specified, this is ignored.}
}
\details{
  It calculates likelihood interval of sd and variance assuming normal distribution in one group. The likelihood interval is asymmetric and there is no standard error in the output.
}
\value{
  \item{PE}{maximum likelihood estimate}
  \item{LL}{lower limit of likelihood interval}
  \item{UL}{upper limit of likelihood interval}
}
\examples{
  x = c(-5.3, -4.5, -1.0, -0.7, 3.7, 3.9, 4.2, 5.5, 6.8, 7.4, 9.3)
  LInormVar(x, k=1/0.15) # Pawitan Ex10-9 p289
  LInormVar(x)
}
\author{Kyun-Seop Bae k@acr.kr}
