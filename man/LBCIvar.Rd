\name{LBCIvar}
\alias{LBCIvar}
\title{Likelihood Based Confidence Interval of sd and variance assuming Norml Distribution}
\description{
  Likelihood based confidence interval of sd and variance assuming normal distribution. It usually shows narrower interval than convenrtional chi-square interval. This uses estimated likelihood, not profile likelihood. 
}
\usage{
  LBCIvar(x, conf.level=0.95)
}
\arguments{
  \item{x}{a vector of observation}
  \item{conf.level}{confidence level}
}
\details{
  It calculates (same height) likelihood based confidence interval of sd and variance assuming normal distribution in one group. The likelihood interval is asymmetric and there is no standard error in the output.
}
\value{
  \item{PE}{point estimation for the proportion}
  \item{LL}{lower limit of likelihood interval}
  \item{UL}{upper limit of likelihood interval}
}
\examples{
  LBCIvar(lh)
  (length(lh) - 1)*var(lh)/qchisq(c(0.975, 0.025), length(lh) - 1)
}
\author{Kyun-Seop Bae k@acr.kr}
