\name{LIvtest}
\alias{LIvtest}
\title{Likelihood Interval of variance and sd assuming Norml Distribution using means and SDs}
\description{
  Likelihood interval of variance and sd assuming normal distribution. This is estimated likelihood interval, not profile likelihood interval.
}
\usage{
  LIvtest(s1, n1, s2, n2, k, conf.level=0.95)
}
\arguments{
  \item{s1}{sample standard deviation of the first group}
  \item{n1}{sample size of the first group}
  \item{s2}{sample standard deviation of the second group}
  \item{n2}{sample size of the second group}
  \item{k}{1/k likelihood interval will be calculated}
  \item{conf.level}{approximate confidence level. If k is specified, this is ignored.}
}
\details{
  It calculates likelihood interval of variance and sd using sufficient statistics. There is no standard error in the output.
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
  LIvtest(10.5, 3529, 8.9, 5190)
  LIvtest(3, 10, 2, 10)
  LIvtest(3, 10, 2, 10, k=15)
}
\author{Kyun-Seop Bae k@acr.kr}
