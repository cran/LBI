\name{LInorm}
\alias{LInorm}
\title{Likelihood Interval of mean, sd and variance assuming Normal Distribution}
\description{
  Likelihood interval of mean and sd assuming normal distribution. This is estimated likelihood interval, not profile likelihood interval.
}
\usage{
  LInorm(x, k, conf.level=0.95, PLOT="", LOCATE=FALSE, Resol=201)
}
\arguments{
  \item{x}{a vector of observation}
  \item{k}{1/k likelihood interval will be calculated}
  \item{conf.level}{approximately corresponding confidence level. If k is specified, this is ignored.}
  \item{PLOT}{"1d" for profile plot or "2d" for contour plot.}
  \item{LOCATE}{use locater. This works only with PLOT="2D" option.}
  \item{Resol}{resolution for plot. This works only with PLOT=TRUE option.}
}
\details{
  It calculates likelihood interval of mean and sd assuming normal distribution in one group. There is no standard error in the output.
}
\value{
  \item{PE}{maximum likelihood estimate}
  \item{LL}{lower limit of likelihood interval}
  \item{UL}{upper limit of likelihood interval}
}
\examples{
  x = c(-5.3, -4.5, -1.0, -0.7, 3.7, 3.9, 4.2, 5.5, 6.8, 7.4, 9.3)
  LInorm(x, k=1/0.15) # Pawitan Ex10-9 p289
  LInorm(x)
  LInorm(x, PLOT="1d")
  LInorm(x, PLOT="2d", LOCATE=TRUE)
}
\author{Kyun-Seop Bae k@acr.kr}
