\name{OneTwo}
\alias{OneTwo}
\title{Likelihood Ratio Test for One group vs Two group gaussian mixture model}
\description{
 With a given vector, it performs likelihood ratio test which model - one or two group - is better.
}
\usage{
  OneTwo(x, alpha=0.05)
}
\arguments{
  \item{x}{a vector of numbers}
  \item{alpha}{alpha value for type I error, significance level}
}
\details{
  It performs likelihood ratio test using both F distribution and Chi-square distribution (by Wilks' theorem).
}
\value{
  \item{Estimate}{n, Mean, SD for each group assumption and prior probability of each group in two group model}
  \item{Delta}{delta number of parameters and log-likelihoods}
  \item{Statistic}{Statistics from both the F distribution and Chi-square distribtuion. Cutoff is in terms of log-likelihood not the statistic.}
}
\examples{
  OneTwo(c(7, 5, 17, 13, 16, 5, 7, 3, 8, 10, 8, 14, 14, 11, 14, 17, 2, 12, 15, 19))
  OneTwo(c(5, 3, 0, 6, 5, 2, 6, 6, 4, 4, 15, 13, 18, 18, 19, 14, 19, 13, 19, 18))
}
\author{Kyun-Seop Bae k@acr.kr}
