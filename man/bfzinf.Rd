% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/bfzinf.R
\name{bfzinf}
\alias{bfzinf}
\title{This function computes log(Bayes factor) under Jefferey's prior}
\usage{
bfzinf(y)
}
\arguments{
\item{y}{a data vector}
}
\value{
\item{lbfNBP}{log(Bayes factor) of NB against Poisson}
\item{lbfZINBNB}{log(Bayes factor) of ZINB against NB} 
\item{lbfZIPP}{log(Bayes factor) of ZIP against Poisson}
\item{lbfZINBZIP}{log(Bayes factor) of ZINB against ZIP}
\item{propzero}{Porportion of Zeros present in the data}
}
\description{
This function computes log(Bayes factor) under Jefferey's prior
}
\details{
NB: negative Binomial, ZINB: zero inflated negative Binomial,
         ZIP: zero inflated Poisson
}
\examples{
y <- rpois(n = 100, lambda = 5)
bfzinf(y)

}
\references{
Pramanik, P., and Maity, A. K. (2023) Bayes Factor of Zero Inflated Models under Jefferey's Prior.
}
