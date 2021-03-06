% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/example2.R
\name{example2}
\alias{example2}
\title{example2: Create the data used in the second setting of the simulation section of Li, Goldberg, \& Zheng}
\usage{
example2(
  seed = 1,
  n = 1000,
  var.est = 0.1,
  thres = 0.001,
  upper.age = 60,
  censoring = "simple"
)
}
\arguments{
\item{seed}{The seed for the simulation.}

\item{n}{The number of observations.}

\item{var.est}{variance of the score.}

\item{thres}{risk threshold for calculating the age to start screening}

\item{upper.age}{the upper time denoted by t in the paper}

\item{censoring}{Either "simple" or "complex". Simple means that the censoring time is independent of the failure time. Complex means that it dependent.}
}
\value{
A data frame with covariate 'Z',
observed time 'time'=Z,
the cometing risk 'status' which is 0 for censoring, 1 for the event of interest, and 2 for the competing risk=status,
the censoring time 'C', event of interest time 'T1', competing risk time 'T2',
and 'tau1' and 'tau2' which are the lower and upper times denoted by S and t in the paper
}
\description{
See details in the simulation section
}
