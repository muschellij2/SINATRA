% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ec_computation.R
\name{integrate_ec_curve}
\alias{integrate_ec_curve}
\title{Integrates an EC curve}
\usage{
integrate_ec_curve(ec_curve)
}
\arguments{
\item{ec_curve}{(matrix) : The Euler Characteristic Curve from \code{compute_standardized_ec_curve} or \code{compute_discrete_ec_curve}.}
}
\value{
ec_curve (matrix): The Integrated Euler Characteristic Curve (SECT)
}
\description{
Computes the numerical integral of an Euler Characteristic Curve curve by computinig the area under the EC curve.
}
