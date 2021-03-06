% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/gp_inference.R
\name{LaplaceApproximation}
\alias{LaplaceApproximation}
\title{Use Laplace Approximation to Approximate mean & covariance}
\usage{
LaplaceApproximation(Kn, class_labels)
}
\value{
params (list): list of the posterior mean and variances.
}
\description{
\code{LaplaceApproximation}
Approximates the latent posterior with a Gaussian distribution; it does so by finding the mode of the posterior, and
using the Hessian (second order Taylor expansion) as an approximation of the covariance. Newton Raphson is used to find
the mode of the posterior.

Pseudocode taken from Rasmussen and Williams, Chapter 3. This function outputs the mean and covariance of the
approximated posterior. To actually generate samples from the latent posterior, generate samples from a multivariate
normal with the parameters returned by this function.


 @param Kn (matrix): the covariance matrix for the GP model
 @param class_labels (vector): +/- 1 values indicating the class labels of the data points
}
