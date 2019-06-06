#' Purpose of simulation" understand impact of different kernels in the power of SINATRA procedure.
#'
#' We are currently using the squared exponential / Gaussian kernel, but there are other ideas:
#' some suggestions say that using polynomial kernels might offer advantages for high dimensional classification problems
#' We could also use a more geometrically inclined kernel like in 'Gaussian Process Landmarking' - at least one that reflects the
#' EC-cone feature structure of our data.
#'
#' Perhaps other covariance functions like the Matern, or OU covariance functions can help.
#'
#' We should also try custom ones.
#'
#'
#' Simulation: measure performance via ROC curves on the easy medium difficult simulations?
