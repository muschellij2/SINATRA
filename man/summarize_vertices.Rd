% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/shape_reconstruction.R
\name{summarize_vertices}
\alias{summarize_vertices}
\title{Reconstruct Vertices in a Cone}
\usage{
summarize_vertices(
  dir,
  complex,
  rate_vals,
  len,
  reduction_operation = intersect,
  threshold,
  cone_size,
  ball = TRUE,
  ball_radius = 1,
  radius = 0
)
}
\arguments{
\item{dir}{(nx3 matrix) : The matrix of directions in the cone that were used to compute the (S/D) EC curve over.}

\item{complex}{(list) : The list containing metadata of the Vertices, Edges, and Faces of the mesh (use process_off_file_v3 to obtain).}

\item{rate_vals}{(vector) : Vector of variable importances for each sub-level set across each direction in a given cone.}

\item{len}{(int) : The number of sub-level sets to compute the (S/D) EC curve on in each direction.}

\item{reduction_operation}{(function) The function to summarize the sets of vertices. We default to intersect, and recommend it.}

\item{threshold}{(float) : The threshold for determining which sub-level sets are used for the reconstruction.}

\item{cone_size}{(int) : The number of directions in each cone.}

\item{ball}{(boolean) : Determining whether or not to compute the (S/D) EC curve on a ball for uniform comparisons.}

\item{ball_radius}{(float) : The radius of the ball to compute the (S/D) EC on; if you want the curve to be computed relative to the shape, don't touch this parameter.}

\item{radius}{(int) : The number of sub-level sets "before" and "after" the selected sub-level sets we want to include.
Setting Radius = 1 is recommened.}

\item{ec_type}{(string) : What type of EC curve to compute. Currently we support ECT (Euler Characteristic Curve), SECT (Smooth Euler Characteristic Curve)}
}
\value{
final_selected_vertices (vector) : Vector of the vertex indices that are reconstructed for a given cone.
}
\description{
Given a cone of directions and a threshold, We take the sub-level sets that are above the given threshold, and for each direction in the cone,
we compute the set of vertices that are associated with the selected sub-level sets. We then intersect the set of
vertices that are associated with the sub-level sets across each direction in the cone; only taking the vertices that are reconstructed in each direction of the cone.
}
