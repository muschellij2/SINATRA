% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/clean_mesh.R
\name{clean_files}
\alias{clean_files}
\title{Scale & Translate Meshes}
\usage{
clean_files(input_dir, output_dir)
}
\arguments{
\item{input_dir}{The input directory of the files to scale & translate.}

\item{output_dir}{The output directory to put the files post scaling and translating.}
}
\value{
None; files written in output path.
}
\description{
Scales the meshes to have unit surface area and translates the meshes to the origin.
}
