% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/config_files.R
\name{wn_config_common}
\alias{wn_config_common}
\title{Configuration options that are common across different WindNinja configuration types}
\usage{
wn_config_common(
  num_threads = max(1, availableCores() - 1),
  elevation,
  vegetation = "grass",
  mesh_choice = "coarse",
  output_formats = "asc",
  ascii_out_resolution = -1L,
  units_ascii_out_resolution = "m",
  shape_out_resolution = -1L,
  units_shape_out_resolution = "m",
  goog_out_resolution = -1L,
  units_goog_out_resolution = "m"
)
}
\arguments{
\item{num_threads}{numeric: number of threads to use}

\item{elevation}{string or Raster: a raster layer, or the filename of an elevation file (*.asc, *.lcp, *.tif, *.img)}

\item{vegetation}{string: dominant type of vegetation ("grass", "brush", "trees")}

\item{mesh_choice}{string: mesh resolution choice ("coarse", "medium", "fine")}

\item{output_formats}{character: one or more of "asc", "kmz", "shp"}

\item{ascii_out_resolution}{numeric: resolution of ascii output files (-1 to use mesh resolution)}

\item{units_ascii_out_resolution}{string: units of ascii output file resolution ("ft", "m")}

\item{shape_out_resolution}{numeric: resolution of shapefile output files (-1 to use mesh resolution)}

\item{units_shape_out_resolution}{string: units of shapefile output file resolution ("ft", "m")}

\item{goog_out_resolution}{numeric: resolution of Google Earth output file (-1 to use mesh resolution)}

\item{units_goog_out_resolution}{string: units of Google Earth resolution ("ft", "m")}
}
\value{
A named list of options
}
\description{
Configuration options that are common across different WindNinja configuration types
}
\examples{
opts <- wn_config_common(elevation = wn_demo_file("missoula_valley_elevation"),
                         ascii_out_resolution = 100)

}
\seealso{
\code{\link{wn_config_domain_average}}
}
