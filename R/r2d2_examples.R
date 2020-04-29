#' Example data for the r2d2 function
#'
#' temperature and precipitation example data to run the r2d2 function
#'
#'
#' @format A list with the following elements:
#' \describe{
#'   \item{bc1d}{A matrix with 589 rows and 10 columns. It contains data previously bias corrected with a univariate BC method.
#'   Those 1d-BC data have inter-site and inter-variable dependences that are to be corrected by the r2d2 method.}
#'   \item{refdata}{A matrix with 589 rows and 10 columns. It contains the reference data to correct the bc1d data.
#'   The r2d2 function resamples rank associations from this dataset, depending on the conditioning dimension values from the bc1d dataset to be corrected.}

#'   For those two matrices, the number of rows corresponds to the number of time-steps.
#'   The number of columns corresponds to the numbers of variables to be corrected.
#'   The first 5 columns correspond to temperatures at the grid points representing respectively Paris, Madrid, Rome, Warsaw and Stockholm.
#'   The following 5 columns correspond to precipitation at the same geographical sites.
#'
#'
#'   The reference dataset comes from \href{http://www.eu-watch.org/data_availability}{WATCH-Forcing-Data-ERA-Interim} (WFDEI).
#'   The bc1d dataset comes from simulations of the \href{http://cmc.ipsl.fr/international-projects/cmip5/}{IPSL-CM5A-LR climate model}
#'   that have been (i) regridded to the same spatial resolution as the reference dataset
#'   and then (ii) bias-corrected with the univariate BC method CDFt (Michelangeli et al, 2009).
#'   The bc1d dataset has been prepared by \href{https://theclimatedatafactory.com/}{The Climate Data Factory}.
#' }
#'
#' @source  Weedon, G.P., Balsamo, G., Bellouin, N., Gomes, S., Best, M.J. and Viterbo, P., 2014.
#' The WFDEI meteorological forcing data set:
#' WATCH Forcing Data methodology applied to ERA-Interim
#' reanalysis data. Water Resources Research, 50, doi:10.1002/2014WR015638
#'
#' @source Michelangeli, P.A., Vrac, M., and Loukos, H. ( 2009),
#' Probabilistic downscaling approaches:
#' Application to wind cumulative distribution functions,
#' Geophys. Res. Lett., 36, L11708, doi:10.1029/2009GL038401.
"r2d2_example"

# load(file = "data/r2d2_example.rdata")

