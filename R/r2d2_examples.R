#' data for the r2d2 example
#'
#' temperature and precipitation data needed to run the r2d2 example
#'
#'
#' @format A list with the following elements:
#' \describe{
#'   \item{bc1d}{A matrix with 589 rows and 8334 columns. It contains the data to be corrected by the r2d2 method.}
#'   \item{refdata}{A matrix with 589 rows and 8334 columns. It contains the reference data. The r2D2 reproduces rank associations observed in this dataset depending on the conditioning dimension values present in in the bc1d dataset to be corrected.}

#'   For those two matrices, the number of rows corresponds to the number of observations or time-steps. The number of columns corresponds to the numbers of variables to be corrected. The first 4167 columns corresponds to temperature variables at different geographical sites. The following 4167 columns corresponds to precipitaion variables for the same site as for the temperature.
#'   The reference datasets come from \href{http://www.eu-watch.org/data_availability}{WATCH-Forcing-Data-ERA-Interim}(WFDEI).
#'   The bc1d datasets to be corrected come from \href{http://cmc.ipsl.fr/international-projects/cmip5/}{the IPSL-CM5A-LR climate model} and have been interpolated to the same grid as WFDEI.
#'   For those bc1d matrices to be corrected by R2D2, note that the marginal distributions of these datasets have already been corrected by a unidimensional bias correction method (\code{\link[CDFt]{CDFt}}) to have the same properties as the reference dataset. The bc1d datasets have been prepared by \href{https://theclimatedatafactory.com/}{The Climate Data Factory}.
#'   \item{datesl}{vector of strings with dates used for the two datasets. For this examples, the datasets only consist of days in January during the 1979-1997 period.}
#'   \item{lon}{numeric vector of length 4167 with the longitudes of the 4167 geopgraphical locations where the simulations of temperarure and precipiation are provided.}
#'   \item{lat}{numeric vector of length 4167 with the latitudes of the 4167 geopgraphical locations where the simulations of temperarure and precipiation are provided.}
#'   \item{icond}{vector of integers containings the indices of the closest grid points to the following cities: Paris, Madrid, Rome, Warsaw, Stockholm. Temperature and precipitation at those five location will be used as conditioning dimension in the r2d2 examples. For instance, to have the longitude of the closest grid point closest to Paris, one could type: \code{r2d2_example$lon[r2d2_example$icond["paris"]]}}
#' }
"r2d2_example"

# load(file = "data/r2d2_example.rdata")


