#' Apply the R2D2 multivariate bias-correction method.
#'
#' \code{r2d2} applies the R2D2 algorithm to correct the dependence structure
#' of the (previously univariately corrected) \code{bc1d} dataset.
#'
#' This function performs the R2D2 rank resampling conditionnaly to the conditioning dimensions,
#' to correct the dependence structure of the \code{P} variables of a dataset \code{bc1d}
#' in terms of rank (i.e. copula). The resulting dependence structure is adjusted with respect
#' to dependence structure of the reference dataset, \code{refdata}.
#'
#' The dataset \code{bc1d} contains data whose \code{P} variables have
#' already been corrected by a univariate (1D) bias-correction method. Note that \code{bc1d} and
#' \code{refdata} must not have any NA values.
#'
#' Ideally, the number of time steps in the \code{bc1d} dataset and in the
#' \code{refdata} reference dataset should be equal (\code{T = T_refdata}).
#' If this is not the case, the ranks of the largest dataset are shrinked
#' so that both datasets have the same maximal rank.
#'
#' For the full reference, see : \emph{Vrac, M.
#' Multivariate bias adjustment of high-dimensional climate simulations:
#' the Rank Resampling for Distributions and Dependences
#' (R2D2) Bias Correction. Hydrol. Earth Syst. Sci., 22, 3175-3196},
#' \url{https://doi.org/10.5194/hess-22-3175-2018}
#'
#' @param refdata A matrix of dimension \code{T_ref x P}
#' that contains the reference data, where \code{T_ref} is the number of
#' time steps and \code{P} is the number of variables.
#' @param bc1d A matrix of dimension \code{T x P}.
#' It contains the results of a 1D-bias
#' correction method applied over the T time steps to the \code{P} variables
#' (separately for each variable).
#' @param icond A vector of integers giving the indices of the variables
#' in \code{refdata} and \code{bc1d}  that are to be used as conditioning
#' dimensions when searching the best analogue in \code{refdata} for the
#' rank association of the conditioning dimension in \code{bc1d}
#' @param lag_search An integer corresponding to the number of time lags to account for
#' when searching in \code{refdata} the best analogue of the conditioning dimension rank
#' association observed in \code{bc1d}. The default value is no lag, i.e., \code{lag_search}=0.
#' @param lag_keep An integer corresponding to the number of time lags to keep in the correction
#' for each best analogue of rank associations found.
#' \code{lag_keep} has to be lower or equal to \code{lag_search}.
#' The default value is no lag, i.e., \code{lag_search}=0.

#' @return A list with the following elements:
#' \itemize{
#'   \item r2d2_bc A matrix of dimension \code{T x P} with
#'   the \code{bc1d} data corrected by the R2D2 bias-correction method
#'   \item visited_time A vector of length \code{T_ref} with the number of times
#'   that the rank of a specific time step in \code{refdata} at time \code{T_ref[i]}
#'   has been resampled to create \code{r2d2_bc}.
#'   \item time_bestanalogue A vector indentifying the time indices of the best rank association analogues found in \code{refdata} for each block of conditioning dimensions.
#'   The time index of the best analogue corresponds to the last time step of the rank association block.
#'   \item dist_bestanalogue A vector with the euclidean distances between the rank associations in \code{bc1d} and the best analogues found in \code{refdata} for each block of conditioning dimensions.
#' }
#' @examples
#' # Reproducing the example provided in Vrac (2018)
#'
#' refdata <- matrix(c(0.3, 0.5, 0.9, 0.8,
#'  1.1, 1.7, 1.2, 1.9,
#'  2.1, 1.8, 3.0, 2.7), ncol = 3, nrow = 4)
#'
#' bc1d <- matrix(c(0.7, 0.5, 0.2, 0.9,
#'  1.3, 1.8, 1.1, 1.4,
#'  1.9, 2.9, 2.0, 2.6), ncol = 3, nrow = 4)
#'
#' # 1 conditioning dimension, 0 lag
#' r2d2(refdata = refdata,
#'  bc1d = bc1d,
#'  icond = 1)
#'
#' r2d2(refdata = refdata,
#'  bc1d = bc1d,
#'  icond = 2)
#'
#' r2d2(refdata = refdata,
#'  bc1d = bc1d,
#'  icond = 3)
#'
#' # 1 conditioning dimension, 1 lag search, 1 lag keep
#' r2d2(refdata = refdata,
#'  bc1d = bc1d,
#'  icond = 1,
#'  lag_search = 1,
#'  lag_keep = 1)
#'
#' # 2 conditioning dimensions, 1 lag search, 1 lag keep
#' r2d2(refdata = refdata,
#'  bc1d = bc1d,
#'  icond = 1:2,
#'  lag_search = 1,
#'  lag_keep = 1)
#'
#'
#' # with climate data
#' data("r2d2_example")
#'
#' # conditioning dimension: temperature in Paris
#' r2d2_correction <-  r2d2(
#'   refdata = r2d2_example$refdata,
#'   bc1d =  r2d2_example$bc1d,
#'   icond =  1,
#'   lag_search = 8, lag_keep = 6 )
#'
#' # conditioning dimension: temperature and precipitation in Paris
#' r2d2_correction <-  r2d2(
#'   refdata = r2d2_example$refdata,
#'   bc1d =  r2d2_example$bc1d,
#'   icond =  c(1, 6),
#'   lag_search = 8, lag_keep = 6 )
#' @export
r2d2 <- function(refdata,
                bc1d,
                icond = c(1),
                lag_search = 0,
                lag_keep = 0) {
  # Looks for blocks of size lag_search+1
  # Keeps only the lag_keep+1 last values of the block
  if (lag_search < lag_keep) {
    stop("lag_search has to be greater or equal to lag_keep")
  }

  P <- ncol(bc1d) # number of var (including icond !!)
  Ntimes_BC <- nrow(bc1d)
  Ntimes_REF <-  nrow(refdata)
  stopifnot(lag_search <= min(Ntimes_BC, Ntimes_REF))
  if (lag_search > min(Ntimes_BC, Ntimes_REF)) {
    stop("lag_search has to be less or equal to the number of observations in refdata or bc1d")
  }

  r2d2_bc <- array(NaN, dim = dim(bc1d)) # (N days x P var)

  sorted_BC <-  apply(bc1d, 2, sort)
  ranks_BC <-  apply(bc1d, 2, rank, ties.method = "min")
  ranks_REF <-  apply(refdata, 2, rank, ties.method = "min")

  # shrink the rank
  if (Ntimes_BC < Ntimes_REF) {
    ranks_REF <-  round((ranks_REF - 1) * (Ntimes_BC - 1) / (Ntimes_REF - 1)) + 1
    cat("SHRINK 1\n")
  }
  if (Ntimes_BC > Ntimes_REF) {
    ranks_BC <-  round((ranks_BC - 1) * (Ntimes_REF - 1) / (Ntimes_BC - 1)) + 1
    ranks_sorted_shrinked <- round(seq.int(0, Ntimes_BC - 1) * (Ntimes_REF - 1) / (Ntimes_BC - 1)) + 1
    iduplicates <- (c(diff(ranks_sorted_shrinked), 1) == 0)
    sorted_BC  <- sorted_BC[!iduplicates, ]
    cat("SHRINK 2\n")
  }

  ranks_conddim_REF <-  ranks_REF[, icond, drop = FALSE]
  ranks_conddim_BC <-  ranks_BC[, icond, drop = FALSE]

  nsearch <-
    ((Ntimes_BC - lag_search - 1) %/% (lag_keep + 1)) +
    (((Ntimes_BC - lag_search - 1) %% (lag_keep + 1)) > 0) +
    1
  time_bestanalogue <- numeric(length = nsearch)
  dist_bestanalogue <- numeric(length = nsearch)
  visited_time <-  numeric(length = Ntimes_REF)

  t <-  0
  for (isearch in seq.int(nsearch)) {
    if (isearch == 1) {
      t <-  t + lag_search + 1
      iblockkeep <- seq.int(t - lag_search, t)
    } else if (isearch == nsearch) {
      iblockkeep <- seq.int(t + 1, Ntimes_BC)
      t <- Ntimes_BC
    } else{
      t <-  t + lag_keep + 1
      iblockkeep <- seq.int(t - lag_keep, t)
    }
    iblocksearch <- seq.int(t - lag_search, t)
    block_conddim_BC <-
      ranks_conddim_BC[iblocksearch, , drop = FALSE]
    bestanalogue <-
      find_bestanalogue_time(ranks_conddim_REF, block_conddim_BC)
    tstar <-  bestanalogue$tstar
    dist <-  bestanalogue$dist
    time_bestanalogue[isearch] <-  tstar
    dist_bestanalogue[isearch] <-  dist
    iblockref <- seq.int(tstar - length(iblockkeep) + 1, tstar)
    for (v in 1:P) {
      rank_REF_for_tstar_in_v <-  ranks_REF[iblockref, v]
      r2d2_bc[iblockkeep, v] <-
        sorted_BC[rank_REF_for_tstar_in_v, v]
    }
    visited_time[iblockref] <-  visited_time[iblockref] + 1
  }

  return(
    list(
      r2d2_bc = r2d2_bc,
      visited_time = visited_time,
      time_bestanalogue = time_bestanalogue,
      dist_bestanalogue = dist_bestanalogue
    )
  )
}



#' Finds the time index of the best analogue of rank association
#'
#' Finds the time index of the best analogue of \code{block_conddim_BC} in \code{conddim_REF}
#' in terms of the Euclidean distance.
#'
#' @param conddim_REF A matrix of dimension \code{T_ref x P_cond}
#' that contains the rank data where the best analogue of \code{block_conddim_BC} is searched.
#' \code{T_ref} is the number of
#' time steps and \code{P_cond} is the number of conditioning dimensions.
#' @param block_conddim_BC A matrix of dimension \code{T_block x P_cond}.
#' It corresponds to the block of rank association for which the best analogue is searched in \code{conddim_REF}.
#' \code{T_block} corresponds to a number of time steps in the block (i.e., its length) and
#' should be lower or equal to \code{T_ref}.
#'
#' @return A list with the following elements:
#' \itemize{
#'   \item tstar The time index of the best analogue of \code{block_conddim_BC} found in \code{conddim_REF}.
#'   This time index corresponds to the last time step of the block of rank association.
#'   \item dist The euclidean distance between \code{block_conddim_BC} and
#'   the best analogue found in \code{conddim_REF}.
#'  }
#' @examples
#' refdata <- matrix(c(0.3, 0.5, 0.9, 0.8,
#'                     1.1, 1.7, 1.2, 1.9,
#'                     2.1, 1.8, 3.0, 2.7), ncol = 3, nrow = 4)
#' bc1d <- matrix(c(0.7, 0.5, 0.2, 0.9,
#'                  1.3, 1.8, 1.1, 1.4,
#'                  1.9, 2.9, 2.0, 2.6), ncol = 3, nrow = 4)
#'
#' ranks_refdata <- apply(refdata, 2, rank)
#' ranks_bc1d <- apply(bc1d, 2, rank)
#'
#' # 1 conditioning dimension, 0 lag
#' icond <- 1
#' block_conddim_BC <- ranks_bc1d[1, icond, drop = FALSE]
#' conddim_REF <- ranks_refdata[, icond, drop = FALSE]
#' find_bestanalogue_time(conddim_REF, block_conddim_BC)
#'
#' # 1 conditioning dimension, 1 lag
#' icond <- 1
#' block_conddim_BC <- ranks_bc1d[1:2, icond, drop = FALSE]
#' conddim_REF <- ranks_refdata[, icond, drop = FALSE]
#' find_bestanalogue_time(conddim_REF, block_conddim_BC)
#'
#' # 2 conditioning dimensions, 1 lag
#' icond <- 1:2
#' block_conddim_BC <- ranks_bc1d[1:2, icond, drop = FALSE]
#' conddim_REF <- ranks_refdata[, icond, drop = FALSE]
#' find_bestanalogue_time(conddim_REF, block_conddim_BC)
#' @export
find_bestanalogue_time = function(conddim_REF,
                                  block_conddim_BC) {
  LAG = nrow(block_conddim_BC) - 1
  distmin = +Inf
  tstar = NaN
  start = 1 + LAG

  for (t in start:nrow(conddim_REF)) {
    DISTa = stats::dist(rbind(as.vector(block_conddim_BC),
                       as.vector(conddim_REF[(t - LAG):t, ]))) # Euclidean distance

    if (DISTa < distmin) {
      distmin = DISTa
      tstar = t
    }
  }

  return(list(tstar = tstar,
              dist = distmin))
}
#### end function find_analogue_time
