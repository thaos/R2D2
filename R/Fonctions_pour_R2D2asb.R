#' Apply the R2D2 multivariate bias-correction method.
#'
#' \code{r2d2} applies the R2D2 shuffling so that the \code{bc1d} dataset possesses
#' the same dependences structure as the reference dataset, \code{refdata}.
#'
#' This function performs the R2D2 shuffling to correct the
#' dependence structure of the \code{P} variables in terms of rank (i.e. copula) of
#' a dataset \code{bc1d} so that it matches the dependence structure observed in the
#' reference dataset, \code{refdata}, conditionnaly to the conditioning dimensions.
#'
#' The dataset \code{bc1d} contains data whose \code{P} variables have
#' already been corrected by a 1D bias-correction method. Note that \code{bc2d} and
#' \code{refdata} must not have any NA values.
#'
#' If the conditionning dimension is univariate, the conditioning dimension is is
#' the only variable for which the time-serie of the rank is garanteed to be
#' unaltered.
#'
#' Ideally, the number of timesteps in the \code{bc1d} dataset and in the
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
#' that contains the data of reference where \code{T_ref} is the number of
#' timesteps and \code{P} is the number of variables.
#' @param bc1d A matrix of dimension \code{T x P}.
#' It contains the results of a 1D-bias
#' correction method applied over the T timesteps to the N variables
#' (separately for each variable).
#' @param icond A vector of integers giving the index of the variables
#' in \code{refdata} and \code{bc1d}.
#' that will be used as the reference dimensionwhen searching in \code{refdata}
#' the best analogue of the conditioning dimension rank association observed in \code{bc1d}
#' @param lag_search An integer denoting the number of time lags to account for
#' when searching in \code{refdata} for the best analogue of the conditioning dimension rank association observed in \code{bc1d}. The default value is no lag.
#' @param lag_keep An integer denoting the number of time lags to keep in the r2d2 correction for each best analogue of rank associations found. \code{lag_keep} has to be less or equal to \code{lag_search. The default value is no lag. \code{lag_search}.

#' @return A list with the following elements:
#' \itemize{
#'   \item r2d2_bc A matrix of dimension \code{T x P} with
#'   the \code{bc1d} data corrected by the R2D2 bias-correction method
#'   \item visited_time A vector of length \code{T_ref} with the number of times
#'   that the rank association observed in \code{refdata} at time \code{T_ref[i]}
#'   has been used in  \code{r2d2_bc}.
#'   \item time_bestanalogue A vector of time indices indentifying the position of the best analogue of rank association found in \code{refdata} for each block of conditioning dimensions. The time index of the best analogue is given by the last time step of the block of rank association.
#'   \item dist_bestanalogue A vector of with the euclidean distances between the rank associations observed in \code{bc1d} and the best analogue found in \code{refdata} for each block of conditioning dimensions.
#' }
#' @examples
#' # Reproducing the example provided in Vrac(2018)
#'
#' refdata <- matrix(c(0.3, 0.5, 0.9, 0.8,
#'  1.1, 1.7, 1.2, 1.9,
#'  2.1, 1.8, 3.0, 2.7), ncol = 3, nrow = 4)
#'
#' bc1d <- matrix(c(0.7, 0.5, 0.2, 0.9,
#'  1.3, 1.8, 1.1, 1.4,
#'  1.9, 2.9, 2.0, 2.6), ncol = 3, nrow = 4)
#'
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
#'  # with climate data from r2d2_examples
#'  data("r2d2_example")
#'  r2d2_correction <- with(r2d2_example,
#'   R2D2(
#'     refdata = refdata,
#'     bc1d = bc1d,
#'     icond = icond,
#'     lag_search = 8,
#'     lag_keep = 6
#'    )
#'  )
#' @export
R2D2 <- function(refdata,
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
    ranks_REF <-  round(ranks_REF * Ntimes_BC / Ntimes_REF)
    cat("SHRINK 1\n")
  }
  if (Ntimes_BC > Ntimes_REF) {
    ranks_BC <-  round(ranks_BC * Ntimes_REF / Ntimes_BC)
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
  visited_time <-  numeric(length = Ntimes_BC)

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



#' Find the time index of the best analogue of rank association
#'
#' Find the time index of the best analogue of \code{block_conddim_BC} in \code{conddim_REF}.
#'
#' @param conddim_REF A matrix of dimension \code{T_ref x P_cond}
#' that contains the data where the best analogue of \code{block_conddim_BC} is searched for
#' \code{T_ref} Is the number of
#' time-steps and \code{P_cond} is the number of conditioning dimensions.
#' @param block_conddim_BC A matrix of dimension \code{T_block x P_cond}.
#' The block of rank association for which
#' the best analogue is searched for in \code{conddim_REF}.
#' \code{T_block} corresponds to a number of time-steps in the block and
#' should be less or equal to \code{T_ref}
#'
#' @return A list with the following elements:
#' \itemize{
#'   \item tstar The time index of the best analogue of \code{block_conddim_BC} found in \code{conddim_REF}. The time index of the best analogue corresponds to the time index of the last time step of the block of rank association.
#'   \item dist The euclidean distance between \code{block_conddim_BC} and
#'   the best analogue found in \code{conddim_REF}.
#'  }
#'  @examples
#'  # Reproducing the example provided in Vrac(2018)
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
#
#' # 2 conditioning dimension, 1 lag
#' icond <- 1:2
#' block_conddim_BC <- ranks_bc1d[1:2, icond, drop = FALSE]
#' conddim_REF <- ranks_refdata[, icond, drop = FALSE]
#' find_bestanalogue_time(conddim_REF, block_conddim_BC)

find_bestanalogue_time = function(conddim_REF,
                                  block_conddim_BC) {
  LAG = nrow(block_conddim_BC) - 1
  distmin = +Inf
  tstar = NaN
  start = 1 + LAG

  for (t in start:nrow(conddim_REF)) {
    DISTa = dist(rbind(as.vector(block_conddim_BC),
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
