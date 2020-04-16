impute_ranks <- function(ranks_refdim, ranks_p, na.last="keep", ties.method = "first"){
  ina <- which(is.na(ranks_p))
  if(length(ina) == 0) return(ranks_p)
  ranks_p_completed <- ranks_p
  ranks_refdim_na <- ranks_refdim[ina] 
  ranks_refdim_nna <- ranks_refdim[-ina] 
  ranks_p_nna <- ranks_p[-ina] 
  distmat <- matrix(ranks_refdim_nna,
                    nrow = length(ranks_refdim_na),
                    ncol = length(ranks_refdim_nna),
                    byrow = TRUE)
  distmat <- abs(ranks_refdim_na - distmat)
  imin <- apply(distmat, 1, which.min)
  ranks_p_completed[ina] <- ranks_p_nna[imin] + 0.5
  ranks_p_completed = rank(ranks_p_completed, na.last = na.last, ties.method = ties.method)
  return(ranks_p_completed)
}

#' Imputation of ranks in the reference dataset rank matrix
#'
#' \code{impute_refdata_ranks} complete the NAs in the rank matrix of the reference
#' dataset \code{refdata_ranks}.
#'
#' This function replaces NAs  in \code{refdata_ranks} the rank matrix of the
#' reference dataset. This is done by taking the rank
#' corresponding to the nearest-neighbor of NA value where the distance 
#' is evaluated with respect to the rank of the reference dimension \code{ranks_refdim}.
#' The ranks are then rescaled so that ranks uniformly go from 1 to \code{T_refdata}, the number of timesteps.
#' When rescaling ranks, when two ranks have the same value, the rank appearing first in the 
#' chronological order is given the smaller rank.
#' For more details, see the small example below.
#'
#' @param ranks_refdim A vector containing the ranks of the reference dimension
#' @param ranks_refdata A matrix of dimension \code{T_refdata x N} 
#' that contains the ranks of the reference data.
#' @param na.last Argument for the function rank used when imputed ranks are rescaled
#' @param ties.method Argument for the function rank used when imputed ranks are rescaled
#'
#' @return A matrix of dimension \code{T_refdata x N}. 
#' It is the matrix of ranks where NAs have been imputed by taking the value
#' of the nearest-neighbor, where nearest-neighbors are defined with respect to 
#' the reference dimension \code{ranks_refdim}
#' 
#' @examples
#' # Example from Vrac(2018)
#'  
#' refdata <- matrix(c(0.3, 0.5, 0.9, 0.8,
#'  1.1, 1.7, 1.2, 1.9,
#'  2.1, 1.8, 3.0, 2.7), ncol = 3, nrow = 4)
#' 
#' ranks_refdata <- apply(refdata, 2, rank, na.last = "keep", ties = "first")
#' ranks_refdata_na <- ranks_refdata
#' 
#' # First NA, at second variable, second line, corresponds to reference dimension rank 2.
#' # Closest rank is 1 or 3 in reference dimension, by default select the lowest rank.
#' # Rank 1 in reference dimension corresponds to rank 1 in variable 2.
#' # NA  at second variable, second line is replaced by 1.5
#' # Ranks for variable two are rescaled (1, 1.5, 2, 4) -->  (1, 2, 3, 4)
#' ranks_refdata_na[2, 2] <- NA
#' 
#' # Second NA, at third variable, third line, corresponds to reference dimension rank 4.
#' # Closest rank is 3 in reference dimension.
#' # Rank 3 in reference dimension corresponds to rank 4 in variable 3.
#' # NA  at third variable, third line is replaced by 4.5
#' # Ranks for variable three are rescaled (2, 1, 4.5, 3) -->  (2, 1, 4, 3)
#' ranks_refdata_na[3, 3] <- NA
#' 
#' cat("Imputation of NA values \n")
#' ranks_refdata_imputed <- impute_refdata_ranks(ranks_refdim = ranks_refdata_na[, 1],
#'   ranks_refdata = ranks_refdata_na,
#'   na.last="keep",
#'   ties.method = "first")
#'   
#' @export
impute_refdata_ranks <- function(ranks_refdim, ranks_refdata, na.last="keep", ties.method = "first"){
  if(any(is.na(ranks_refdim))) stop("NAs in rank_refdim")
  ranks_refdata_modified <- apply(ranks_refdata,
                                  2,
                                  impute_ranks,
                                  ranks_refdim = ranks_refdim,
                                  na.last="keep",
                                  ties.method = "first")
  return(ranks_refdata_modified)
}


r2d2_shuffle_1ref <- function(iref, ranks_refdata, ranks_bc1d, sorted_bc1d){
  np <- ncol(ranks_bc1d)
  nt <- nrow(ranks_bc1d)
  order_iref <-  order(ranks_refdata[,iref]) # = num de times (not rank values)
  ranks_refdata_ordered <-  ranks_refdata[order_iref, ]
  # needed in case of shrinking
  iduplicate <- which(diff(ranks_refdata_ordered[, iref]) == 0)
  if(length(iduplicate) > 0){
    ranks_refdata_ordered <- ranks_refdata_ordered[-iduplicate, ] # matrix of desired ranks
  }
  rneeded <-  ranks_refdata_ordered[ranks_bc1d[, iref],]
  # indexes to find the values corresponding to the desired ranks in the 1d-bias-corrected data
  shuffle_idx <- matrix(c(rneeded, rep(seq.int(ncol(rneeded)), rep(nt, np))),
                        ncol = 2)
  bc1d_shuffled <- array(sorted_bc1d[shuffle_idx], dim = dim(rneeded))
  return(bc1d_shuffled)
}

#' Apply the R2D2 multivariate bias-correction method.
#'
#' \code{r2d2} applies the R2D2 shuffling so that the \code{bc1d} dataset possesses
#' at each timestep the same intervariable dependences as
#' the reference dataset, \code{refdata}.
#'
#' This function performs the R2D2 shuffling to correct the 
#' dependence structure of the N variables in terms of rank (i.e. cppula) of
#' a dataset \code{bc1d} so that at each timestep
#' it matches the dependence structure observed in the
#' reference dataset, \code{refdata}.
#'
#' The dataset \code{bc1d} contains data whose \code{N} variables have
#' already been corrected by a 1D bias-correction method. Note that \code{bc1d}
#' must not have any NA values.
#'
#' By applying the R2D2 method several times with different \code{iref} values,
#' R2D2 can provide up to \code{N} different corrections that respect the rank
#' dependence structure of the reference dataset. These \code{N} corrections differ
#' only based on the choice of a reference variable in the dataset \code{bc1d} which
#' dictates the temporal evolution of the variable. The reference variable is
#' the only variable for which the time-serie of the rank is garanteed to be
#' unaltered.
#'
#' Ideally, the number of timesteps in the \code{bc1d} dataset and in the
#' \code{refdata} reference dataset should be equal (\code{T = T_refdata}).
#' If this is not the case, the ranks of the largest dataset are shrinked
#' so that both datasets have the same maximal rank.
#'
#' The functions only works if the reference dimension in the \code{refdata} 
#' reference dataset does not have any NA values. If for other variables
#' of the \code{refdata} reference dataset , NA values are present,
#' the ranks for those NA values are imputed with 
#' the function \code{\link{impute_refdata_ranks}}.
#'
#' For the full reference, see : \emph{Vrac, M.
#' Multivariate bias adjustment of high-dimensional climate simulations:
#' the Rank Resampling for Distributions and Dependences
#' (R2D2) Bias Correction. Hydrol. Earth Syst. Sci., 22, 3175-3196},
#' \url{https://doi.org/10.5194/hess-22-3175-2018}
#'
#' @param refdata A matrix of dimension \code{T_ref x N} 
#' that contains the data of reference where \code{T_ref} is the number of 
#' timesteps and \code{N} is the number of variables.
#' @param bc1d A matrix of dimension \code{T x N}.
#' It contains the results of a 1D-bias
#' correction method applied over the T timesteps to the N variables 
#' (separately for each variable). 
#' @param iref An integer giving the index of the variable
#' that will be used as the reference dimension (iref must be between 1 and N).
#'
#' @return A matrix of dimension \code{T x N}. 
#' It contains one possible bias-correction
#' dataset that respects the dependence structure of the N variables
#' in terms of ranks (i.e. the copula) exhibited in the reference dataset.
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
#'  iref = 1)
#'
#' r2d2(refdata = refdata,
#'  bc1d = bc1d,
#'  iref = 2)
#'  
#' r2d2(refdata = refdata,
#'  bc1d = bc1d,
#'  iref = 3)
#' @export
r2d2 <- function(refdata, bc1d, iref = 1){
  
  # Checks
  if(ncol(bc1d) != ncol(refdata)) stop("ncol(ranks_bc1d) != ncol(ranks_refdata)")
  ranks_refdata <- apply(refdata, 2, rank, na.last = "keep", ties = "first")
  if(any(is.na(ranks_refdata[, iref]))) stop("NAs in the reference dimension")
  if(any(is.na(bc1d))) stop("NAs in the bc1d dataset")
  
  # Imputation for NA
  if(any(is.na(ranks_refdata))){
    warning("NA values in refdata \n => imputing ranks !!!")
    ranks_refdata <- impute_refdata_ranks(ranks_refdim = ranks_refdata[, iref],
                                          ranks_refdata = ranks_refdata)
  }
  
  # Rank shrinking
  ranks_bc1d <- apply(bc1d, 2, rank, na.last = "keep", ties = "first")
  if(nrow(ranks_refdata) > nrow(ranks_bc1d)){
    warning("nrow(ranks_refda) > nrow(ranks_bc1d) \n => shrinking applied !!!")
    ranks_refdata <- round((ranks_refdata - 1) * (nrow(ranks_bc1d) - 1) / (nrow(ranks_refdata) - 1)) + 1
  }
  if(nrow(ranks_bc1d) > nrow(ranks_refdata)){
    warning("nrow(ranks_bc1d) > nrow(ranks_refdata) \n => shrinking applied !!!")
    rank_bc1d <- round((ranks_bc1d - 1) * (nrow(ranks_refdata) - 1) / (nrow(ranks_bc1d) - 1)) + 1
  }
  
  sorted_bc1d <- apply(bc1d, 2, sort)
  bc_r2d2d <- r2d2_shuffle_1ref(iref = iref,
                                   ranks_refdata = ranks_refdata,
                                   ranks_bc1d = ranks_bc1d,
                                   sorted_bc1d = sorted_bc1d)
  dimnames(bc_r2d2d) <- dimnames(ranks_refdata)
  bc_r2d2d
}
