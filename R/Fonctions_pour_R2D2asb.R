
# Contains 2 functions:
# - R2D2_asb : performing R2D2 with analogs-based conditioning with semi-block approach
# - find_analogue_time : called by R2D2asb

#####################################################

R2D2 = function(
  REF, BC,
  icond = c(1), lag_search = 2, lag_keep = 1
){

  # Looks for blocks of size lag_search+1
  # Keeps only the lag_keep+1 last values of the block

  if(lag_search < lag_keep){
    stop("lag_search < lag_keep => Stop\n")
  }

  # Nicond = length(icond)

  # dim(REF) Ntimes x NVar
  P = ncol(BC) # number of var (including icond !!)
  Ntimes_BC = Max_rank_BC = nrow(BC)
  Ntimes_REF = Max_rank_REF = nrow(REF)

  BC_corrdep = array(NaN, dim = dim(BC)) # (N days x p var)

  TS_rank_BC = array(NaN, dim = dim(BC)) # N x p
  TS_rank_REF = array(NaN, dim = dim(REF)) # N x p
  for(v in 1:P){
    TS_rank_BC[, v] = rank(BC[, v], ties.method = "min")
    TS_rank_REF[, v] = rank(REF[, v], ties.method = "min")
  }

  #############
  # SHRINK THE RANK VECTORS SHRINK THE RANK VECTORS SHRINK THE RANK VECTORS SHRINK THE RANK VECTORS SHRINK THE RANK VECTORS
  if(Ntimes_BC < Ntimes_REF){
    TS_rank_REF = round(TS_rank_REF * Max_rank_BC / Max_rank_REF)
    cat("SHRINK 1\n")
  }
  if(Ntimes_BC > Ntimes_REF){
    TS_rank_BC = round(TS_rank_BC * Max_rank_REF / Max_rank_BC)
    cat("SHRINK 2\n")
  }
  # SHRINK THE RANK VECTORS SHRINK THE RANK VECTORS SHRINK THE RANK VECTORS SHRINK THE RANK VECTORS SHRINK THE RANK VECTORS
  #############


  Sorted_BC = array(NaN, dim = c(Ntimes_BC, P))
  for(v in 1:P){
    Sorted_BC[, v] = sort(BC[, v])
  }

  #####################
  TS_rank_REF_icond = TS_rank_REF[, icond, drop = FALSE]
  TS_rank_BC_icond = TS_rank_BC[, icond, drop = FALSE]


  ##########################

  Visited_t = numeric(length = Ntimes_BC)
  TSTAR = c() ; DIST = c()

  t = 1+lag_search

  while(t <= Ntimes_BC){
    MultiVec_rank_BC_lag_t = TS_rank_BC_icond[(t-lag_search) : t,, drop = FALSE]   # t x r = 1 x 800 for 400_2_0
    FF = find_analogue_time(TS_rank_REF_icond, MultiVec_rank_BC_lag_t)
    tstar = FF$tstar
    dist = FF$dist
    TSTAR = c(TSTAR, tstar)
    DIST = c(DIST, dist)

    if(t == (1 + lag_search)){
      for(v in 1:P){
        rank_REF_for_tstar_in_v = TS_rank_REF[(tstar - lag_search) : tstar, v]
        BC_corrdep[(t-lag_search) : t, v] = Sorted_BC[rank_REF_for_tstar_in_v, v]
      }
      Visited_t[(tstar-lag_search) : tstar] = Visited_t[(tstar - lag_search) : tstar] +1
    }

    if(t > (1 + lag_search)){
      for(v in 1:P){
        rank_REF_for_tstar_in_v = TS_rank_REF[(tstar - lag_keep) : tstar, v]
        BC_corrdep[(t-lag_keep) : t, v] = Sorted_BC[rank_REF_for_tstar_in_v, v]
      }
      Visited_t[(tstar-lag_keep) : tstar] = Visited_t[(tstar -lag_keep) : tstar] +1
    }

    t = t + lag_keep + 1
  }

  ###################
  if( (t - 1 - lag_keep) < Ntimes_BC){
    Diff_lag = Ntimes_BC - (t - 1 - lag_keep)
    MultiVec_rank_BC_lag_t = TS_rank_BC_icond[(Ntimes_BC - lag_search) : Ntimes_BC,, drop = FALSE]

    FF = find_analogue_time(TS_rank_REF_icond, MultiVec_rank_BC_lag_t)
    tstar = FF$tstar
    dist = FF$dist
    TSTAR = c(TSTAR, tstar)
    DIST = c(DIST, dist)

    for(v in 1:P){
      rank_REF_for_tstar_in_v = TS_rank_REF[(tstar - Diff_lag + 1) : tstar, v]
      BC_corrdep[(Ntimes_BC - Diff_lag + 1) : Ntimes_BC, v] = Sorted_BC[rank_REF_for_tstar_in_v, v]
    }
    Visited_t[(tstar-Diff_lag+1):tstar] = Visited_t[(tstar-Diff_lag+1):tstar] + 1
  }

  ##########################

  return(
    list(
      BC_corrdep = BC_corrdep,
      Visited_t = Visited_t,
      TSTAR = TSTAR,
      DIST = DIST
    )
  )
}
############# end function R2D2asb


###############################################
# This function is called by R2D2asb

find_analogue_time = function(TS_rank_REF_icond, MultiVec_rank_BC_lag_t){

  LAG = nrow(MultiVec_rank_BC_lag_t) - 1
  distmin = +Inf
  tstar = NaN
  start = 1 + LAG

  for(t in start : nrow(TS_rank_REF_icond)){

    DISTa = dist(
      rbind(
        as.vector(MultiVec_rank_BC_lag_t),
        as.vector(TS_rank_REF_icond[(t-LAG) : t,])
      )
    ) # Euclidean distance

    if(DISTa < distmin){
      distmin = DISTa
      tstar = t
    }
  }

  return(
    list(
      tstar = tstar,
      dist = distmin
    )
  )
}
#### end function find_analogue_time



