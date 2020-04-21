
# Contains 2 functions:
# - R2D2_asb : performing R2D2 with analogs-based conditioning with semi-block approach
# - Find_t_in_MULTIVEC_rank_REF_with_closest_Vec_rank_REF_lag_with_dist : called by R2D2asb

#####################################################

R2D2_asb = function(REF, BC,
                    iref=c(1), lag_all=0, keep_lag=0){

  # Looks for blocks of size lag_all+1
  # Keeps only the keep_lag+1 last values of the block

  if(lag_all<keep_lag){
    stop("lag_all < keep_lag => Stop\n")
  }

  Niref = length(iref)

  # dim(REF) Ntimes x NVar
  P = dim(BC)[2] # number of var (including iref !!)
  Ntimes_BC = dim(BC)[1]
  Ntimes_Ref = dim(REF)[1]

  BC_corrdep = array(NaN, dim=c(Ntimes_BC,P)) # (N days x p var)


  TS_rank_BC = array(NaN, dim=dim(BC))    # N x p
  Max_rank_BC_period = dim(TS_rank_BC)[1]

  TS_rank_REF = array(NaN, dim=dim(REF))  # N x p
  Max_rank_REF = dim(TS_rank_REF)[1]

  for(v in 1:P){
    TS_rank_BC[,v] = rank(BC[,v], ties.method = "min")
    TS_rank_REF[,v] = rank(REF[,v], ties.method = "min")
  }

  #############
  # SHRINK THE RANK VECTORS SHRINK THE RANK VECTORS SHRINK THE RANK VECTORS SHRINK THE RANK VECTORS SHRINK THE RANK VECTORS
  if(Ntimes_BC < Ntimes_Ref){
    TS_rank_REF = round(TS_rank_REF * (Max_rank_BC_period / Max_rank_REF))
    cat("SHRINK 1\n")
  }
  if(Ntimes_BC > Ntimes_Ref){
    TS_rank_BC = round(TS_rank_BC * (Max_rank_REF / Max_rank_BC_period))
    cat("SHRINK 2\n")
  }
  # SHRINK THE RANK VECTORS SHRINK THE RANK VECTORS SHRINK THE RANK VECTORS SHRINK THE RANK VECTORS SHRINK THE RANK VECTORS
  #############


  Sorted_BC = array(NaN, dim=c(Ntimes_BC,P))
  for(v in 1:P){
    Sorted_BC[,v] = sort(BC[,v])
  }

  #####################

  TS_rank_REF_iref = array(NaN, dim=c(dim(TS_rank_REF)[1], Niref)) # 1:Ntimes x 800 for 400_2_0
  for(r in 1:Niref){
    TS_rank_REF_iref[,r] = TS_rank_REF[,iref[r]]
  }

  TS_rank_BC_iref = array(NaN, dim=c(dim(TS_rank_BC)[1], Niref))
  for(r in 1:Niref){
    TS_rank_BC_iref[,r] = TS_rank_BC[,iref[r]]
  }

  ##########################

  #Visited_t = array(0, dim=Ntimes_BC)
  Visited_t = numeric(length = Ntimes_BC)

  TSTAR = c() ; DIST = c()

  t = 1+lag_all


  while(t<=Ntimes_BC){
    MultiVec_rank_BC_lag_t = array(NaN, dim=c((lag_all+1),Niref))
    for(r in 1:Niref){
      MultiVec_rank_BC_lag_t[,r] = TS_rank_BC_iref[(t-lag_all):t,r]   # t x r = 1 x 800 for 400_2_0
    }
    FF = Find_t_in_MULTIVEC_rank_REF_with_closest_Vec_rank_REF_lag_with_dist(TS_rank_REF_iref, MultiVec_rank_BC_lag_t)
    tstar = FF$tstar
    dist = FF$dist
    TSTAR = c(TSTAR,tstar)
    DIST = c(DIST, dist)

    if(t==(1+lag_all)){
      for(v in 1:P){
        rank_REF_for_tstar_in_v = TS_rank_REF[(tstar-lag_all):tstar,v]
        BC_corrdep[(t-lag_all):t,v] = Sorted_BC[rank_REF_for_tstar_in_v,v]
      }
      Visited_t[(tstar-lag_all):tstar] = Visited_t[(tstar-lag_all):tstar] +1
    }

    if(t>(1+lag_all)){
      for(v in 1:P){
        rank_REF_for_tstar_in_v = TS_rank_REF[(tstar-keep_lag):tstar,v]
        BC_corrdep[(t-keep_lag):t,v] = Sorted_BC[rank_REF_for_tstar_in_v,v]
      }
      Visited_t[(tstar-keep_lag):tstar] = Visited_t[(tstar-keep_lag):tstar] +1
    }

    t = t + keep_lag + 1
  }

  ###################
  if((t-1-keep_lag)<Ntimes_BC){
    Diff_lag = Ntimes_BC - (t-1-keep_lag)
    MultiVec_rank_BC_lag_t = array(NaN, dim=c(length(c((Ntimes_BC-lag_all):Ntimes_BC)) ,Niref))
    for(r in 1:Niref){
      MultiVec_rank_BC_lag_t[,r] = TS_rank_BC_iref[(Ntimes_BC-lag_all):Ntimes_BC,r]
    }

    FF = Find_t_in_MULTIVEC_rank_REF_with_closest_Vec_rank_REF_lag_with_dist(TS_rank_REF_iref, MultiVec_rank_BC_lag_t)
    tstar = FF$tstar
    dist = FF$dist
    TSTAR = c(TSTAR,tstar)
    DIST = c(DIST, dist)

    for(v in 1:P){
      rank_REF_for_tstar_in_v = TS_rank_REF[(tstar-Diff_lag+1):tstar,v]
      BC_corrdep[(Ntimes_BC-Diff_lag+1):Ntimes_BC,v] = Sorted_BC[rank_REF_for_tstar_in_v,v]
    }
    Visited_t[(tstar-Diff_lag+1):tstar] = Visited_t[(tstar-Diff_lag+1):tstar] + 1
  }

  ##########################

  return(list(BC_corrdep=BC_corrdep, Visited_t=Visited_t, TSTAR=TSTAR, DIST=DIST))
}
############# end function R2D2asb


###############################################
# This function is called by R2D2asb

Find_t_in_MULTIVEC_rank_REF_with_closest_Vec_rank_REF_lag_with_dist = function(TS_rank_REF_iref, MultiVec_rank_BC_lag_t){

  LAG = dim(MultiVec_rank_BC_lag_t)[1] - 1
  distmin = +Inf
  tstar = NaN
  start = 1 + LAG

  for(t in start:(dim(TS_rank_REF_iref)[1])){

    DISTa = dist(rbind(as.vector(MultiVec_rank_BC_lag_t), as.vector(TS_rank_REF_iref[(t-LAG):t,]))) # Euclidean distance

    if(DISTa<distmin){
      distmin = DISTa
      tstar = t
    }
  }

  return(list(tstar=tstar, dist=distmin))
}
#### end function Find_t_in_MULTIVEC_rank_REF_with_closest_Vec_rank_REF_lag_with_dist



