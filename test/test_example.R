devtools::load_all()
data("r2d2_example")
system.time(r2d2_vrac <-
  with(r2d2_example, R2D2_asb(REF = refdata, BC = bc1d)))
system.time(r2d2_sthao <-
  with(r2d2_example, R2D2(refdata = refdata, bc1d = bc1d)))
all.equal(r2d2_vrac, r2d2_sthao)

system.time(r2d2_vrac <-
              with(r2d2_example, R2D2_asb(
                REF = refdata, BC = bc1d, iref = icond
              )))
system.time(r2d2_sthao <-
              with(r2d2_example, R2D2(
                refdata = refdata, bc1d = bc1d, icond = icond
              )))
all.equal(r2d2_vrac, r2d2_sthao)

system.time(r2d2_vrac <-
              with(r2d2_example, R2D2_asb(
                REF = refdata, BC = bc1d, iref = icond,
                lag_all = 18, keep_lag = 18
              )))
system.time(r2d2_sthao <-
              with(r2d2_example, R2D2(
                refdata = refdata, bc1d = bc1d, icond = icond,
                lag_search = 18, lag_keep = 18
              )))
all.equal(r2d2_vrac, r2d2_sthao)

system.time(r2d2_vrac <-
              with(r2d2_example, R2D2_asb(
                REF = refdata, BC = bc1d, iref = icond,
                lag_all = 8, keep_lag = 6
              )))
system.time(r2d2_sthao <-
              with(r2d2_example, R2D2(
                refdata = refdata, bc1d = bc1d, icond = icond,
                lag_search = 8, lag_keep = 6
              )))
all.equal(r2d2_vrac, r2d2_sthao)

# Reproducing the example provided in Vrac(2018)

refdata <- matrix(c(0.3, 0.5, 0.9, 0.8,
                    1.1, 1.7, 1.2, 1.9,
                    2.1, 1.8, 3.0, 2.7), ncol = 3, nrow = 4)

bc1d <- matrix(c(0.7, 0.5, 0.2, 0.9,
                 1.3, 1.8, 1.1, 1.4,
                 1.9, 2.9, 2.0, 2.6), ncol = 3, nrow = 4)

R2D2(refdata = refdata,
     bc1d = bc1d,
     icond = 1)

R2D2(refdata = refdata,
     bc1d = bc1d,
     icond = 2)

R2D2(refdata = refdata,
     bc1d = bc1d,
     icond = 3)


refdata <- matrix(c(0.3, 0.5, 0.9, 0.8,
                    1.1, 1.7, 1.2, 1.9,
                    2.1, 1.8, 3.0, 2.7), ncol = 3, nrow = 4)
bc1d <- matrix(c(0.7, 0.5, 0.2, 0.9,
                 1.3, 1.8, 1.1, 1.4,
                 1.9, 2.9, 2.0, 2.6), ncol = 3, nrow = 4)

ranks_refdata <- apply(refdata, 2, rank)
ranks_bc1d <- apply(bc1d, 2, rank)

# 1 conditioning dimension, 0 lag
icond <- 1
block_conddim_BC <- ranks_bc1d[1, icond, drop = FALSE]
conddim_REF <- ranks_refdata[, icond, drop = FALSE]
find_bestanalogue_time(conddim_REF, block_conddim_BC)

# 1 conditioning dimension, 1 lag
icond <- 1
block_conddim_BC <- ranks_bc1d[1:2, icond, drop = FALSE]
conddim_REF <- ranks_refdata[, icond, drop = FALSE]
find_bestanalogue_time(conddim_REF, block_conddim_BC)

# 2 conditioning dimension, 1 lag
icond <- 1:2
block_conddim_BC <- ranks_bc1d[1:2, icond, drop = FALSE]
conddim_REF <- ranks_refdata[, icond, drop = FALSE]
find_bestanalogue_time(conddim_REF, block_conddim_BC)