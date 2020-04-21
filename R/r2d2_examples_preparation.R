#  icond_df <- data.frame(
#    "city" = c("paris", "madrid", "rome", "warsaw", "stockholm"),
#    "lon" = c(2.35, -3.70, 12.49, 21.01, 18.06),
#    "lat" = c(48.85, 40.41, 41.90, 52.22, 59.32)
#  )
#
#  load(file = "data/Example_JAN.RData")
#
#  icond <- integer(nrow(icond_df))
#  for(irow in seq.int(nrow(icond_df))){
#    icond[irow] <- which.min(
#      apply(
#        sweep(
#          cbind("lon" = LON_for_comp, "lat" = LAT_for_comp),
#         MARGIN = 2,
#          STATS = unlist(icond_df[irow, c("lon", "lat")]),
#          FUN = "-"
#        ),
#        MARGIN = 1,
#        FUN = function(x) sum(x^2)
#      )
#    )
#  }
#  icond
#
#  icond_df <- cbind(icond_df, "icond" = icond)
#  names(icond) <- icond_df$city
#
#  r2d2_example <- list(
#    bc1d = BC1d_JAN_79_97,
#   # bc1d_val = BC1d_JAN_98_16,
#    refdata= Ref_JAN_79_97,
#   # ref_val = Ref_JAN_98_16,
#    dates = sprintf("%04d-%02d-%02d", DATES_JAN_79_97$Y, DATES_JAN_79_97$m, DATES_JAN_79_97$d),
#   # dates_val = sprintf("%04d-%02d-%02d", DATES_JAN_98_16$Y, DATES_JAN_98_16$m, DATES_JAN_98_16$d),
#    lon = LON_for_comp,
#    lat = LAT_for_comp,
#    icond = icond
#  )
# # save(r2d2_example, file = "data/r2d2_example.rdata")
# usethis::use_data(r2d2_example, overwrite = TRUE)
#
