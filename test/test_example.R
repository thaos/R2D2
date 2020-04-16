devtools::load_all()
data("r2D2_example")
r2d2_vrac <- with(r2d2_example, R2D2_asb(REF = ref_cal, BC = bc1d_cal))
r2d2_sthao <- with(r2d2_example, R2D2(REF = ref_cal, BC = bc1d_cal))
all.equal(r2d2_vrac, r2d2_sthao)

r2d2_vrac <- with(r2d2_example, R2D2_asb(REF = ref_cal, BC = bc1d_cal, iref = icond))
r2d2_sthao <- with(r2d2_example, R2D2(REF = ref_cal, BC = bc1d_cal, icond = iconc))
all.equal(r2d2_vrac, r2d2_sthao)