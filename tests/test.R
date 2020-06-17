library(R2D2)

mat <- matrix(
  rep(1:10, 2),
  ncol  = 2 
)
print(mat)

mat2 <- matrix(
  rep(1:15, 2),
  ncol  = 2 
)
print(mat2)

r2d2_bc = r2d2(mat, mat, 1, 0, 0)
stopifnot(
  all.equal(mat, r2d2_bc[[4]])
)
stopifnot(
  all.equal(r2d2_bc[[5]], rep(1, nrow(mat)))
)

r2d2_bc = r2d2(mat, mat, 2, 0, 0)
stopifnot(
  all.equal(mat, r2d2_bc[[4]])
)
stopifnot(
  all.equal(r2d2_bc[[5]], rep(1, nrow(mat)))
)

r2d2_bc = r2d2(mat, mat2, 1, 0, 0)
stopifnot(
  all.equal(
    matrix(rep(r2d2_bc[[3]], rep(r2d2_bc[[5]], 2)), ncol = 2),
    r2d2_bc[[4]]
  )
)

r2d2_bc = r2d2(mat2, mat, 1, 0, 0)
stopifnot(
  all.equal(mat, r2d2_bc[[4]])
)

