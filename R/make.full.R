.make.full <-
function(X) {
        svd.X <- svd(X)
        r <- max(which(svd.X$d > 1e-08))
        return(as.matrix(svd.X$u[, 1:r]))
          }
