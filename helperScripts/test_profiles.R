
library("Matrix")

rankcheck = function(Z, tol=1e-10) {
  #### identifiability conditions in Park, E. S., Spiegelman, C. H., and Henry, R. C. (2002), “Bilinear Estimation of Pollution Source Profiles and Amounts by Using Multivariate Receptor Models,” Environmetrics, 13, 775–798. [68]
  # 1. At least K-1 zeros in each column of Lam
  # 2. Rank of each sub-matrix at least K-1
  # 3. ColSums of Lam = 1
  ####
  
  (K = ncol(Z))

  for ( k in 1:K ) {
    zindx = which(Z[,k] <= tol)
    if (length(zindx) > 0 && K > 2) {
      rankinfo = rankMatrix(Z[zindx, -k])
      cat("k =", k, "\n")
      cat("number of zero entries: ", length(zindx), "\n")
      print(rankinfo)
      cat("Rank of submatrix is", ifelse(rankinfo >= (K-1), "INDEED", "NOT"), "sufficient \n\n")
    } else if (length(zindx) > 0 && K == 2) {
      cat("Truth K=2. Other column entries for column", k, "zeros: ", (Z[zindx, -k, drop=FALSE]), "\n" )
    }
  }
  
  return(NULL)
}
