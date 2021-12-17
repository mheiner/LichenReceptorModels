
library("Matrix")

#### identifiability conditions in Park, E. S., Spiegelman, C. H., and Henry, R. C. (2002), “Bilinear Estimation of Pollution Source Profiles and Amounts by Using Multivariate Receptor Models,” Environmetrics, 13, 775–798. [68]
# 1. At least K-1 zeros in each column of Lam
# 2. Rank of each sub-matrix at least K-1
# 3. ColSums of Lam = 1
####

Z = Lam_all # from model spec
Z = Lam_orig
dim(Z)
(K = ncol(Z))
Z[Z==1e-9] = 0
Z

for ( k in 1:K ) {
  zindx = which(Z[,k] == 0)
  if (length(zindx) > 0) {
    rankinfo = rankMatrix(Z[zindx, -k])
    cat("k =", k, "\n")
    print(rankinfo)
    cat("Rank of submatrix is", ifelse(rankinfo >= (K-1), "INDEED", "NOT"), " sufficient \n\n")
  }
}

library("lattice")
levelplot(1*(Z > 0.0))

