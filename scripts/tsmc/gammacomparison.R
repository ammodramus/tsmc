dat.c = read.table("gammaout",header=F)
dat.mma = read.csv("gamma_mma.csv",header=F)

max(abs(c(as.matrix(dat.c)-as.matrix(dat.mma))))
# [1] 3.736536e-07
# working
