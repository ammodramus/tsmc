dat.c = read.table("backwardout",header=F)
dat.mma = read.csv("backward_mma.csv",header=F)

max(abs(c(as.matrix(dat.c)-as.matrix(dat.mma))))
# [1] 3.633357e-07
# okay, good
