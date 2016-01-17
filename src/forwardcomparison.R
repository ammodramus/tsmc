dat.c = read.table("forwardout",header=F)
dim(dat.c)

dat.mma = read.csv("forward_mma.csv",header=F)

max(abs(c(as.matrix(dat.c)-as.matrix(dat.mma))))
# [1] 3.739775e-07
# good
