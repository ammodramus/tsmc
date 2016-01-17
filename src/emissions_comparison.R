dat.mma = as.matrix(read.csv("emissions_mma_output.csv",header=F))
dat.c = as.matrix(read.table("emissions_c_output.csv",header=F))

hist(abs(c(dat.mma-dat.c)))
max(abs(c(dat.mma-dat.c)))

# good, checks out
