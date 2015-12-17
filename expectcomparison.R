dat.c = read.table("expectout",header=F)
dat.mma = read.csv("expect_mma.csv",header=F)
dat.c

which(as.matrix(dat.c) == max(as.matrix(dat.c)), arr.ind=T)

dat.c[21,21]

dat.mma[21,21]

hist(abs(c(as.matrix(dat.c)-as.matrix(dat.mma))))
which(abs(as.matrix(dat.c)-as.matrix(dat.mma)) > 1e-6, arr.ind = T)

# one entry that is off by 1e-6, others at 1e-7 or less
# (index is 21, 21, which corresponds to (i,j) = (1,10), (k,l) = (1,10)
# hmm...

# test likelihood, see if matches
