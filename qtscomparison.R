get_rowcol_index = function(i, j, n){
    idx = i*n-i*(i-1)/2+j+1  # + 1 for 1-based indexing in R
    return(idx)
}

n = 10
numStates = (n+1)*(n+2)/2
s3idx = numeric(numStates)
s2idx = numeric(numStates)
for(i in 0:n){
    for(j in i:n){
        idx = get_rowcol_index(i, j, n)
        s3idx[idx] = i
        s2idx[idx] = j
    }
}
states = data.frame(i = s3idx, j = s2idx, ijidx = 1:numStates)

dat.tsmc = as.matrix(read.csv("qtsprobs",header=F))
dat.mma = as.matrix(read.csv("mathematica_qtsprobs.csv",header=F))

hist(c(dat.tsmc)-c(dat.mma))

summary(c(dat.tsmc)-c(dat.mma))
states[which(abs(dat.tsmc-dat.mma) > 0.1, arr.ind=T),c(1,2)]
states
dim(which(abs(dat.tsmc-dat.mma) > 0.1, arr.ind=T))

badIdxs = which((dat.tsmc-dat.mma) > 0.1, arr.ind=T)
badi = sapply(badIdxs[,1], function(x) states[x,1])
badj = sapply(badIdxs[,1], function(x) states[x,2])
badk = sapply(badIdxs[,2], function(x) states[x,1])
badl = sapply(badIdxs[,2], function(x) states[x,2])
badOnes = data.frame(i=badi, j = badj, k = badk, l = badl)
badOnes
