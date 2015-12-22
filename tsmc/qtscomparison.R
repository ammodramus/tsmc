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
dat.int = as.matrix(read.csv("mathematica_intqtsprobs.csv",header=F))

head(dat.tsmc)
entry.tsmc = function(i,j,k,l){
    ijidx = get_rowcol_index(i,j,10)
    klidx = get_rowcol_index(k,l,10)
    return(dat.tsmc[ijidx,klidx])
}
entry.mma = function(i,j,k,l){
    ijidx = get_rowcol_index(i,j,10)
    klidx = get_rowcol_index(k,l,10)
    return(dat.mma[ijidx,klidx])
}

diffs = abs(c(dat.tsmc)-c(dat.mma))
hist(abs(c(dat.tsmc)-c(dat.mma)))

# all working!
summary(abs(c(dat.tsmc)-c(dat.mma)))
summary(abs(c(dat.tsmc)-c(dat.int)))

################################
# working with ms transitions
################################

ms.tallies = as.matrix(read.csv("qtstallies",header=F))
ms.tallies2 = as.matrix(read.csv("qtstallies_mspms_2",header=F))
qtsprobs = as.matrix(read.csv("qtstestout",header=F))
dim(ms.tallies)
hist(abs(c(ms.tallies) - c(qtsprobs)), breaks = 100)
ms.tallies[1:4,1:4]
qtsprobs[1:4,1:4]

sum((c(ms.tallies) > 0) & (c(qtsprobs) == 0))

which((ms.tallies > 0) & (qtsprobs == 0), arr.ind = T)

which((ms.tallies2 > 0) & (qtsprobs == 0), arr.ind = T)

states[states$ijidx %in% c(46,30),]

# one explanation:
# sites hit twice by recombination
