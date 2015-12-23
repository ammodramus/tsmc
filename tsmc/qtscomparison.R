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
qtsprobs = as.matrix(read.csv("qtstestout",header=F))
hist(abs(c(ms.tallies) - c(qtsprobs)), breaks = 100)
ms.tallies[1:10,1:5]
qtsprobs[1:10,1:5]

#cutoff = 0.08
#mstal = c(ms.tallies)[which(abs(c(ms.tallies) - c(qtsprobs)) > cutoff)]
#qtsprob = c(qtsprobs)[which(abs(c(ms.tallies) - c(qtsprobs)) > cutoff)]
#which(abs(ms.tallies - qtsprobs) > cutoff, arr.ind=T)
#badIjs = which(abs(ms.tallies - qtsprobs) > cutoff, arr.ind=T)[,1]
#badKls = which(abs(ms.tallies - qtsprobs) > cutoff, arr.ind=T)[,2]
#badIs = sapply(badIjs, function(x) states$i[states$ijidx==x])
#badJs = sapply(badIjs, function(x) states$j[states$ijidx==x])
#badKs = sapply(badKls, function(x) states$i[states$ijidx==x])
#badLs = sapply(badKls, function(x) states$j[states$ijidx==x])

#bad.dat = data.frame(ms = mstal, qts = qtsprob, ij = badIjs, kl = badKls,
#                     i = badIs, j = badJs, k = badKs, l = badLs)
#bad.dat

#summary(abs(c(ms.tallies) - c(qtsprobs)), breaks = 100)

# one explanation:
# sites hit twice by recombination

#hist(rowSums(qtsprobs)-rowSums(ms.tallies))

