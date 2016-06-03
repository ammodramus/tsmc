dat = read.table("pistest",header=F)
names(dat) = c("i","j", "ti", "tip1", "tj", "tjp1", "prob")

numReps = 10000000
coalt3 = rexp(numReps, rate = 3)
coalt2 = coalt3 + rexp(numReps, rate = 1)

counts = rep(0, nrow(dat))
for(i in 1:nrow(dat)){
    ti = dat$ti[i]
    tip1 = dat$tip1[i]
    tj = dat$tj[i]
    tjp1 = dat$tjp1[i]
    counts[i] = sum(coalt3 >= ti & coalt3 < tip1 & coalt2 >= tj & coalt2 < tjp1)
}

plot(dat$prob, counts/sum(counts))
#text(dat$prob, counts/sum(counts), labels = sapply(1:nrow(dat), function(x) sprintf("(%i,%i)", dat$i[x], dat$j[x])))
abline(0,1)
# spot on

dat2 = read.table("testexpectations",header=F)
names(dat2) = c("i","j", "ti", "tip1", "tj", "tjp1", "et3", "et2")
head(dat2)
dat2$et3
dat2$et2

meanst3 = numeric(nrow(dat2))
meanst2 = numeric(nrow(dat2))
counts = rep(0, nrow(dat))
for(i in 1:nrow(dat)){
    ti = dat$ti[i]
    tip1 = dat$tip1[i]
    tj = dat$tj[i]
    tjp1 = dat$tjp1[i]
    indices = which(coalt3 >= ti & coalt3 < tip1 & coalt2 >= tj & coalt2 < tjp1)
    meanst3[i] = mean(coalt3[indices])
    meanst2[i] = mean(coalt2[indices])
}

# again, spot-on now
plot(dat2$et3, meanst3)
#text(dat2$et3, meanst3, labels = sapply(1:nrow(dat), function(x) sprintf("(%i,%i)", dat$i[x], dat$j[x])))
abline(0,1)
plot(dat2$et2, meanst2)
abline(0,1)
# both spot-on again (after correcting rate 2 -> rate 1)

# simulations where lambda = 1.0 for [0,0.45198)
# and 2.0 for [0.45198,infinity)

changepoint = 0.451985

numReps = 1000000
coalt3.2 = rexp(numReps, rate = 3)
coalt2.2 = coalt3.2 + rexp(numReps, rate = 1)

bothlessindices = which(coalt3.2 < changepoint & coalt2.2 < changepoint)
firstlessindices = which(coalt3.2 < changepoint & coalt2.2 >= changepoint)
neitherlessindices = which(coalt3.2 >= changepoint & coalt2.2 >= changepoint)

numFirstLess = length(firstlessindices)
coalt2.2[firstlessindices] = changepoint + rexp(numFirstLess, rate = 1/2.0)
numNeitherless = length(neitherlessindices)
coalt3.2[neitherlessindices] = changepoint + rexp(numNeitherless, rate = 3/2.0)
coalt2.2[neitherlessindices] = coalt3.2[neitherlessindices] + 
    rexp(numNeitherless, rate = 1/2.0)

dat3 = read.table("testpis2",header=F)
names(dat3) = c("i","j", "ti", "tip1", "tj", "tjp1", "prob")

counts3 = rep(0, nrow(dat3))
for(i in 1:nrow(dat)){
    ti = dat3$ti[i]
    tip1 = dat3$tip1[i]
    tj = dat3$tj[i]
    tjp1 = dat3$tjp1[i]
    counts3[i] = sum(coalt3.2 >= ti & coalt3.2 < tip1 & coalt2.2 >= tj & coalt2.2 < tjp1)
}

plot(dat3$prob, counts3/sum(counts3))
abline(0,1)
# spot on (?)

dat4 = read.table("testexpectations2",header=F)
names(dat4) = c("i","j", "ti", "tip1", "tj", "tjp1", "et3", "et2")
dat4$et3
dat4$et2

meanst3 = numeric(nrow(dat4))
meanst2 = numeric(nrow(dat4))
counts = rep(0, nrow(dat4))
for(i in 1:nrow(dat4)){
    ti = dat4$ti[i]
    tip1 = dat4$tip1[i]
    tj = dat4$tj[i]
    tjp1 = dat4$tjp1[i]
    indices = which(coalt3.2 >= ti & coalt3.2 < tip1 & coalt2.2 >= tj & coalt2.2 < tjp1)
    meanst3[i] = mean(coalt3.2[indices])
    meanst2[i] = mean(coalt2.2[indices])
}

plot(dat4$et3, meanst3)
text(dat4$et3, meanst3, labels = sapply(1:nrow(dat), function(x) sprintf("(%i,%i)", dat4$i[x], dat4$j[x])))
abline(0,1)
# spot-on
plot(dat4$et2, meanst2)
text(dat4$et2, meanst2, labels = sapply(1:nrow(dat), function(x) sprintf("(%i,%i)", dat4$i[x], dat4$j[x])))
abline(0,1)
# spot-on again

######################################################
# Now extracting qts transitions from ms simulations #
######################################################

sims = as.matrix(read.csv("qtssims.txt",header=F))
probs = as.matrix(read.csv("qtsprobs.txt",header=F))

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

states = data.frame(i = s3idx, j = s2idx, ijidx = 1:66)

rowSums(probs)
which(rowSums(probs) >= 1)
states[which(rowSums(probs)>=1),c(1,2)]

plot(c(sims), c(probs), asp = 1)
abline(0,1)

par(mfrow=c(1,2))
hist(rowSums(probs))
hist(rowSums(sims))


