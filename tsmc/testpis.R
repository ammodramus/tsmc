dat = read.table("pistest",header=F)
names(dat) = c("i","j", "ti", "tip1", "tj", "tjp1", "prob")

numReps = 10000000
coalt3 = rexp(numReps, rate = 3)
coalt2 = coalt3 + rexp(numReps, rate = 2)

counts = rep(0, nrow(dat))
for(i in 1:nrow(dat)){
    ti = dat$ti[i]
    tip1 = dat$tip1[i]
    tj = dat$tj[i]
    tjp1 = dat$tjp1[i]
    counts[i] = sum(coalt3 >= ti & coalt3 < tip1 & coalt2 >= tj & coalt2 < tjp1)
}

plot(dat$prob, counts/sum(counts))
# spot on
