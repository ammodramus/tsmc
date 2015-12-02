n = 8
numStates = (n+1)*(n+2)/2

row_index = function(i, M){
    ii = M(M+1)/2-1-i
    K = floor((sqrt(8*ii+1)-1)/2)
    return(M-1-K)
}

column_index = function(i, M){
    ii = M(M+1)/2-1-i
    K = floor((sqrt(8*ii+1)-1)/2)
    return(i - M(M+1)/2 + (K+1)(K+2)/2)
}

get_idx = function(i,j, M){
    return(i*M-i*(i-1)/2+j+1) # 1-based indexing
}

plot_done = function(done){
    plot(0:nrow(done), 0:ncol(done), type = "n", xaxt = "n", yaxt = "n", xlab = "", ylab = "", asp = 1)
    for(i in 1:nrow(done)){
        for(j in 1:ncol(done)){
            rect(i,j,i+1,j+1, col = done[i,j])
        }
    }
}
    
done = matrix(0, numStates, numStates)

indices = matrix(,numStates,2)
for(i in 0:n){
    for(j in i:n){
        idx = get_idx(i,j,n)
        indices[idx,1] = i
        indices[idx,2] = j
    }
}

for(i in 0:n){
    for(j in i:n){
        for(k in 0:n){    
            for(l in k:n){
                rowIdx = get_idx(i,j,n)
                colIdx = get_idx(k,l,n)
                if((i!=k)&&(j!=l))
                    done[rowIdx,colIdx] = 1 # (this may be overridden by another below)
                if(i==k&&j==l)
                    done[rowIdx,colIdx] = 1
                if(i==k&&k<=j&&j<l)
                    done[rowIdx,colIdx] = 2
                if(i==k&&k<l&&l<j)
                    done[rowIdx,colIdx] = 3
                if(k<i&&i==l&&l<j)
                    done[rowIdx,colIdx] = 4
                if(k<i&&i<=j&&j==l)
                    done[rowIdx,colIdx] = 5
                if(i<k&&k<j&&j==l)
                    done[rowIdx,colIdx] = 6
                if(i<j&&j==k&&k<l)
                    done[rowIdx,colIdx] = 7
                # special
                if(i==k&&k==l&&l<j)
                    done[rowIdx,colIdx] = 8
                if(i<k&&k==l&&l==j)
                    done[rowIdx,colIdx] = 9
            }
        }
    }
}

plot_done(done)

which(done == 0, arr.ind=T)
# (none)
