
#######  Test bash/dash smoothing on a data   ###################

dat <- get(load("reads_all_1_93297593_93307481.Robj"))
adipose_data <- dat[[1]][[2]]
dim(adipose_data)
plot(adipose_data[1,])

require(Rcpp)
require(inline)

interleave=function(x,y){
  return(as.vector(rbind(x,y)))
}

rshift = function(x){L=length(x); return(c(x[L],x[-L]))}
lshift = function(x){return(c(x[-1],x[1]))}

sig <- adipose_data[1,]

ParentTItable=function(sig){
  n = length(sig)
  J = log2(n)

  # Create decomposition table of signal, using pairwise sums,
  # keeping just the values that are *not* redundant under the
  # shift-invariant scheme.  This is very similar to TI-tables
  # in Donoho and Coifman's TI-denoising framework.
  dmat = matrix(0, nrow=J+1, ncol=n)
  dmat[1,] = sig
  #dmat[1,] = as.matrix(sig)
  dmat2 = matrix(0, nrow=J, ncol=2*n) #the parent table

  for(D in 0:(J-1)){
    nD = 2^(J-D);
    nDo2 = nD/2;
    twonD = 2*nD;
    for(l in 0:(2^D-1)){
      ind = (l*nD+1):((l+1)*nD)
      ind2 = (l*twonD+1):((l+1)*twonD)
      x = dmat[D+1,ind]
      lsumx = x[seq(from=1,to=nD-1, by=2)] + x[seq(from=2,to=nD,by=2)]
      rx = rshift(x);
      rsumx = rx[seq(from=1,to=nD-1, by=2)] + rx[seq(from=2,to=nD,by=2)]
      dmat[D+2,ind] = c(lsumx,rsumx)
      dmat2[D+1,ind2] = c(x,rx)
    }
  }
  return(list(TItable=dmat,parent=dmat2))
}


