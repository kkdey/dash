

##################  dirichlet likelihood matrix   ######################

logfac = function(x){
  if(x < 1 && x > 0 || x < 0){
    stop("x cannot be less than 0 or a fraction")
  }else if(x < 10){
    out <- log(factorial(x))
  }else if (x == 0){
    out <- 0
  }else{
    out <- sum(log(1:x))
  }
  return(out)
}


concentration <- c(10^5, 100, 50, 20, 10, 5, 2, 1, 0.5, 02, 0.1, 0.0001)
usr_mean <- c(4, 2, 2, 1)
usr_mean_prop <- usr_mean/sum(usr_mean)
alpha_mat <- t(sapply(concentration,function(x) return(x*usr_mean_prop)))
prior <- c(10, 5, rep(1, (dim(alpha_mat)[1]-2)))
k <- 1

xmat <- rbind(c(5, 0, 2, 0),
              c(1, 1, 0, 1),
              c(100, 100, 50, 100),
              c(20, 50, 100, 10),
              c(10, 10, 200, 20),
              c(50, 54, 58, 53),
              c(1,1,1,3),
              c(2, 4, 1, 1))

xmat2 <- t(apply(xmat, 1, function(x) return(x/sum(x))))

if(length(usr_mean) != dim(xmat)[2]){
  stop("the length of the dirichlet mode must match with the number of columns in PWM matrix")
}

matrix_lik <- matrix(0, dim(xmat)[1], dim(alpha_mat)[1])

for(n in 1:dim(xmat)[1]){
  for(k in 2:dim(alpha_mat)[1]){
    x <- xmat[n,]
    numero <- sum(x)*beta(sum(alpha_mat[k,]), sum(x))
    lognumero <- log(sum(x)) - LaplacesDemon::ddirichlet(rep(1,2), alpha = c(sum(alpha_mat[k,]), sum(x)), log=TRUE)
    if(numero == 0){
      matrix_lik[n, k] <- 0
    }else{
      index1 <- which(x > 0)
      deno <- prod(x[index1] * beta(alpha_mat[k, index1], x[index1]))
      matrix_lik[n,k] <- numero/deno
    }
  }
  matrix_lik[n,1] <- exp(logfac(n) - sum(sapply(x, function(y) return(logfac(y)))) + sum(x*log(alpha_mat[k,]/sum(alpha_mat[k,]))))
}

ddir <- function (x, alpha, log = FALSE)
{
  if (missing(x))
    stop("x is a required argument.")
  if (missing(alpha))
    stop("alpha is a required argument.")
  if (!is.matrix(x))
    x <- rbind(x)
  if (!is.matrix(alpha))
    alpha <- matrix(alpha, nrow(x), length(alpha), byrow = TRUE)
  if (any(rowSums(x) != 1))
    x/rowSums(x)
  if (any(x < 0))
    stop("x must be non-negative.")
  if (any(alpha <= 0))
    stop("alpha must be positive.")
  dens <- as.vector(lgamma(rowSums(alpha)) - rowSums(lgamma(alpha)) +
                      rowSums((alpha - 1) * log(x)))
  if (log == FALSE)
    dens <- exp(dens)
  return(dens)
}

out <- mixEM(matrix_lik, prior)


##################  the mixture proportions  ########################

out$pihat

##################  posterior weights on grid  ###########################

pi_complete <- rep(1, dim(matrix_lik)[1]) %*% t(out$pihat)
matrix_lik_adj <- matrix_lik*pi_complete
posterior_weights <- t(apply(matrix_lik_adj, 1, function(x) return(x/sum(x))))
posterior_weights

##################  posterior mean  ##################################

posmean_comp <- array(0, c(dim(xmat)[1], dim(xmat)[2], dim(alpha_mat)[1]))
for(n in 1:dim(xmat)[1]){
  for(k in 1:dim(alpha_mat)[1]){
   temp <-  xmat[n,]+ alpha_mat[k,]
   posmean_comp[n,,k] <- temp/sum(temp)
  }
}

posmean <- matrix(0, dim(xmat)[1], dim(xmat)[2])

for(n in 1:dim(xmat)[1]){
  posmean[n,] <- posmean_comp[n,,]%*%posterior_weights[n,]
}

##################  covariance and correlation structure  ##################

poscov_comp <- array(0, c(dim(xmat)[1], dim(xmat)[2], dim(xmat)[2], dim(alpha_mat)[1]))
for(n in 1:dim(xmat)[1]){
  for(k in 1:dim(alpha_mat)[1]){
    temp <-  xmat[n,]+ alpha_mat[k,]
    alpha_0 <- sum(temp)
    posvar <- (alpha_0* temp)/((alpha_0)^2*(alpha_0 + 1))
    poscov_comp[n,,,k] <- diag(posvar) -(temp %*% t(temp))/((alpha_0)^2*(alpha_0 + 1))
  }
}

poscov <- array(0, c(dim(xmat)[2], dim(xmat)[2], dim(xmat)[1]))
poscor <- array(0, c(dim(xmat)[2], dim(xmat)[2], dim(xmat)[1]))

for(n in 1:dim(xmat)[1]){
  dsum <- matrix(0, dim(xmat)[2], dim(xmat)[2])
  for(k in 1:dim(alpha_mat)[1]){
    dsum <- dsum + poscov_comp[n,,,k]*posterior_weights[n,k]
  }
  poscov[,,n] <- dsum
  poscor[,,n] <- cov2cor(dsum)
}

#################   false discovery & enruchment rates  ##############################

lfsr <- posterior_weights[,1]
lfdr <- (posterior_weights[,1] + posterior_weights[,2])

enrichment_prob <- rowSums(posterior_weights[, which(concentration < 1)])

#########   entropy of the posterior distribution   #######################

posentropy_comp <- matrix(0, dim(xmat)[1], dim(alpha_mat)[1])

for(n in 1:dim(xmat)[1]){
  for(k in 1:dim(alpha_mat)[1]){
    temp <- xmat[n,]+ alpha_mat[k,]
    first_term <- -LaplacesDemon::ddirichlet(rep(1,dim(xmat)[2]), alpha = temp, log=TRUE)
    alpha_0 <- sum(temp)
    sec_term <- - (length(temp) - alpha_0)*digamma(alpha_0) - sum((temp-1)*digamma(temp))
    posentropy_comp[n,k] <- first_term + sec_term
  }
}

posentropy <- array(0, dim(xmat)[1])

for(n in 1:dim(xmat)[1]){
  posentropy[n] <- sum(posentropy_comp[n,]*posterior_weights[n,])
}



