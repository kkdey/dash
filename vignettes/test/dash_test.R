

##############  dash test  ######################

xmat <- rbind(c(5, 0, 2, 0),
              c(1, 1, 0, 1),
              c(100, 100, 50, 100),
              c(20, 50, 100, 10),
              c(10, 10, 200, 20),
              c(50, 54, 58, 53),
              c(1,1,1,3),
              c(2, 4, 1, 1))


out <- dash(xmat, optmethod = "w_mixEM", verbose=TRUE)

xmat2 <- t(apply(xmat, 1, function(x) return(x/sum(x))))

plot(out$posterior_weights[1,], pch=20, col="black", cex=1.5,
     ylab = "posterior weight", xlab = "component")

out$posmean
xmat2





comp_data <- xmat
concentration = NULL
mode = NULL
weight = list("center" = 10, "null" = 10, "corner" = 10)
def_positions = list("center" = Inf, "null" = 1, "corner" = 0.5)
optmethod = "mixEM"
verbose = FALSE
fdr_bound = 50
sample_weights = NULL
pi_init = NULL
control = list()
