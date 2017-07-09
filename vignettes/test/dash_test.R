

##############  dash test  ######################

xmat <- rbind(c(5, 0, 2, 0),
              c(1, 1, 0, 1),
              c(100, 100, 50, 100),
              c(20, 50, 100, 10),
              c(10, 10, 200, 20),
              c(50, 54, 58, 53),
              c(1,1,1,3),
              c(2, 4, 1, 1))


out <- dash(xmat, optmethod = "mixEM", verbose=TRUE, bf=TRUE)

x1 <- rbind(c(1,	0, 	0, 	0, 	0, 	0, 	2, 	0, 	0, 	0, 	1, 	0, 	0, 	1),
            c(2,	0,	3,	5,	3,	0,	1,	0,	0,	4,	1,	4,	1,	1),
            c(1,	5,	1,	0,	0,	5,	2,	5,	5,	1,	2,	1,	4,	0),
            c(1,	0,	1,	0,	2,	0,	0,	0,	0,	0,	1,	0,	0,	3))

xmat <- t(x1)
ll <- dash(xmat, optmethod = "mixEM", bf=FALSE)




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

LaplacesDemon::ddirichlet(rep(1,2), alpha = c(sum(alpha_mat[k,]), sum(x)), log=TRUE)
