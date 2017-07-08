
##########  dirichlet model plots  ###################

install.packages("Compositional")

library(Compositional)

diri.contour(c(0.4, 0.4, 0.2), n=100)

x <- as.matrix( iris[, 1:3] )
x <- x / rowSums(x)

xmat2 <- xmat/rowSums(xmat)
xmat3 <- posmean
par(mfrow=c(1,1))
diri.contour( a = c(1.5, 1.5, 1.5))
diri.contour( a = c(50, 30, 10), x = xmat2[6,])
diri.contour( a = c(50, 30, 10), x = posmean[6,])

diri.contour( a = c(10, 10, 10), x = xmat2[6,])
diri.contour( a = c(10, 10, 10), x = posmean[6,])

diri.contour( a = c(5, 5, 5), x = xmat2)
diri.contour( a = c(5, 5, 5), x = posmean)

diri.contour( a = c(10, 10, 10))
diri.contour( a = c(1, 1, 1))
diri.contour( a = c(0.1, 0.1, 0.1))
diri.contour( a = c(5, 5, 5))

