

#############  dash_test hocomoco examplee  ##################

x1 <- rbind(c(1,	0, 	0, 	0, 	0, 	0, 	2, 	0, 	0, 	0, 	1, 	0, 	0, 	1),
            c(2,	0,	3,	5,	3,	0,	1,	0,	0,	4,	1,	4,	1,	1),
            c(1,	5,	1,	0,	0,	5,	2,	5,	5,	1,	2,	1,	4,	0),
            c(1,	0,	1,	0,	2,	0,	0,	0,	0,	0,	1,	0,	0,	3))

xmat <- t(x1)

comp_data <- xmat
concentration = NULL
mode = NULL
weight = list("center" = 10, "null" = 1, "corner" = 1)
def_positions = list("center" = Inf, "null" = 1, "corner" = 0.5)
optmethod = c("mixEM", "w_mixEM", "mixIP")
sample_weights = NULL
verbose = FALSE
fdr_bound = 50
bf = TRUE
pi_init = NULL
control = list()
reportcov = FALSE
