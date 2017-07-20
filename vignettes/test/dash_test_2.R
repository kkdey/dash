
#devtools::install_github("kkdey/singleCellRNASeqHumanLengESC", force=TRUE)

library(singleCellRNASeqMouseDeng2014)
deng.counts <- Biobase::exprs(Deng2014MouseESC)
deng.meta_data <- Biobase::pData(Deng2014MouseESC)
deng.gene_names <- rownames(deng.counts)


comp_data <- t(deng.counts)[1:4,]+1
system.time(out <- dash(comp_data = comp_data,
            optmethod = "mixEM",
            mode = colMeans(comp_data),
            def_positions = list("center" = Inf, "null" = 1, "corner" = 1),
            concentration = c(Inf, 100, 50, 20, 10, 5, 2),
            weight = list("center" = 100, "null" = 1, "corner" = 1),
            bf = TRUE,
            verbose=TRUE))




x1 <- rep(1, 10000)
x1[floor(seq(1, 10000, length.out = 5000))] <- 0
x1[10] <- 100
x2 <- rep(1, 10000)

mat <- rbind(x1, x2)

system.time(out <- dash(comp_data = mat,
                        optmethod = "mixEM",
                        mode = x2,
                        def_positions = list("center" = Inf, "null" = 1, "corner" = 1),
                        concentration = c(Inf, 100, 50, 20, 10, 5, 2, 1),
                        weight = list("center" = 100, "null" = 1, "corner" = 1),
                        bf = TRUE,
                        verbose=TRUE))



comp_data <- t(deng.counts)
concentration = NULL
mode = colMeans(t(deng.counts))
weight = list("center" = 10, "null" = 10, "corner" = 10)
def_positions = list("center" = Inf, "null" = 1, "corner" = 0.5)
optmethod = "mixEM"
verbose = FALSE
fdr_bound = 50
sample_weights = NULL
pi_init = NULL
control = list()
bf <- TRUE
