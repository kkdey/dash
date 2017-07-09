
#devtools::install_github("kkdey/singleCellRNASeqHumanLengESC", force=TRUE)

library(singleCellRNASeqMouseDeng2014)
deng.counts <- Biobase::exprs(Deng2014MouseESC)
deng.meta_data <- Biobase::pData(Deng2014MouseESC)
deng.gene_names <- rownames(deng.counts)


comp_data <- t(deng.counts)
system.time(out <- dash(comp_data = t(deng.counts),
            optmethod = "mixEM",
            mode = colMeans(comp_data),
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
