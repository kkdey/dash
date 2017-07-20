

###############  dash_test_3  (Himalyan bird abundance)  ###################

library(devtools)

#install_github("kkdey/ecostructure")
library(ecostructure)

data <- get(load(system.file("extdata", "HimalayanBirdsData.rda",
                             package = "ecostructure")))
taxonomic_counts <- t(exprs(data))

rowSums(taxonomic_counts)

taxonomic_counts_1 <- taxonomic_counts +1

m1 <- colMeans(taxonomic_counts_1)
m2 <- colMeans(taxonomic_counts)

system.time(out <- dash(comp_data = taxonomic_counts_1,
                        optmethod = "mixEM",
                        mode = m1,
                        def_positions = list("center" = Inf, "null" = 1, "corner" = 1),
                        concentration = c(Inf, 100, 50, 20, 10, 5, 2, 1),
                        weight = list("center" = 10, "null" = 1, "corner" = 1),
                        bf = TRUE,
                        verbose=TRUE))


