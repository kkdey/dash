

###############  dash_test_3  (Himalyan bird abundance)  ###################

library(devtools)

#install_github("kkdey/ecostructure")
library(ecostructure)

data <- get(load(system.file("extdata", "HimalayanBirdsData.rda",
                             package = "ecostructure")))
taxonomic_counts <- t(exprs(data))

rowSums(taxonomic_counts)

taxonomic_counts_1 <- taxonomic_counts +1

m1 <- colMeans(taxonomic_counts)

system.time(out <- dash(comp_data = taxonomic_counts_1,
                        optmethod = "mixEM",
                        mode = m1,
                        def_positions = list("center" = Inf, "null" = 1, "corner" = 1),
                        concentration = c(Inf, 100, 50, 20, 10, 5, 2, 1),
                        weight = list("center" = 10, "null" = 1, "corner" = 1),
                        bf = TRUE,
                        verbose=TRUE))

plot(out$posterior_weights[,1], rowSums(taxonomic_counts))

mod_counts <- round((rowSums(taxonomic_counts) %*% t(rep(1, dim(taxonomic_counts)[2])))*out$posmean)

grid_metadata <- pData(phenoData(data))
head(grid_metadata)


elevation_metadata=grid_metadata$Elevation;
east_west_dir = grid_metadata$WorE;
gom_fit <- CountClust::FitGoM(taxonomic_counts, K=2:4, tol=0.1)

elevation_metadata=grid_metadata$Elevation;
east_west_dir = grid_metadata$WorE;
omega <- gom_fit[[2]]$omega

BlockStructure(omega, blocker_metadata = east_west_dir,
               order_metadata = elevation_metadata,
               yaxis_label = "Elevation",
               levels_decreasing = FALSE)






rep_counts <- do.call(rbind, replicate(10, taxonomic_counts[c(29, 35),], simplify=FALSE))
zero_counts <- matrix(1, 100, dim(rep_counts)[2])
comb_counts <- rbind(rep_counts, zero_counts)
m1 <- colMeans(comb_counts)
system.time(out <- dash(comp_data = comb_counts,
                        optmethod = "mixEM",
                        mode = m1,
                        def_positions = list("center" = Inf, "null" = 1, "corner" = 1),
                        concentration = c(Inf, 100, 50, 20, 10, 5, 2, 1),
                        weight = list("center" = 100, "null" = 1, "corner" = 1),
                        bf = TRUE,
                        verbose=TRUE))

mat <- Matrix::rsparsematrix(nrow = 500, ncol = 100, density = 0.3)
mat2 <- as.matrix(mat)
mat2[mat2 != 0] = 1
m1 <- colMeans(as.matrix(mat2))
system.time(out <- dash(comp_data = mat2,
                        optmethod = "mixEM",
                        mode = m1,
                        def_positions = list("center" = Inf, "null" = 1, "corner" = 1),
                        concentration = c(Inf, 100, 50, 20, 10, 5, 2, 1),
                        weight = list("center" = 100, "null" = 1, "corner" = 1),
                        bf = TRUE,
                        verbose=TRUE))

