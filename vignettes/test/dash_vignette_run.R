

################  vignette run dash  ###########################

xmat <- cbind(c(5, 0, 2, 0),
              c(1, 1, 0, 1),
              c(100, 100, 50, 100),
              c(20, 50, 100, 10),
              c(10, 10, 200, 20),
              c(50, 54, 58, 53),
              c(1,1,1,3),
              c(2, 4, 1, 1))
rownames(xmat) <- c("A", "C", "G", "T")
colnames(xmat) <- paste0("pos-", 1:dim(xmat)[2])
xmat_norm <- apply(xmat, 2, function(x) return(x/sum(x)))

xmat

out <- dash(xmat, optmethod = "mixEM", verbose=FALSE, bf=TRUE)

grid.newpage()
layout.rows <- 1
layout.cols <- 2
top.vp <- viewport(layout=grid.layout(layout.rows, layout.cols,
                                      widths=unit(rep(6,layout.cols), rep("null", 2)),
                                      heights=unit(c(20,50), rep("lines", 2))))

plot_reg <- vpList()
l <- 1
for(i in 1:layout.rows){
  for(j in 1:layout.cols){
    plot_reg[[l]] <- viewport(layout.pos.col = j, layout.pos.row = i, name = paste0("plotlogo", l))
    l <- l+1
  }
}


plot_tree <- vpTree(top.vp, plot_reg)

color_profile = list("type" = "per_row",
                     "col" = RColorBrewer::brewer.pal(4,name ="Spectral"))

pushViewport(plot_tree)
seekViewport(paste0("plotlogo", 1))
logomaker(xmat_norm,color_profile = color_profile,
          frame_width = 1,
          pop_name = "pre dash PWM",
          newpage = F)

seekViewport(paste0('plotlogo',2))
logomaker(out$posmean,color_profile = color_profile,
          frame_width = 1,
          pop_name = "post dash PWM",
          newpage = F)
