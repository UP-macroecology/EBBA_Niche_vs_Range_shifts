### Run Ecospat Analyses with full climatic background to quantify overall climate shift between time periods ###

# File pathways made relative to '.Data' need to be checked

# load packages
library(sp)
library(rgdal)
require(terra)
require(raster)
require(ggmap)
require(maptools)
require(rgeos)
require(ecospat)
require(scales)
library(ade4)

#load bioclimate raster stacks
Bioclim80s <- Raster("./Data/climate_inputs/Bioclim80s.tif")
Bioclim10s <- Raster("./Data/climate_inputs/Bioclim10s.tif")

#convert and remove NAs
Background80s <- as.matrix(Bioclim80s)
Background80s <- Background80s[!rowSums(!is.finite(Background80s)),]

Background10s <- as.matrix(Bioclim10s)
Background10s <- Background10s[!rowSums(!is.finite(Background10s)),]

#calculate the PCA scores
pca.cal <- dudi.pca(rbind(Background10s, Background80s), scannf= FALSE, nf= 2)

historic.env.scores <- suprow(pca.cal, Background80s[,1:19])$li
recent.env.scores <- suprow(pca.cal, Background10s[,1:19])$li
climate.allcell.scores <- pca.cal$li

#Calculate the occurrence density grids
grid.clim.historic <- ecospat.grid.clim.dyn(climate.allcell.scores, historic.env.scores, historic.env.scores, R = 100)

grid.clim.recent <- ecospat.grid.clim.dyn(climate.allcell.scores, recent.env.scores, recent.env.scores, R = 100)

# Calculate niche overlap-
# Compute Schnoer´s D, index of niche overlap
D.overlap <- ecospat.niche.overlap(grid.clim.historic, grid.clim.recent, cor = TRUE)$D 
paste("Schoners D=", D.overlap)

# Niche similarity test
#It is recommended to use at least 1000 replications for the similarity test.
sim.test <- ecospat.niche.similarity.test(grid.clim.historic,grid.clim.recent, rep = 1500, overlap.alternative = "higher", expansion.alternative = "lower", stability.alternative = "higher", unfilling.alternative = "lower", rand.type = 2)

saveRDS(sim.test, file = "./Data/Europe_ebba_species/European_Bioclimate_Niche_similarity_test.Rds")

#readRDS("./Data/Europe_ebba_species/European_Bioclimate_Niche_similarity_test.Rds")

# plot similarity test
tiff(filename = paste("./plots/European_Bioclimate_Niche_Similarity_test.tiff",sep=""))

ecospat.plot.overlap.test(sim.test, "D", "Niche Similarity")

dev.off()
# delimiting niche categories and quantifying niche dynamics in analogue climates
full.euro.niche.dyn<- ecospat.niche.dyn.index(grid.clim.historic, grid.clim.recent, intersection = NA)#run with 0=only overlapping enviros, NA= all available enviros

paste("Euopean climate")
full.euro.niche.dyn$dynamic.index.w # only overlapping:unfilling = 0; stability = 1; expansion = 0. BUT, for all climates:unfilling = 0.00409; stability = 0.999; expansion = 0.00120.

# visualizing niche categories, niche dynamics and climate analogy between ranges
tiff(filename = paste("./plots/European_Bioclimate_enviro_space.tiff",sep=""),width=1200,height=1200,res=180)

ecospat.plot.niche.dyn(grid.clim.historic, grid.clim.recent, quant = 0.25, interest = 2, title = "Niche overlap", name.axis1 = "PC1", name.axis2 = "PC2",
                       col.exp="green", col.unf = "red", col.stab = "blue", colZ1 = "blue", colZ2="green")

ecospat.shift.centroids(historic.env.scores, recent.env.scores, climate.allcell.scores, climate.allcell.scores)

#legend(3, 7.5, legend = c("expansion", "stability", "unfilling"), col = c("green", "blue", "red"), pch = 15, cex = 1.5, pt.cex = 2)

dev.off()

##Extract niche-space centroids##
full.euro.niche.dyn$historic.niche.centroid.x <- mean(historic.env.scores$Axis1)
full.euro.niche.dyn$historic.niche.centroid.y <- mean(historic.env.scores$Axis2)

full.euro.niche.dyn$recent.niche.centroid.x <- mean(recent.env.scores$Axis1)
full.euro.niche.dyn$recent.niche.centroid.y <- mean(recent.env.scores$Axis2)

##Calculate magnitute of niche centroid shift with some phythagoras##
full.euro.niche.dyn$niche.centroid.shift <- sqrt(((full.euro.niche.dyn$recent.niche.centroid.x - full.euro.niche.dyn$historic.niche.centroid.x)^2) + ((full.euro.niche.dyn$recent.niche.centroid.y - full.euro.niche.dyn$historic.niche.centroid.y)^2))

print(rbind(c("Historic bioclimate centroid =", paste(full.euro.niche.dyn$historic.niche.centroid.x, full.euro.niche.dyn$historic.niche.centroid.y, sep=",")), c("Recnent bioclimate centroid =", paste(full.euro.niche.dyn$recent.niche.centroid.x, full.euro.niche.dyn$recent.niche.centroid.y,sep=",")), c("Bioclimate centroid shift =", full.euro.niche.dyn$niche.centroid.shift)))


