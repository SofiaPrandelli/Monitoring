# install.packages("raster")
library(raster)

setwd("/Users/sofiaprandelli/lab")

snt <- brick("snt_r10.tif")    # we brick the image Sentinel
plot(snt)

# B1 blue
# B2 green
# B3 red
# b4 NIR 

# R3 G2 B1 
plotRGB(snt, 3,2,1, stretch="lin")    # we can observe a very complex system:from mountains to fields
plotRGB(snt,4,3,2, stretch="lin")     # VEGETATION IS REFLECTING IN THE NIR, so if we put NIR at the top, we'll see a lot of red

pairs(snt)


#### PCA analysis

library(RStoolbox) # this is for PCA

sntpca <- rasterPCA(snt)
sntpca

summary(sntpca$model)
# 70% of information 
plot(sntpca$map)

plotRGB(sntpca$map, 1, 2, 3, stretch="lin")

# set the moving window 

window <- matrix(1, nrow = 5, ncol = 5)
window


sd_snt <- focal(sntpca$map$PC1, w=window, fun=sd) 

# focal function: Calculate focal ("moving window") values for the neighborhood of focal cells
#                 using a matrix of weights, perhaps in combination with a function.

cl <- colorRampPalette(c('dark blue','green','orange','red'))(100) 
plot(sd_snt, col=cl, main=diversity)



####################### 27/05 - Cladonia example

# Cladonia stellaris calaita: index of clean air; it it has an important role in regulating nutrient cycling and soil microbiological communities.
# FOCAL ON Cladonia
# we'll make a windows and we'll calculate variability on the image

library(raster) 
library(RStoolbox)
setwd("/Users/sofiaprandelli/lab") 
# to import the file (we can use 2 functions: Raster(it imports one single band) or Brick(it can imports several layers at the time))
# here we'll use BRICK since we have several layers 

clad <- brick("cladonia_stellaris_calaita.JPG")

# now we use the focal fn on top of the image 
# first, we should decide the window (frazione dell'immagine, in pixel) then calculate the standard deviation and report the SD in evry pixel  

window <- matrix(1, nrow=3, ncol=3) # number 1 doesnt matter for the calculation
window # matrix of 3 x 3 pixels 

# Focal function: it caluclates values for the neighborhood of focal cells

### PCA analysis for Cladonia 
cladpca <- rasterPCA(clad)
cladpca # to seeall the datas that are putput of this function 

summary(cladpcs$model)
#the proportion of variance is 0.98, that means 98% of component 1 -> l'immagine è fatta nella parte visibile dello schermo 

# now we see how much information is explained by the PCA
plotRGB(cladpca$map, 1, 2, 3, stretch="lin")

# set the moving window
window <- matrix(1, nrow = 5, ncol = 5)
window

# focal function: calculate values for neighborhood of focal cells
sd_clad <- focal(sntpca$map$PC1, w=window, fun=sd) 

# aggregate function 
PC1_agg <- aggregate(cladpca$map$PC1, fact=10)
sd_clad_agg <- focal(PC1_agg, w=window, fun=sd)


par(mfrow=c(1,2)) # to see different graphs in the same image 
cl <- colorRampPalette(c('yellow','violet','black'))(100)
plot(sd_clad, col=cl)
plot(sd_clad_agg, col=cl) # Cladonia set aggregated 

#plot the calculation 
par(mfrow=c(1,2)) 
cl <- colorRampPalette(c('yellow','violet','black'))(100) 
plotRGB(clad, 1,2,3, stretch="lin")
plot(sd_clad, col=cl)

dev.off()
q()
