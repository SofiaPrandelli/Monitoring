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
plot(sd_snt, col=cl)






