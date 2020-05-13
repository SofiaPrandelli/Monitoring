# R_code_faPAR.r
# how to look at chemical cycling from satellites 

# install.packages("raster")

 setwd("/Users/sofiaprandelli/lab")

library(raster)
library(rasterVis)
library(rasterdiv)

plot(copNDVI9)
copNDVI <- reclassify(copNDVI, cbind(253:255, NA))

levelplot(copNDVI)
faPAR10 <- raster("faPAR10.tif")

levelplot(faPAR10)
|
# WHAT HAPPEN IN THIS CASE? differences with the previous graph? there's a big rate of photosyntesis in the euator band
# while on top of the northern part the NDVI it's not so high --> most of the light is going through plants and down to the soil so not used for photosyntesis
# also in conifer forest: light is not so used like in tropical forest (where all the light is udes) 
# DIFFERENT CAPABILITY OF TAKING LIGHT 
#carbon uptake of plants--> very high in tropical forests 

#save the plot as pdf
pdf("copNDVI.pdf")
levelpot(copNDVI)
dev.off()


pdf("faPAR.pdf")
levelplot(faPAR10)
dev.off()                     # to safe plots as PDF from R



# estimate where is the relation between NDVI and faPAR 

####################################################### day 2

setwd("C:/lab/") 
load("faPAR.RData")
library(raster)
library(rasterdiv)
library (rasterVis)

#the original faPAR from Copernicus is  2GB
# let's see how much space is needed for an 8-bit set

writeRaster( copNDVI, "copNDVI.tif")
# 5.3 MB

# to make the level plot of the faPAR
levelplot(faPAR10) #faPAR = fraction of the solar radiation absorbed by live leaves 

