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

##### regression model between faPAR and NDVI : relationship between the 2 variables
erosion <- c(12, 14, 16, 24, 26, 40, 50, 67) #example of values of erosion in a certain area
hm <- c(30, 100, 150, 200, 260, 340, 460, 600) #heavy metals in this area
# now we do a plot between EROSION AND HEAVY METALS

plot(erosion, hm, col="red", pch=19, xlab="erosion", ylab="heavy metals")

# HOW MUCH the 2 variables are related to each others? --> LINEAR MODEL 

model1 <- lm(hm ~ erosion)
summary(model1)
abline(model1)


####### faPAR vs NDVI model 

library(raster)
library(rasterdiv)
install.packages("sf")
library(sf)

setwd("/Users/sofiaprandelli/lab")
faPAR10 <- raster("faPAR10.tif")

plot(faPAR10)
plot(copNDVI)
copNDVI <- reclassify(copNDVI, cbind(253:255, NA), right=TRUE)


library(sf) # to call st_* functions
random.points <- function(x,n)   # x= raster file; n = points
{
lin <- rasterToContour(is.na(x))
pol <- as(st_union(st_polygonize(st_as_sf(lin))), 'Spatial') # st_union to dissolve geometries
pts <- spsample(pol[1,], n, type = 'random')
}


pts <- random.points(faPAR10,1000)
copNDVIp <- extract(copNDVI, pts)
faPAR10p <- extract(faPAR10,pts)


# photosyntesis vs biomass
model2 <- lm(faPAR10p ~ copNDVIp)

plot(copNDVIp, faPAR10p, col="green", xlab="biomass", ylab="photosynthesis")
abline(model2, col="red")

