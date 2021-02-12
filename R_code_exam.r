# R CODE EXAM 

# 1. R First code
# 2. R code Multipanel
# 3. R code Spatial
# 4. R code Point pattern
# 5. R code for Multivariate analysis
# 6. R code DVI deforestation 
# 7. R code Remote sensing
# 8. R code Ecosystem functions
# 9. R code pca Remote sensing
# 10. R code faPAR
# 11. R code Radiance
# 12. R code EBV
# 13. R code Cladonia example 
# 14. R code Snow
# 15. R code NO2
# 16. R code Interpolation
# 17. R code sdm
# 18. R code Final project


################################################################## 1. R_code_first.r
install.packages("sp")

library(sp)
data(meuse)

#Let's see how the meuse dataset is structured:
meuse 

#let's look at the first row of the set
head(meuse)

#let's plot two variables
#let's see if the zinc concentration is relate to that of copper
attach(meuse)
plot(zinc,copper)
plot(zinc,copper,col="green")
plot(zinc,copper,col="green",pch=19)
plot(zinc,copper,col="green",pch=19,cex=2)
plot(zinc,copper,col="orange",pch=6,cex=1)


################################################################# 2. R_code_multipanel.r
# we'll use the install package fn
install.packages("sp")
install.packages("GGally") #this is used for the function ggpairs()  
                           # GGally is an extension of ggplot2; it adds functions to reduce the complexity of combining geometric objects with transformed data

library(sp)     # require(sp) will also do the job
library(GGally) # used for the other packages installed:
               
data(meuse) #there is a dataset available named meuse
attach(meuse)

# Excercise: see the names of the variables and plot cadmium versus zinc  
# there are two ways to see the names of the variable:

# first one: 
names(meuse)

# second one: 
head(meuse)
plot(cadmium,zinc,pch=15,col="red",cex=2)

#EXCERCISE: make all the possible pairwise plots of the dataset

#plot(x,cadmium)
#plot(x,zinc)
#plot(...)
#plot is not a good idea! so we use: 

pairs(meuse)
# with the last fn, we did a MULTIPANEL, that are all the graphs in one figure 
# in case you receive the error "the size is too large", you have to reshape with the mouse the graph window

pairs(~ cadmium + copper + lead + zinc, data=meuse) 

# grouping variabiles 
# ~ called tilde --> blocnum alt+126
# comma separates different arguments

pairs(meuse[,3:6])

#EXERCISE: prettify this graph 
pairs(meuse[,3:6],pch=12,col="green",cex=1.5)

#SURPRISEEE--> go up at the beginning of the code and istall another package
#GGally package will prettify the graph
#install the package and use the library(GGally)


ggpairs(meuse[,3:6]) # it's going to be a pairs but by using ggally
                     # reading of the graph: cadium and copper have very high correlation (0<correlation<1)

 
################################################################## 3. R_code_spatial.r
# R code for spatial view of points 

install.packages("sp")
library(sp) 

data(meuse) 

head(meuse)

# coordinates: a simple manner to do it is to use the fn "coordinates"
coordinates(meuse) = ~X+Y          # ~ è tilde

plot(meuse) #it will show the points in a graph 

spplot(meuse, "zinc")  # declare which dataset we want to use and add a variable 
# the zinc is concentrated in the yellow part (upper part of the graph), numbers in the graph are the measurments (amount of the element)


# exercise: plot the spatial amount of copper 
spplot(meuse, "copper", main="Copper concentration") # main = title

# exercise: bubble copper in red 
bubble(meuse, "copper", main="Copper concentrtion"), col="red")


#### importing new data

# download covid_agg.csv from our teaching site and build a folder called lab into dataC
#put the covid_agg.csv file into the folder lab

# setting the working directory: lab
# Mac users
setwd("/Users/yourname/lab/")

covid <- read.table("covid_agg.csv", head=TRUE) # second argument means that the first line isn't data, but is just the title 
                                                # head=TRUE or head=T
head(covid)

attach(covid)
plot(country,cases)
#plot(covid$country, covid$cases)

#this graph doesn't show all the countries, so we have to change type of graph
plot(country,cases, las=0) #las=0 parallel labels to axes
plot(country,cases, las=1) #las=1 horizontal labels to axes
plot(country,cases, las=2) #las=2 perpendicular labels to axes
plot(country,cases, las=3) #las=3 vertical labels to axes

#decrease the size of the axes label with cex.axis=0.5 
plot(country,cases, las=3, cex.axis=0.5) 

#plot spatially with ggplot 
install.packages("ggplot2")
library(ggplot2)
data(mpg)
head(mpg)

#save project



################################################################## 4. R_code_point_pattern_analyses.r
# POINT PATTERN ANALYSIS: DENSITY MAP --> package spatstat: statystical analysis of spatial point pattern

install.packages("spatstat")
library(spatstat)

attach(covid)
head(covid)

# give a name to the object we are goingo to prepare
covids <- ppp(lon,lat,c(-180,180), c(-90,90))  
# ppp(x,y, range of longiture, range of latitude) -> point pattern dataset in the two-dimensional plane

# without attaching the covid set
# the ppp commang becomes: covids <-ppp(covid$lon, covid$lat,c(-180,180), c(-90,90))  

#######FUNCTION DENSITY
d<-density(covids)

plot(d) #to show the density graph
points(covids)



#--- 08/04/2020

setws("C:/lab/") #windows
load("point_pattern.RData")
ls()
#covids:point pattern
#d: density map
library(spatstat) 
install.packages("rgdal")
library(rgdald)

plot(d) 
points(covids)

#let's input vector lines (xoyo, x1y1, x2y2)
#import coastline
coastlines<- readOGR("ne_10m_coastline.shp") #OGR is written in capital letter!

plot(d)
points(covids)
plot(coastline, add=T)

# Change the colour and make the graph beautiful
cl <- colorRampPalette (c("yellow","orange","red"))(100) #color scheme: series of colors
                                                         #100 colors between yellow to red (range of colors)
plot(d, col=clr, main="Densities of covid-19") #title of the map
plot(d,col=cl)
points (covids)
plot(coastlines, add=T)

# EXERCISE: new colour ramp palette
cl <- colorRampPalette (c("blue","green","yellow"))(100)

plot(d,col=cl)
points (covids)
plot(coastlines, add=T)

# export graph
pdf("covid_density.pdf")  # or png ()

# we need to copy all the functions used to define the plot
cl <- colorRampPalette (c("yellow","orange","red"))(100)
plot(d, col=clr, main="Densities of covid-19")
plot(d,col=cl)
points (covids)
plot(coastlines, add=T)
dev.off () #dev=device 


################################################################## 5. R_code_multivariate.r
# R code for multivariate analysis 


install.packages("vegan") # community ecology package
library(vegan)

setwd("/Users/sofiaprandelli/lab")

#import the two tables from IOL: biomes and biomes_types (from Chrome, not Safari) 
biomes <- read.table("biomes.csv", header=T, sep=",") #in the biomes there is an header (so T=true) values are separated by coma
head(biomes) # or view(biomes), biomes


# Multivariate analysis: how the species are related each others?
# DEtrended CORrrespondence ANAlysis = DECORANA
multivar <- decorana(biomes)
plot(multivar)

# analysis of the graph: 
# red colubus, giant orb, tree fern and raflesa are related each otehr --> they are in the same part of the graph
# the same occurs for another part of the graph, for ex. bufo, fox, squirrel, alnus, mosses

biomes_types <- read.table("biomes_types.csv", head=T, sep=",")
head(biomes_types)

attach(biomes_types) 

# we are going to draw an ellipse that connect all the points
# 4 different biomes, so 4 different colors or we can write col=c("green","blue","red","black")
# kind= type of graph --> hull is a convex shape and "e" is for ellipse

ordiellipse(multivar, type, col=1:4, kind = "ehull", lwd=3)

ordispider(multivar, type, col=1:4, label=T)


################################################################## 6. R code DVI deforestation 
# DVI on remote sensing data and deforestation
setwd("/Users/sofiaprandelli/lab")

load("rs.RData")
ls()  # to have a list of all these data into this set

library(raster)
p224r63_1988 <- brick("p224r63_1988_masked.grd")     # DIFFERENT IMAGE FROM THE LAST EXCERCISE 

plot(p224r63_1988)


# plotRGB
#band of Landsat
#B1: blue
#B2: green
#B3:red
#B4:NIR

# EXERCISE: plot in visible RGB 321 both images (1988 and 2011)
par(mfrow=c(2,1))

plotRGB(p224r63_2011, r=3, g=2, b=1, stretch="Lin") # STRETCHING THE DATA: we'll have small changing in the colors..more colors? 
plotRGB(p224r63_1988, r=3, g=2, b=1, stretch="Lin")

# EXERCISE: plot in false color RGB 432 both images

par(mfrow=c(2,1))

plotRGB(p224r63_2011, r=4, g=3, b=2, stretch="Lin")
plotRGB(p224r63_1988, r=4, g=3, b=2, stretch="Lin")


# enhace the noise of the images!
par(mfrow=c(2,1))

plotRGB(p224r63_2011, r=4, g=3, b=2, stretch="hist")
plotRGB(p224r63_1988, r=4, g=3, b=2, stretch="hist")

#DVI= NIR - RED --> stressed plants have very low value of difference vegetation index

# plotRGB
#band of Landsat
#B1: blue
#B2: green
#B3:red: b3_sre
#B4:NIR: b4_sre

dvi2011 <- p224r63_2011$B4_sre - p224r63_2011$B3_sre
cl <- colorRampPalette(c('"yellow"','light blue','lightpink4'))(100) 
plot(dvi2011)


dvi2011 <- p224r63_2011$B4_sre - p224r63_2011$B3_sre
cl <- colorRampPalette(c('darkorchid3','light blue','lightpink4'))(100) 
plot(dvi2011,col=cl)


#EXERCISE: dvi for 1988
dvi1988 <- p224r63_1988$B4_sre - p224r63_1988$B3_sre
cl <- colorRampPalette(c('darkorchid3','light blue','lightpink4'))(100) 
plot(dvi1988,col=cl)


 # DIFFERENCE BETWEEN DVI 1988 AND DVI 2011
diff <- dvi2011 - dvi1988 
plot(diff)

# changing the grain of the image!
# resampling
p224r63_2011res <- aggregate(p224r63_2011, fact=10)
p224r63_2011res100 <- aggregate(p224r63_2011, fact=100)
 
par(mfrow=c(3,1))
plotRGB(p224r63_2011, r=4, g=3, b=2, stretch="Lin")
plotRGB(p224r63_2011res, r=4, g=3, b=2, stretch="Lin")
plotRGB(p224r63_2011res100, r=4, g=3, b=2, stretch="Lin")

p224r63_2011 # per vedere tutte le caratteristiche di un'immagine: es. riga 3--> risoluzione 30metriX30metri
# capire cosa è utile e cosa no dal punto di vista ecologico
# la risoluzione può essere anche data da un problema economico 
# spostamento di una popolazione--> very high resolution data
# monitoring forest changes and not interest in small changes --> medium resolution data
# global changes analysis --> pixel of 1 km can be good 
# high amount of datas can produce a lot of noise, so high resolution is not always useful!




################################################################## 7. R_code_RS.r
# Codice R PER ANALISI DI IMMAGINI SATELLITARI 
# R code remote sensing: remote sensing permits us to see new things into landscape, planets, 
#                        than our eyes can't see

setwd("/Users/sofiaprandelli/lab")
install.packages("raster")
library(raster)

#importiamo l'immagine in R
p224r63_2011 <- brick("p224r63_2011_masked.grd") #la parte iniziale è il nome che diamo all'immagine
plot(p224r63_2011) #grafico

#let's change the color palette: first, introduce the color "cl"
cl <- colorRampPalette(c('black','grey','light grey'))(100)

#EXERCISE: plot the image with the new color ramp palette
plot(p224r63_2011, col=cl)

#let's make a different plot of the image: in order to make a multiframe of different plots --> fn par()
par(mfrow=c(2,2)) # 2x2=4 images

#first image
clb <- colorRampPalette(c('dark blue','blue','light blue'))(100)
plot(p224r63_2011$B1_sre, col=clb)      #dollar $ to LINK different parts: in this case to link the bands

# just do the same for the green band B2_sre
clg <- colorRampPalette(c('dark green','green','light green'))(100)
plot(p224r63_2011$B2_sre, col=clg)

# now red band: B3 RED
clr <- colorRampPalette(c('dark red','red','pink'))(100)
plot(p224r63_2011$B3_sre, col=clr)

#B4 Near Infra Red= NIR
clr <- colorRampPalette(c('red','orange','yellow'))(100)
plot(p224r63_2011$B4_sre, col=clr)

# cambiamo la disposizione delle immagini e mettiamole in un'unica colonna

par(mfrow=c(4,1))
# B1: blue
clb <- colorRampPalette(c('dark blue','blue','light blue'))(100) 
plot(p224r63_2011$B1_sre, col=clb)


# B2 green
# Exercise: do the same for the green band B2_sre
clg <- colorRampPalette(c('dark green','green','light green'))(100) 
plot(p224r63_2011$B2_sre, col=clg)


# B3 red
clr <- colorRampPalette(c('dark red','red','pink'))(100) 
plot(p224r63_2011$B3_sre, col=clr)


# B4 NIR
cln <- colorRampPalette(c('red','orange','yellow'))(100) 
plot(p224r63_2011$B4_sre, col=cln)

#plotRGB : association between the 3 components (red blue and green) 
#          and exacting linearly stretching the colors

#band of Landsat
#B1: blue
#B2: green
#B3:red
#B4:NIR

dev.off() #serve per chiudere i grafici precedenti

plotRGB(p224r63_2011, r=3, g=2, b=1, stretch="Lin")

#now we change the components in NIR, red, green --> vegetation is becoming red (forest is really red) 
#     white and grey parts are the opened areas --> forest in danger by human ctivity 
plotRGB(p224r63_2011, r=4, g=3, b=2, stretch="Lin") 

#EXERCISE: NIR on top of the RGB (Band number 4) 
plotRGB(p224r63_2011, r=3, g=4, b=2, stretch="Lin") #it can also measure the humidity in the forest, 
# this combination allows to see differences between main vegetation and nude soil (agricoltural fields)
#                    (we have put the near IR into the green component)

plotRGB(p224r63_2011, r=3, g=2, b=4, stretch="Lin") 


################################################################## 8. R_code_ecosystem_functions.r
# R CODE TI VIEW BIOMASS OVER THE WORLD AND CALCULATE CHANGES IN ECOSYSTEM FUNCTIONS
# ENERGY
# CHEMICAL CYCLING
# PROXIES

# Global trends in biodiversity and ecosystem services from 1900 to 2050: https://www.biorxiv.org/content/10.1101/2020.04.14.031716v1.full

install.packages("rasterdiv") #RASTER DIVERSITY
install.packages("rasterVis") #raster VISUALIZATION 

library(rasterdiv)
library(rasterVis)

data(copNDVI) # NDVI is the same of DVI but instead of doing NR-ReD

plot(copNDVI) # --> a world map opens, we now remove values from 253 to 255 with the argument cbind

copNDVI <- reclassify(copNDVI, cbind(253:255, NA)) #so it's to remove data not useful for us
levelplot(copNDVI)
levelplot #in the final graph, it shows the biomass related to the plants of the last 30 years at global scale 
          # powerful fn --> with just a fn we have a big output 
          # this map is based on pixels of 8kmX8km
          # how much this ecosyst. fn like energy flows and chemical cycling are 
          
          
#let's change the resolution: 
copNDVI10 <- aggregate(copNDVI, fact=10)    # 10kmX10km pixels --> smooth effect from the original one
levelplot(copNDVI10)

copNDVI100 <- aggregate(copNDVI, fact=100) #100km#100km
levelplot(copNDVI100)

#######################################################

setwd("/Users/sofiaprandelli/lab")

library(raster)

defor1 <- brick("defor1_.jpg")
defor2 <- brick("defor2_.jpg")    #to import both images

#band1: NIR
#band2: RED
#band3: GREEN 

plotRGB(defor1, r=1, g=2, b=3, stretch="Lin")
# we see in the image, forest before the agricolture destroyed them 

plotRGB(defor2, r=1, g=2, b=3, stretch="Lin")
#energy flow and cheical cycles destroyed from agricolture

par(mfrow=c(1,2))
plotRGB(defor1, r=1, g=2, b=3, stretch="Lin")
plotRGB(defor2, r=1, g=2, b=3, stretch="Lin")
#comparison between the 2 images: before and after

dvi1 <- defor1$defor1_.1 - defor1$defor1_.2       #dvi number 1


#defor2
#band1:
#
#
dvi2 <- defor2$defor2_.1 - defor2$defor2_.2



################################################################## 9. R code_pca_remote_sensing
# install both packages raster and RStoolbox
library(raster) 
library(RStoolbox) 

setwd("/Users/sofiaprandelli/lab") 

p224r63_2011 <- brick("p224r63_2011_masked.grd")

#b1 blue
#b2 green
#b3 red
#b4 NIR
#b5 SWIR
#b6 THERMAL IR
#b7 SWIR
#b8 PANCHROMATIC

# RGB:
plot(p224r63_2011, r=5, g=4, b=3, stretch="Lin")
ggRGB(p224r63_2011, 5, 4, 3)


# EXERCIE: DO the same with the 1988 image!
p224r63_1988 <- brick("p224r63_1988_masked.grd")
plotRGB(p224r63_1988, r=5, g=4, b=3, stretch="Lin")


par(mfrow=c(1,2))
plotRGB(p224r63_1988, r=5, g=4, b=3, stretch="Lin")
plotRGB(p224r63_2011, r=5, g=4, b=3, stretch="Lin")


names(p224r63_2011)
plot(p224r63_2011$B1_sre, p224r63_2011$B3_sre)
plot(p224r63_2011$B1_sre, p224r63_2011$B3_sre)


# PCA!
# DECREASE THE RESOLUTION 
p224r63_2011_res <- aggregate(p224r63_2011, fact=10)

# library (RStoolbox) is now needed!
p224r63_2011_pca <- rasterPCA(p224r63_2011_res)

plot(p224r63_2011_pca$map)


# LET'S MAKE a new color ramp different from the basic one 

cl <- colorRampPalette(c('dark grey','grey','light grey'))(100) # 
plot(p224r63_2011_pca$map, col=cl)


summary(p224r63_2011_pca$model)
# PC1 99.83% of the whole variation 

pairs(p224r63_2011)
plotRGB(p224r63_2011_pca$map, r=1, g=2, b=3, stretch="Lin")

#1988
p224r63_1988_res <- aggregate(p224r63_1988, fact=10)
p224r63_1988_pca <- rasterPCA(p224r63_1988_res) 
plot(p224r63_1988_pca$map, col=cl)

summary(p224r63_1988_pca$model)
pairs(p224r63_1988)

# NOW we make a difference between the 2 imagesand then we plot the difference
difpca <- p224r63_2011_pca$map - p224r63_1988_pca$map
plot(difpca)
cldif <- colorRampPalette(c('blue','black','yellow'))(100)
plot(difpca$PC1,col=cldif)


################################################################## 10. R code_faPAR
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



################################################################## 11. R_code_radiance.r

library(raster)
library(rasterVis)
library(rasterdiv)

toy <- raster(ncol=2, nrow=2, xmn=1, xmx=2, ymn=1, ymx=2)

# now we are going to put some values in the raster toy 
values(toy) <- c(1.13,1.44,1.55,3.4)

# now we put the datas into EACH PIXEL
plot(toy)
text(toy, digits=2)

toy2bits <- stretch(toy,minv=0,maxv=3)
storage.mode(toy2bits[]) = "integer"
plot(toy2bits)
text(toy2bits, digits=2)

toy4bits <- stretch(toy,minv=0,maxv=15)
storage.mode(toy4bits[]) = "integer"
plot(toy4bits)
text(toy4bits, digits=2)

toy8bits <- stretch(toy,minv=0,maxv=255)
storage.mode(toy8bits[]) = "integer"

plot(toy8bits)
text(toy8bits, digits=2)

par(mfrow=c(1,4))

plot(toy)
text(toy, digits=2) 

plot(toy2bits)
text(toy2bits, digits=2)

plot(toy4bits)
text(toy4bits, digits=2)

plot(toy8bits)
text(toy8bits, digits=2)



################################################################## 12. R_code_EBVs.r

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


#### PCA analysis

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


################################################################## 13.R_code_cladonia.r

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



################################################################## 14.R_code_snow.r
# snow as indicator

setwd("/Users/sofiaprandelli/lab/")

install.packages("ncdf4")
library(ncdf4)
library(raster)

snowmay <- raster("c_gls_SCE500_202005180000_CEURO_MODIS_V1.0.1.nc") #giving the name "snow may"
# warning message: in general these images are covering all the world, but into copernicus we are downloading just a part of the image.
#                  It says we are just using a part of the reference system, so the software cananot process the part of the image we didnt download

# snow: associated with blue color, so we choose a blue palette 
cl <- colorRampPalette(c('darkblue','blue','light blue'))(100)

#EXERCISE: plot the snow cover with the color ramp palette
plot(snowmay, col=cl)

# now, we import all the snow tif files in snow folder in lab folder

setwd("/Users/sofiaprandelli/lab/snow/")

snow2000 <- raster("snow2000r.tif") 
snow2005 <- raster("snow2005r.tif") 
snow2010 <- raster("snow2010r.tif") 
snow2015 <- raster("snow2015r.tif") 
snow2020 <- raster("snow2020r.tif") 

par(mfrow=c(2,3))   # to import and plot all of the datas
plot(snow2000, col=cl)
plot(snow2005, col=cl)
plot(snow2010, col=cl)
plot(snow2015, col=cl)
plot(snow2020, col=cl)
# we can observe how the amount of snow decreases during the years 

############## another (FAST) method to import and plot in R all the DATAS
# we are going to apply the Lapply function: returns a list of the same length as X, 
# each element of which is the result of applying FUN to the corresponding element of X.

rlist <- list.files(pattern="snow")     # pattern is the argument 
rlist   # to have the output of this function 

# now we can apply the fn Raster to all the list of these 5 files, just one time
import <- lapply(rlist, raster)

# RASTER STACK: the function stack takes different layers and put them in the same stack
snow.multitemp <- stack(import) 
snow.multitemp

# we can simply plot the snow.multitemp in just one line 
plot(snow.multitemp, col=cl)

# SUMMARY: we did the list, then we apply the raster function to every single layer(file) using lapply fn, thenk we stack on different layers all together

# LET'S SEE how to PREDICT how the state of snow will be in the future
# continuing the regression line (negative slope) we can have the next points 
# in this case we have a quite straight line, but depending on the phenomena we have different fn

#GO TO IOL and save the link prediction.r in the folder snow in lab

source("prediction.r")    #source function: it causes R to accept its input from the named file/URL/connection or expressions directly

plot(predicted.snow.2025.norm, col=cl)


############################# Lesson 3/05 - day 2nd

setwd("/Users/sofiaprandelli/lab/snow")

# Exercise: import all together the snow cover images
library(raster)

rlist <- list.files(pattern="snow")   # as we did previously 
rlist

import <- lapply(rlist, raster)
snow.multitemp <- stack(import)
plot(snow.multitemp, col=cl)

cl <- colorRampPalette(c('darkblue','blue','light blue'))(100) 
plot(snow.multitemp, col=cl)

# if you have saved the document on lab --> load("prediction.RData")
# if not: download into snow folder the file prediction outputFile on IOL:
prediction <- raster("predicted.2025.norm.tif")

plot(prediction, col=cl)   # let's plot the prediction using the color we choose before (blue ramp palette) 

# let's export the output: for example if we made the calculation and we want to send the output to a collegue

writeRaster(prediction, "final.tif") #the fn creates a data, not just exporting a graph
                                     # it's the opposite of raster fn: whiteR let us to export from R to the folder
# final stack 
final.stack <- stack(snow.multitemp, prediction)
# 1) at the end of the process, we had the input layers as a stack (pila) wich for us is snow.multitemp
# 2) then we made the map "prediction"
# 3) now we made a final stack composed by the first stack and the prediction: we have 2000,2005,2010,2015,2020 + 2025 
final.stack <- stack(snow.multitemp, prediction)
plot(final.stack, col=cl)

# let's  EXPORT THE R GRAPH AS PDF in the lab folder!
pdf("my_final_exciting_graph.pdf")
plot(final.stack, col=cl)
dev.off()



################################################################## 15. R_code_NO2
# monitoring the trend of no2 in atmosphere during quarantine

setwd("/Users/sofiaprandelli/lab/no2")
library(raster)

# we are going to apply the LAPPLY function: returns a list of the same length as X, 
# each element of which is the result of applying FUN to the corresponding element of X.


# EXERCISE: import all of the no2 data in R BY THE LAPPLY FUNCITON
rlist <- list.files(pattern="EN")
rlist

import <- lapply(rlist, raster)
EN <- stack(import)
cl <- colorRampPalette(c('red','orange','yellow'))(100)
plot(EN, col=cl)


# january and march
par(mfrow=c(1,2)) # with Graphical parameter mfrow we put multiple graphs in a single plot
#let's now plot the first and the last files
plot(EN$EN_0001, col=cl)
plot(EN$EN_0013, col=cl)
# we can see the higher value of no2 in the first image and the lower value: reduction in northern Italy 

dev.off() #to cancel what we did before
# RGB space
plotRGB(EN, r=1, g=7, b=13, stretch="lin") # stretch is to stretch the colors and see them better
                                           # 1st image of the stack is 0001, 7th is 0007, the last one(13rd) is 0013
# with plotRGB, we have put 3 different layers (ANY LAYER) in the RGB space, having the layers contemporary 
# we can manage different images in the same graph
# between Lombardy and Liguria, there is a very big red square --> very first spread of NO2 at the beginning (red in january) 
# red parts are those with higher values, in january, while higher values in blue means in march-->
# we dont have blue peaks


# difference map between the 2 situations 
dif <- EN$EN_0013 - EN$EN_0001
cld <- colorRampPalette(c('blue', 'white', 'red')) (100)
plot(dif, col=cld)
# this graph is just making the difference between the first image and the last image so it doesnt consider the "real" output
# ex: in march there was the lockdown but still using heating


# quantitative estimate of the decrase of NO2
# BOX PLOT FUNCTION: statystical function--> Produce box-and-whisker plot(s) of the given (grouped) values
dev.off()  #to cancel what we did before

boxplot(EN)
boxplot(EN, outline=F) # with the outliens removed 
boxplot(EN,outline=F, horizontal=T) # putting orizonally instead of vertically  
boxplot(EN,outline=F, horizontal=T, axes=T)


# plot 
plot(EN$EN_0001, EN$EN_0013)
 
# we can compare each pixels of 2 stituations plotting the 2 images in order to see if NO2 decreases or increases in each pixel
plot(EN$EN_0001, EN$EN_0013) 
abline(0,1,col="red")   # most of the point are under the line, that means that the NO2 decreases



################################################################## 16. R_code_interpolation.r
# measuring datas where they have not been misurated in a field 

setwd("/Users/sofiaprandelli/lab")

#install.packages(spatstat")
library(spatstat) #

inp <- read.table("dati_plot55_LAST3.csv", sep=";", head=T)
head(inp)

attach(inp)
# instead of writing " plot(inp$X, inp$Y) " we can just write: 
plot(X,Y)

summary(inp)
inppp <- ppp(x=X, y=Y, c(716000,718000),c(4859000,4861000))
marks(inppp) <- Canopy.cov  # spatstat fn is ready to make the estimate: we have points, marks and labels for every point

canopy <- Smooth(inppp)
plot(canopy)
points(inppp, col="green")

marks(inppp) <- cop.lich.mean
lichs <- Smooth(inppp)
plot(lichs)
points(inppp)

par(mfrow=c(1,2))
plot(canopy)
points(inppp)
plot(lichs)
points(inppp)

plot(Canopy.cov, cop.lich.mean, col="red", pch=19, cex=2)
# making comparison between canopy cover anc lichens cover, we see they are inversely proportional


###################### excercise n° 2: Psammophilus vegetation 
# estimates dataset from a sand dunes environment 
# relaionship between vegtation and different conditions of the soil 
#we'll have an idea about how much previous organisms lived in the are casue we'll visualize the actual amount of carbon 

inp.psam <- read.table("dati_psammofile.csv", sep=";" , head=T)
attach(inp.psam)

head(inp.psam)
plot(E,N)
inp.psam.ppp <- ppp(x=E,y=N,c(356450,372240),c(5059800,5064150))
marks(inp.psam.ppp) <- C_org

C <- Smooth(inp.psam.ppp)   # C = Carbon 

plot(C)
plot(C)
points(inp.psam.ppp)

# we can see there are some small powerful predictions, and all the points are grouped in just some areas
# rather than having a clamping effect, we'll make a graph with mean values
# we'll change our view seeing just big points were the little pints are grouped 



################################################################## 17. R_code_SDM.r
# species distribution modeling SDM 
# examples of how to model distribution 
# SDM package: for developing sdm using individual and community-based approaches, generate ensembles of models,
#              evaluate the models and predict species potential distributions in space and time

install.packages("sdm") 
library(sdm) 
library(raster) # for the ecological variables that can be used to predict species distribution 
library(rgdal) # for the species: Geospatial Data Abstraction Library

file <- system.file("external/species.shp", package="sdm") # to make the import of species data
# let's convert the file to a real shape file:
species <- shapefile(file)

plot(species[species$Occurrence == 1,],col='blue',pch=16) # graph where the condition is that a specie is present (occurrence==1)
points(species[species$Occurrence == 0,],col='red',pch=16) # now we add also the absences (red points) 

# now we can import the predictors (ECOLOGICAL VARIABLES) which are going to hepl us predicting the expansion of species
# ECOLOGICAL VARIABLES indexes = elevation, precipitation, temperature, vegetation 
# for this, we'll use a list of files
path <- system.file("external", package="sdm")

lst <- list.files(path=path,pattern=asc$,full.names = T) 
lst

preds <- stack(lst)

cl <- colorRampPalette(c('blue','orange','red','yellow')) (100)
plot(preds, col=cl)

plot(preds$elevation, col=cl)
points(species[species$Occurrence == 1,], pch=16)

plot(preds$temperature, col=cl)
points(species[species$Occurrence == 1,], pch=16)

plot(preds$precipitation, col=cl)
points(species[species$Occurrence == 1,], pch=16)

plot(preds$vegetation, col=cl)
points(species[species$Occurrence == 1,], pch=16)

# model
# now we put all the infos together and make the model
# we explain what to do to the software using the fn sdmData

d <- sdmData(train=species, predictors=preds)

m1 <- sdm(Occurrence ~ elevation + precipitation + temperature + vegetation, data=d, methods="glm")

p1 <- predict(m1, newdata=preds)

# final plot of the prediction 
plot(p1, col=cl) 
points(species[species$Occurrence == 1,], pch=16)

# we can see the highest probability of finding species depends on our 4 variables

s1 <- stack(preds,p1)
plot(s1, col=cl)



################################################################## 18. R_code_final_project
# LNU Lightning Complex fires - California

# I used Onda Dias for the database and Sentinel2 for the analysis

setwd("/Users/sofiaprandelli/lab/project")
library(raster) 

    # trying to convert all the jp2 files in Tiff format: 
    library(rgdal)
    gdal_translate("T10SEH_20200807T184919_B02_20m.jp2", "T10SEH_20200807T184919_B02_20m.tif")
    band1 <- readGDAL("T10SEH_20200807T184919_B02_20m.tif")
    gdal_translate("T10SEH_20200807T184919_B03_20m.jp2", "T10SEH_20200807T184919_B03_20m.tif")
    band2 <- readGDAL("T10SEH_20200807T184919_B03_20m.tif")
    gdal_translate("T10SEH_20200807T184919_B04_20m.jp2", "T10SEH_20200807T184919_B04_20m.tif")
    band3 <- readGDAL("T10SEH_20200807T184919_B04_20m.tif")
    gdal_translate("T10SEH_20200807T184919_B8A_20m.jp2", "T10SEH_20200807T184919_B8A_20m.tif")
    band4 <- readGDAL("T10SEH_20200807T184919_B8A_20m.tif")
    gdal_translate("T10SEH_20200807T184919_B11_20m.jp2", "T10SEH_20200807T184919_B11_20m.tif")
    band5 <- readGDAL("T10SEH_20200807T184919_B11_20m.tif")
    # It doesn't work, I try to convert the files with QGis, using the translate function (in Raster): it works!

rlist20200807 <- list.files(pattern="20200807")
rlist20200807
# B02 Blue -> Band 1
# BO3 Green -> Band 2
# B04 Red -> Band 3
# B11 SWIR -> Band 4
# B08A Vegetation Red Edge -> Band 5 

# I can start the analyses because the resolution of all the bands is the same (20 m)

####### Beginning of August 2020: normal situation and vegetation cover
# Applying the raster function to every single layer using Lapply function: fast method to import and plot in R all the data
import20200807 <- lapply(rlist20200807,raster)

beforeLNU <- stack(import20200807)
plot(beforeLNU)

# Showing the park in human eye colors
plotRGB(beforeLNU, r=3, g=2, b=1, stretch="lin")

# Vegetation Red Edge analysis: vegetation underlined in red
plotRGB(beforeLNU, r=5, g=3, b=2, stretch="lin")

###### Beginning of the fires: August 17th 2020
rlist20200822 <- list.files(pattern="20200822")
rlist20200822
# B02 Blue -> Band 1
# BO3 Green -> Band 2
# B04 Red -> Band 3
# B11 SWIR -> Band 4
# B08A Vegetation Red Edge -> Band 5

# Applying the raster function to every single layer using Lapply function
import20200822 <- lapply(rlist20200822,raster)

# I can make the analyses because resolutions of all the bands are the same (20 m)
august17 <- stack(import20200822)
plot(august17)

# Showing the park in human eye colors
plotRGB(august17, r=3, g=2, b=1, stretch="lin")

# Vegetation Red Edge analysis: vegetation underlined in red
plotRGB(august17, r=5, g=3, b=2, stretch="lin")

###### End of the fires: October 11th 2020
      # In September, fire activity decreased significantly within the complex
      # By mid-September, only the Hennessey and Walbridge Fires were still burning
      # On October 2nd, CAL FIRE reported that the entire complex had been extinguished
rlist20201011 <- list.files(pattern="20201011")
rlist20201011
# B02 Blue -> Band 1
# BO3 Green -> Band 2
# B04 Red -> Band 3
# B11 SWIR -> Band 4
# B08A Vegetation Red Edge -> Band 5

# Applying the raster function to every single layer using Lapply function
import20201011 <- lapply(rlist20201011,raster)

# I can make the analyses because resolutions of all the bands are the same (20 m)
october11 <- stack(import20201011)
plot(october11)

# Showing the park in human eye colors
plotRGB(october11, r=3, g=2, b=1, stretch="lin")

# Vegetation Red Edge analysis: vegetation underlined in red
plotRGB(october11, r=5, g=3, b=2, stretch="lin")

###### January 9th 2021: after 3 months the end of LNU
rlist20210109 <- list.files(pattern="20210109")
rlist20210109
# B02 Blue -> Band 1
# BO3 Green -> Band 2
# B04 Red -> Band 3
# B11 SWIR -> Band 4
# B08A Vegetation Red Edge -> Band 5

# Applying the raster function to every single layer using Lapply function
import20210109 <- lapply(rlist20210109,raster)

# I can make the analyses because resolutions of all the bands are the same (20 m)
january9 <- stack(import20210109)
plot(january9)

# Showing the park in human eye colors
plotRGB(january9, r=3, g=2, b=1, stretch="lin")

# Vegetation Red Edge analysis: vegetation underlined in red
plotRGB(january9, r=5, g=3, b=2, stretch="lin")


# Now I can compare all the pictures to see how much the vegetation changed:                  

# Human eye analysis showing the differences from August 2020 to January 2021
# Par function: in order to make a multiframe of different plots
par(mfrow=c(1,4))
plotRGB(beforeLNU, r=3, g=2, b=1, stretch="lin", main="07/08/2020", axes = TRUE)
plotRGB(august17, r=3, g=2, b=1, stretch="lin", main="17/08/2020", axes = TRUE)
plotRGB(october11, r=3, g=2, b=1, stretch="lin", main="11/10/2020", axes = TRUE)
plotRGB(january9, r=3, g=2, b=1, stretch="lin", main="09/01/2021", axes = TRUE)

# Vegetation Red Edge analysis showing the differences from August 2020 to January 2021        
par(mfrow=c(1,4)) 
plotRGB(beforeLNU, r=5, g=3, b=2, stretch="lin", main="07/08/2020", axes = TRUE)
plotRGB(august17, r=5, g=3, b=2, stretch="lin", main="17/08/2020", axes = TRUE)
plotRGB(october11, r=5, g=3, b=2, stretch="lin", main="11/10/2020", axes = TRUE)
plotRGB(january9, r=5, g=3, b=2, stretch="lin", main="09/01/2021", axes = TRUE) 

# Burnt area analysis showing the differences from August 2020 to January 2021
     # vegetation in green and burnt area in natural colors
plotRGB(beforeLNU, r=4, g=5, b=3, stretch="lin", main="Burnt area 07/08/2020", axes = TRUE)
plotRGB(august17, r=4, g=5, b=3, stretch="lin", main="Burnt area 17/08/2020", axes = TRUE)
plotRGB(october11, r=4, g=5, b=3, stretch="lin", main="Burnt area 11/10/2020", axes = TRUE)
plotRGB(january9, r=4, g=5, b=3, stretch="lin", main="Burnt area 09/01/2021", axes = TRUE)

# Differences between shortly after the end of the fires and after 3 months the end of the fires:
    # Human eye analysis           
    par(mfrow=c(1,2))                    
    plotRGB(october11, r=3, g=2, b=1, stretch="lin", main="11/10/2020", axes = TRUE)
    plotRGB(january9, r=3, g=2, b=1, stretch="lin", main="09/01/2021", axes = TRUE)

    # Vegetation Red Edge analysis
    par(mfrow=c(1,2))
    plotRGB(october11, r=5, g=3, b=2, stretch="lin", main="11/10/2020", axes = TRUE)
    plotRGB(january9, r=5, g=3, b=2, stretch="lin", main="09/01/2021", axes = TRUE)


###### NDVI calculation: Normalized Difference Vegetation Index
       # NIR - RED / NIR + RED 
       # Using Vegetation Red Edge instead of Infra Red

ndvibeforeLNU <- (beforeLNU$T10SEH_20200807T184919_B8A_20m - beforeLNU$T10SEH_20200807T184919_B04_20m) / 
    (beforeLNU$T10SEH_20200807T184919_B8A_20m + beforeLNU$T10SEH_20200807T184919_B04_20m)
ndviAugust <- (august17$T10SEH_20200822T184921_B8A_20m - august17$T10SEH_20200822T184921_B04_20m) / 
    (august17$T10SEH_20200822T184921_B8A_20m + august17$T10SEH_20200822T184921_B04_20m)
ndviOctober <- (october11$T10SEH_20201011T185321_B8A_20m - october11$T10SEH_20201011T185321_B04_20m) / 
    (october11$T10SEH_20201011T185321_B8A_20m + october11$T10SEH_20201011T185321_B04_20m)
ndviJanuary <- (january9$T10SEH_20210109T185751_B8A_20m - january9$T10SEH_20210109T185751_B04_20m) / 
    (january9$T10SEH_20210109T185751_B8A_20m + january9$T10SEH_20210109T185751_B04_20m)

# Plotting NDVI
clNDVI = colorRampPalette(c("darkblue","yellow","red","black"))(100) 
par(mfrow=c(1,4))
plot(ndvibeforeLNU, col = clNDVI, main = "07/08/2020")
plot(ndviAugust, col = clNDVI, main = "22/08/2020")
plot(ndviOctober, col = clNDVI, main = "11/10/2020")
plot(ndviJanuary, col = clNDVI, main = "09/01/2021")

clNDVI = colorRampPalette(c("blue", "white", "red"))(256)
par(mfrow=c(1,4))
plot(ndvibeforeLNU, col = clNDVI, main = "07/08/2020")
plot(ndviAugust, col = clNDVI, main = "22/08/2020")
plot(ndviOctober, col = clNDVI, main = "11/10/2020")
plot(ndviJanuary, col = clNDVI, main = "09/01/2021")

