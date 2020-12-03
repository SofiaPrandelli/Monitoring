########## R_code_spatial.r
# R code for spatial view of points 

install.packages("sp")
library(sp) 

data(meuse) 

head(meuse)

# coordinates: a simple manner to do it is to use the fn "coordinates"
coordinates(meuse) = ~X+Y          # ~ --> tilde

plot(meuse) #it will show the points in a graph 

spplot(meuse, "zinc")

# exercise: plot the spatial amount of copper 
spplot(meuse, "copper", main="Copper concentration")

# exercise: bubble copper in red 
bubble(meuse, "copper", main="Copper concentrtion"), col="red")

#### importing new data

# download covid_agg.csv from our teaching site and build a folder called lab into dataC
#put the covid_agg.csv file into the folder lab

# setting the working directory: lab
# Mac users
setwd("/Users/yourname/lab/")

covid <- read.table("covid_agg.csv", head=TRUE)



########## R_code_snow.r
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

# let's  EXPORT THE R GRAPH AS PDF --> you'll find it in lab folder! ˆIMPORTANT!
pdf("my_final_exciting_graph.pdf")
plot(final.stack, col=cl)
dev.off()




########## R_code_sdm.r

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



########## R_code_rs.r
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


#################


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
#capire cosa è utile e cosa no dal punto di vista ecologico
# la risoluzione può essere anche data da un problema economico 
# spostamento di una popolazione--> very high resolution data
# monitoring forest changes and not instereste in small changes --> medium resolution data
# global changes analysis --> pixel of 1 km can be good 
# high amount of datas can produce a lot of noise, so high resolution is not always useful!














