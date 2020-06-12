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



