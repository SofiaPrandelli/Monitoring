# R CODE TI VIEW BIOMASS OVER THE WORLD AND CALCULATE CHANGES IN ECOSYSTEM FUNCTIONS
# ENERGY
# CHEMICAL CYCLING
# PROXIES

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




