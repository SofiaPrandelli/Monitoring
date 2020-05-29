# R_code_snow.r
# snow as indicator

setwd("/Users/sofiaprandelli/lab/")

install.packages("ncdf4")
library(ncdf4)
library(raster)

snowmay <- raster("c_gls_SCE500_202005180000_CEURO_MODIS_V1.0.1.nc") #giving the name "snow may"
# warning message: in general these images are covering all the world, but into copernicus we are downloading just a part of the image.
#                 It says we are just using a part of the reference system, so the software cananot process the part of the image we didnt download

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






