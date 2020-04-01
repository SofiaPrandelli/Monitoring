# R code for spatial view of points 

install.packages("sp")
library(sp) 

data(meuse) 

head(meuse)

# coordinates: a simple manner to do it is to use the fn "coordinates"
coordinates(meuse) = ~X+Y          # ~ è tilde

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
