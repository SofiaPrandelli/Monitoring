install.packages("sp")

library(sp)
data(meuse)

#Let's see how the meuse dataset is structured:
meuse 

#let'ss look at the first row of the set
head(meuse)

#let's plot two variables
#let's see if the zinc concentration is realte to that of copper
attach(meuse)
plot(zinc,copper)
plot(zinc,copper,col="green")
plot(zinc,copper,col="green",pch=19)
plot(zinc,copper,col="green",pch=19,cex=2)
plot(zinc,copper,col="orange",pch=6,cex=1)
