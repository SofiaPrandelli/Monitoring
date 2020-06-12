# R_code_interpolation.r
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






