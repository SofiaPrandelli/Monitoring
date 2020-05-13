# R_code_radiance.r

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
