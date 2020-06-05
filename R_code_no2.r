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
cl <- colorRampPalette(c('red','orange','yellow'))(100) #
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



