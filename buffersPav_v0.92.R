############################################################################
#   A script to calculate forest cover (and the cover of other classes)    #
#    and some landscape metrics in multiple buffers around multiple sites. #
#   by Pavel Dodonov - pdodonov@gmail.com                                  #
#   Spatial Ecology Lab, Institute of Biology, Federal University          #
#    of Bahia                                                              #
#   No rights reserved                                                     #
############################################################################

### Load the required packages

library(maptools)
library(rgdal)
library(raster)
library(landscapemetrics) # not needed for forest cover, needed for edge density and number of patches
library(rgeos)
library(igraph)


### The follwoing objects are used for the function's input:
### A land use map, which has to be a RasterLayer object, with different numbers representing different land use classes;
### An object showing the location of the sites, which has to be a SpatialPointsDataFrame;
### Note that both objects must use the same coordinate system. It should be a metric system, i.e. UTM (Albers should work as well).

### The calculate.cover function has the following arguments:
### map - raster object with land uses
### sites - SpatialPointsDataFrame object with the study sites
### buffers - size (in map units, therefore meters) of the buffers to be calculated
### classes=list() - list of the land uses of interest and their codes, e.g. list(forest=11, open=c(31,35))
### if classes=="all", the frequency of all classes is returned
### site.name.slot - name of the part of the sites object that contains site identification (usually "ID")
### print.progress=T - should the progress be printed as the loop works? Might not work in RStudio.

### The calculate.config function calculates two landscape metrics (edge density and number of patches) for binary landscapes and uses the landscapemetrics package. It has the following arguments:
### map - raster object with land uses
### sites - SpatialPointsDataFrame object with the study sites
### buffers - size (in map units, therefore meters) of the buffers to be calculated
### patch.class - the class number (in the binary map) corresponding to habitat or patches
### site.name.slot - name of the part of the sites object that contains site identification (usually "ID")
### print.progress=T - should the progress be printed as the loop works? Might not work in RStudio.
### calculate.ED - should edge density be calculated? TRUE or FALSE
### calculate.NP - should number of patches be calculated? TRUE or FALSE

### Load the new functions:


calculate.cover <- function(map, sites, buffers, classes="all", site.name.slot, print.progress=T)
{
	# This function calcul	ates the cover of different classes in the landscape
	Nsites <- length(sites)
	Nbuffers <- length(buffers)
	site.names <- as.character(sites@data[site.name.slot][[1]])
	foo <- matrix(nrow=Nsites, ncol=Nbuffers)
	rownames(foo) <- site.names
	colnames(foo) <- paste("b",buffers,sep="_")
	if(length(classes)==1) {
		if(classes=="all") {
			classes <- list()
			temp <- unique(map)
			for(i in 1:length(temp)){
				classes[[i]] <- temp[i]
			}
			names(classes) <- temp
		}
	}
	result <- list()
	for(k in 1:length(classes)) result[[k]] <- foo
	names(result) <- names(classes)
	result[[length(classes)+1]] <- foo
	names(result)[length(classes)+1] <- "NAvaluesInBuffer"
	for(i in 1:Nsites) {
		name.now <- site.names[i]
		site.now <- sites[sites@data[site.name.slot][[1]]==name.now,]
		for(j in 1:Nbuffers) {
			if(print.progress) print(paste("site", i, "of", Nsites, "and buffer", j, "of", Nbuffers))
			buffer.now <- buffer(site.now, width=buffers[j])
			map.now <- crop(map, extent(buffer.now))
			result[["NAvaluesInBuffer"]][i,j] <- sum(is.na(unlist(map.now@data@values)))
			raster.now <- mask(x=map.now, mask=buffer.now, inverse=F)
			freq.now <- freq(raster.now)
			for(k in 1:length(classes)) {
				result[[k]][i,j] <- sum(freq.now[,2][freq.now[,1] %in% classes[[k]]])/sum(freq.now[,2][!is.na(freq.now[,1])])
			}
		}
	}
	return(result)
}



calculate.config <- function(map, sites, buffers, patch.class, site.name.slot, print.progress=T, calculate.ED=T, calculate.NP=T)
{
	Nsites <- length(sites)
	Nbuffers <- length(buffers)
	site.names <- as.character(sites@data[site.name.slot][[1]])
	foo <- matrix(nrow=Nsites, ncol=Nbuffers)
	rownames(foo) <- site.names
	colnames(foo) <- paste("b",buffers,sep="_")
	result <- list()
	if(calculate.ED + calculate.NP == 0) {
		print("It makes no sense to run this function with neither edge density nor patch number; please rethink your choices in life!")
		} else {
		for(k in 1:sum(c(calculate.ED,calculate.NP))) result[[k]] <- foo
		names(result) <- c("EdgeDensity", "NumberOfPatches")[c(calculate.ED, calculate.NP)]
		result[[k+1]] <- foo
		names(result)[k+1] <- "NAvaluesInBuffer"
		for(i in 1:Nsites) {
			name.now <- site.names[i]
			site.now <- sites[sites@data[site.name.slot][[1]]==name.now,]
			for(j in 1:Nbuffers) {
				if(print.progress) print(paste("site", i, "of", Nsites, "and buffer", j, "of", Nbuffers))
				buffer.now <- buffer(site.now, width=buffers[j])
				map.now <- crop(map, extent(buffer.now))
				result[["NAvaluesInBuffer"]][i,j] <- sum(is.na(unlist(map.now@data@values)))
				raster.now <- mask(x=map.now, mask=buffer.now, inverse=F)
				if(calculate.ED) {
					ED.now <- lsm_l_ed(raster.now)
					result[["EdgeDensity"]][i,j] <- ED.now$value
				}
				if(calculate.NP) {
					NP.now <- lsm_c_np(raster.now)
					foo <- ifelse(any(NP.now$class==patch.class), NP.now$value[NP.now$class==patch.class], 0)
					result[["NumberOfPatches"]][i,j] <- foo
				}
			}
		}
	return(result)
	}
}


coverlist2df <- function(x) {
	for(i in 1:length(x)){
		names(landuse.real)[i]
		foo <- landuse.real[[i]]
		foo <- data.frame(foo)
		bar <- reshape(data=foo, varying=colnames(foo), v.names="cover", timevar="buffer", times=colnames(foo), direction="long", ids=rownames(foo))
		rownames(bar) <- 1:nrow(bar)
		bar$class <- names(x)[i]
		if (i == 1) foobar <- bar else foobar <- rbind(foobar, bar)
	}
	foobar$buffer <- gsub("b_","",foobar$buffer)
	foobar$buffer <- as.numeric(foobar$buffer)
	foobar <- data.frame(foobar)
	if(all(colnames(foobar) %in% c("id","buffer","class","cover"))) foobar <- foobar[,c("id","buffer","class","cover")]
	return(foobar)
}




### Quick start for experienced users (just replace with your files):

setwd("/home/pavel/Profissional/script_buffers/")
map.ras <- raster("Belmonte_raster.tiff")
values(map.ras)[values(map.ras)==0] <- NA
sites.real <- readOGR("Pontos_Morcegos.shp")
sites.real <- sites.real[sites.real$regiao=="sul",]

### To calculate cover of all classes and transform object into dataframe:

landuse.real <- calculate.cover(map=map.ras, sites=sites.real, buffers=c(42, 666, 2020), classes="all", site.name.slot="codigo_sit", print.progress=T)

### This results in a list; to convert to dataframe:

landuse.real.df <- coverlist2df(landuse.real)
View(landuse.real.df)

### If you want to export the results as text file:

write.table(landuse.real.df, file="landuse_real.txt", quote=F, sep="\t")



### In case you don't want all classes, but only, say, forest (classes 11, 12, 13) and something that I don't remember what it is (classes 32 and 33):

landuse.real.2 <- calculate.cover(map=map.ras, sites=sites.real, buffers=c(42, 666, 2020), classes=list(forest=c(11,12,13), something=c(32,32)), site.name.slot="codigo_sit", print.progress=T)

### The result will be a list, with each element corresponding to a class, rows as sites and columns as buffer sizes:
landuse.real.2

### In case you prefer a data.frame, use the function as above:
landuse.real.2.df <- coverlist2df(landuse.real.2)


### To calculate edge density and number of patches:

### First we need to reclassify the raster into a binary one:
reclass <- matrix(ncol=3, nrow=3) # reclassification matrix
reclass[1,] <- c(0,2.10,1) # classes smaller than 11 will be "1"
reclass[2,] <- c(10.9,13.1,2) # classes 11 to 13 (forest) will be "2"
reclass[3,] <- c(13.9,99,1) # all greater than 11 will be "2"

map.bin <- reclassify(map.ras, rcl=reclass)


### And now use this binary raster for the calculations; the output will be a list

config.real <- calculate.config(map=map.bin, sites=sites.real, buffers=c(42,666,2020), patch.class=2, site.name.slot="codigo_sit", print.progress=T, calculate.ED=T, calculate.NP=T)
config.real

### If you want as data.frame, this might work:
config.real.df <- coverlist2df(config.real) # Here "class" is actually "metric".






##### Below are detailed instructions for people familar with R but not familiar with GIS-stuff in R.



### Load the data; here is an example that worked for me

setwd("/home/pavel/Profissional/script_buffers/")

### If your map is in shapefile, you'll have to rasterize it.

map.shp <- readOGR("Sisbiota_Sul_2011-07.shp")
map.shp

### To rasterize, we need to replace land use names with classes. I'll use the classe2 column because it has no special characters and special characters, therefore no problems with encoding.
### There's probably a simpler way to make the replacement, but I'm feeling lazy to look for it :-)

bar <- map.shp@data$classe2 ### extract the classes; just to facilitate coding.

unique(bar)

### This shows the values in classe2: forest    pasture   water     eucalypt  cabruca   rural     urban grassland road 

foo <- numeric(length(bar)) # an empty vector with number of elements equal to the number of features in the shapefile.
foo[bar == "forest"] <- 11 # I chose to give the number 11 to all forest
foo[bar == "pasture"] <- 35
foo[bar == "water"] <- 91
foo[bar == "eucalypt"] <- 25
foo[bar == "cabruca"] <- 21
foo[bar == "rural"] <- 29
foo[bar == "urban"] <- 41
foo[bar == "road"] <- 45
foo[bar == "grassland"] <- 31
### These numbers are arbitrary but sort of increasing in order of decreasing matrix permeability.
any(foo==0) # If this states FALSE, it means that all values have been succesfully replaced.

### To insert this as a new column in the map file:

map.shp@data$classeNum <- foo

mapColors <- topo.colors(length(unique(foo)))[as.numeric(as.factor(foo))]

plot(map.shp, col=mapColors, lty=0)

### A map just for fun

### Now we need rasterize it. The finer the resolution, the longer it takes. So if you just want to make a test, I suggest using a really coarse resolution, such as resolution=100. For the true calculations use the resolution that your study requires (finer than 1 m would probably be overkill).

map.ras <- raster(ext=extent(map.shp), resolution=10, crs=crs(map.shp))
map.ras <- rasterize(x=map.shp, y=map.ras, field="classeNum")

### You migh get some error like "Error in x$.self$finalize() : attempt to apply non-function"  If you ever find out what it means, let me know!

plot(map.ras)

### Or, if your map is already a raster, you can just open it:

map.ras2 <- raster("Belmonte_raster.tiff")

### You can see the objects' properties by writing their names:
map.shp
map.ras
map.ras2

### Important detail: all projectios must be in meters. This may be shown as "+units=m" in the crs field of these objects. If your projection is not in meters, fix it :-)

### It may happen that some value (such as 0) is actually missing data. In this case, do the following:

values(map.ras2)[values(map.ras2)==0] <- NA

### Now we need points around which we'll calculate the landscape metrics. You can either open s shapefile or make, for example, a grid of points. Here we'll use both methods.

### Reading a file
sites.real <- readOGR("Pontos_Morcegos.shp")
sites.real <- sites.real[sites.real$regiao=="sul",]

sites.real

### It should have the same crs as the map. If it doesn't, you should probably do something about it.


### Making a grid with 666 points. Thanks to Ivana Cardoso for help with this part!
contorno <- as(extent(map.ras), "SpatialPolygons")
sites.grid <- makegrid(contorno, n=666)
sites.grid <- SpatialPointsDataFrame(sites.grid, data=data.frame(ID=1:nrow(sites.grid)))
proj4string(sites.grid) <- CRS("+init=epsg:32724")

### And now we can run the the functions! Let's say we want to calculate the amount of forest (class 11), open areas (classes 31 and 35) and eucalypt (class 25) for buffers with widths of 42, 666 and 2020 m around our sampling points.
### In the sites.real object, the sites are identified by the "codigo_sit" column in the @data slot.
### In the sites.grid object, the sites are identified by "ID".

landuse.real <- calculate.cover(map=map.ras, sites=sites.real, buffers=c(42, 666, 2020), classes=list(forest=11, open=c(31,35), eucalypt=25), site.name.slot="codigo_sit", print.progress=T)


landuse.real <- calculate.cover(map=map.ras2, sites=sites.real, buffers=c(42, 666, 2020), classes="all", site.name.slot="codigo_sit", print.progress=T)

landuse.real # Please note that there are some NAs in the largest buffer - perhaps these areas have not been mapped.

landuse.grid <- calculate.cover(map=map.ras, sites=sites.grid, buffers=c(42, 666, 2020), classes=list(forest=11, open=c(31,35), eucalypt=25), site.name.slot="ID", print.progress=T)

landuse.grid # Here many buffers include areas that had not been mapped. Consider removing them, as, without a full mapping, land use and cover in these areas cannot be thrutfully assessed.



### To calculate edge density (ED) and number of patches (NP), I recommend using a binary raster so that it's clear what exactly are patches and what exactly are edges. Here we'll consider forest (class 11 in our example) as patches.

reclass <- matrix(ncol=3, nrow=3) # reclassification matrix
reclass[1,] <- c(0,2.10,1) # classes smaller than 11 will be "1"
reclass[2,] <- c(10.9,11.1,2) # class 11 (forest) will be "2"
reclass[3,] <- c(11.9,99,1) # all greater than 11 will be "2"

map.bin <- reclassify(map.ras, rcl=reclass)


config.real <- calculate.config(map=map.bin, sites=sites.real, buffers=c(42,666,2020), patch.class=2, site.name.slot="codigo_sit", print.progress=T, calculate.ED=T, calculate.NP=T)
config.real

config.grid <- calculate.config(map=map.bin, sites=sites.grid, buffers=c(42,666,2020), patch.class=2, site.name.slot="ID", print.progress=T, calculate.ED=T, calculate.NP=T) 
config.grid


### So, if this works, now you have a hopefully user-friendly way to rapidly calculate some basic landscape metrics. :-)


### Comparing with scale_sample from landscapemetrics:

config.real.2 <- scale_sample(landscape=map.bin, y=sites.real, shape="circle", size=666, max_size=2020, progress=T, what = c("lsm_l_ed", "lsm_c_np"))
config.real.2
View(config.real.2)

### The results are similar but shown in a different (more flexible) way. In addition, the buffer sizes appear to be calculated automatically; but I haven't looked for other ways of doing this using landscapemetrics.


