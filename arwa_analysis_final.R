#analysis of arctic warbler and other Alaska breeding bird geolocator data
#conducted by Evan Adams
#based on an analysis framework by Michael Hallworth

library(GeoLocTools)
#setupGeolocation()

#library(slider)
library(zoo)
library(FLightR)
library(BAStag)
library(SGAT)
library(GeoLight)
library(TwGeos)
library(probGLS)
library(ggmap)
library(lubridate)
library(fields)
library(stringr)
library(reshape2)
library(MASS)
library(maptools)
library(raster)
library(sp)
library(rgeos)
library(geosphere)
library(rgdal)

data(wrld_simpl)

#homemade twilightedit fn to remove outlier twilights
twilightEdit2 <- function (twilights, offset = 17, window = 4, outlier.mins = 45, 
                           stationary.mins = 15, zlim = c(0, 64), plot = T) 
{
  day <- twilights$Twilight
  hour <- hourOffset(as.hour(twilights$Twilight), offset)
  sunr <- which(twilights$Rise)
  suns <- which(!twilights$Rise)
  fnc <- function(x) {
    if(nrow(x) > 2){
      ind0 <- abs(x[(window/2) + 1, 1] - x[-((window/2) + 1), 
                                           1]) > median(diff(as.numeric(day[suns]))) * ((window/2) + 
                                                                                          1)
      # if (any(ind0)) 
      #   x <- x[-((window/2) + 1), ][-which(ind0), ]
      if (nrow(x) < window/2) {
        out <- cbind(x[(window/2 + 1), 1], FALSE, FALSE, 
                     x[(window/2 + 1), 1])
      } else {
        diffr <- abs(x[(window/2) + 1, 2] - median(x[-((window/2) + 
                                                         1), 2])) * 60
        if (diffr >= outlier.mins & all(dist(x[-((window/2) + 
                                                 1), 2]) * 60 <= stationary.mins)) {
          out <- cbind(x[(window/2) + 1, 1] + (median(x[c(window/2, 
                                                          (window/2) + 2), 2]) - x[(window/2) + 1, 2]) * 
                         60 * 60, FALSE, TRUE, x[(window/2) + 1, 1])
        }
        if (diffr >= outlier.mins & !all(dist(x[-((window/2) + 
                                                  1), 2]) * 60 <= stationary.mins)) {
          out <- cbind(x[(window/2) + 1, 1], TRUE, FALSE, 
                       x[(window/2) + 1, 1])
        }
        if (diffr < outlier.mins) 
          out <- cbind(x[(window/2) + 1, 1], FALSE, FALSE, 
                       x[(window/2) + 1, 1])
      }
    } else{
      
      out <- c(NA,NA,NA,NA)
    }
    return(out)
  }
  sunrT <- rollapply(cbind(day, hour)[sunr, ], width = window + 
                       1, FUN = fnc, fill = FALSE, by.column = F, partial = FALSE)
  sunsT <- rollapply(cbind(day, hour)[suns, ], width = window + 
                       1, FUN = fnc, fill = FALSE, by.column = F, partial = FALSE)
  sunrT[which(sunrT[, 1] == 0), c(1, 4)] <- day[sunr][which(sunrT[, 
                                                                  1] == 0)]
  sunsT[which(sunsT[, 1] == 0), c(1, 4)] <- day[suns][which(sunsT[, 
                                                                  1] == 0)]
  out <- data.frame(Twilight = as.POSIXct(c(sunrT[, 1], sunsT[, 
                                                              1]), origin = "1970-01-01", tz = "GMT"), Rise = c(rep(TRUE, 
                                                                                                                    nrow(sunrT)), rep(FALSE, nrow(sunsT))), Deleted = ifelse(c(sunrT[, 
                                                                                                                                                                                     2], sunsT[, 2]) == 1, TRUE, FALSE), Edited = ifelse(c(sunrT[, 
                                                                                                                                                                                                                                                 3], sunsT[, 3]) == 1, TRUE, FALSE), Twilight0 = as.POSIXct(c(sunrT[, 
                                                                                                                                                                                                                                                                                                                    4], sunsT[, 4]), origin = "1970-01-01", tz = "GMT"))
  out <- out[order(out[, 1]), ]
  rownames(out) <- 1:nrow(out)
  if (plot) {
    day0 <- out$Twilight0
    hour0 <- hourOffset(as.hour(out$Twilight0), offset)
    day <- out$Twilight
    hour <- hourOffset(as.hour(out$Twilight), offset)
    plot(day0, hour0, type = "n", xlab = "Date", ylab = "Hour")
    points(day[!out$Deleted], hour[!out$Deleted], pch = 16, 
           cex = 0.5, col = ifelse(out$Rise[!out$Deleted], "firebrick", 
                                   "cornflowerblue"))
    arrows(day0[out$Edited], hour0[out$Edited], day[out$Edited], 
           hour[out$Edited], length = 0.1)
    points(day0[out$Deleted | out$Edited], hour0[out$Deleted | 
                                                   out$Edited], pch = 16, col = "grey50")
    points(day[out$Edited], hour[out$Edited], pch = 16, col = ifelse(out$Rise[out$Edited], 
                                                                     "firebrick", "cornflowerblue"))
    points(day0[out$Deleted], hour0[out$Deleted], pch = "X")
  }
  out
}

##############ESTIMATING DENALI SONGBIRD POSITIONS USING GEOLOCATORS######################

#pull in the recovery data

recovs <- read.csv('geolocator data/tag_recov_update.csv')

#pull in the range maps

WGS84<-"+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"

arwa <- raster::shapefile('P:/Denali NP_critical connections/geolocator data/range_maps/arctic_warbler_range.shp')
blpw <- raster::shapefile('P:/Denali NP_critical connections/geolocator data/range_maps/blackpoll_warbler_range.shp')
fosp <- raster::shapefile('P:/Denali NP_critical connections/geolocator data/range_maps/fox_sparrow_range.shp')
gcth <- raster::shapefile('P:/Denali NP_critical connections/geolocator data/range_maps/graycheeked_thrush_range.shp')
heth <- raster::shapefile('P:/Denali NP_critical connections/geolocator data/range_maps/hermit_thrush_range.shp')
swth <- raster::shapefile('P:/Denali NP_critical connections/geolocator data/range_maps/swainsons_thrush_range.shp')
wiwa <- raster::shapefile('P:/Denali NP_critical connections/geolocator data/range_maps/wilsons_warbler_range.shp')

crs(arwa)<-WGS84
crs(blpw)<-WGS84
crs(fosp)<-WGS84
crs(gcth)<-WGS84
crs(swth)<-WGS84
crs(wiwa)<-WGS84

#look at the original BAStag data to make sure that they are all fine

dir <- 'P:/Denali NP_critical connections/geolocator data/processed data'

files <- list.files(dir)
files <- files[grep('.Rdata', files)]

tag.name <- str_split_fixed(files, '_', n = 5)

rawdat.list <- list()

for(i in 1:length(files)){
  
  rawdat.list[[i]] <- mget(load(paste0(dir, '/', files[i])))
  
}#i

#stitch together the data when there are multiple files due to processing issues


ids <- data.frame(ID = c('1381-88032', '1381-88010', '1381-88006', '1381-88017', '1381-88024',
                         '1381-88009', '1381-88004', '1381-88004', '1381-88004', '1381-88004',
                         '1381-88004','1381-88004', '1381-88004', '1381-88004', '1381-88022',
                         '1381-88029', '1381-88066', '1381-88066', '1381-88090', '1760-53551',
                         '2790-30314', '1780-53921', '1760-53502_2', '1760-53502_2', '1760-53518_2',
                         '1760-53518_2', '1780-53957', '1780-53957', '2850-21906', '2850-21906',
                         '1760-53526', '1760-53518_1', '1760-53502_1', '1760-53562', '1381-88074',
                         '1381-88049', '1760-53542', '0991-58917', '1760-53520', '2790-30313',
                         '1760-53506', '1760-53590', '1381-88005', '1381-88020', '0991-58902',
                         '0991-58926', '1381-88052'))

## note from Emily: 1760-53502 & 1760-53518 both were BLPW that were initially tagged in 2016 (geo recovered in 2017; tagged in 2018; geo recovered in 2019) 

recov.id <- merge(ids, recovs, by.x = 'ID', by.y = 'ID', all.x = TRUE)

unique.id <- unique(ids$ID)

unique.recov <- merge(data.frame(unique.id), recovs, by.x = 'unique.id', by.y = 'ID', all.x = TRUE)

unique.recov$DDate.UTC <- as.POSIXct(unique.recov$Start.Date, '%m/%d/%Y', tz = "US/Alaska")



#combine data sets for birds with more than one file

which(ids$ID %in% unique.recov$unique.id[4])

#1381-88004 (4)  7  8  9 10 11 12 13 14

rd4 <- list(twl = rbind(rawdat.list[[7]]$twl, rawdat.list[[8]]$twl, rawdat.list[[9]]$twl, rawdat.list[[10]]$twl, rawdat.list[[11]]$twl, rawdat.list[[12]]$twl, rawdat.list[[13]]$twl, rawdat.list[[14]]$twl))

# 1381-88066 (17) 17 18

rd17 <- list(twl = rbind(rawdat.list[[17]]$twl, rawdat.list[[18]]$twl))

# 1760-53502_2 (21)  23 24

rd21 <- list(twl = rbind(rawdat.list[[23]]$twl, rawdat.list[[24]]$twl))

#1760-53518_2 (24) 25 26

rd24 <- list(twl = rbind(rawdat.list[[25]]$twl, rawdat.list[[26]]$twl))

#1780-53957 (32) 27 28

rd32 <- list(twl = rbind(rawdat.list[[27]]$twl, rawdat.list[[28]]$twl))

#2850-21906 (35) 29 30

rd35 <- list(twl = rbind(rawdat.list[[29]]$twl, rawdat.list[[30]]$twl))

rawdat.final <- list()

to.replace <- c(4, 17, 21, 24, 32, 35)

for(i in which(!(1:nrow(unique.recov) %in% to.replace))){
  
  rawdat.final[[i]] <- rawdat.list[[which(ids$ID %in% unique.recov$unique.id[i])]]
  
}

rawdat.final[[4]] <- rd4
rawdat.final[[17]] <- rd17
rawdat.final[[21]] <- rd21
rawdat.final[[24]] <- rd24
rawdat.final[[32]] <- rd32
rawdat.final[[35]] <- rd35

twl <- rawdat.final

#use SGAT to build movement models

#first make sure that all the twilight data use the appropriate interval (should 120 s)

for(i in 1:length(twl)){
  dat <- twl[[i]]$twl
  twl[[i]]<-twilightAdjust(twilights = dat, interval=60)
}

#remove obvious outliers where twilight times shifted from the baseline considerably by more than 60 mins

for(i in 1:length(twl)){
if(nrow(twl[[i]]) > 20){
twl[[i]] <- twilightEdit2(twl[[i]], outlier.mins = 45, window = 4)}
}


#create a list of calibration dates

cal.per <- list()
for(i in 1:length(twl)){
  cal.per[[i]] <- data.frame(
    calibration.start = as.POSIXct(str_sub(unique.recov$DDate.UTC[i] + 1*(60*60*24), end = -1)),
    calibration.stop = as.POSIXct(paste0(unique.recov$D.Year[[i]], '-07-31')),
    recap.date = as.POSIXct(unique.recov$R.Date_L[i], format = '%m/%d/%Y'),
    lat = unique.recov$D.Lat[i], lon = unique.recov$D.Lon[i])
}

#adjust this for one bird that was caught very late and we need extra time to calibrate the date

cal.per[[31]]$calibration.stop <- as.POSIXct('2018-08-10')

#extract twilight data from the calibration dates

cal.data <- list()
for(i in 1:length(twl)){
  cal.data[[i]]<-subset(twl[[i]],twl[[i]]$Twilight>=cal.per[[i]]$calibration.start & twl[[i]]$Twilight<=cal.per[[i]]$calibration.stop)
}

#now calibrate the twilight data 

sun<-z<-zenith0<-zenith1<-twl_t<-twl_dev<-alpha<-fitml<- list()

for(i in 1:length(twl)){
  if(nrow(twl[[i]]) > 50){
    sun[[i]]<- solar(cal.data[[i]][, 1])
    z[[i]]<- refracted(zenith(sun[[i]], unique.recov$D.Lon[i], unique.recov$D.Lat[i]))
    twl_t[[i]]   <- twilight(cal.data[[i]][,1], unique.recov$D.Lon[i], unique.recov$D.Lat[i], rise = cal.data[[i]][,2], zenith = max(z[[i]])+0.1)
    twl_dev[[i]] <- ifelse(cal.data[[i]]$Rise, as.numeric(difftime(cal.data[[i]][,1], twl_t[[i]], units = "mins")),
                           as.numeric(difftime(twl_t[[i]], cal.data[[i]][,1], units = "mins")))
    twl_dev[[i]] <- ifelse(twl_dev[[i]] < 0, NA, twl_dev[[i]])
    fitml[[i]] <- fitdistr(twl_dev[[i]][which(is.na(twl_dev[[i]]) == FALSE)], "log-Normal")
    alpha[[i]] <- c(fitml[[i]]$estimate[1], fitml[[i]]$estimate[2])

  }
}

#31 has a bad calibration so let's replace it with a conspecific (25)

z[[31]] <- z[[25]]

#10 has a bad calibration, let's use the calibration from 13(another GCTH) as our best approximation

z[[10]] <- z[[13]]

#create preliminary location data

d.twl <- path <- na.path <- list()
zenith0 <- zenith1 <- z.adj <- rep(NA, length(twl))

for(i in 1:length(twl)){
  if(nrow(twl[[i]]) > 50){
    
    zenith0[i] <-quantile(z[[i]],prob=0.5)
    zenith1[i] <-quantile(z[[i]],prob=0.95)
    
    d.twl[[i]]<-subset(twl[[i]],twl[[i]]$Twilight>=cal.per[[i]]$calibration.start & !Deleted)
    path[[i]] <- thresholdPath(d.twl[[i]]$Twilight, d.twl[[i]]$Rise, zenith=zenith0[i],tol=0.15)
    na.path[[i]] <- thresholdLocation(d.twl[[i]]$Twilight, d.twl[[i]]$Rise, zenith=zenith0[i],tol=0.15)
    
  }
  
  
}


for(i in 1:length(twl)){
  
  if(nrow(twl[[i]]) > 50){
    png(file = paste0('P:/Denali NP_critical connections/geolocator data/prelim_maps_tol15/', unique.recov$unique.id[i],'prelim_map.png'))
    par(mfrow = c(1, 2), mar = c(2,4,1,1)+0.1)
    
    plot(path[[i]]$time, path[[i]]$x[, 2], type = "b", pch = 16, cex = 0.5, ylab = "Lat", xlab = '')
    abline(h = cal.per[[i]]$lat)
    abline(v = as.POSIXct("2011-09-23"),col="red",lty=2,lwd=1.5)
    abline(v = as.POSIXct("2012-03-20"),col="red",lty=2,lwd=1.5)
    
    plot(path[[i]]$x, type = "n")
    plot(wrld_simpl, add = T, col = "grey95")
    box()
    lines(path[[i]]$x, col = "blue")
    points(path[[i]]$x, pch = 16, cex = 0.5, col = "blue")
    dev.off()
  }
  
}


x0<-z0<-list()
for(i in 1:length(twl)){
  if(nrow(twl[[i]]) > 50){
    # Take the location estimates created above
    x0[[i]]<- path[[i]]$x
    # the model also needs the mid-points - generate those here
    z0[[i]]<- trackMidpts(x0[[i]])
  }
}

#setting initial values for the beta estimates in the movement model
beta <- c(0.7, 0.08)

#fix the values where we know the location of the individual (i.e., the breeding grounds right after capture)
fixedx<-list()
for(i in 1:length(twl)){
  if(nrow(twl[[i]]) > 50){
    fixedx[[i]]<- rep(F, nrow(x0[[i]]))
    fixedx[[i]][1:5] <- T
    #fixedx[[i]][(nrow(x0[[i]])-5):nrow(x0[[i]])] <-T
    
    x0[[i]][fixedx[[i]], 1] <- unique.recov$D.Lon[i]
    x0[[i]][fixedx[[i]], 2] <- unique.recov$D.Lat[i]
    
    z0[[i]] <- trackMidpts(x0[[i]]) # update z0 positions
  }
}

#now we can run an MCMC anlaysis to describe these movement behaviors

# set plot parameters

xlim.spp <- ylim.spp <- list()

flr.spp <- data.frame(Species = unique(recov.id$Species),
                      sppnum = 1:length(unique(recov.id$Species)),
                      gridx1 = c(-165, -165, -165, -165, -165, -130, -165),
                      gridx2 = c(-60, -40, -50, -50, -40, 70, -50),
                      gridy1 = c(10, -40, -15, 0, -20, -20, 0),
                      gridy2 = c(70, 70, 70, 70, 70, 70, 70),
                      M.mean = c(300, 300, 300, 300, 300, 300, 300),
                      M.sd = c(500, 500, 500, 500, 500, 500, 500),
                      dir1 = c(135, 135, 135, 135, 135, 225, 135),
                      dir2 = c(315, 315, 315, 315, 315, 45, 315))

for(i in 1:length(twl)){
  if(nrow(twl[[i]]) > 50){
    spp <- flr.spp[which(flr.spp$Species == unique.recov$Species[i]), ]
    
    xlim.spp[[i]] <- c(spp$gridx1, spp$gridx2)
    ylim.spp[[i]] <- c(spp$gridy1, spp$gridy2)
    
  }
}

opar <- par(mar=c(2,2,2,2)+0.1, mfrow=c(1,2))


#using a spatial mask to make overland locations more likley in the model output

earthseaMask <- function(xlim, ylim, n = 2, pacific=FALSE) {
  
  if (pacific) { wrld_simpl <- nowrapRecenter(wrld_simpl, avoidGEOS = TRUE)}
  
  # create empty raster with desired resolution
  r = raster(nrows = n * diff(ylim), ncols = n * diff(xlim), xmn = xlim[1],
             xmx = xlim[2], ymn = ylim[1], ymx = ylim[2], crs = proj4string(wrld_simpl))
  
  # create a raster for the stationary period, in this case by giving land a value of 1 and sea NA
  mask = cover(rasterize(elide(wrld_simpl, shift = c(-360, 0)), r, 1, silent = TRUE),
               rasterize(wrld_simpl, r, 1, silent = TRUE), 
               rasterize(elide(wrld_simpl,shift = c(360, 0)), r, 1, silent = TRUE))
  
  xbin = seq(xmin(mask),xmax(mask),length=ncol(mask)+1)
  ybin = seq(ymin(mask),ymax(mask),length=nrow(mask)+1)
  
  function(p) mask[cbind(.bincode(p[,2],ybin),.bincode(p[,1],xbin))]
}

xlim <- range(c(min(unlist(xlim.spp)), max(unlist(xlim.spp)) + c(-5,5)))
ylim <- range(c(min(unlist(ylim.spp)), max(unlist(ylim.spp)) + c(-5,5)))

mask <- earthseaMask(xlim, ylim, n = 1)

## Define the log prior for x and z
log.prior <- function(p) {
  f <- mask(p)
  ifelse(f | is.na(f), log(2), log(1))
}
 

x0<-z0<-model<-list()
for(i in (1:length(twl))[-27]){ 
  if(nrow(twl[[i]]) > 50){
    # Take the location estimates created above
    x0[[i]]<- path[[i]]$x
    # the model also needs the mid-points - generate those here
    z0[[i]]<- trackMidpts(x0[[i]])  
    # Define the threshold model - similar to above #
    model[[i]] <- thresholdModel(d.twl[[i]]$Twilight,d.twl[[i]]$Rise,
                                 twilight.model="ModifiedLogNormal",
                                 alpha=alpha[[i]],beta=beta,
                                 # Here is where we will set the constraints for land/spp dist
                                 logp.x = log.prior, logp.z = log.prior, 
                                 x0=x0[[i]],z0=z0[[i]],zenith=zenith1[i],fixedx=fixedx[[i]])
    
    # Here you need to set the first few locations as "known locations"-fixed locations. These are periods  
    # when you know the bird was at a specific location - capture site - when first deployed and captured. 
    model[[i]]$fixedx<-c(model[[i]]$fixedx,rep(FALSE,(dim(model[[i]]$x0)[1]-length(model[[i]]$fixedx))))
    #model[[i]]$fixedx[(length(model[[i]]$fixedx)-3):length(model[[i]]$fixedx)]<-TRUE
  }
}

proposal.x<-proposal.z<-list()
for(i in (1:length(twl))[-27]){
  if(nrow(twl[[i]]) > 50){
    proposal.x[[i]] <- mvnorm(S=diag(c(0.005,0.005)),n=nlocation(x0[[i]]))
    proposal.z[[i]] <- mvnorm(S=diag(c(0.005,0.005)),n=nlocation(z0[[i]]))
  }
}

fit <- list()
for(i in (1:length(twl))[-27]){
  if(nrow(twl[[i]]) > 50){
    #grab the appropriate species distribution map
    spp <- flr.spp[which(flr.spp$Species == unique.recov$Species[i]), ]
    #is.dist <- is.distspp[[spp$sppnum]]
    
    #fit the model
    fit[[i]] <- estelleMetropolis(model[[i]],proposal.x[[i]],proposal.z[[i]],iters=200,thin=25,chains=3)
  }
}

# Gather the data needed to make the plot 
s<-fixedz<-dt<-im<-list()

for(i in (1:length(twl))[-27]){
  if(nrow(twl[[i]]) > 50){
    png(file = paste0('P:/Denali NP_critical connections/geolocator data/movemod_rr_maps2_tol15/', unique.recov$unique.id[i],'prelim_map.png'))
    xlim <- xlim.spp[[i]]
    ylim <- ylim.spp[[i]]
    s[[i]] <- locationSummary(fit[[i]]$x,time=model[[i]]$time,collapse=F)
    
    fixedz[[i]] <- fixedx[[i]][-length(fixedx[[i]])] > 0 & fixedx[[i]][-length(fixedx[[i]])]==fixedx[[i]][-1]
    dt[[i]] <- ifelse(fixedz[[i]],0,model[[i]]$dt)
    im[[i]] <- locationImage(fit[[i]]$z,xlim=xlim,ylim=ylim,nx=4*diff(xlim),ny=4*diff(ylim),
                             weight=dt[[i]][1:dim(fit[[i]]$z[[1]])[1]],collapse=TRUE)
    
    # Generate the plot #
    plot(wrld_simpl,col= "grey90",border="grey10" , xlim = range(s[[i]][[1]][,"Lon.mean"]), 
         ylim = range(s[[i]][[1]][,"Lat.mean"]))
    plot(elide(wrld_simpl,shift=c(360,0)),xlim=xlim,ylim=ylim,add=TRUE,col= "grey90",border="grey10")
    image(im[[i]]$x,im[[i]]$y,im[[i]]$W,xlab="",ylab="",cex.axis=0.7, add = T, 
          col = c("transparent", rev(topo.colors(200))))
    plot(wrld_simpl, add = TRUE)
    for(k in 1:2) {
      lines(s[[i]][[k]][,"Lon.mean"],s[[i]][[k]][,"Lat.mean"],
            col=rgb(t(col2rgb(c("darkblue", "darkgreen")))/255,alpha=0.4))
    }
    abline(h=unique.recov$D.Lat[i],v=unique.recov$D.Lon[i])
    box()
    dev.off()
  }
}


#fine tuning the model
x0<-z0<-model<-list()
for(i in (1:length(twl))[-27]){
  if(nrow(twl[[i]]) > 50){
    x0[[i]] <- chainLast(fit[[i]]$x)
    z0[[i]] <- chainLast(fit[[i]]$z)
    # Define the threshold model - slimilar to above #
    model[[i]] <- thresholdModel(d.twl[[i]]$Twilight,d.twl[[i]]$Rise,
                                 twilight.model="ModifiedLogNormal",
                                 alpha=alpha[[i]],beta=beta,
                                 # Here is where we will set the constraints for land/spp dist
                                 logp.x = log.prior, logp.z = log.prior, 
                                 x0=x0[[i]],z0=z0[[i]],zenith=zenith1[i],fixedx=fixedx[[i]])
    
    # Here you need to set the first few locations as "known locations"-fixed locations. These are periods  
    # when you know the bird was at a specific location - capture site - when first deployed and captured. 
    model[[i]]$fixedx<-c(model[[i]]$fixedx,rep(FALSE,(dim(model[[i]]$x0[[2]])[1]-length(model[[i]]$fixedx))))
    #model[[i]]$fixedx[(length(model[[i]]$fixedx)-3):length(model[[i]]$fixedx)]<-TRUE
  }
}

proposal.x<-proposal.z<-list()
for(i in (1:length(twl))[-27]){
  if(nrow(twl[[i]]) > 50){
    #The observed data were parsed into transects then truncated to those less than 500m from the transect line. ){
    proposal.x[[i]] <- mvnorm(S=diag(c(0.005,0.005)),n=nlocation(x0[[i]]))
    proposal.z[[i]] <- mvnorm(S=diag(c(0.005,0.005)),n=nlocation(z0[[i]]))
  }
}

fit<-list()
for(i in (1:length(twl))[-27]){
  if(nrow(twl[[i]]) > 50){
    fit[[i]] <- estelleMetropolis(model[[i]],proposal.x[[i]],proposal.z[[i]],iters=5000,thin=25,chains=3)
  }
}

# Gather the data needed to make the plot 
s<-fixedz<-dt<-im<-list()

#plotting
for(i in (1:length(twl))[-27]){
  if(nrow(twl[[i]]) > 50){
    png(file = paste0('P:/Denali NP_critical connections/geolocator data/movemod_rr_maps2_tol15/', unique.recov$unique.id[i],'refined_map.png'))
    xlim <- xlim.spp[[i]]
    ylim <- ylim.spp[[i]]
    s[[i]] <- locationSummary(fit[[i]]$x,time=model[[i]]$time,collapse=F)
    
    fixedz[[i]] <- fixedx[[i]][-length(fixedx[[i]])] > 0 & fixedx[[i]][-length(fixedx[[i]])]==fixedx[[i]][-1]
    dt[[i]] <- ifelse(fixedz[[i]],0,model[[i]]$dt)
    im[[i]] <- locationImage(fit[[i]]$z,xlim=xlim,ylim=ylim,nx=4*diff(xlim),ny=4*diff(ylim),
                             weight=dt[[i]][1:dim(fit[[i]]$z[[1]])[1]],collapse=TRUE)
    
    # Generate the plot #
    plot(wrld_simpl,col= "grey90",border="grey10" , xlim = range(s[[i]][[1]][,"Lon.mean"]), 
         ylim = range(s[[i]][[1]][,"Lat.mean"]))
    plot(elide(wrld_simpl,shift=c(360,0)),xlim=xlim,ylim=ylim,add=TRUE,col= "grey90",border="grey10")
    image(im[[i]]$x,im[[i]]$y,im[[i]]$W,xlab="",ylab="",cex.axis=0.7, add = T, 
          col = c("transparent", rev(topo.colors(200))))
    plot(wrld_simpl, add = TRUE)
    for(k in 1:2) {
      lines(s[[i]][[k]][,"Lon.mean"],s[[i]][[k]][,"Lat.mean"],
            col=rgb(t(col2rgb(c("darkblue", "darkgreen")))/255,alpha=0.4))
    }
    abline(h=unique.recov$D.Lat[i],v=unique.recov$D.Lon[i])
    box()
    dev.off()
  }
}

#another5K iterations

x0<-z0<-model<-list()
for(i in (1:length(twl))[-27]){
  if(nrow(twl[[i]]) > 50){
    x0[[i]] <- chainLast(fit[[i]]$x)
    z0[[i]] <- chainLast(fit[[i]]$z)
    # Define the threshold model - slimilar to above #
    model[[i]] <- thresholdModel(d.twl[[i]]$Twilight,d.twl[[i]]$Rise,
                                 twilight.model="ModifiedLogNormal",
                                 alpha=alpha[[i]],beta=beta,
                                 # Here is where we will set the constraints for land/spp dist
                                 logp.x = log.prior, logp.z = log.prior, 
                                 x0=x0[[i]],z0=z0[[i]],zenith=zenith1[i],fixedx=fixedx[[i]])
    
    # Here you need to set the first few locations as "known locations"-fixed locations. These are periods  
    # when you know the bird was at a specific location - capture site - when first deployed and captured. 
    model[[i]]$fixedx<-c(model[[i]]$fixedx,rep(FALSE,(dim(model[[i]]$x0[[2]])[1]-length(model[[i]]$fixedx))))
    #model[[i]]$fixedx[(length(model[[i]]$fixedx)-3):length(model[[i]]$fixedx)]<-TRUE
  }
}

proposal.x<-proposal.z<-list()
for(i in (1:length(twl))[-27]){
  if(nrow(twl[[i]]) > 50){
    #The observed data were parsed into transects then truncated to those less than 500m from the transect line. ){
    proposal.x[[i]] <- mvnorm(S=diag(c(0.005,0.005)),n=nlocation(x0[[i]]))
    proposal.z[[i]] <- mvnorm(S=diag(c(0.005,0.005)),n=nlocation(z0[[i]]))
  }
}

fit<-list()
for(i in (1:length(twl))[-27]){
  if(nrow(twl[[i]]) > 50){
    fit[[i]] <- estelleMetropolis(model[[i]],proposal.x[[i]],proposal.z[[i]],iters=5000,thin=25,chains=3)
  }
}

#assessing model convergence

for(i in (1:length(twl))[-27]){
  if(nrow(twl[[i]]) > 50){
  opar <- par(mfrow = c(2, 1), mar = c(3, 5, 2, 1) + 0.1)
  matplot(t(fit[[i]]$x[[1]][!fixedx[[i]], 1, ]), type = "l", lty = 1, col = "dodgerblue", ylab = "Lon")
  matplot(t(fit[[i]]$x[[1]][!fixedx[[i]], 2, ]), type = "l", lty = 1, col = "firebrick", ylab = "Lat")
  par(opar)
  }
}

#appears to have good mixing though it can be hard to tell for some of these individuals

# Gather the data needed to make the plot 
s<-s.sum<-fixedz<-dt<-im<-r<-list()

xlimplot.spp <- xlim.spp
ylimplot.spp <- ylim.spp

#adjust ARWA ranges to deal with dateline issues

xlimplot.spp[[25]] <- c(-270, -130)
xlimplot.spp[[31]] <- c(-270, -130)
ylimplot.spp[[25]] <- c(-40, 80)
ylimplot.spp[[31]] <- c(-40, 80)

#plotting
for(i in (1:length(twl))[-27]){
  if(nrow(twl[[i]]) > 50){
    png(file = paste0('P:/Denali NP_critical connections/geolocator data/movemod_rr_maps2_tol15/', unique.recov$unique.id[i],'final_map.png'))
    xlim <- xlimplot.spp[[i]]
    ylim <- ylimplot.spp[[i]]
    s[[i]] <- locationSummary(fit[[i]]$x,time=model[[i]]$time,collapse=F)
    s.sum[[i]] <- locationSummary(fit[[i]]$x,time=model[[i]]$time,collapse=T)
    
    fixedz[[i]] <- fixedx[[i]][-length(fixedx[[i]])] > 0 & fixedx[[i]][-length(fixedx[[i]])]==fixedx[[i]][-1]
    dt[[i]] <- ifelse(fixedz[[i]],0,model[[i]]$dt)
    im[[i]] <- locationImage(fit[[i]]$z,xlim=xlim,ylim=ylim,nx=4*diff(xlim),ny=4*diff(ylim),
                             weight=dt[[i]][1:dim(fit[[i]]$z[[1]])[1]],collapse=TRUE)
    
    # Generate the plot #
    plot(wrld_simpl,col= "grey90",border="grey10" , xlim = range(s[[i]][[1]][,"Lon.mean"]), 
         ylim = range(s[[i]][[1]][,"Lat.mean"]))
    plot(elide(wrld_simpl,shift=c(360,0)),xlim=xlim,ylim=ylim,add=TRUE,col= "grey90",border="grey10")
    image(im[[i]]$x,im[[i]]$y,im[[i]]$W,xlab="",ylab="",cex.axis=0.7, add = T, 
          col = c("transparent", rev(topo.colors(200))))
    plot(wrld_simpl, add = TRUE)
    for(k in 1:3) {
      lines(s[[i]][[k]][,"Lon.mean"],s[[i]][[k]][,"Lat.mean"],
            col=rgb(t(col2rgb(c("darkblue", "darkgreen")))/255,alpha=0.4))
    }
    abline(h=unique.recov$D.Lat[i],v=unique.recov$D.Lon[i])
    box()
    dev.off()
    dput(s[[i]], file = paste0('P:/Denali NP_critical connections/geolocator data/final_results_output_tol15/', unique.recov$Species[i], unique.recov$unique.id[i],'_locs.out'))
    dput(s.sum[[i]], file = paste0('P:/Denali NP_critical connections/geolocator data/final_results_output_tol15/', unique.recov$Species[i], unique.recov$unique.id[i],'_locsum.out'))
    dput(im[[i]], file = paste0('P:/Denali NP_critical connections/geolocator data/final_results_output_tol15/', unique.recov$Species[i], unique.recov$unique.id[i],'_image.out'))
    #dput(fit[[i]], file = paste0('P:/Denali NP_critical connections/geolocator data/final_results_output_tol15/', unique.recov$Species[i], unique.recov$unique.id[i],'_fit.out'))
    
    # #convert to a shapefile
    # tmp.sp <- SpatialPointsDataFrame(cbind(s.sum[[i]]$Lon.mean, s.sum[[i]]$Lat.mean), s.sum[[i]])
    # writeOGR(tmp.sp, dsn = 'P:/Denali NP_critical connections/geolocator data/final_results_shapefiles_tol15', layer = paste0(unique.recov$Species[i], unique.recov$unique.id[i], '_modeled_position'), driver = 'ESRI Shapefile', overwrite = TRUE)
    # 
    #convert to a shapefile
    tmp.sp <- SpatialPointsDataFrame(cbind(s.sum[[i]]$Lon.mean, s.sum[[i]]$Lat.mean), s.sum[[i]], proj4string = CRS("+init=epsg:4326"),)
    tmp.sp$Month <- month(tmp.sp$Time)
    tmp.sp$Year <- year(tmp.sp$Time)
    writeOGR(tmp.sp, dsn = 'P:/Denali NP_critical connections/geolocator data/final_results_shapefiles_tol15', layer = paste0(unique.recov$Species[i], unique.recov$unique.id[i], '_modeled_position'), driver = 'ESRI Shapefile', overwrite = TRUE)
    
    #convert an image file into a raster
    
    w <- matrix(im[[i]]$W, ncol = length(im[[i]]$x) - 1, nrow = length(im[[i]]$y) - 1, byrow = TRUE)
    w <- w[nrow(w):1, ]
    
    r[[i]] <- raster(w, xmn = xlimplot.spp[[i]][1], xmx = xlimplot.spp[[i]][2],
                     ymn = ylimplot.spp[[i]][1], ymx = ylimplot.spp[[i]][2])
    
    crs(r[[i]]) <- CRS("+init=epsg:4326")
    
    writeRaster(r[[i]], paste0('P:/Denali NP_critical connections/geolocator data/final_results_raster_tol15/', unique.recov$Species[i], unique.recov$unique.id[i], '_modeled_space_use_new'), format = 'GTiff', prj = TRUE,  overwrite = TRUE)
    
    
  }
}


save.image('denali_geo_112020_tol15_update.RData')
