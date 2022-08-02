#results summary
require(data.table)
require(tidyverse)
require(lubridate)
require(changepoint)
require(amt)
require(sp)
require(rgdal)
require(sf)
require(SGAT)
require(scales)
require(grid)
require(raster)

arwa1 <- raster('P:/Denali NP_critical connections/geolocator data/arwa_revisions/ARWA1760-53520_modeled_space_use_new.tif')
arwa2 <- raster('P:/Denali NP_critical connections/geolocator data/arwa_revisions/ARWA1780-53921_modeled_space_use_new.tif')

plot(arwa1)
plot(arwa2)



#load the modeled data

load('denali_geo_060321_arwa_revis.RData')

arwa.im1 <- list(im[[25]], im[[31]])

#make some summary figures of lat/lon over time by species

track.sum <- rbindlist(s.sum, idcol = TRUE)

#combine with the unique recoveries data table

unique.recov$.id <- as.numeric(row.names(unique.recov))
track.join <- left_join(track.sum, unique.recov, by = '.id')
track.join$yday <- yday(track.join$Time)
track.join$diffdeploy <- track.join$Time - track.join$DDate.UTC

track.join$ref.date <- as.POSIXct(paste0(6, '/', 1, '/', track.join$D.Year), format = '%m/%d/%Y')

track.join$diff_ref <- track.join$Time - track.join$ref.date

#plot shifts in latitude by individual and grouped by species

p <- ggplot(track.join, aes(x = diff_ref/24, y = Lat.mean, ymin = `Lat.2.5%`, ymax =`Lat.97.5%`))
p <- p + geom_line(aes(color = band..)) + geom_ribbon(aes(fill = band..)) + facet_grid(~Species)
p <- p + theme_bw() + ylab('Mean Latitude') + xlab('Days Since June 1') + guides(color = FALSE, fill = FALSE)
p <- p + theme(axis.text.x = element_text(angle = 0))

plot(p)


#make a table of arrival/departure times from breeding/wintering grounds

#lets us changepoint analysis to help find those transitions

cps.lat <- list()

for(i in 1:length(s.sum)){
  
  if(length(s.sum[[i]]$Lat.mean) > 0)
  cps.lat[[i]] <- cpt.mean(s.sum[[i]]$Lat.mean, method = 'SegNeigh', Q = 5, pen.value = 0.9, penalty = 'Manual', test.stat = 'CUSUM')
  
}

cpts(cps.lat[[25]])
plot(cps.lat[[25]])

cpts(cps.lat[[31]])
plot(cps.lat[[31]])

cps.lon <- list()

for(i in 1:length(s.sum)){
  
  if(length(s.sum[[i]]$Lon.mean) > 0)
    cps.lon[[i]] <- cpt.mean(s.sum[[i]]$Lon.mean, method = 'SegNeigh', Q = 5, pen.value = 0.9, penalty = 'Manual', test.stat = 'CUSUM')
  
}

cpts(cps.lon[[25]])
plot(cps.lon[[25]])

cpts(cps.lon[[31]])
plot(cps.lon[[31]])

#use the change points to determine when they are leaving the breeding/wintering grounds

#look up the lat changepoint dates

lat.cds <- list()

for(i in 1:length(cps.lat)){

  if(length(s.sum[[i]]$Lon.mean) > 0)  
  lat.cds[[i]] <- s.sum[[i]]$Time[cpts(cps.lat[[i]])]
  
}

plot(cps.lat[[25]])
plot(cps.lat[[31]])


#write the .csv file of unique recovs

#I'm going to have to compile the results by hand so I can filter out spurious crap still

write.csv(unique.recov, 'denali_recovs.csv')

#okay, pulling this all together by hand...

#look at denali_recovs_with_arrival_departure_dates if you want to see what this looks like

sea.chg <- read.csv('denali_recovs_with_arrival_departure_dates.csv')

sea.chg$BG_Dep <- as.POSIXct(as.character(sea.chg$BG_Dep))
sea.chg$BG_Dep2 <- as.POSIXct(as.character(sea.chg$BG_Dep2))
sea.chg$WG_Dep <- as.POSIXct(as.character(sea.chg$WG_Dep))
sea.chg$WG_Dep2 <- as.POSIXct(as.character(sea.chg$WG_Dep2))

sea.chg$BG_Dep_yday <- yday(sea.chg$BG_Dep)
sea.chg$WG_Dep_yday <- yday(sea.chg$WG_Dep)

#the relationships between breeding ground and wintering ground departure and species 

p <- ggplot(sea.chg, aes(x = Species))
p <- p + geom_point(aes(y = sea.chg$BG_Dep_yday), color = 'Dark Orange', size = 6)
p <- p + geom_point(aes(y = sea.chg$WG_Dep_yday), color = 'Dark Green', size = 6)
p <- p + theme_bw() + ylab('Ordinal Date of Departure') + coord_flip()

plot(p)

p <- ggplot(sea.chg, aes(xmin = WG_Dep_yday, xmax = BG_Dep_yday, y = Species))
p <- p + geom_()
p <- p + theme_bw() + ylab('Ordinal Date of Departure') 

plot(p)


#let's look at how daily movement distance and timing

xy <- s.sum[[25]][, c(4, 9)]
proj4string(xy) <- CRS('+init=epsg:4326')

  
  '+proj=aea +lat_1=55 +lat_2=65 +lat_0=50 +lon_0=-154 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m +no_defs'

track1 <- make_track(s.sum[[25]], `Lon.50%`, `Lat.50%`, Time)
step1 <- step_lengths(track1)
direl1 <- direction_rel(track1)
direabs1 <- direction_abs(track1)

track2 <- make_track(s.sum[[31]], `Lon.50%`, `Lat.50%`, Time)
step2 <- step_lengths(track2)
direl2 <- direction_rel(track2)
direabs2 <- direction_abs(track2)


out <- matrix(NA, ncol = 3, nrow = length(ids))

for(i in 1:length(ids)){
  track1 <- make_track(dat1@data[which(dat1@data$bird_id == ids[i]),], x, y, datetime)
  out[i, 1] <- mean(step_lengths(track1), na.rm = TRUE)
  out[i, 2] <- mean(direction_rel(track1), na.rm = TRUE)
  out[i, 3] <- mean(direction_abs(track1, full_circle = TRUE), na.rm = TRUE)
}


#convert points to tracklines

arwa1.pts <- s.sum[[25]]
arwa2.pts <- s.sum[[31]]

WGS84_lon_wrap <- "+proj=longlat +datum=WGS84 +lon_wrap=180"

arwa1.pts$Lon.mean2 <- ifelse(arwa1.pts$Lon.mean < -180, 180 - abs(arwa1.pts$Lon.mean + 180), arwa1.pts$Lon.mean)
arwa2.pts$Lon.mean2 <- ifelse(arwa2.pts$Lon.mean < -180, 180 - abs(arwa2.pts$Lon.mean + 180), arwa2.pts$Lon.mean)

arwa1 <- as.matrix(arwa1.pts[, c('Lon.mean2', 'Lat.mean')])
arwa1_lines <- st_sfc(st_linestring(arwa1), crs=WGS84_lon_wrap) %>%
  sf::st_wrap_dateline(options = c("WRAPDATELINE=YES", "DATELINEOFFSET=180"))
plot(arwa1_lines)
write_sf(arwa1_lines, file.path('P:/Denali NP_critical connections/geolocator data/arwa_revisions/', "arwa1_lines.shp"))

arwa2 <- as.matrix(arwa2.pts[, c('Lon.mean2', 'Lat.mean')])
arwa2_lines <- st_sfc(st_linestring(arwa2), crs=WGS84_lon_wrap) %>%
  sf::st_wrap_dateline(options = c("WRAPDATELINE=YES", "DATELINEOFFSET=180"))
plot(arwa2_lines)
write_sf(arwa2_lines, file.path('P:/Denali NP_critical connections/geolocator data/arwa_revisions/', "arwa2_lines.shp"))


plot(arwa1[,2])

#creating season specific slices from each space use raster

startDate <- "2016-11-15"
endDate   <- "2017-04-10"

rast.month <- slices(mcmc = fit[[25]], breaks = 'month', grid = r[[25]])

rast.mig <- slice(rast.month, c(4))

plot(rast.mig)

writeRaster(rast.mig, paste0('P:/Denali NP_critical connections/geolocator data/arwa_revisions/', unique.recov$Species[25], unique.recov$unique.id[25], '_modeled_space_use_mig'), format = 'GTiff', prj = TRUE,  overwrite = TRUE)

rast.win <- slice(rast.month, c(5, 6, 7, 8))

plot(rast.win)

writeRaster(rast.win, paste0('P:/Denali NP_critical connections/geolocator data/arwa_revisions/', unique.recov$Species[25], unique.recov$unique.id[25], '_modeled_space_use_win'), format = 'GTiff', prj = TRUE,  overwrite = TRUE)

rast.month <- slices(mcmc = fit[[31]], breaks = 'month', grid = r[[31]])

rast.mig <- slice(rast.month, c(2))

plot(rast.mig)

writeRaster(rast.mig, paste0('P:/Denali NP_critical connections/geolocator data/arwa_revisions/', unique.recov$Species[31], unique.recov$unique.id[31], '_modeled_space_use_mig'), format = 'GTiff', prj = TRUE,  overwrite = TRUE)

rast.win <- slice(rast.month, c(4,5,6,7))

plot(rast.win)

writeRaster(rast.win, paste0('P:/Denali NP_critical connections/geolocator data/arwa_revisions/', unique.recov$Species[31], unique.recov$unique.id[31], '_modeled_space_use_win'), format = 'GTiff', prj = TRUE,  overwrite = TRUE)

#making lat/lon change figures
s.sum[[25]]$otime <- s.sum[[25]]$Time + years(2)

#lon
p <- ggplot(s.sum[[25]], aes(x = otime, y = `Lon.50%`, ymin = `Lon.2.5%`, ymax = `Lon.97.5%`))
p <- p + geom_line() + geom_ribbon(alpha = 0.3, fill = 'green') + geom_hline(yintercept = -180, linetype = 'dashed')
p1 <- p + theme_bw() + xlab('') + ylab('Longitude') + scale_x_datetime(limits = c(min(s.sum[[25]]$otime), max(s.sum[[31]]$Time)), labels = date_format('%b'))

print(p1)

#ggsave(paste0('arwa_',  unique.recov$unique.id[25], 'long_change.png'), device = 'png', dpi = 300)

p <- ggplot(s.sum[[31]], aes(x = Time, y = `Lon.50%`, ymin = `Lon.2.5%`, ymax = `Lon.97.5%`))
p <- p + geom_line() + geom_ribbon(alpha = 0.3, fill = 'dark blue') + geom_hline(yintercept = -180, linetype = 'dashed')
p2 <- p + theme_bw() + xlab('Month') + ylab('Longitude') + scale_x_datetime(limits = c(min(s.sum[[25]]$otime), max(s.sum[[31]]$Time)), labels = date_format('%b'))

print(p2)

#ggsave(paste0('arwa_',  unique.recov$unique.id[31], 'long_change.png'), device = 'png', dpi = 300)

#put them on the same scale

pdf('arwa_lon_joint.pdf')
grid.newpage()
grid.draw(rbind(ggplotGrob(p1), ggplotGrob(p2), size = 'last'))

dev.off()

#lat
p <- ggplot(s.sum[[25]], aes(x = otime, y = `Lat.50%`, ymin = `Lat.2.5%`, ymax = `Lat.97.5%`))
p <- p + geom_line() + geom_ribbon(alpha = 0.3, fill = 'green') + geom_hline(yintercept = 0, linetype = 'dashed')
p <- p + geom_rect(aes(xmin = s.sum[[25]]$otime[145], xmax = s.sum[[25]]$otime[191], ymin = -Inf, ymax = Inf), fill = 'yellow', alpha = 0.007)
p <- p + geom_rect(aes(xmin = s.sum[[25]]$otime[501], xmax = s.sum[[25]]$otime[546], ymin = -Inf, ymax = Inf), fill = 'yellow', alpha = 0.007)
p1 <- p + theme_bw() + xlab('') + ylab('Latitude') + scale_x_datetime(limits = c(min(s.sum[[25]]$otime), max(s.sum[[31]]$Time)), labels = date_format('%b'))

print(p1)

#ggsave(paste0('arwa_',  unique.recov$unique.id[25], 'lat_change.png'), device = 'png', dpi = 300)

#old code for the line method
# p <- ggplot(s.sum[[31]], aes(x = Time, y = `Lat.50%`, ymin = `Lat.2.5%`, ymax = `Lat.97.5%`))
# p <- p + geom_line() + geom_ribbon(alpha = 0.3, fill = 'dark blue') + geom_hline(yintercept = 0, linetype = 'dashed')
# p <- p + geom_vline(xintercept = as.numeric(s.sum[[25]]$otime[145]), col = 'yellow', size = 2)
# p <- p + geom_vline(xintercept = as.numeric(s.sum[[25]]$otime[191]), col = 'yellow', size = 2)
# p <- p + geom_vline(xintercept = as.numeric(s.sum[[25]]$otime[501]), col = 'yellow', size = 2)
# p <- p + geom_vline(xintercept = as.numeric(s.sum[[25]]$otime[546]), col = 'yellow', size = 2)
# p2 <- p + theme_bw() + xlab('Month') + ylab('Latitude') + scale_x_datetime(limits = c(min(s.sum[[25]]$otime), max(s.sum[[31]]$Time)), labels = date_format('%b'))
# 
# print(p2)

p <- ggplot(s.sum[[31]], aes(x = Time, y = `Lat.50%`, ymin = `Lat.2.5%`, ymax = `Lat.97.5%`))
p <- p + geom_line() + geom_ribbon(alpha = 0.3, fill = 'dark blue') + geom_hline(yintercept = 0, linetype = 'dashed')
p <- p + geom_rect(aes(xmin = s.sum[[25]]$otime[145], xmax = s.sum[[25]]$otime[191], ymin = -Inf, ymax = Inf), fill = 'yellow', alpha = 0.007)
p <- p + geom_rect(aes(xmin = s.sum[[25]]$otime[501], xmax = s.sum[[25]]$otime[546], ymin = -Inf, ymax = Inf), fill = 'yellow', alpha = 0.007)
p2 <- p + theme_bw() + xlab('Month') + ylab('Latitude') + scale_x_datetime(limits = c(min(s.sum[[25]]$otime), max(s.sum[[31]]$Time)), labels = date_format('%b'))

print(p2)

#ggsave(paste0('arwa_',  unique.recov$unique.id[31], 'lat_change.png'), device = 'png', dpi = 300)

pdf('arwa_lat_joint.pdf')
grid.newpage()
grid.draw(rbind(ggplotGrob(p1), ggplotGrob(p2), size = 'last'))

dev.off()




#note: need to add period of no lat data to these graphs still
