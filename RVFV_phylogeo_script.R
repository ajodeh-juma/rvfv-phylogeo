#!/usr/bin/env Rscript

#'@author John Juma
#'
#'
# clear the working environment and free up memory
rm(list=ls(all.names = TRUE))
gc()


# set options
r = getOption("repos")
r["CRAN"] = "http://cran.us.r-project.org"
options(repos = r)
options(timeout = 1000)

###########################################
######## load/install packages ############
###########################################

library(countrycode)
library(colorspace)
library(data.table)
library(diagram)
library(ggplot2)
library(ggplotify)
library(ggthemes)
library(lemon)
library(lubridate)
library(patchwork)
library(rgdal)
library(rstudioapi)
library(sf)
library(spData)
library(stringr)
library(terra)
library(scales)
library(seraphim)
library(showtext)
# additional fonts
# check the current search path for fonts
font_paths()
# List available font files in the search path
font_files()  
# Add desired font
font_add("Arial Narrow", "/System/Library/Fonts/Supplemental/Arial Narrow.ttf")
# list loaded fonts
font_families()

# get path of the script
if (rstudioapi::isAvailable()){
  if (require('rstudioapi') != TRUE){
    install.packages('rstudioapi', repos = r)
  } else {
    library(rstudioapi)
  }
  wd <- dirname(getActiveDocumentContext()$path)
} else {
  wd <- getwd()
}


setwd(wd)

# imports
source("mccTreeExtractions.r")
source("mccTree_R.r")

# directories
figuresDirectory <- file.path(wd, "figures")
if (!dir.exists(figuresDirectory)) dir.create(figuresDirectory, recursive = TRUE)

segments <- c("L", "M", "nss")
nberOfExtractionFiles <- 1000
mostRecentSamplingDates <- c(decimal_date(ymd("2022-09-12")), 
                             decimal_date(ymd("2022-09-12")), 
                             decimal_date(ymd("2022-08-25")))


# 1. Preparing all the different environmental factors (rasters) to test

template = raster("All_rasters_1/Croplands.asc")
template[!is.na(template[])] = 1; template[is.na(template[])] = 0
template = raster::aggregate(template, 10, fun=mean)
template[template[]<0.5] = NA; template[template[]>=0.5] = 0
rasterNames = list.files("All_rasters_1")
for (i in 1:length(rasterNames))
{
  r1 = raster(paste0("All_rasters_1/",rasterNames[i]))
  r2 = raster::aggregate(r1, 10, fun=mean); r2[is.na(template[])] = NA
  writeRaster(r2, paste0("All_rasters_2/",rasterNames[i]), overwrite=T)
}


# 2. Extracting spatio-temporal information embedded in 1000 posterior trees

for (h in 1:length(segments))
{
  localTreesDirectory <- paste0("RVFV-", segments[h], "_exts")
  burnIn <- 1001 # 10% of sampled trees
  randomSampling <- FALSE
  nberOfTreesToSample <- nberOfExtractionFiles
  mostRecentSamplingDatum <- mostRecentSamplingDates[h]
  coordinateAttributeName <- "location"
  nberOfCores <- 10
  
  if (dir.exists(localTreesDirectory)) {
    fns <- file.path(localTreesDirectory, list.files(path = localTreesDirectory,
                                                     pattern = "^TreeExtractions"))
    if (length(fns) != nberOfTreesToSample) {
      for (fn in fns) {
        file.remove(fn)
      }
      allTrees <- scan(paste0(segments[h],"-lineage-C.trees"), 
                       what="", sep="\n", 
                       quiet=TRUE)
      treeExtractions(localTreesDirectory, allTrees, burnIn, randomSampling, 
                      nberOfTreesToSample, mostRecentSamplingDatum, 
                      coordinateAttributeName, nberOfCores)
    }
  } else {
      allTrees <- scan(paste0(segments[h],"-lineage-C.trees"), 
                       what="", sep="\n", 
                       quiet=TRUE)
      treeExtractions(localTreesDirectory, allTrees, burnIn, randomSampling, 
                      nberOfTreesToSample, mostRecentSamplingDatum, 
                      coordinateAttributeName, nberOfCores)
    }
}

# area of study
xmin <- -17.520278 # 16.45148
xmax <- 63.47519
ymin <- -34.65302
ymax <- 37.34698

e_Figure <- extent(xmin, xmax, ymin, ymax)
country.codes <- c("AGO", "BDI", "BEN", "BFA", "BWA", "CAF", "CIV", "CMR", 
                   "COD", "COG", "CPV", "DJI", "DZA", "EGY", "ERI", "ESH",
                   "ETH", "GAB", "GHA", "GIN", "GMB", "GNB", "GNQ", "KEN",
                   "LBR", "LBY", "LSO", "MAR", "MDG", "MLI", "MOZ", "MRT",
                   "MUS", "MWI", "MYT", "NAM", "NER", "NGA", "RWA", "SDN",
                   "SEN", "SLE", "SOM", "SSD", "SWZ", "SYC", "TCD", "TGO",
                   "TUN", "TZA", "UGA", "ZAF", "ZMB", "ZWE", "SAU", "YEM",
                   "OMN")

country.names <- countrycode(country.codes,
                             origin = 'iso3c',
                             destination = 'country.name')
country.names <- append(country.names, "Somaliland")
spdf <- spData::world %>% 
  dplyr::filter((continent == 'Africa' ) | 
                  (continent == 'Asia' & name_long %in% 
                     c('Saudi Arabia', 'Yemen', 'Oman')))
sp <- as(st_geometry(spdf), "Spatial")
sp@proj4string <- CRS("+init=epsg:4326")

# subset on countries to label
spdf.df <- spData::world %>% 
  dplyr::filter(name_long %in% c("Tanzania", "Kenya", "Sudan", "South Africa", 
                                 "Rwanda", "Burundi", "Madagascar", 
                                 "Saudi Arabia", "Yemen", "Uganda", "Ethiopia", 
                                 "Zimbabwe"))
spdf.sp <- as(st_geometry(spdf.df), "Spatial")
spdf.sp@proj4string <- CRS("+init=epsg:4326")

gadmDirectory <- file.path(wd)
if (!dir.exists(gadmDirectory)) dir.create(gadmDirectory, recursive = TRUE)
borders <- terra::crop(geodata::gadm(country = country.codes, level = 0, 
                                     path = gadmDirectory), e_Figure)
writeVector(borders, file.path(gadmDirectory, "StudyArea.shp"), overwrite=TRUE)

# study area
# shp <- shapefile("StudyArea.shp")
# study_area <- crop(shp, extent(e_Figure))

# features
background <- terra::crop(raster(file.path(wd, "Natural_Earth/Gray_background.tif")), e_Figure)
background[background[]==106] <- NA
lakes <- terra::crop(terra::vect(x=file.path(wd, "Natural_Earth/"), layer="Natural_Earth_lakes"), e_Figure)
lakes <- sf::st_as_sf(lakes)
coastlines <- terra::crop(terra::vect(x=file.path(wd, "Natural_Earth/"), layer="Coastline_borders"), e_Figure)
coastlines <- sf::st_as_sf(coastlines)

pols_list <- list()
cols_list <- list()
end_years_list <- list()
start_year_list <- list()
mcc_list <- list()

for (h in 1:length(segments))
{
  # 5.1. Extracting spatio-temporal information embedded in the MCC tree
  
  mcc_tre <- readAnnotatedNexus(paste0(segments[h],"-lineage-C.mcc.tree"))
  mcc_tab <- mccTree_R(mcc_tre, mostRecentSamplingDates[h])
  write.csv(mcc_tab, paste0("RVFV-", segments[h],".csv"), row.names=F, quote=F)
  
  # 5.2. Estimating the HPD region for each time slice
  percentage <- 80
  precision <- 2
  prob = percentage/100
  
  
  minYears <- rep(NA, nberOfExtractionFiles)
  maxYears <- rep(NA, nberOfExtractionFiles)
  points <- c()
  for (i in 1:nberOfExtractionFiles)
  {
    localTreesDirectory <- paste0("RVFV-", segments[h], "_exts")
    tab <- read.csv(paste0(localTreesDirectory,"/TreeExtractions_",i,".csv"), head=T)
    minYears[i] <- min(tab[,"startYear"])
  }
  startDatum <- HDInterval::hdi(minYears)[1] # to get the lower bound of the 95% HPD interval for the starting year to consider

  polygons <- suppressWarnings(spreadGraphic2(localTreesDirectory, 
                                              nberOfExtractionFiles, 
                                              prob, startDatum, precision))

  # 5.3. Defining the different colour scales to use
  # cols <- gsub("FF","",rev(viridis::mako(151))[1:101])
  cols = gsub("FF","",rev(viridis::magma(151))[1:101])
  #cols <- colorRampPalette(brewer.pal(11,'RdYlGn'))(131)[21:121]
  minYear <- startDatum
  maxYear <- mostRecentSamplingDates[h]
  mcc <- read.csv(paste0("RVFV-",segments[h],".csv"), head=T)
  mcc <- mcc[order(mcc[,"endYear"],mcc[,"startYear"],decreasing=F),]
  startYear_index <- (((min(mcc[,"startYear"])-minYear)/(maxYear-minYear))*100)+1
  startYear_colour <- cols[startYear_index]
  endYears_indices <- (((mcc[,"endYear"]-minYear)/(maxYear-minYear))*100)+1
  endYears_indices[which(endYears_indices<1)] = NA
  endYears_colours <- cols[endYears_indices]
  polygons_colours <- rep(NA, length(polygons))
  
  pols <- list()
  for (i in 1:length(polygons))
  {
    date <- as.numeric(names(polygons[[i]]))
    polygon_index <- round((((date-minYear)/(maxYear-minYear))*100)+1)
    if (polygon_index >= 1)
    {
      polygons_colours[i] <- paste0(cols[polygon_index], "20")
      polygons[[i]]@proj4string <- CRS("+init=epsg:4326")
      pols[[i]] <- polygons[[i]]
      pols[[i]]$year <- date
    }
  }
  mcc_list[[h]] <- mcc
  cols_list[[h]] <- cols
  end_years_list[[h]] <- endYears_colours
  start_year_list[[h]] <- startYear_colour
  
  joined.pols <- do.call(bind, pols)
  joined.pols.df <- st_as_sf(joined.pols)
  pols_list[[h]] <- joined.pols.df

  # 5.4. Co-plotting the HPD regions and MCC tree
  
  pdf(file.path(figuresDirectory, paste0("RVFV-", segments[h], ".pdf")),
      width=7.5, height=6.0)
  par(mfrow=c(1,1), mar=c(0,0,0,0), oma=c(1.5,2.7,1.0,0), mgp=c(0,0.4,0), lwd=0.2, bty="o")
  plot(background, main="",  col="#EEEEE0", cex.main=1, bty="n", box=F, axes=F,
       legend=F, colNA="#FFFFFF")
  plot(sp, asp=1, add=T, lwd=0.5, border="#000000")
  text(spdf.sp, labels = spdf.df$name_long, col="#000000", cex=0.7)
  xs <- seq(-17, 63, b=20)
  axis(3, lwd=2,at=xs,lab=parse(text=degreeLabelsEW(xs)),tck=0.01,
       col="#FFFFFF", mgp=c(0,-1.3,0))
  ys <- pretty(par()$usr[3:4])
  axis(2, las=1,tck=0.01,lwd=2,at=ys,lab=parse(text=degreeLabelsNS(ys)),
       tck=0.01, col="#FFFFFF", mgp=c(0,-2.0,0))
  for (i in 1:length(polygons))
  {
    plot(polygons[[i]], axes=F, col=polygons_colours[i], add=T, border=NA)
  }
  for (i in 1:dim(mcc)[1])
  {
    curvedarrow(cbind(mcc[i,"startLon"],mcc[i,"startLat"]), cbind(mcc[i,"endLon"],mcc[i,"endLat"]), arr.length=0,
                arr.width=0, lwd=0.3, lty=1, lcol="#4D4D4D", arr.col=NA, arr.pos=F, curve=0.1, dr=NA, endhead=F)
  }
  for (i in dim(mcc)[1]:1)
  {
    if (i == 1) # for the root (only correct if the first line of the CSV file corresponds to one of the two most ancestral branches)
    {
      points(mcc[i,"startLon"], mcc[i,"startLat"], pch=16, col=startYear_colour, cex=0.7)
      points(mcc[i,"startLon"], mcc[i,"startLat"], pch=1, col="#4D4D4D", lwd=0.2, cex=0.7)
    }
    if (mcc[i,"node2"]%in%mcc[,"node1"]) # for internal nodes (dots)
    {
      points(mcc[i,"endLon"], mcc[i,"endLat"], pch=16, col=endYears_colours[i], cex=0.7)
      points(mcc[i,"endLon"], mcc[i,"endLat"], pch=1, col="#4D4D4D", lwd=0.2, cex=0.7)
    }	else	{	# for the tip nodes (squares)
      points(mcc[i,"endLon"], mcc[i,"endLat"], pch=15, col=endYears_colours[i], cex=0.6)
      points(mcc[i,"endLon"], mcc[i,"endLat"], pch=0, col="#4D4D4D", lwd=0.2, cex=0.6)
    }
  }
  rast <- raster(matrix(nrow=1, ncol=2))
  rast[1] <- startDatum
  rast[2] <- max(mcc[,"endYear"])

  plot(rast, legend.only=T, add=T, col=cols, legend.width=0.5,
       legend.shrink=0.3, smallplot=c(0.88,0.89,0.18,0.82), alpha=1.0,
       legend.args=list(text="", cex=0.9, line=0.5, col="#4D4D4D"),
       axis.args=list(col="#4D4D4D", cex.axis=0.6, col.axis="#4D4D4D",
                      lwd=0, lwd.tick=0.5, tck=-0.05, col.tick="#4D4D4D",
                      line=0, mgp=c(0,0.4,0)))
  dev.off()
}


# plot with ggplot2
plot_dispersal <- function(mcc, spdf, spdf.df, joined.pols.df,
                           startYearCol, endYearsCols, cols){
  p <- ggplot() +
    geom_sf(data = spdf, linewidth=0.2, fill="#EEEEE0", colour="#000000") +
    geom_sf(data = st_as_sf(lakes), colour="#D0D0D0") +
    geom_text(
      data = spdf.df,
      aes(label = name_long, geometry = geom),
      stat = "sf_coordinates",
      colour = "#000000", 
      fontface = "bold") +
    geom_sf(data = joined.pols.df,
            lwd = 0,
            aes(fill=year),
            alpha=0.2) +
    scale_fill_gradientn(
      name = "Year",
      colours = cols,
      breaks = seq(as.numeric(format(date_decimal(min(joined.pols.df$year)), "%Y")),
                   as.numeric(format(date_decimal(max(joined.pols.df$year)), "%Y")),
                   by=10)
    ) +
    geom_point(mcc, mapping = aes(x = endLon, y = endLat, colour=endYearsCols),
               size=1.0, shape=15, show.legend = F) +
    geom_point(mcc, mapping = aes(x = endLon, y = endLat, colour="#4D4D4D"),
               size=1.0, shape=1, show.legend = F) +
    geom_point(mcc[1,], mapping = aes(x = startLon, y = startLat, colour=startYearCol),
               size=1.0, shape=16, show.legend = F) +
    geom_point(mcc[1,], mapping = aes(x = startLon, y = startLat, colour="#4D4D4D"),
               size=1.0, shape=1, show.legend = F) +
    geom_curve(
      aes(x = startLon, y = startLat, xend = endLon, yend = endLat),
      arrow = arrow(
        length = unit(0.008, "npc"),
        type="closed"
      ),
      colour = endYearsCols,
      size = 0.1,
      angle = 90,
      data = mcc,
      curvature = 0.5,
      linewidth = 0.2,
      linetype = "solid"
    ) +
    ggtitle("", subtitle = "RVFV lineage C dispersal graph") +
    theme_map() +
    coord_sf(xlim = c(xmin, xmax), ylim = c(ymin, ymax), expand = FALSE) +
    xlab("Longitude") + 
    ylab("Latitude") +
    theme(
      text = element_text(family = "Arial Narrow"),
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      axis.text = element_blank(),
      axis.line.x = element_blank(),
      axis.line.y = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      plot.tag = element_text(size = 30, face = "bold", colour = "#000000"),
      legend.position = c(0.02, 0.015),
      legend.background = element_rect(color = "#FFFFFF",
                                       fill = "transparent",
                                       linetype = "blank"),
      legend.title = element_text(size = 18, face = "bold", colour = "#000000"),
      legend.text =  element_text(size = 13, face = "bold", colour = "#000000"),
      legend.key.height= unit(1.5, 'cm'),
      legend.key.width= unit(0.4, 'cm'),
      plot.subtitle=element_text(size=22, face = "bold", colour = "#000000"))
  return(p)
}


l.plot <- plot_dispersal(mcc = mcc_list[[1]],
                         spdf = spdf, 
                         spdf.df = spdf.df,
                         joined.pols.df = pols_list[[1]],
                         startYearCol = start_year_list[[1]],
                         endYearsCols = end_years_list[[1]],
                         cols = cols_list[[1]])
l.plot
m.plot <- plot_dispersal(mcc = mcc_list[[2]],
                         spdf = spdf, 
                         spdf.df = spdf.df,
                         joined.pols.df = pols_list[[2]],
                         startYearCol = start_year_list[[2]],
                         endYearsCols = end_years_list[[2]],
                         cols = cols_list[[2]])
m.plot

nss.plot <- plot_dispersal(mcc = mcc_list[[3]],
                           spdf = spdf,
                           spdf.df = spdf.df,
                           joined.pols.df = pols_list[[3]],
                           startYearCol = start_year_list[[3]],
                           endYearsCols = end_years_list[[3]],
                           cols = cols_list[[3]])
nss.plot




for (h in 1:length(segments))
{
  # . Estimate of several dispersal/epidemiological statistics
  
  nberOfExtractionFiles <- 1000
  timeSlices <- 200
  onlyTipBranches <- FALSE
  showingPlots <- FALSE
  nberOfCores <- 10
  slidingWindow <- 1
  outputName <- paste0("RVFV-",segments[h])
  
  localTreesDirectory <- paste0("RVFV-", segments[h], "_exts")
  
  stats <- file.path(wd, paste0(outputName,"_estimated_dispersal_statistics.txt"))
  if (!file.exists(stats)) {
    spreadStatistics(localTreesDirectory, nberOfExtractionFiles, timeSlices,
                     onlyTipBranches, showingPlots, outputName,
                     nberOfCores, slidingWindow)
  } else {
    stats = read.table(paste0(outputName, "_estimated_dispersal_statistics.txt"), header=T)
    v <- stats[,"mean_branch_dispersal_velocity"]
    D <- stats[,"original_diffusion_coefficient"]
    D <- D/365
    cat(paste0("Mean branch dispersal velocity:\t", round(mean(v),1)," km/year,\t 95% HPD = [",round(quantile(v,0.025),1),",",round(quantile(v,0.975),1),"]"), "\n")
    cat(paste0("Original diffusion coefficient:\t", round(mean(D),1)," km2/day,\t 95% HPD = [",round(quantile(D,0.025),1),",",round(quantile(D,0.975),1),"]"), "\n")

    v <- stats[,"weighted_branch_dispersal_velocity"]
    D <- stats[,"weighted_diffusion_coefficient"]
    cat(paste0("Weighted branch dispersal velocity:\t", round(mean(v),1)," km/year, \t 95% HPD = [",round(quantile(v,0.025),1),",",round(quantile(v,0.975),1),"]"), "\n")
    cat(paste0("Weighted diffusion coefficient:\t", round(mean(D))," km2/year, \t 95% HPD = [",round(quantile(D,0.025)),",",round(quantile(D,0.975)),"]"), "\n")
    
    wdlc = stats[,"weighted_branch_dispersal_velocity"]
    median = round(median(wdlc),1)
    HPD = round(HDInterval::hdi(wdlc)[1:2],1)
    cat(segments[h],": median WLDV = ",median,", 95% HPD = [",HPD[1],", ",HPD[2],"]","\n",sep="")
    
    
    wdc = stats[,"weighted_diffusion_coefficient"]
    median = round(median(wdc),0)
    HPD = round(HDInterval::hdi(wdc)[1:2],0)
    cat(segments[h],": median WDC = ",median,", 95% HPD = [",HPD[1],", ",HPD[2],"]","\n",sep="")
  }
}

# [1] "RVFV-L"
# Estimation of summary statistics
# Median value of mean branch dispersal velocity = 582.9766
# 95% HPD = [207.4312, 2040.729]
# Median value of weighted branch dispersal velocity = 132.5512
# 95% HPD = [96.25616, 168.6628]
# Median value of original diffusion coefficient = 220347
# 95% HPD = [65621.11, 956223.1]
# Median value of weighted diffusion coefficient = 55959.9
# 95% HPD = [40955.2, 75493.94]
# Building wavefront distance evolution graphs
# Building branch dispersal velocity evolution graphs
# [1] "RVFV-M"
# Estimation of summary statistics
# Median value of mean branch dispersal velocity = 554.4614
# 95% HPD = [246.8895, 1659.01]
# Median value of weighted branch dispersal velocity = 141.7737
# 95% HPD = [93.10701, 182.7575]
# Median value of original diffusion coefficient = 143192.1
# 95% HPD = [45570.89, 652299.3]
# Median value of weighted diffusion coefficient = 54644.68
# 95% HPD = [35615.68, 77580.15]
# Building wavefront distance evolution graphs
# Building branch dispersal velocity evolution graphs
# [1] "RVFV-nss"
# Estimation of summary statistics
# Median value of mean branch dispersal velocity = 401.3954
# 95% HPD = [173.3204, 980.1356]
# Median value of weighted branch dispersal velocity = 177.6489
# 95% HPD = [135.7197, 216.4322]
# Median value of original diffusion coefficient = 113701.7
# 95% HPD = [47657.36, 342418.3]
# Median value of weighted diffusion coefficient = 82247.92
# 95% HPD = [59839.77, 108917.1]
# Building wavefront distance evolution graphs
# Building branch dispersal velocity evolution graphs


# plot dispersal stats
plot_dispersal <- function(tab1, tab2, y, xlab, ylab, subtitle, col){
  
  df1 <- fread(tab1)
  df2 <- fread(tab2)
  
  df <- dplyr::left_join(df1, df2, by="time") %>% 
    dplyr::rename(lower_hpd="95%HPD_lower_value",
                  upper_hpd="95%HPD_higher_value")
  
  # convert to class Date
  df$date <- as.Date(as.character(df$time), "%Y")
  
  p <- ggplot(data = df, mapping = aes(x=date)) +
    geom_line(aes(y = .data[[y]]), colour = col, show.legend = T, alpha = 0.8) +
    geom_ribbon(aes(ymin = lower_hpd, ymax = upper_hpd), fill = col, alpha = 0.2) +
    scale_x_date(date_labels = "%Y", date_breaks = "8 year",
                 limits = c(as.Date(as.character(minYear), "%Y"),
                            as.Date(as.character(maxYear), "%Y"))) +
    ggtitle("", subtitle = subtitle) +
    labs(y = ylab, x = xlab) +
    theme_classic() +
    theme(
      text = element_text(family = "Arial Narrow"),
      plot.subtitle=element_text(size=16, face = "bold", colour = "#000000"),
      axis.title = element_text(size=15, face = "bold", colour = "#000000"),
      axis.text = element_text(size = 14, face = "bold", colour = "#000000"),
      axis.line.x = element_line(colour = 'black', linewidth = 0.5, linetype='solid'),
      axis.line.y = element_line(colour = 'black', linewidth = 0.5, linetype='solid')
    )
}

wavefront.plots <- list()
mean.velocity.plots <- list()
weighted.velocity.plots <- list()

for (h in 1:length(segments))
{
  # . Plot dispersal stats 
  tab1 <- file.path(wd, paste0(outputName, "_mean_patristic_wavefront_distance.txt"))
  tab2 <- file.path(wd, paste0(outputName, "_95%HPD_patristic_wavefront_distance.txt"))
  
  p2 <- plot_dispersal(tab1 = tab1,
                       tab2 = tab2,
                       y = "distance",
                       ylab = "Wavefront distance (km)",
                       xlab = "Time (year)",
                       subtitle = "Wavefront distance from epidemic origin",
                       col = "#9F2B68FF")
  p2 <- p2 + theme(
    axis.title = element_text(size = 15, face = "bold", colour = "#000000"),
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold", colour = "#000000"),
    # axis.title.x = element_blank(),
    plot.subtitle=element_text(size=16, face = "bold", colour = "#000000"))
  p2
  
  wavefront.plots[[h]] <- p2
  
  
  tab1 <- file.path(wd, paste0(outputName, "_mean_mean_branch_dispersal_velocity.txt"))
  tab2 <- file.path(wd, paste0(outputName, "_95%HPD_mean_branch_dispersal_velocity.txt"))
  
  p3 <- plot_dispersal(tab1 = tab1,
                       tab2 = tab2,
                       y = "velocity",
                       ylab = "Velocity (km/year)",
                       xlab = "Time (year)",
                       subtitle = "Evolution of mean branch dispersal velocity",
                       col = "#0000FF")
  p3 <- p3 + theme(
    axis.title = element_text(size = 15, face = "bold", colour = "#000000"),
    plot.title = element_text(hjust = 0.5, vjust = -35, size = 16, face = "bold", colour = "#000000"),
    # axis.title.x = element_blank(),
    plot.subtitle=element_text(size=16, face = "bold", colour = "#000000"))
  p3
  
  mean.velocity.plots[[h]] <- p3
  
  
  tab1 <- file.path(wd, paste0(outputName, "_mean_weighted_branch_dispersal_velocity.txt"))
  tab2 <- file.path(wd, paste0(outputName, "_95%HPD_weighted_branch_dispersal_velocity.txt"))
  
  p4 <- plot_dispersal(tab1 = tab1,
                       tab2 = tab2,
                       y = "velocity",
                       ylab = "Velocity (km/year)",
                       xlab = "Time (year)",
                       subtitle = "Weighted branch dispersal velocity",
                       col = "#0000FF")
  
  p4 <- p4 + theme(
    axis.title = element_text(size = 15, face = "bold", colour = "#000000"),
    plot.title = element_text(hjust = 0.5, vjust = -35, size = 16, face = "bold", colour = "#000000"),
    # axis.title.x = element_blank(),
    plot.subtitle=element_text(size=16, face = "bold", colour = "#000000"))
  p4
  
  weighted.velocity.plots[[h]] <- p4
}



# . plot skygrid population history of the virus
plot_pop <- function(file, age, lower_tmrca, upper_tmrca, youngest_date, 
                     subtitle, col, segment) {
  
  # Read model output from BEAST
  df <- read.csv(file, skip = 1, header = T, sep = "\t")
  df <- df %>% dplyr::rename(time=1,
                             mean=2,
                             median=3,
                             upper=4,
                             lower=5)
  base_n <- dim(df)[1]
  BaseDateVector <- rep(youngest_date, base_n) # base number from BEAST
  print(BaseDateVector)
  decimalDates <- (BaseDateVector - df$time)
  calendarDates = format(date_decimal(df$time), "%Y") # convert to calendar dates
  
  # add calender and decimal dates difference to dataframe
  dfCalendar = data.frame(df, calendarDates)
  dfCalendar <- data.frame(dfCalendar, decimalDates)
  
  # convert to class Date
  dfCalendar$calendarDates <- as.Date(dfCalendar$calendarDates, "%Y")
  
  # define TMRCA from BEAST tracer summary
  tmrca <- as.Date(age, "%Y")
  minTMRCA <- as.Date(lower_tmrca, "%Y")
  maxTMRCA <- as.Date(upper_tmrca, "%Y")
  
  tmrca_year <- as.numeric(format(as.POSIXct(tmrca, format = "%Y-%m-%d %H:%M:%S"), format = "%Y"))
  min_tmrca_year <- as.numeric(format(as.POSIXct(minTMRCA, format = "%Y-%m-%d %H:%M:%S"), format = "%Y"))
  max_tmrca_year <- as.numeric(format(as.POSIXct(maxTMRCA, format = "%Y-%m-%d %H:%M:%S"), format = "%Y"))
  
  
  # define range of median values
  min_lower <- min(dfCalendar$lower)
  max_upper <- max(dfCalendar$upper)
  
  print(dfCalendar)
  
  # plot
  
  # adjust these values according to your data
  n = length(seq(from=log(min_lower), to=log(max_upper)))
  by = 2
  breaks = seq(from=log(min_lower), to=log(max_upper), by=by)
  
  if (segment == 'L'){
    labels = c(0.1, 1, 10, 100, 1000, 10000)
  }
  if (segment == 'M'){
    labels = c(0.01, 0.1, 1, 10, 100, 1000, 10000)
  }
  if (segment == 'nss'){
    labels = c(0.1, 1, 10)
  }
  
  
  p <- ggplot(data = dfCalendar, mapping = aes(x=calendarDates)) +
    geom_line(aes(y = log(mean)), colour = col , alpha = 0.6) +
    geom_line(aes(y = log(median)), colour = col, alpha = 0.8, linetype="dashed") +
    geom_ribbon(aes(ymin = log(lower), ymax = log(upper)), fill = col, alpha = 0.2) +
    scale_x_date(date_labels = "%Y", date_breaks = "6 year") +
    scale_y_continuous(breaks = seq(from=log(min_lower), to=log(max_upper), by=by),
                       labels = labels) +
    labs(y=expression(paste(italic("N"[e]),tau," ","(years)")), x = "Year") +
    ggtitle("", subtitle = subtitle) +
    theme_classic() +
    theme(
      text = element_text(family = "Arial Narrow"),
      axis.title = element_text(size = 15, face = "bold", colour = "#000000"),
      axis.text = element_text(size = 14, face = "bold", colour = "#000000"),
      axis.line.x = element_line(colour = "#000000", linewidth =0.5, linetype='solid'),
      axis.line.y = element_line(colour = "#000000", linewidth =0.5, linetype='solid'),
      plot.subtitle=element_text(size=16, face = "bold", colour = "#000000")
    )
}

pop.L <- plot_pop(file = file.path(wd, "RVFV-L-skygrid.tsv"),
                  age = "1965.705",
                  lower_tmrca = "1952.6239",
                  upper_tmrca = "1974.9542",
                  youngest_date = 2022.695,
                  col = "#FFA500",
                  subtitle = "RVFV demographic history",
                  segment = "L")

pop.M <- plot_pop(file = file.path(wd, "RVFV-M-skygrid.tsv"),
                  age = "1965.8001",
                  lower_tmrca = "1949.1966",
                  upper_tmrca = "1975.7566",
                  youngest_date = 2022.695,
                  col = "#FFA500",
                  subtitle = "RVFV demographic history",
                  segment = "M")

pop.nss <- plot_pop(file = file.path(wd, "RVFV-nss-skygrid.tsv"),
                    age = "1971.9732",
                    lower_tmrca = "1962.758",
                    upper_tmrca = "1975.9996",
                    youngest_date = maxYear,
                    col = "#FFA500",
                    subtitle = "RVFV demographic history",
                    segment = "nss")



# plot diffusion coefficients
datalist <- NULL
for (i in 1:length(segments)){
  outputName <- paste0("RVFV-",segments[i])
  tab <- read.table(paste0(outputName, "_estimated_dispersal_statistics.txt"), header=T)
  tab <- tab[,c("weighted_branch_dispersal_velocity", "weighted_diffusion_coefficient")]
  tab <- tab %>% dplyr::rename("WLDC" = "weighted_branch_dispersal_velocity",
                               "WDC"="weighted_diffusion_coefficient")
  tab$Segment <- segments[i]
  print(head(tab, 10))
  datalist[[i]] <- tab
}

data = do.call(rbind, datalist)
long.data <- reshape2::melt(data)


plot_diffusion_coefficient <- function(data, column, subtitle, col){
  column <- ensym(column)
  median = round(median(data$value),1)
  HPD = round(HDInterval::hdi(data$value)[1:2],1)
  p <- ggplot(data, aes(x=value, fill = !!column)) +
    scale_fill_manual(name = "Segment", values = col) +
    geom_density(alpha=0.1, show.legend = T) +
    geom_vline(data = data,
               aes(xintercept = HPD[1]), colour = "grey50",
               linewidth=0.3, show.legend = FALSE, linetype="dashed") +
    geom_vline(data = data,
               aes(xintercept = HPD[2]), colour = "grey50",
               linewidth=0.3, show.legend = FALSE, linetype="dashed") +
    geom_vline(data = data,
               aes(xintercept = median), colour = "grey40",
               linewidth=0.3, show.legend = FALSE) +
    ggtitle("", subtitle = subtitle) +
    labs(y="Density", x="WDC") +
    theme_classic() +
    theme(
      text = element_text(family = "Arial Narrow"),
      axis.title = element_text(size=15, face = "bold", colour = "#000000"),
      axis.text = element_text(size=14, face = "bold", colour = "#000000"),
      axis.line.x = element_line(colour = "#000000", linewidth =0.5, linetype='solid'),
      axis.line.y = element_line(colour = "#000000", linewidth =0.5, linetype='solid'),
      plot.subtitle=element_text(size=16, face = "bold", colour = "#000000"),
      legend.position = "none",
      legend.background = element_blank(),
      legend.title = element_text(size = 14, face = "bold", colour = "#000000"),
      legend.text =  element_text(size = 12, face = "bold", colour = "#000000")
    )
  print(p)
  return(p)
}
dataWDClarge <- long.data[long.data$variable=='WDC' & long.data$Segment=='L',]
WDCsL <- plot_diffusion_coefficient(data = dataWDClarge, 
                           column = "variable", 
                           subtitle = "Weighted Diffusion Coefficient",
                           col = "yellow4")

dataWDCmedium <- long.data[long.data$variable=='WDC' & long.data$Segment=='M',]
WDCsM <- plot_diffusion_coefficient(data = dataWDCmedium, 
                           column = "variable", 
                           subtitle = "Weighted Diffusion Coefficient",
                           col = "yellow4")
dataWDCnss <- long.data[long.data$variable=='WDC' & long.data$Segment=='nss',]
WDCsNSS <- plot_diffusion_coefficient(data = dataWDCnss, 
                           column = "variable", 
                           subtitle = "Weighted Diffusion Coefficient",
                           col = "yellow4")

plots.L <- (l.plot | 
  ((wavefront.plots[[1]] / weighted.velocity.plots[[1]]) | 
     (WDCsL / pop.L))) +
  plot_layout(widths = c(24, 24),
              heights = c(9, 9),
              byrow = T) +
  plot_annotation(tag_levels = "A") & 
  theme(plot.tag = element_text(size = 30, face = "bold"))
plots.L

ggsave(file.path(figuresDirectory, "FigureS6.pdf"), plots.L,
       width = 50, height = 25, units = "cm", limitsize = FALSE,
       dpi = 300, bg="white", device = cairo_pdf)

ggsave(file.path(figuresDirectory, "FigureS6.png"), plots.L,
       width = 50, height = 25, units = "cm", limitsize = FALSE,
       dpi = 300, bg="white", device = "png")



plots.M <- (m.plot | 
  ((wavefront.plots[[2]] / weighted.velocity.plots[[2]]) |
     (WDCsM / pop.M))) + 
  plot_layout(widths = c(24, 24), 
              heights = c(9, 9), 
              byrow = T) +
  plot_annotation(tag_levels = "A") & 
  theme(plot.tag = element_text(size = 30, face = "bold"))

ggsave(file.path(figuresDirectory, "Figure3.pdf"), plots.M,
       width = 50, height = 25, units = "cm", limitsize = FALSE,
       dpi = 300, bg="white", device = cairo_pdf)

ggsave(file.path(figuresDirectory, "Figure3.png"), plots.M,
       width = 50, height = 25, units = "cm", limitsize = FALSE,
       dpi = 300, bg="white", device = "png")



plots.nss <- (nss.plot | 
              ((wavefront.plots[[3]] / weighted.velocity.plots[[3]]) |
              (WDCsNSS / pop.nss))) + 
  plot_layout(widths = c(24, 24), 
              heights = c(9, 9), 
              byrow = T) +
  plot_annotation(tag_levels = "A") & 
  theme(plot.tag = element_text(size = 30, face = "bold"))
plots.nss


ggsave(file.path(figuresDirectory, "FigureS7.pdf"), plots.nss,
       width = 50, height = 25, units = "cm", limitsize = FALSE,
       dpi = 300, bg="white", device = cairo_pdf)

ggsave(file.path(figuresDirectory, "FigureS7.png"), plots.nss,
       width = 50, height = 25, units = "cm", limitsize = FALSE,
       dpi = 300, bg="white", device = "png")