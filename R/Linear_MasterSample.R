########################
#Master Sample code for linear features
#Paul van Dam-Bates & Ollie Gansell
########################
#Take in a polygon and generate points from the fw master sample for linear features
#Step 1: Determine Island that the polygon falls in
#Step 2: Use random seed from that island to start Halton Sequence
#Step 3: Output number of points required clipped for that region.

#' @import dplyr
#' @import rgeos
#' @import tidyverse
#' @import Rcpp
#' @import raster
#' @import stars
NULL

## usethis namespace: start
#' @useDynLib NZMastersample, .registration = TRUE
#' @importFrom Rcpp sourceCpp
## usethis namespace: end
NULL

#Rcpp::sourceCpp("R/BASFunctions.cpp")

# shp <-sample.frame
#
# stratum <-sample.frame$Stratum
#
# island <-"North"

# n = 10
#
# oversamples = 1

#' @name linear_masterSample
#' @title Generate sample points on linear features in New Zealand using a HIP based master sample
#' @description Generates HIP sample points in a specified sample frame based on New Zealand linear features mastersample. Users need to specify 'island' the sample frame is in and how many points are required.
#' @export
linear_masterSample <-function(island = "North", shp, stratum, oversamples = 1, n = 10){

  if (!island %in% c("South", "North"))
    return("Define the island please.")

  nztm <- getProj(island)


  if (class(shp)[1] == "sf")
  {
    shp <-  as_Spatial(shp)
  }

  if (proj4string(shp) != nztm)
    shp <- spTransform(shp, nztm)

  stratum <- unique(stratum)

  stratum <-as.character(stratum)

  #pts.sample = SpatialPointsDataFrame(data.frame(x = 0, y = 0),data.frame(x = 0, y = 0))[-1,]
  #proj4string(pts.sample) <-getProj(island = island)
  pts.sample <-NULL

for(j in unique(stratum))
{

  i = island

  shp <-st_as_sf(shp)

  dat1 <- filter(shp,Stratum == j)

  dat1 <-as_Spatial(dat1)

  dat1$present <- 1

  #What bounding box and seed?
  bb <- getBB(i)
  seed.list <- getFWSeed(i)
  hip.seed <- getStratumPermutation(island = i, stratum = j)

  # Choose optimal size to rasterize with the Halton Frame.
  J <- round(log((bb[,2] - bb[,1])/c(100,100))/log(c(2,3)),0)
  base <- c(2,3)
  r <- raster(extent(as.matrix(bb)), nrow=base[2]^J[2], ncol=base[1]^J[1])
  values(r) <- 0
  r <-raster::rasterize(dat1[,"present"],r,field = "present", background =-1,fun='last')

  halt.info <- data.table(raster::xyFromCell(r, cell = c(1:ncell(r))))
  setnames(halt.info, c("X", "Y"))
  halt.info[, "Present" := values(r)]	# Couldn't figure out an easier way to extract values.
  halt.bank <- halt.info[Present != -1]

  # Draw HIP Partitioning:
  index <- getIndividualBoxIndices(input = halt.bank, J = J, bb, base = c(2,3), seed = seed.list$seed, s1 = seed.list$s1, s2 = seed.list$s2)
  dat <- data.table(X = halt.bank$X, Y = halt.bank$Y, index = index)

  hip <- getHipSample(X = halt.bank$X, Y = halt.bank$Y, index = index, N = 400, bb = bb,  seed = seed.list$seed, hipS1 = hip.seed$s1, hipS2 = hip.seed$s2)
  smp <- getOverSamples(hip, n, oversamples)
  smp.pts1 <- getPolyRast(dat1, smp, bb, J)
  smp.pts1$Stratum <- j
  smp.pts1$Island <- i
  #smp.pts1$SiteOrder <- smp.pts1$sampleRank # I'm not sure how you're using this column but be careful with it....
  smp.pts1$Stratum <-j
  smp.pts1 <-st_as_sf(smp.pts1, coords = c("X","Y"))
  pts.sample <-rbind(pts.sample,smp.pts1)
  print(paste("Stratum=", j))
}
  return(pts.sample)
}


# linear_masterSample(island = "North", sample.frame,
#                        stratum = sample.frame$Stratum, oversamples = 1, n = 1)
