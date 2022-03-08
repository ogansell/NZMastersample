#' @import sf
#' @import data.table

# Code specifically for assigning Halton Index to points.
# This is equal sized boxes code!
getIndividualBoxIndices <- function(input, J = NULL, bb, base = c(2,3), seed = c(587, 8750, 3061), s1 = 0:1, s2 = 0:2)
{
	scale.bas <- bb[,2] - bb[,1]
	shift.bas <- bb[,1]
	
	if(is.null(J)) J <- round(log((bb[,2] - bb[,1])/100)/log(base), 0)

	B <- prod(base^J)
	
	if(any(class(input) == "sf"))
	{
		dat <- data.table(st_coordinates(input))
	} else if(ncol(coordinates(input[1,])) == 2){
		# Assume sp input:
		dat <- data.table(coordinates(input))
		setnames(dat, c("X", "Y"))
	}else{
		dat <- data.table(input$X, input$Y)
		setnames(dat, c("X", "Y"))
	}
	Bx <- base[1]^J[1]
	By <- base[2]^J[2]
	# Shift points to the lower box coordinate.
	dat <- dat[,c("x", "y") := list((X - shift.bas[1])/scale.bas[1], (Y - shift.bas[2])/scale.bas[2])]
	# Put it in 0-1 terms, with an adjustment to deal with "zeros" and machine precision.
	dat <- dat[,c("Ax", "Ay") := list(floor((x + 2*.Machine$double.eps)*Bx), 
		floor((y + 2*.Machine$double.eps)*By))]
		
	# This code seems to be doing something wrong. Switched to a more explicit merge above instead...
	dat[, c("Ax.new", "Ay.new") := list(permutation_B(Ji = J[1], bi = 2, si = s1)[Ax+1], permutation_B(Ji = J[2], bi = 3, si = s2)[Ay+1])] # Add one since indices include 0

	haltonIndex <- SolveCongruence(cbind(dat$Ax.new, dat$Ay.new), base = c(2,3), J = J)

	# Adjust everything for the Master Sample random seed.
	a1 <- seed[1:2] %% base^J
	boxInit <- SolveCongruence(matrix(a1, ncol = 2), base = c(2,3), J = J)
	
	# Adjusted index:
	haltonIndex <- ifelse(haltonIndex < boxInit, B + (haltonIndex - boxInit), haltonIndex - boxInit)
	# Return the Halton Index for all "pts" in dat that are passed to this function.
	return(haltonIndex)
}

# Subset Code:
subset_dist <- function(n, xi, index, slice = 2)
{
	K <- n %% slice

	# If they don't have an index then randomly remove points...
	if(is.null(index))
	{
		index <- sample.int(n, n, replace = FALSE)
	}
	
	out <- vector("integer", length = n)

	if(K == 0){
		out <- (1:n > floor(n/slice))*1 + (1:n > floor(2*n/slice))*(slice == 3)
		return(out)
	}else if(K == 1)
	{
		max.i <- which.max(index)
		out[-max.i] <- (1:(n-1) > floor((n - 1)/slice))*1 + (1:(n - 1) > floor(2*(n - 1)/slice))*(slice == 3)
		out[max.i] <- -1
		return(out)
	}else if(K == 2){
		maxi.1 <- max(index)
		maxi.2 <- max(index[index != maxi.1])
		rm.max <- which(index %in% c(maxi.1, maxi.2))
		out[-(rm.max)] <- (1:(n - 2) > floor((n - 2)/slice))*1 + (1:(n - 2) > floor(2*(n - 2)/slice))*(slice == 3)
		out[rm.max] <- -1
	}
	return(out)
}

	
# Permutation:
permutation_B <- function(Ji = 3, bi = 2, si = c(1, 0))
{
	if(Ji == 1) return(si)
	#Permute Base bi
	I <- si	
	for(k in 1:(Ji - 1))
	{
		v <- c()
		for(j in 1:bi^k)
		{
			v <- c(v, rep(I[j], bi) + bi^k*si)
		}
	I <- v
	}
	return(I)
}


##################################
# MAIN MASTER SAMPLE CODE:
##################################
#' @export
getHipSample <- function(X, Y, index = NULL, N = NULL, bb,  base = c(2,3), seed = c(587, 8750), quiet = TRUE, Ps1 = 0:1, Ps2 = 0:2, hipS1 = 0:1, hipS2 = 0:2)
{
	if(is.null(N))
	{
		N <- 50
		print("Making 50 sites.")
	}
	
	scale.bas <- bb[,2] - bb[,1]
	shift.bas <- bb[,1]
	 
	dat <- data.table(X, Y)
	if(!is.null(index)) dat$index <- index

	# Shift points to the lower box coordinate.
	dat <- dat[,c("x", "y") := list((X - shift.bas[1])/scale.bas[1], (Y - shift.bas[2])/scale.bas[2])]
	
	dat[, ID := 1:.N]
	
	# Choosing boxes slightly as close to 100 m^2 as possible but maybe a
	# bit on the bigger side.
	J.index <- ceiling(log((bb[,2] - bb[,1])/100)/log(base))

	#Produces individual indices
	if(is.null(index)){
	# Put it in 0-1 terms, with an adjustment to deal with "zeros" and machine precision.
	Bx <- base[1]^J.index[1]
	By <- base[2]^J.index[2]

	dat <- dat[,c("Ix", "Iy") := list(floor((x + 2*.Machine$double.eps)*Bx), 
		floor((y + 2*.Machine$double.eps)*By))]

		# Permute the boxes for the new method of adding "extra" randomness.
		dat[, c("Ix.new", "Iy.new") := list(permutation_B(Ji = J[1], bi = 2, si = Ps1)[Ix+1], permutation_B(Ji = J[2], bi = 3, si = Ps2)[Iy+1])] # Add one since indices include 0

		haltonIndex <- SolveCongruence(cbind(dat$Ix.new, dat$Iy.new), base = c(2,3), J = J.index)

		dat[, index := haltonIndex]

		# Adjust everything for the Master Sample random seed.
		a1 <- seed %% base^J.index
		boxInit <- SolveCongruence(matrix(a1, ncol = 2), base = c(2,3), J = J.index)
		dat[, index := ifelse(index < boxInit, prod(base^J.index) + (index - boxInit), index - boxInit)]
	}
	
	dat[, c("Ax", "Ay") := list(0, 0)]
	
	# Before we do HIP I want to subset for all boxes that might have more than one observation...
	# The reason to do this is that we are considering the box itself as the sample unit to "partition" and not
	# the actual point. We can add it back later...
	dat.rm <- dat[duplicated(index)]
	dat <- dat[!duplicated(index)]
	dat <- doHIP(dat, n = N, quiet = quiet, base = base)	# This is the iterative while loop now in a pretty function :)

	J <- attributes(dat)$J
	B <- prod(base^J)
	
	# HIP specific permutatations that are random each time
	# s1 <- sample(0:1, 2, replace = FALSE)
	# s2 <- sample(0:2, 3, replace = FALSE)
	s1 <- hipS1
	s2 <- hipS2
	
	# Add permutation from HIP Paper to add more "randomness" to the seed selection.
	dat[, c("Ax.new", "Ay.new") := list(permutation_B(Ji = J[1], bi = 2, si = s1)[Ax+1], permutation_B(Ji = J[2], bi = 3, si = s2)[Ay+1])] # Add one since indices include 0
	
	# Use the new permutation to solve the congruence.
	dat[, "HIPIndex" := SolveCongruence(cbind(Ax.new, Ay.new), base = c(2,3), J = J)]
	
	#Find initial box.
	# boxSeed <- seed %% base^J
	# boxSeed <- SolveCongruence(matrix(boxSeed, ncol = 2, nrow = 1), base = c(2,3), J = J)
	boxSeed <- dat[which.min(index),]$HIPIndex
	dat[, "HIPIndex" := ifelse(HIPIndex < boxSeed, B + (HIPIndex - boxSeed), HIPIndex - boxSeed)]
			
	tmp <- dat[, c("index", colnames(dat)[!names(dat) %in% names(dat.rm)]), with = FALSE]
	
	dat.tmp <- merge(dat.rm, tmp, by = "index", all.x = TRUE, all.y = FALSE)
	dat <- rbind(dat, dat.tmp)
	setorder(dat, ID)
	setattr(dat, "J", J)
	setattr(dat, "B", B)
	setattr(dat, "Permutation1", s1)
	setattr(dat, "Permutation2", s2)
	
	# Return the Halton Index for all "pts" in dat that are passed to this function.
	return(dat)
}

#' @export
doHIP <- function(dat2, n = NULL, quiet = quiet, base = base)
{
	if(is.null(n)) n = nrow(dat2)*2

	grp.size <- dat2[,.N, by = c("Ax", "Ay")]

	i <- 0
	j <- 0
	while(all(grp.size$N > 1) & prod(base^c(i,j)) < n)
	{
		setorder(dat2, x)
		
		dat2[, "Ax2" := subset_dist(.N, x, index, 2), by = c("Ax","Ay")]
		dat2 <- dat2[Ax2 != -1]		
		dat2[, "Ax" := Ax*2 + Ax2]

		grp.size <- dat2[,.N, by = c("Ax", "Ay")]
		i <- i + 1
		if(all(grp.size$N > 2) & base[1]^i > base[2]^j){	# Don't divide by 3 every time so we can achieve reasonably square divisions. Need to play with this.
			setorder(dat2, y)
			dat2[, "Ay2" := subset_dist(.N, y, index, 3), by = c("Ax","Ay")]
			dat2 <- dat2[Ay2 != -1]
			dat2[, "Ay" := Ay*3 + Ay2]			
			j <- j + 1
		}
		grp.size <- dat2[,.N, by = c("Ax", "Ay")]
		
		if(quiet == FALSE){
			print(paste0("J1 = ", i, " and J2 = ", j))
			print(min(grp.size$N))
		}
	}

	J <- c(i,j)
	setattr(dat2, "J", J)
	return(dat2)
}

#' @export
# Get Correct Number of Samples:
getSamples <- function(dat, n)
{

	B <- max(dat$HIPIndex) + 1 
	dat[ , samplei := (rank(index) - 1)*B + HIPIndex, by = "HIPIndex"]
	return(dat[samplei < n, ])
}

#' @export
# Get Correct Number of Samples:
getOverSamples <- function(dat, n, ni)
{

	B <- max(dat$HIPIndex) + 1 
	dat[ , sampleRank := (rank(index) - 1), by = "HIPIndex"]
	return(dat[HIPIndex < n & sampleRank < ni, ])
}

#' @export
# Subset Code:
subset_boxes <- function(subdat, boxes, slice = 2)
{
	bb.new <- list()
	for(i in sort(unique(subdat$Ax)))
	{
		for(j in sort(unique(subdat$Ay)))
		{
			bb.i <- boxes[[paste(c(i, j), collapse = " ")]]
			tmp <- subdat[Ax == i & Ay == j]
			tmp <- setorder(tmp, -index)			
			if(slice == 2)
			{
				miss <- nrow(tmp) %% 2

				if(miss >= 1 ) tmp <- tmp[-(1:miss)]
				md <- median(tmp$x)
				
				bb1 <- bb.i
				bb2 <- bb.i
				bb1[1,2] <- md
				bb2[1,1] <- md
				bb.new[[paste(i*2, j, collapse = " ")]] <- bb1
				bb.new[[paste(i*2 + 1, j, collapse = " ")]] <- bb2
			}
			if(slice == 3)
			{
				miss <- nrow(tmp) %% 3
				if(miss >= 1 ) tmp <- tmp[-(1:miss)]
				md <- quantile(tmp$y, c(1/3, 2/3))
				
				bb1 <- bb.i
				bb2 <- bb.i
				bb3 <- bb.i
				
				bb1[2,2] <- md[1]
				bb2[2,1] <- md[1]
				bb2[2,2] <- md[2]
				bb3[2,1] <- md[2]

				bb.new[[paste(i, j*3, collapse = " ")]] <- bb1
				bb.new[[paste(i, j*3 + 1, collapse = " ")]] <- bb2
				bb.new[[paste(i, j*3 + 2, collapse = " ")]] <- bb3
				
			}
		}
	}
	return(bb.new)
}

# This function returns a bunch of polygons that represent the HIP splits of the population.
# It really just does HIP again and could be combined with doHIP in the future with just the "boxes" 
# added and create Polygon = TRUE or something.
#' @export
getHIPBoxes <- function(hip, bb, n)
{
	boxes <- list("0 0" = bb)

	dat2 <- hip[,.(x = X, y = Y, index, Ax = 0, Ay = 0)]
	
	grp.size <- dat2[,.N, by = c("Ax", "Ay")]

	i <- 0
	j <- 0
	while(all(grp.size$N > 1) & length(boxes) < n)
	{
		setorder(dat2, x)
		boxes <- subset_boxes(dat2, boxes, slice = 2)
		
		dat2[, "Ax2" := subset_dist(.N, x, index, 2), by = c("Ax","Ay")]
		dat2 <- dat2[Ax2 != -1]	
		dat2[, "Ax" := Ax*2 + Ax2]

		grp.size <- dat2[,.N, by = c("Ax", "Ay")]
		i <- i + 1
		if(all(grp.size$N > 2) & base[1]^i > base[2]^j & length(boxes) < n){	# Don't divide by 3 every time so we can achieve reasonably square divisions. Need to play with this.		
			setorder(dat2, y)
			boxes <- subset_boxes(dat2, boxes, slice = 3)

			dat2[, "Ay2" := subset_dist(.N, y, index, 3), by = c("Ax","Ay")]
			dat2 <- dat2[Ay2 != -1]
			dat2[, "Ay" := Ay*3 + Ay2]			
			j <- j + 1
		}
		grp.size <- dat2[,.N, by = c("Ax", "Ay")]
		
	}

	polys <- lapply(1:length(boxes), FUN = function(x){as(extent(as.matrix(boxes[[x]])), 'SpatialPolygons')})
	for(i in 1:length(polys))
	{
		polys[[i]] <- polys[[i]]@polygons[[1]]
		slot(polys[[i]], "ID") <- names(boxes)[i]
	}
	
	polys <- SpatialPolygonsDataFrame(SpatialPolygons(polys), data.frame(ID = names(boxes), row.names = names(boxes)))
	return(polys)
}

#' @export
# Create the Sample Raster:
getPolyRast <- function(shp, smp, bb, J = c(5,3), base = c(2,3))
{
	scale.bas <- bb[,2] - bb[,1]
	shift.bas <- bb[,1]
	box.size <- scale.bas/base^J

	pts <- list()
	for(i in 1:nrow(smp))
	{
		# Find the box, clip the shape and get the dang sample!
		x <- c(smp[i]$X - box.size[1]/2, smp[i]$X + box.size[1]/2)
		y <- c(smp[i]$Y - box.size[2]/2, smp[i]$Y + box.size[2]/2)		
		b_poly <- as(extent(c(x,y)), "SpatialPolygons")
		if(any(class(shp) %in% "sf"))
		{
			rv <- st_intersection(shp, st_set_crs(st_as_sf(b_poly), st_crs(shp)))	
			rv <- as_Spatial(rv)
		}else{
			proj4string(b_poly) <- proj4string(shp)
			rv <- gIntersection(shp, b_poly, byid = T)
		}
		pts[[i]] <- sampRiv(n = 1, x = rv, seed = smp[i]$Ax.new)
	}
	pts <- do.call("rbind", pts)
	pts <- SpatialPointsDataFrame(pts, data = smp)
	return(pts)
}

#' @export
# Great now how to do we turn a sample into an actual point on a river feature???!!!!
# Sampling a linear feature:
sampRiv <- function(n = 10, x, seed = 0)
{
	cc <- coordinates(x)
	cc <- do.call("rbind", cc)
	cc.mat <- as.matrix(do.call("rbind", cc))
	lengths = LineLength(cc.mat, longlat = FALSE, sum = FALSE)
	csl = c(0, cumsum(lengths))
	maxl = csl[length(csl)]

	pts = HaltonSeq(seed, 2, n)* maxl
	int = findInterval(pts, csl, all.inside = TRUE)
	where = (pts - csl[int])/diff(csl)[int]
	xy = cc.mat[int, , drop = FALSE] + where * (cc.mat[int + 1, , drop = FALSE] - cc.mat[int, , drop = FALSE])
	SpatialPoints(xy, proj4string = CRS(proj4string(x)))
}

#########################
# Test Spatial Balance:
#########################
#' @export
getBalance <- function(hip, n)
{
	smp <- getSamples(hip, n)
	smp.pts <- getPolyRast(banks, smp, bb, J)
	pts <- SpatialPoints(cbind(hip$X, hip$Y)) 
	vp <- SDraw::voronoi.polygons(smp.pts)
	proj4string(pts) <- proj4string(vp)
	Nvp <- over(vp, pts, fn = length)
	return(var(Nvp))
}

#' @export
# Get samples with extra per hip box.
getOverSamples <- function(dat, n, ni)
{
  
  B <- max(dat$HIPIndex) + 1 
  dat[ , sampleRank := (rank(index) - 1), by = "HIPIndex"]
  return(dat[HIPIndex < n & sampleRank < ni, ])
}

#' @export	
getFWSeed <- function(island = "South")
{
	if(island == "South"){
		seed <-  c(7815, 699)
		s1 <- c(0,1)
		s2 <- c(2,0,1)		
	}
	if(island == "North"){
		seed <- c(601, 5024)
		s1 <- c(0,1)
		s2 <- c(2,1,0)
	}
	return(list(seed = seed, s1 = s1, s2 = s2))
}

#' @export
#Halton Sequence:
getStratumPermutation <-function(island = "South", stratum = "1")
{
  if(island == "South"){
    if(stratum == as.character(stratum[1])) return(list(s1 = c(1,0), s2 = c(2,1,0)))
    if(stratum == as.character(stratum[2])) return(list(s1 = c(1,0), s2 = c(0,2,1)))
    if(stratum == as.character(stratum[3])) return(list(s1 = c(1,0), s2 = c(2,0,1)))
    if(stratum == as.character(stratum[4])) return(list(s1 = c(1,0), s2 = c(2,0,1)))
  }
  if(island == "North"){
    if(stratum == as.character(stratum[1])) return(list(s1 = c(0,1), s2 = c(1,0,2)))
    if(stratum == as.character(stratum[2])) return(list(s1 = c(0,1), s2 = c(0,2,1)))
    if(stratum == as.character(stratum[3])) return(list(s1 = c(0,1), s2 = c(0,1,2)))	
    if(stratum == as.character(stratum[4])) return(list(s1 = c(0,1), s2 = c(0,1,2)))
  }
  return("ERROR")
}

