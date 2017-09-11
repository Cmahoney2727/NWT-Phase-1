# Predict area wide vegetation heights
#
# Author: craig.mahoney
# Created: 25/04/2014
# Last Revised: 31/07/2014
#
#	source("E:/NWT/ModelFramework/Scripts/R/PredictionskNN.R")
##########################################################################
rm(list = objects())

require(raster)
require(rgdal)
require(matlab)
require(itertools)
require(parallel)
require(foreach)
require(doSNOW)
require(kknn)

#####FUNCTIONS
kNNsdev <- function(X, Predictions, Training, Name = NULL){
					SDev <- rep(NA, length(Predictions[[X]]))
				for(i in 1:nrow(Predictions[X, "C"][[1]]))
					SDev[i] <- sd(Training[Predictions[X, "C"][[1]][i, ], Name])
			return(SDev)	
}

#####PREAMBLE
	ttt <- Sys.time()
	
Filters <- list.files("E:/NWT/NWT_Phase1/ModelFramework/Observation_Layers/", "Level_")[-5]
	
	Pth <- "E:/NWT/ModelFramework/Data_Layers/Composite/CompositeStack.tif"
		Rstack <- stack()
	for(bands in 1:8)
		Rstack <- stack(Rstack, raster(Pth, band = bands))
		names(Rstack) <- c("Band3", "Band4", "Band5", "CDEM", "CMI", "CTI", "EOSD", "SMI")
		
		Mulitplier <- 1500
		xdim <- ceiling( (extent(Rstack)@xmax - extent(Rstack)@xmin) / (res(Rstack)[1] * Mulitplier) )# - 1
		ydim <- ceiling( (extent(Rstack)@ymax - extent(Rstack)@ymin) / (res(Rstack)[1] * Mulitplier) )# - 1
	
	Grd <- GridTopology(c(extent(Rstack)@xmin, extent(Rstack)@ymin) + (0.5 * (res(Rstack) * Mulitplier)), res(Rstack) * Mulitplier, c(xdim, ydim))
	Sg <- SpatialGrid(Grd, proj4string = CRS(projection(Rstack)))
	Sp <- as(Sg, "SpatialPolygons")
	
		image(Rstack[["CTI"]], xlim = c(extent(Sg)@xmin, extent(Sg)@xmax), ylim = c(extent(Sg)@ymin, extent(Sg)@ymax))
		plot(Sp, add = TRUE)	
	
for(fff in 4){	#1:length(Filters)
for(Att in c("SHt"))	#c("SHt", "CC")
	Opth <- paste0("E:/NWT/NWT_Phase1/ModelFramework/kNN/Outputs/", Filters[fff])
	dir.create(Opth)
	dir.create(paste0(Opth, "/", Att))
	dir.create(paste0(Opth, "/", Att, "/Tiles"))
	dir.create(paste0(Opth, "/", Att, "/SDev"))

		Training <- get(load(paste("E:/NWT/ModelFramework/kNN/TrainedModels/", Filters[fff], "/", Att, "/TrainingData_", Att, ".dat", sep = "")))
		OptInfo <- get(load(paste("E:/NWT/ModelFramework/kNN/TrainedModels/", Filters[fff], "/", Att, "/Optimized_kNN.out", sep = "")))
	
	for(i in 38){	#1:length(Sp@polygons)
		P <- SpatialPolygons( list( Polygons( list( Polygon(Sp@polygons[[i]]@Polygons[[1]]@coords) ), ID = i) ), proj4string = CRS(projection(Rstack)))
		system.time(T <- crop(Rstack, P))
		system.time(Zmat <- getValues(T))
		Zmat <- round(Zmat, 4)
		#Remove non-vegetation land cover classes
			Zmat[Zmat[, "EOSD"] == 0, ] <- NA
			Zmat[Zmat[, "EOSD"] == 11, ] <- NA
			Zmat[Zmat[, "EOSD"] == 12, ] <- NA
			Zmat[Zmat[, "EOSD"] == 20, ] <- NA
			Zmat[Zmat[, "EOSD"] == 31, ] <- NA
			Zmat[Zmat[, "EOSD"] == 32, ] <- NA
			Zmat[Zmat[, "EOSD"] == 33, ] <- NA
			Zmat[Zmat[, "EOSD"] == 40, ] <- NA
			Zmat[Zmat[, "EOSD"] == 51, ] <- NA
			Zmat[Zmat[, "EOSD"] == 52, ] <- NA
			Zmat[Zmat[, "EOSD"] == 100, ] <- NA
			
				F <- function(X, Zmat){which(is.na(Zmat[, X]))}
			Rmvs <- lapply(1:8, F, Zmat)
				Rmvs <- unique(unlist(Rmvs))
		
			Ncores <- 5
			cl <- makeCluster(Ncores)
			registerDoSNOW(cl)	
		system.time(predkNN <- foreach(d = isplitRows(Zmat[-Rmvs, names(Training)[-1]], chunks = Ncores), .combine = rbind, .packages = c("kknn")) %dopar% { kknn(as.formula(paste("ShtGLAS ~", paste(names(Training)[-1], collapse = " + "))), train = Training, test = data.frame(d), k = OptInfo$k, kernel = OptInfo$kernal) })	#OptInfo$predictors
			print("Finished kNN predictions...")
		
		system.time(SD <- parLapply(cl, x = 1:Ncores, fun = kNNsdev, Predictions = predkNN, Training, Name = Att)	)
			stopCluster(cl)
		#SD <- sapply(X = 1:Ncores, kNNsdev, Predictions = predkNN, Name = "ShtGLAS")
		SDev <- rep(NA, (ncol(T) * nrow(T)) - length(Rmvs))
		for(i in 1:length(SD)){
			if(i == 1)
				i1 <- i
				i2 <- i1 + (length(SD[[i]]) - 1)
				SDev[i1:i2] <- SD[[i]]
				i1 <- i2 + 1
		}
		#Sdev	
		Smat <- rep(-1, ncol(T) * nrow(T))
			rrr <- seq(1, ncol(T) * nrow(T), 1)[-Rmvs]
		Smat[rrr] <- SDev	
		
		sp <- SpatialPointsDataFrame((expand.grid(x = seq(extent(T)@xmin, extent(T)@xmax - res(T)[1], res(T)[1]) + (0.5 * res(T)[1]), y = rev(seq(extent(T)@ymin, extent(T)@ymax - res(T)[2], res(T)[2]) + (0.5 * res(T)[1])))), data = data.frame(SDev = Smat), proj4string = CRS(projection(Rstack)))
			gridded(sp) <- TRUE
				print("Saving SDev...")
				writeGDAL(sp, fname = paste0(Opth, "/", Att, "/SDev/S_", i, ".tif"), drivername = "GTiff", setStatistics = TRUE, type = "Float32", mvFlag = -1)
		#Predictions	
		Omat <- rep(-1, ncol(T) * nrow(T))
			rrr <- seq(1, ncol(T) * nrow(T), 1)[-Rmvs]
		Omat[rrr] <- unlist(predkNN[, "fitted.values"])
				
		sp <- SpatialPointsDataFrame((expand.grid(x = seq(extent(T)@xmin, extent(T)@xmax - res(T)[1], res(T)[1]) + (0.5 * res(T)[1]), y = rev(seq(extent(T)@ymin, extent(T)@ymax - res(T)[2], res(T)[2]) + (0.5 * res(T)[1])))), data = data.frame(kNN = Omat), proj4string = CRS(projection(Rstack)))
			gridded(sp) <- TRUE
				print("Saving predictions...")
				writeGDAL(sp, fname = paste0(Opth, "/", Att, "/Tiles/T_", i, ".tif"), drivername = "GTiff", setStatistics = TRUE, type = "Float32", mvFlag = -1)
		
	}
	
}#End Att
}#End fff
