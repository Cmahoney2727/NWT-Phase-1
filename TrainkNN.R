# Train kNN models
#
# Author: craig.mahoney
# Created: 25/04/2014
# Last Revised: 31/07/2014
#
#	source("E:/NWT/ModelFramework/Scripts/R/TrainkNN.R")
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

#####PREAMBLE
	ttt <- Sys.time()

Pth <- "E:/NWT/NWT_Phase1/ModelFramework/Observation_Layers/"
	Filters <- list.files(Pth, "Level_")[-5]
	Var <- c("SHt", "CC")

#####START
for(ddd in 1:length(Filters)){	#1:length(Filters)
	Gpth <- paste(Pth, Filters[ddd], "/", sep = "")
	Gdata <- get(load(paste(Gpth, list.files(Gpth, ".out"), sep = "")))
for(vvv in 1:length(Var)){
		#Out Directory
		OutDir <- paste("E:/NWT/ModelFramework/kNN/TrainedModels/", Filters[ddd], "/", Var[vvv], "/", sep = "")
		dir.create(paste("E:/NWT/ModelFramework/kNN/TrainedModels/", Filters[ddd], "/", sep = ""))
		dir.create(OutDir)
#Prepare for CC predictions		
	if(Var[vvv] == "CC"){
			ULz <- unique(Gdata[, "CDM"])
			Lz <- matrix(NA, ncol = 1, nrow = nrow(Gdata), dimnames = list(NULL, c("CC")))
		for(cc in 1:4){
				iii <- which(Gdata[, "CDM"] == ULz[cc])
			if((ULz[cc] == "CON") || (ULz[cc] == "WET")){
				Lz[iii, ] <- Gdata[iii, "LzCON"]
			}else if(ULz[cc] == "MIX"){
				Lz[iii, ] <- Gdata[iii, "LzMIX"]
			}else if(ULz[cc] == "DEC"){
				Lz[iii, ] <- Gdata[iii, "LzDEC"]
			}
		}
		if(any(Lz < 0))
			Lz[which(Lz < 0)] <- 0
	}
	
#####kNN
	DataLayers <- list.files("E:/NWT/ModelFramework/Data_Layers", ".tif$")
	DataLayers <- DataLayers[c(7, 1, 2, 3, 4, 5, 6, 8)]
	count <- 0
	
	if(Var[vvv] == "SHt"){
		Training <- data.frame(ShtGLAS = Gdata[, "ShtGLAS"], round(Gdata[, substr(DataLayers, 1, nchar(DataLayers) - nchar(".tif"))], 3))
			save(Training, file = paste(OutDir, "TrainingData_", Var[vvv], ".dat", sep = ""), compress = TRUE)

		TunekNN <- train.kknn(as.formula("ShtGLAS ~ EOSD + CDEM + CMI + SMI"), data = Training, kmax = 10, kernel = c("epanechnikov", "triangular", "biweight", "triweight", "cos", "inv", "gaussian", "optimal"))
	
		kNN <- list(Name = "ShtGLAS", k = TunekNN$best.parameters$k, kernel = TunekNN$best.parameters$kernel, predictors = c("EOSD", "CDEM", "CMI", "SMI"))
			save(kNN, file = paste0(OutDir, "Optimized_kNN.out"), compress = TRUE)
	}else if(Var[vvv] == "CC"){
		Training <- data.frame(CclGLAS = Gdata[, "CclGLAS"], Gdata[, substr(DataLayers, 1, nchar(DataLayers) - nchar(".tif"))])
			save(Training, file = paste(OutDir, "TrainingData_", Var[vvv], ".dat", sep = ""), compress = TRUE)
			
		TunekNN <- train.kknn(as.formula("CclGLAS ~ EOSD + Band3 + Band4 + Band5 + CDEM + CMI + SMI"), data = Training, kmax = 25, kernel = c("epanechnikov", "triangular", "biweight", "triweight", "cos", "inv", "gaussian", "optimal"))
	
		kNN <- list(Name = "CclGLAS", k = TunekNN$best.parameters$k, kernel = TunekNN$best.parameters$kernel, predictors = c("EOSD", "CDEM", "CMI", "SMI"))
			save(kNN, file = paste0(OutDir, "Optimized_kNN.out"), compress = TRUE)
	}
	
}#End vvv
}#End ddd

require(kknn)
require(astro)
require(TeachingDemos)

Spth <- "E:/NWT/NWT_Phase1/ModelFramework/kNN/TrainedModels/Level_4/CC/TrainingData_CC.dat"
	Training <- get(load(Spth))

	Cmb <- list()
for(j in 2:length(names(Training)[-1]))	
	Cmb[[j]] <- combn(names(Training)[-1], j)
	Cmb <- Cmb[-1]
	
	ttt <- Sys.time()
	Bmse <- list()
for(i in 1:length(Cmb)){

	print(paste(i, "of", length(Cmb)))
	
		V <- NULL
	for(j in 1:ncol(Cmb[[i]])){
		tmp <- Cmb[[i]][, j]
			Pred <- paste(tmp, collapse = " + ")
			Formula <- paste("CclGLAS ~", Pred)
	
	set.seed(12345)
	system.time(TunekNN <- train.kknn(as.formula(Formula), data = Training, kmax = 25, kernel = c("epanechnikov", "triangular", "biweight", "triweight", "cos", "inv", "gaussian", "optimal")))
	
	V <- c(V, min(TunekNN$MEAN.SQU))
		#plot(TunekNN)
	}
	Bmse[[i]] <- V
}
	print(abs(ttt - Sys.time()))
	
	Mcomb <- NULL
for(k in 1:length(Bmse))
	Mcomb <- c(Mcomb, min(Bmse[[k]]))
	
	iii <- which(Mcomb == min(Mcomb))
	jjj <- which(Bmse[[iii]] == min(Bmse[[iii]]))
	
	
		tmp <- Cmb[[iii]][, jjj]
			Pred <- paste(tmp, collapse = " + ")
			Formula <- paste("CclGLAS ~", Pred)
	
	set.seed(12345)
	system.time(TunekNN <- train.kknn(as.formula(Formula), data = Training, kmax = 25, kernel = c("epanechnikov", "triangular", "biweight", "triweight", "cos", "inv", "gaussian", "optimal")))
	qqq <- which(TunekNN$MEAN.SQU == min(TunekNN$MEAN.SQU)) / 25
	Dec <- (qqq %% 1)
	Int <- qqq - Dec
	
	Col <- Int + 1
	Row <- 25 * Dec
	
#####Plot data
		Names <- c("Epanechinkov", "Triangular", "Biweight", "Triweight", "Cosine", "Inverse", "Gaussian", "Optimal")
		cols <- c("azure4", "cadetblue3", "aquamarine1", "chocolate1", "coral3", "brown1", "darkorchid2", "darkblue")
		ltys <- c(1, 2, 3, 1, 2, 3, 1, 2)

jpeg("E:/NWT/NWT_Phase1/ModelFramework/kNN/TrainedModels/Optimized_kNN_Cclab.jpg", width = 1024, height = 1024, res = 300, pointsize = 12, bg = "white")	
par(mar = c(3.5, 3.5, 1, 1), mfrow = c(2, 1))
	plot(Mcomb, type = "b", xlab = "", ylab = "", xaxt = "n", yaxt = "n", ylim = c(42, 56), bty = "l", lwd = 2, pch = 6, cex = 0.75)
	axis(1, at = seq(1, 7, 1), labels = seq(2, 8, 1))
	title(xlab = "# predictors in combination", line = 2.25)
	axis(2, at = seq(42, 56, 7), labels = seq(42, 56, 7), cex.axis = 1.0, tick = TRUE)
	shadowtext(length(Mcomb), 42, "a", "black", "white", font = 2, adj = c(1, 0))
	
	plot(Bmse[[iii]], type = "l", xlab = "", ylab = "", lwd = 2, col = "grey50", bty = "l", yaxt = "n")
	axis(2, at = seq(44, 50, 3), labels = seq(44, 50, 3))
	title(xlab = "Index of 7 predictor combinations", line = 2.25)
	shadowtext(length(Bmse[[iii]]), min(Bmse[[iii]]), "b", "black", "white", font = 2, adj = c(1, 0))
	
	mtext(side = 2, "Mean squared error (m)", outer = TRUE, line = -1)
dev.off()
	
	
jpeg("E:/NWT/NWT_Phase1/ModelFramework/kNN/TrainedModels/Optimized_kNN_Cclc.jpg", width = 1024, height = 1024, res = 300, pointsize = 12, bg = "white")
	par(mar = c(3.5, 3.5, 1, 1))
	plot(seq(1, 25, 1), TunekNN$MEAN.SQU[, 1], type = "l", xlim = c(1, 25), ylim = range(TunekNN$MEAN.SQU), col = cols[1], lty = ltys[1], lwd = 2, bty = "l", xlab = "", ylab = "")
	for(q in 2:8)
		points(seq(1, 25, 1), type = "l", TunekNN$MEAN.SQU[, q], col = cols[q], lty = ltys[q], lwd = 2)
	points(Row + 25, TunekNN$MEAN.SQU[Row + 25, Col], pch = 6, cex = 0.5, lwd = 3)
	
	title(xlab = "k", ylab = "Mean squared error (m)", line = 2.25)
	#legend(25, max(TunekNN$MEAN.SQU), Names, col = cols, lty = ltys, lwd = rep(2, 8), xjust = 1, cex = 0.8)
	shadowtext(25, max(TunekNN$MEAN.SQU), "c", col = "black", "white", font = 2)
dev.off()
