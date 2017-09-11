
rm(list = objects())

library(foreign)
library(rgdal)
library(raster)
library(maptools)
library(hydroGOF)
library(sp)
library(astro)

SS2 <- function (x, y = NULL, nbin = 128, bandwidth, colramp = colorRampPalette(c("white", blues9)), nrpoints = 0, ret.selection = FALSE, pch = ".", cex = 1, col = "black", transformation = function(x) x^0.25, postPlotHook = box(bty = "l"), xlab = NULL, ylab = NULL, xlim, ylim, xaxs = par("xaxs"), yaxs = par("yaxs"), ...){
    if (!is.numeric(nrpoints) || nrpoints < 0 || length(nrpoints) != 
        1) 
        stop("'nrpoints' should be numeric scalar with value >= 0.")
    nrpoints <- round(nrpoints)
    ret.selection <- ret.selection && nrpoints > 0
    xlabel <- if (!missing(x)) 
        deparse(substitute(x))
    ylabel <- if (!missing(y)) 
        deparse(substitute(y))
    xy <- xy.coords(x, y, xlabel, ylabel)
    xlab <- if (is.null(xlab)) 
        xy$xlab
    else xlab
    ylab <- if (is.null(ylab)) 
        xy$ylab
    else ylab
    x <- cbind(xy$x, xy$y)[I <- is.finite(xy$x) & is.finite(xy$y), 
        , drop = FALSE]
    if (ret.selection) 
        iS <- which(I)
    if (!missing(xlim)) {
        stopifnot(is.numeric(xlim), length(xlim) == 2, is.finite(xlim))
        x <- x[I <- min(xlim) <= x[, 1] & x[, 1] <= max(xlim), 
            , drop = FALSE]
        if (ret.selection) 
            iS <- iS[I]
    }
    else {
        xlim <- range(x[, 1])
    }
    if (!missing(ylim)) {
        stopifnot(is.numeric(ylim), length(ylim) == 2, is.finite(ylim))
        x <- x[I <- min(ylim) <= x[, 2] & x[, 2] <= max(ylim), 
            , drop = FALSE]
        if (ret.selection) 
            iS <- iS[I]
    }
    else {
        ylim <- range(x[, 2])
    }
    map <- grDevices:::.smoothScatterCalcDensity(x, nbin, bandwidth)
    xm <- map$x1
    ym <- map$x2
    dens <- map$fhat
    dens[] <- transformation(dens)
    image(xm, ym, z = dens, col = colramp(256), xlab = xlab, 
        ylab = ylab, xlim = xlim, ylim = ylim, xaxs = xaxs, yaxs = yaxs, bty = "n",
        ...)
    if (!is.null(postPlotHook)) 
        postPlotHook()
    if (nrpoints > 0) {
        nrpoints <- min(nrow(x), ceiling(nrpoints))
        stopifnot((nx <- length(xm)) == nrow(dens), (ny <- length(ym)) == 
            ncol(dens))
        ixm <- 1L + as.integer((nx - 1) * (x[, 1] - xm[1])/(xm[nx] - 
            xm[1]))
        iym <- 1L + as.integer((ny - 1) * (x[, 2] - ym[1])/(ym[ny] - 
            ym[1]))
        sel <- order(dens[cbind(ixm, iym)])[seq_len(nrpoints)]
        x <- x[sel, , drop = FALSE]
        points(x, pch = pch, cex = cex, col = col)
        if (ret.selection) 
            iS[sel]
    }
}

P <- "E:/NWT/NWT_Phase1/ModelFramework/Validation/BT_ALS/ValidationSubsetSlope"
Sfiles <- list.files(P, ".shp$")
Snames <- substr(Sfiles, 1, nchar(Sfiles) - 4)

	M <- NULL
for(i in 1:length(Snames)){
	print(paste("i =", i))
	tmp <- readOGR(P, Snames[i])

	M <- rbind(M, tmp@data)
}
	colnames(M) <- names(tmp)
	Sdiff <- M[, "ALSsht"] - M[, "L4_SH"]
	Cdiff <- M[, "ALScc"] - M[, "L4_CC"]
	M <- cbind(M, Sdiff, Cdiff)
		Rmv <- which(is.na(M[, "Name"]))
	if(length(Rmv) > 0)
		M <- M[-Rmv, ]
	#Round slope
	M[, "Slope"] <- round(M[, "Slope"] / 5) * 5
		
#####STATISTICS#####
#All
	MSdiff <- mean(Sdiff, na.rm = TRUE)
	MCdiff <- mean(Cdiff, na.rm = TRUE)
		print(c(length(Sdiff), MSdiff, MCdiff))
#CDM
	#CON
		iii <- which(M[, "CDM"] == "CON")
	MSdiff <- mean(Sdiff[iii], na.rm = TRUE)
	MCdiff <- mean(Cdiff[iii], na.rm = TRUE)
		print(c(length(iii), MSdiff, MCdiff))
	#DEC
		iii <- which(M[, "CDM"] == "DEC")
	MSdiff <- mean(Sdiff[iii], na.rm = TRUE)
	MCdiff <- mean(Cdiff[iii], na.rm = TRUE)
		print(c(length(iii), MSdiff, MCdiff))
	#MIX
		iii <- which(M[, "CDM"] == "MIX")
	MSdiff <- mean(Sdiff[iii], na.rm = TRUE)
	MCdiff <- mean(Cdiff[iii], na.rm = TRUE)
		print(c(length(iii), MSdiff, MCdiff))
#ECOREGION
	Ueco <- unique(M[, "Name"])
	for(i in 1:length(Ueco)){
			iii <- which(M[, "Name"] == Ueco[i])
		MSdiff <- mean(Sdiff[iii], na.rm = TRUE)
		MCdiff <- mean(Cdiff[iii], na.rm = TRUE)
			#print(data.frame(EcoR = Ueco[i], N = length(iii), MSdiff, MCdiff))
			
			T <- table(M[iii, "CDM"])
			Pdm <- round(100 * ( sum(c(T[2], T[3])) / sum(c(T[1], T[2], T[3])) ))
			Pslp <- ( length(which(M[iii, "Slope"] > 0)) / nrow(M[iii,]) ) * 100
			print(as.character(Ueco[i]))
			print(Pdm)
			print(Pslp)
	}
#SLOPE
	Uslope <- unique(M[, "Slope"])
	for(i in 1:length(Uslope)){
			iii <- which(M[, "Slope"] == Uslope[i])
		MSdiff <- mean(Sdiff[iii], na.rm = TRUE)
		MCdiff <- mean(Cdiff[iii], na.rm = TRUE)
			print(data.frame(EcoR = Ueco[i], N = length(iii), MSdiff, MCdiff))
	}

#####BOXPLOTS CDM#####		
	cols <- c("azure2", "cadetblue3", "aquamarine3", "chocolate1", "coral3", "brown1", "darkorchid2")
#Stand height
jpeg("E:/NWT/NWT_Phase1/ModelFramework/Validation/BT_ALS/FigureXa.jpg", width = 1024, height = 1024, res = 300, pointsize = 12, bg = "white")
	par(mar = c(3.5, 3.5, 1, 1))
	B <- boxplot(Sdiff ~ CDM, data = M, cex = 0.2, pch = 2, varwidth = TRUE, col = cols[c(1, 3, 5)], xlab = "", ylab = "", frame = FALSE, outline = FALSE)
	title(xlab = "EOSD vegetation type", ylab = "Stand height (LiDAR plots - kNN) (m)", line = 2.5)
	box(bty = "l")
	shadowtext(3.5, min(B$stats), "a", "black", "white", font = 2)
dev.off()
#Crown closure
jpeg("E:/NWT/NWT_Phase1/ModelFramework/Validation/BT_ALS/FigureXb.jpg", width = 1024, height = 1024, res = 300, pointsize = 12, bg = "white")
	par(mar = c(3.5, 3.5, 1, 1))
	B <- boxplot(Cdiff ~ CDM, data = M, cex = 0.2, pch = 2, varwidth = TRUE, col = cols[c(1, 3, 5)], frame = FALSE, outline = FALSE)
	title(xlab = "EOSD vegetation type", ylab = "Crown closure (LiDAR plots - kNN) (%)", line = 2.5)
	box(bty = "l")
	shadowtext(3.5, min(B$stats), "b", "black", "white", font = 2)
dev.off()

#####BOXPLOTS ECOREGION#####
	Names <- c("Great Slave Lake Plain", "Hay River Lowland", "Northern AB Lowland", "Sibbeston Lake Plain", "Hylan Highland", "Muskwa Plain", "Horn Plateau")
#Stand height
jpeg("E:/NWT/NWT_Phase1/ModelFramework/Validation/BT_ALS/FigureYa.jpg", width = 1024, height = 1260, res = 300, pointsize = 12, bg = "white")
	par(mar = c(7.5, 3.5, 1, 1))
	B <- boxplot(Sdiff ~ Name, data = M, cex = 0.2, pch = 2, varwidth = TRUE, col = cols, xlab = "", ylab = "", frame = FALSE, outline = FALSE, xaxt = "n")
	title(xlab = "Ecoregion", line = 6.5)
	title(ylab = "Stand height (LiDAR plots - kNN) (m)", line = 2.5)
	axis(1, at = seq(1, 7, 1), labels = Names, las = 3, cex.axis = 0.6)
	box(bty = "l")
	shadowtext(7.5, min(B$stats), "c", "black", "white", font = 2)
dev.off()
#Crown closure
jpeg("E:/NWT/NWT_Phase1/ModelFramework/Validation/BT_ALS/FigureYb.jpg", width = 1024, height = 1260, res = 300, pointsize = 12, bg = "white")
	par(mar = c(7.5, 3.5, 1, 1))
	B <- boxplot(Cdiff ~ Name, data = M, cex = 0.2, pch = 2, varwidth = TRUE, col = cols, xlab = "", ylab = "", frame = FALSE, outline = FALSE, xaxt = "n")
	title(xlab = "Ecoregion", line = 6.5)
	title(ylab = "Crown closure (LiDAR plots - kNN) (%)", line = 2.5)
	axis(1, at = seq(1, 7, 1), labels = Names, las = 3, cex.axis = 0.6)
	box(bty = "l")
	shadowtext(7.5, min(B$stats), "d", "black", "white", font = 2)
dev.off()

#####BOXPLOTS SLOPE#####
#Stand height
jpeg("E:/NWT/NWT_Phase1/ModelFramework/Validation/BT_ALS/FigureZa.jpg", width = 1024, height = 1024, res = 300, pointsize = 12, bg = "white")
	par(mar = c(3.5, 3.5, 1, 1))
	B <- boxplot(Sdiff ~ Slope, data = M, cex = 0.2, pch = 2, varwidth = TRUE, col = cols, xlab = "", ylab = "", frame = FALSE, outline = FALSE, xaxt = "n")
	title(xlab = substitute(paste("Slope (", degree, ")", sep = "")), ylab = "Stand height (LiDAR plots - kNN) (m)", line = 2.5)
	axis(1, at = seq(1, 8, 1), labels = unique(M[, "Slope"]), cex.axis = 0.85)
	box(bty = "l")
	shadowtext(8.5, min(B$stats), "c", "black", "white", font = 2)
dev.off()
#Crown closure
jpeg("E:/NWT/NWT_Phase1/ModelFramework/Validation/BT_ALS/FigureZb.jpg", width = 1024, height = 1024, res = 300, pointsize = 12, bg = "white")
	par(mar = c(3.5, 3.5, 1, 1))
	B <- boxplot(Cdiff ~ Slope, data = M, cex = 0.2, pch = 2, varwidth = TRUE, col = cols, xlab = "", ylab = "", frame = FALSE, outline = FALSE, xaxt = "n")
	title(xlab = substitute(paste("Slope (", degree, ")", sep = "")), ylab = "Crown closure (LiDAR plots - kNN) (%)", line = 2.5)
	axis(1, at = seq(1, 8, 1), labels = unique(M[, "Slope"]), cex.axis = 0.85)
	box(bty = "l")
	shadowtext(8.5, min(B$stats), "d", "black", "white", font = 2)
dev.off()
	
jpeg("E:/NWT/NWT_Phase1/ModelFramework/Validation/BT_ALS/Figure8a.jpg", width = 1024, height = 1024, res = 300, pointsize = 12, bg = "white")
par(mar = c(3.5, 3.5, 1, 1))
#Stand height histogram
	hist(M[, "ALSsht"], seq(0, 30, 2), main = "", xlab = "", ylab = "", ylim = c(0, 3500))
	hist(M[, "L4_SH"], seq(0, 30, 2), add = TRUE, lty = 2, lwd = 2, border = "grey80")
	title(xlab = "Stand height (m)", ylab = "Frequency", line = 2.5)
	box(bty = "l")
	shadowtext(30, 0, "a", "black", "white", font = 2)
	legend(30, 3500, c("LiDAR plots", "kNN"), lty = c(1, 2), lwd = c(1, 2), col = c("black", "grey80"), pch = c(26, 26), xjust = 1)
dev.off()

jpeg("E:/NWT/NWT_Phase1/ModelFramework/Validation/BT_ALS/Figure8b.jpg", width = 1024, height = 1024, res = 300, pointsize = 12, bg = "white")
par(mar = c(3.5, 3.5, 1, 1))
#Crown closure histogram
	hist(M[, "ALScc"], seq(10, 80, 5), main = "", xlab = "", ylab = "", ylim = c(0, 3500))
	hist(M[, "L4_CC"], seq(10, 80, 5), add = TRUE, lty = 2, lwd = 2, border = "grey80")
	title(xlab = "Crown closure (%)", ylab = "Frequency", line = 2.5)
	box(bty = "l")
	shadowtext(80, 0, "b", "black", "white", font = 2)
dev.off()

jpeg("E:/NWT/NWT_Phase1/ModelFramework/Validation/BT_ALS/Figure8c.jpg", width = 1024, height = 1024, res = 300, pointsize = 12, bg = "white")
par(mar = c(3.5, 3.5, 1, 1))
#Stand height density scatter
	SS2(M[, "L4_SH"], M[, "L4_SH_SD"], nrpoints = 0, xlab = "", ylab = "", xlim = c(5, 25), ylim = c(0, 12), postPlotHook = NULL)
	title(xlab = "kNN stand height (m)", ylab = "Uncertainty (m)", line = 2.5)
	box(bty = "l")
	shadowtext(25, 0, "c", "black", "white", font = 2)
dev.off()

jpeg("E:/NWT/NWT_Phase1/ModelFramework/Validation/BT_ALS/Figure8d.jpg", width = 1024, height = 1024, res = 300, pointsize = 12, bg = "white")
par(mar = c(3.5, 3.5, 1, 1))
#Crown closure density scatter
	SS2(M[, "L4_CC"], M[, "L4_CC_SD"], nrpoints = 0, xlab = "", ylab = "", xlim = c(30, 60), ylim = c(2, 16), postPlotHook = NULL)
	title(xlab = "kNN crown closure (%)", ylab = "Uncertainty (%)", line = 2.5)
	box(bty = "l")
	shadowtext(60, 2, "d", "black", "white", font = 2)
dev.off()
	