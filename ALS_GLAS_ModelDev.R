
data <- read.table("E:/NWT/NWT_Phase1/Data/ALS_GLAS_CrossValidation/CC.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE)
data <- data[-27, ]
data[, "ScaledGLz"] <- data[, "ScaledGLz"] * log(10)

	DepV <- "FScc"
		Plotname <- "a) Stand height"
		Metrics <- c("ScaledGLz")
		Mnames <- c("Lz")


	Stats <- data.frame(matrix(NA, ncol = 9, nrow = length(Metrics), dimnames = list(NULL, c("Metric", "N", "a", "b", "Rsq", "RMSE", "Mpe", "Min", "Max"))))
for(m in 1:length(Metrics)){
	#Define training dataset
		Ytrain <- data[, DepV]
		Xtrain <- data[, Metrics[m]]
	#Regress training data to find model coefficients
	if(DepV == "ALSht")
		TrainReg <- lm(Ytrain ~ Xtrain)
	if(DepV == "FScc")
		TrainReg <- nls(Ytrain ~ a * Xtrain ^ b, start = list(a = 64.63, b = 0.25))
	#Define testing dataset
		Ytest <- data[, DepV]
	if(DepV == "ALSht")
		Xtest <- coef(TrainReg)[2] * data[, Metrics[m]] + coef(TrainReg)[1]	
	if(DepV == "FScc")
		Xtest <- coef(TrainReg)[1] * data[, Metrics[m]] ^ coef(TrainReg)[2]	
	#Calculate statistics for testing data
	if(DepV == "ALSht"){
		a <- round(coef(TrainReg)[2], 2)
		b <- round(coef(TrainReg)[1], 2)
	}else if(DepV == "FScc"){
		a <- round(coef(TrainReg)[1], 2)
		b <- round(coef(TrainReg)[2], 2)
	}
	Tdata <- data.frame(X = Xtrain, Y = Ytrain)
	Rsq <- round(summary(lm(Ytest ~ Xtest))$adj.r.squared, 2)
	Rmse <- round(rmse(Xtest, Ytest), 2)
	Mpe <- round(mean(cvTool(call("lm", formula = Y ~ X), data = Tdata, x = Tdata[, "X"], y = Tdata[, "Y"], folds = cvFolds(nrow(Tdata), K = 5, R = 65))), 1)
	Min <- round(min(Xtest), 1)
	Max <- round(max(Xtest), 1)
		Stats[m, ] <- c(Metrics[m], nrow(data), a, b, Rsq, Rmse, Mpe, Min, Max)
}
	print(Stats)
