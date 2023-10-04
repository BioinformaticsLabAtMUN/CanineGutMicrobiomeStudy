library(caret)
library(randomForest)


sampleda <- read.csv(".Data/samplevariables.csv", header=TRUE, sep=",", row.names=1, comment.char = "", stringsAsFactors = T)
ilr <-  read.table("Data/ILR.tsv", header=TRUE, sep="\t", row.names=1, comment.char = "")

abundance <- read.table("Data/SeqTab_NoChim_SamplesInColumns_Taxa.tsv", header = T, sep = "\t", stringsAsFactors = T)

head(sampleda)
head(ilr)

setdiff(row.names(sampleda), row.names(ilr))
setdiff(row.names(ilr), row.names(sampleda))

length(intersect(colnames(abundance), row.names(sampleda)))

#Remove sOTUS without genus
abundance <- abundance[!is.na(abundance[,"Genus"]),]
#sum abundance same genus
abundanceG <- data.frame(do.call("rbind", by(abundance[, grepl("Sample", colnames(abundance))], abundance$Genus, colSums)))

#make samples rows
abundanceG <- t(abundanceG)

dim(abundanceG) #48 samples, 83 genus
head(abundanceG)

#---------------------------------  


performGridSearch <- function(D, s = 123, reps = 10, filter = T){
   Variables <- c("Aggression", "Anxiety")

   fits <- list()
   importance <- list()
   trainingDataF <- list()


	set.seed(s)                         
	for (var in Variables){
		trainingData <- data.frame(D, class = sampleda[row.names(D),var])
		if (filter){
			RF <- randomForest(class ~ ., data = trainingData, importance = TRUE, keep.forest = FALSE)
			imp <- importance(RF, type = 1)
			importance[[var]] <- imp[order(imp, decreasing = T),]
			Features2keep <- names(importance[[var]])[importance[[var]] > 0]
			print(length(Features2keep))
			trainingData <- trainingData[,c(Features2keep, "class")]
		}
		trainingDataF[[var]] <- trainingData
	
		fitControl <- trainControl(## 10-fold CV
                           method = "repeatedcv",
                           number = 10,
                           repeats = reps,
                           ## Estimate class probabilities
                           classProbs = TRUE,
                           ##Keep predictions to be able to draw ROC and PRC
                           savePredictions = T,
                           ## Evaluate performance using 
                           ## the following function
                           summaryFunction = prSummary)

		##### Random Forest                          
		rfGrid <-  expand.grid(mtry = c(round(log2(ncol(trainingData) -1),0), round(sqrt(ncol(trainingData)-1),0), 
						       floor(ncol(trainingData)/5), floor(ncol(trainingData)/4),
                                floor(ncol(trainingData)/3),floor(ncol(trainingData)/2)))
		rfFit <- train(class ~ ., data = trainingData, method = "rf", ntree = 1000, trControl = fitControl, tuneGrid = rfGrid, metric = "F")
    		fits[[var]] <- rfFit
	}
	
	res <- list(fits = fits, imp = importance, data = trainingDataF)
	return(res)
}

resILR <- performGridSearch(ilr, reps = 10)
resILR_NF <- performGridSearch(ilr, reps = 10, filter = FALSE)


resAbun <- performGridSearch(abundanceG, reps = 10)
resAbun_NF <- performGridSearch(abundanceG, reps = 10, filter = FALSE)



#### Variable Importance
### Take hyperparameters of best model found by grid search

RF <- randomForest(class ~ ., data = resILR_NF$data[["Aggression"]], ntree = 1000, mtry = resILR_NF$fits$Aggression$bestTune$mtry, importance = TRUE, keep.forest = FALSE)
imp <- importance(RF, type = 1)
resILR_NF$imp[["Aggression"]] <- imp[order(imp, decreasing = T),]

RF <- randomForest(class ~ ., data = resILR_NF$data[["Anxiety"]], ntree = 1000, mtry = resILR_NF$fits$Anxiety$bestTune$mtry, importance = TRUE, keep.forest = FALSE)
imp <- importance(RF, type = 1)
resILR_NF$imp[["Anxiety"]] <- imp[order(imp, decreasing = T),]



RF <- randomForest(class ~ ., data = resAbun_NF$data[["Aggression"]], ntree = 1000, mtry = resAbun_NF$fits$Aggression$bestTune$mtry, importance = TRUE,  keep.forest = FALSE)
imp <- importance(RF, type = 1)
resAbun_NF$imp[["Aggression"]] <- imp[order(imp, decreasing = T),]

RF <- randomForest(class ~ ., data = resAbun_NF$data[["Anxiety"]], ntree = 1000, mtry = resAbun_NF$fits$Anxiety$bestTune$mtry, importance = TRUE,  keep.forest = FALSE)
imp <- importance(RF, type = 1)
resAbun_NF$imp[["Anxiety"]] <- imp[order(imp, decreasing = T),]


#Variable importance final filtered model

RF <- randomForest(class ~ ., data = resILR$data[["Aggression"]], ntree = 1000, mtry = resILR$fits$Aggression$bestTune$mtry, importance = TRUE, keep.forest = FALSE)
imp <- importance(RF, type = 1)
resILR$impF <- list()
resILR$impF[["Aggression"]] <- imp[order(imp, decreasing = T),]

RF <- randomForest(class ~ ., data = resILR$data[["Anxiety"]], ntree = 1000, mtry = resILR$fits$Anxiety$bestTune$mtry, importance = TRUE, keep.forest = FALSE)
imp <- importance(RF, type = 1)
resILR$impF[["Anxiety"]] <- imp[order(imp, decreasing = T),]



RF <- randomForest(class ~ ., data = resAbun$data[["Aggression"]], ntree = 1000, mtry = resAbun$fits$Aggression$bestTune$mtry, importance = TRUE, keep.forest = FALSE)
imp <- importance(RF, type = 1)
resAbun$impF <- list()
resAbun$impF[["Aggression"]] <- imp[order(imp, decreasing = T),]


RF <- randomForest(class ~ ., data = resAbun$data[["Anxiety"]], ntree = 1000, mtry = resAbun$fits$Anxiety$bestTune$mtry, importance = TRUE, keep.forest = FALSE)
imp <- importance(RF, type = 1)
resAbun$impF[["Anxiety"]] <- imp[order(imp, decreasing = T),]



#-------------------------------------
# Compare performance models: Using ILR with/without prefiltering  and using abundance with/without prefiltering
# Without prefiltering are indicated with NF

### For Aggression
resamps <- resamples(list(
                          ILR = resILR$fits[["Aggression"]],
                          ILR_NF = resILR_NF$fits[["Aggression"]],
                          Abun = resAbun$fits[["Aggression"]],
                          Abun_NF = resAbun_NF$fits[["Aggression"]]
                          ))
resamps
summary(resamps)


trellis.par.set(caretTheme())
dotplot(resamps, metric = "AUC", main = "Classification performance (AUPRC) - Aggression")

### For Anxiety

resamps <- resamples(list(
                          ILR = resILR$fits[["Anxiety"]],
                          ILR_NF = resILR_NF$fits[["Anxiety"]],
                          Abun = resAbun$fits[["Anxiety"]],
                          Abun_NF = resAbun_NF$fits[["Anxiety"]]
                          ))
resamps
summary(resamps)


trellis.par.set(caretTheme())
dotplot(resamps, metric = "AUC",  main = "Classification performance (AUPRC) - Anxiety")


###### Confusion matrices of models constructed filtering attributes based on decrease in accuracy

confusionMatrix(resILR$fits[["Aggression"]], "average")
 
confusionMatrix(resILR$fits[["Anxiety"]], "average")


confusionMatrix(resAbun$fits[["Aggression"]], "average")



confusionMatrix(resAbun$fits[["Anxiety"]], "average")



### See individual PR curves - Use this to find the AUPRC and AUC-ROC
library(MLeval)
### CHANGE for a different model and different behavioural group
x <- evalm(resAbun$fits[["Anxiety"]], positive = "High", plots="pr")
plot(x$proc)

#####################
### Generate figures for manuscript

#1) Get the predicted probabilities of high aggression/anxiety from CV
#2) Calculate PR values
#3) Get range
#4) Create ribbon plot

evaluatePredictions <- function(predD, testLabels){
	res <- list()
	res$predD <- predD
	require(ROCR)
	require(PRROC)
	#Evaluate predictions
	res$pred <- prediction(res$predD, as.logical(testLabels))
	res$PR <- performance(res$pred, measure = "prec", x.measure = "rec")
    res$pr <- pr.curve(scores.class0 = res$predD[testLabels == 1], scores.class1 = res$predD[testLabels == 0], curve = T, max.compute=T,min.compute=T,rand.compute=T)
    return(res)
}



getPR_CV <- function(cvResults, m) {
	print(m)
	df_PR <- list()
	best<- cvResults$bestTune$mtry
	print(paste("mtry ", best, sep ="")
	reps <- c(paste("Rep0", 1:9, sep=""), "Rep10")
	cvPred <- cvResults$pred[cvResults$pred$mtry == best,c("rowIndex", "obs", "High", "mtry", "Resample")]
	commonRecall <- c()
    for (i in reps){
		perfFold <- evaluatePredictions(cvPred[grep(i, cvPred$Resample),"High"], ifelse(cvPred[grep(i, cvPred$Resample),"obs"] == "High", 1, 0))
		print(perfFold$pr)
		df_PR[[i]] <- unlist(attr(perfFold$PR, "y.values")) #precision
		names(df_PR[[i]]) <- unlist(attr(perfFold$PR, "x.values")) #recall
		df_PR[[i]][1] <- 1
		if (i == "Rep01") {
			commonRecall <- names(df_PR[[i]])
		} else {
			commonRecall <- intersect(commonRecall, names(df_PR[[i]]))
		}
	}
	df_tmp <- sapply(df_PR, function(l){l[commonRecall]})
	df_range <- data.frame(recall = as.numeric(commonRecall), min = apply(df_tmp,1, min), max=apply(df_tmp,1, max), mean = rowMeans(df_tmp), method = rep(m, length(commonRecall)))
	row.names(df_range) <- NULL
	return(df_range)
}

PR_CV_aggression <- getPR_CV(resILR$fits$Aggression, "ILR + FS")
PR_CV_aggression <- rbind(PR_CV_aggression, getPR_CV(resILR_NF$fits$Aggression, "ILR"))
PR_CV_aggression <- rbind(PR_CV_aggression, getPR_CV(resAbun$fits$Aggression, "Abundance + FS"))
PR_CV_aggression <- rbind(PR_CV_aggression, getPR_CV(resAbun_NF$fits$Aggression, "Abundance"))

library(ggplot2)
library(RColorBrewer)

colores <- brewer.pal(4, "Pastel2")
#display.brewer.pal(4, "Pastel2")

table(sampleda$Aggression) #to find location default performance

ggplot(PR_CV_aggression, aes(x = recall)) + geom_ribbon(stat = "identity", aes(ymin = min, ymax = max, fill = method, group = method), alpha = 0.5) + ylim(0,1) + labs(title ="Precision-recall curves for models predicting higher aggression", y = "Precision", x = "Recall") + geom_line(aes(x=recall, y=mean, group = method, color = method), linewidth = 1, show.legend = F) + guides(fill = guide_legend(title = "")) + theme(legend.position = "top") + scale_fill_manual(values = colores) + scale_color_manual(values = c("darkseagreen", "orange", "slategray", "hotpink")) + geom_hline(yintercept = round(25/48,2), color = "grey50", linewidth = 0.25)



PR_CV_anxiety <- getPR_CV(resILR$fits$Anxiety, "ILR + FS")
PR_CV_anxiety <- rbind(PR_CV_anxiety , getPR_CV(resILR_NF$fits$Anxiety, "ILR"))
PR_CV_anxiety <- rbind(PR_CV_anxiety , getPR_CV(resAbun$fits$Anxiety, "Abundance + FS"))
PR_CV_anxiety <- rbind(PR_CV_anxiety , getPR_CV(resAbun_NF$fits$Anxiety, "Abundance"))

table(sampleda$Anxiety)


ggplot(PR_CV_anxiety, aes(x = recall)) + geom_ribbon(stat = "identity", aes(ymin = min, ymax = max, fill = method, group = method), alpha = 0.5) + ylim(0,1) + labs(title ="Precision-Recall curves for models predicting higher anxiety", y = "Precision", x = "Recall") + geom_line(aes(x=recall, y=mean, group = method, color = method), linewidth = 1, show.legend = F) + guides(fill = guide_legend(title = "")) + theme(legend.position = "top") + scale_fill_manual(values = colores) + scale_color_manual(values = c("darkseagreen", "orange", "slategray", "hotpink")) + geom_hline(yintercept = round(23/48,2), color = "grey50", linewidth = 0.25)


#Plots variable importance
plot(resILR$impF[["Aggression"]], ylab = "Decrease in Accuracy", xlab = "Features", pch = 20, col = ifelse(resILR$impF[["Aggression"]] > 2.23, "firebrick", "black"), main = "Importance of 176 features in the model ILR+FS for aggression")
text(x= c(1,3,5,8,10,14), y=resILR$impF[["Aggression"]][c(1,3,5,8,10,14)], labels = c(names(resILR$impF[["Aggression"]][1]),  paste(names(resILR$impF[["Aggression"]])[2:4], collapse = ", "), paste(names(resILR$impF[["Aggression"]])[5:7], collapse = ", "),  paste(names(resILR$impF[["Aggression"]])[8:9], collapse = ", "), paste(names(resILR$impF[["Aggression"]])[10:11], collapse = ", "), paste(names(resILR$impF[["Aggression"]])[12:14], collapse = ", ") ), cex = 0.7, pos = 4, offset = 0.5)
abline(h=0, col = "grey")

plot(resAbun$impF[["Aggression"]], ylab = "Decrease in Accuracy", xlab = "Features", pch = 20, col = ifelse(resAbun$impF[["Aggression"]] > 3, "firebrick", "black"), main = "Importance of 23 features in the model Abundance+FS for aggression")
text(x= c(1,2,3), y=resAbun$impF[["Aggression"]][c(1,2,3)], labels = c(names(resAbun$impF[["Aggression"]][1:2]),paste(names(resAbun$impF[["Aggression"]])[3:5], collapse = ", ") ), cex = 0.7, pos = 4, offset = 0.5)
abline(h=0, col = "grey")

plot(resILR$impF[["Anxiety"]], ylab = "Decrease in Accuracy", xlab = "Features", pch = 20, col = ifelse(resILR$impF[["Anxiety"]] > 2.7, "firebrick", "black"), main = "Importance of 171 features in the model ILR+FS for anxiety")
text(x= c(1,2,4,5,7, 10), y= resILR$impF[["Anxiety"]][c(1,2,4,5,7, 10)], labels = c(names(resILR$impF[["Anxiety"]][1:2]), paste(names(resILR$impF[["Anxiety"]])[3:4], collapse = ", "), names(resILR$impF[["Anxiety"]][5]), paste(names(resILR$impF[["Anxiety"]])[6:9], collapse = ", "), paste(names(resILR$impF[["Anxiety"]])[10:12], collapse = ", ")), cex = 0.7, pos = 4, offset = 0.5)
abline(h=0, col = "grey")


plot(resAbun$impF[["Anxiety"]], ylab = "Decrease in Accuracy", xlab = "Features", pch = 20, col = ifelse(resAbun$impF[["Anxiety"]] > 4, "firebrick", "black"), main = "Importance of 23 features in the model Abundance+FS for anxiety")
text(x= c(1,2,4,5), y=resAbun$impF[["Anxiety"]][c(1,2,4,5)], labels = c(names(resAbun$impF[["Anxiety"]][1:2]),paste(names(resAbun$impF[["Anxiety"]])[3:4], collapse = ", "), names(resAbun$impF[["Anxiety"]][5]) ), cex = 0.7, pos = 4, offset = 0.5)
abline(h=0, col = "grey")

length(resAbun$impF[["Anxiety"]])

