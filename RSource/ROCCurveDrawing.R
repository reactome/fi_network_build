draw.thresholds <- function(points) {
	for (i in points) {
		fpr <- rocData$False_Positive_Rate[i];
		tpr <- rocData$True_Positive_Rate[i];
		x <- c(0.0, fpr, fpr);
		y <- c(tpr, tpr, 0.0);
		lines(x, y, type="l");
	}
}

calculate.AUC <- function(rocData) {
	# Need to use ROC package from bioconductor
	order <- order(rocData$False_Positive_Rate);
	auc <- trapezint(rocData$False_Positive_Rate[order], rocData$True_Positive_Rate[order], 0, 1);
	print(paste("AUC: ", auc, sep=""));
}

#fileName <- "/Users/wgm/v3/ROC_100_072809.txt";
fileName <- "/Users/gwu/Documents/EclipseWorkspace/caBigR3/results/v4/ROC_100_031512.txt";
rocData <- read.table(fileName, sep="\t", header=TRUE);
## Draw x-y plot
#plot(rocData$False_Positive_Rate, rocData$True_Positive_Rate, type="l", main="NBC ROC Curve", xlab="False Positive Rate", ylab="True Positive Rate", xlim=c(0.00185, 0.05), ylim=c(0.035, 0.965));
plot(rocData$False_Positive_Rate, rocData$True_Positive_Rate, type="l", main="NBC ROC Curve", xlab="False Positive Rate", ylab="True Positive Rate", xlim=c(0.0, 0.05), ylim=c(0.035, 0.5));
## Draw three lines with cutoff values 0.25 (26th), 0.50 (51th), 0.75(76th)
#draw.thresholds(c(26, 51, 76));
