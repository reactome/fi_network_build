# This script is used to do survival analysis for calling from a servlet that is used by the Cytoscape FI plug-in.
# This script should be invoked by using Rscript shell command

# Here is the requirements of the file formats to be used in the survival analysis
# score.file.name: it should be a matrix with columns are modules and rows are samples as following:
# Sample\tModule1\tModule2...
# clin.file.name: it should contain at least two columns: OSEVENT and OSDURATION. The unit of duration
# actually is not specified. The user should take care of it. The sheet should be something like the 
# following:
# Sample\tOSEVENT\tOSDURATION
# It is not required the samples in the two files should be the same. This R script will take care of it.

# Load the surv library
require(survival);
# The function is used to do Kaplan-Meier analysis. A module.index, 0-based, must be provided
# for Kaplan-Meier analysis.
kaplan.meyer.analysis <- function(score.file.name, 
		                          clin.file.name,
								  module.index, # Kaplan-Meier analysis can be used for one specific module only
								  plot.file.name = NULL,
								  format = c("PDF", "PNG")) { 
	if (!is.null(plot.file.name)) {
		format <- match.arg(format);
	}
	# load the score matrix
	sample.scores <- read.delim(score.file.name, header=TRUE, sep="\t", check.names = FALSE);
	# load clinical information
	sample.clin <- load.survival.info(sample.scores, clin.file.name);
	surv.tmp <- Surv(sample.clin$survival.time, sample.clin$censoring.status);
	scores <- sample.scores[, as.numeric(module.index) + 2]; # plus 2 to make it consistent with other Java code
	# Check if the scores are binary or not
	values <- unique(scores);
	size.values <- length(values);
	if (size.values > 5) { # This is arbitrary: if there is only 5 values, use these values
		# Use the median to divide the samples into two groups
		median.tmp <- median(scores, na.rm = TRUE);
		groups <- ifelse(scores < median.tmp, 1, 2);
		scores <- groups;
		size.values <- 2;
	}
	survdiff.tmp <- survdiff(surv.tmp ~ scores);
	print(survdiff.tmp);
	# Output to a file
	if (!is.null(plot.file.name)) {
#		print(paste("Format: ", format, sep=""));
		if (format == "PDF") {
			pdf(plot.file.name);
		} else if (format == "PNG") {
			png(plot.file.name);
		}
		plot(survfit(surv.tmp ~ scores), 
				col=2 : (size.values + 2), 
				xlab="Time", 
				ylab="Probability");
		labels <- paste("Group ", 1:size.values, sep="");
		legend("topright", labels, lty=1, col=2 : (size.values + 2), cex=0.5);
		junk <- dev.off(); # Reset to the original device
	}
}

# This function is used to do coxph analysis
coxph.analysis <- function(score.file.name, clin.file.name, module.index = NULL) {
	# load the score matrix. Using "check.names=FALSE" so that original names should be used.
	sample.scores <- read.delim(score.file.name, header=TRUE, sep="\t", check.names = FALSE); 
	# load clinical information
	sample.clin <- load.survival.info(sample.scores, clin.file.name);
	surv.tmp <- Surv(sample.clin$survival.time, sample.clin$censoring.status);
	module.names <- names(sample.scores);
	if (is.null(module.index)) {
#		cat("Univariate CoxPH:\n");
		cat("Module", "Coefficient", "P-value", "\n", sep="\t");
		for (i in 2 : length(sample.scores)) {
			scores <- sample.scores[, i];
			coxph.tmp <- coxph(surv.tmp ~ scores);
			coxph.sum <- summary(coxph.tmp);
			cat(module.names[i], coxph.sum$coef[[1]], coxph.sum$coef[[5]], "\n", sep="\t");
		}
#		cat("\nMulti-variate CoxPH:\n");
#		col.length <- length(sample.scores);
##		print(col.length);
#		scores <- sample.scores[, 2:col.length];
#		coxph.tmp <- coxph(surv.tmp ~ ., data = scores);
#		print(summary(coxph.tmp));
	} else {
		# Do a single module analysis
		cat("Module index: ", module.index, "\n", sep="");
		scores <- sample.scores[, as.numeric(module.index) + 2];
		coxph.tmp <- coxph(surv.tmp ~ scores);
		print(summary(coxph.tmp));
	}
}

# This help function is used to make sure samples in two data frames are the sample
load.survival.info <- function(sample.scores, clin.file.name) {
	# load clinical information
	sample.clin <- read.delim(clin.file.name, header=TRUE, sep="\t", check.names = FALSE);
	sample.names <- sample.clin[, 1]; # The first column should be sample names
	# Using OSDURATION should get columns labeled something like OSDURATION_MONTHS too.
	# So it is based on stringStarts match. This is a good sign so that units can be attached
	# to the column definition.
	survival.time <- sample.clin$OSDURATION;
	survival.status <- sample.clin$OSEVENT;
	which <- match(sample.scores[, 1], sample.names);
	time <- survival.time[which];
	status <- survival.status[which];
	return (list(survival.time = time, censoring.status = status));
}

# Get parameters from the command
args <- commandArgs(TRUE);
# two arguemtns should be provided
if (length(args) < 3) {
	stop("Three arguments (score.file.name, clin.file.name, and method) should be provided in order to run this script!");
}

if (args[3] == 'coxph') {
	if (length(args) > 3) {
		# Module index is 0 based as been used in other Java code.
		coxph.analysis(args[1], args[2], module.index = args[4]);
	} else {
		coxph.analysis(args[1], args[2]);
	}
} else if (args[3] == 'kaplan-meier') { # This is very important: else must be in the same line as }. 
	                                    # Otherwise, a syntical error will be generated.
	if (length(args) < 4) {
		stop("For Kaplan-Meier analysis, a zero-based module index must be provided!");
	}
	# With Kaplan-Meier analysis, a survival plot should be generated.
	if (length(args) == 6) {
		kaplan.meyer.analysis(args[1], args[2], args[4], args[5], args[6]);
	} else {
		kaplan.meyer.analysis(args[1], args[2], args[4]);
	}
}