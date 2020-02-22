# This R script is used to calculate matrix exponential
# expm should take care of Matrix. 
#library("Matrix")
# Required for write.matrx
library("MASS")
library("expm")

calculate.matrix.exponential <- function(in.file.name,
									     out.file.name) {
	# load matrix into R
	in.matrix <- read.delim(in.file.name, header=FALSE, sep="\t");
	print("Loaded the data");
	in.matrix <- as.matrix(in.matrix);
	print("Converted as a matrix");
	out.matrix <- expm(in.matrix);
	print("Done expm");
	# Output matrix
	write.matrix(out.matrix, file = out.file.name, sep="\t");
	print("Exported the matrix to a file");
}

# dir.name <- "/Users/gwu/Documents/temp"
# in.file.name <- paste(dir.name, "SimpleMatrix.txt", sep="/")
# out.file.name <- paste(dir.name, "SimpleMatrixResult.txt", sep="/")

# dir.name <- "/Users/gwu/Documents/EclipseWorkspace/FINetworkBuild/results/2016"
# in.file.name <- paste(dir.name, "HotNet_L_matrix_2016.txt", sep="/")
# out.file.name <- paste(dir.name, "HeatKernel_HotNet_time_01_2016_040417.txt", sep="/")

# dir.name <- "/Users/wug/git/FINetworkBuild/results/2018"
# in.file.name <- paste(dir.name, "HotNet_L_matrix_2018.txt", sep="/")
# out.file.name <- paste(dir.name, "HeatKernel_HotNet_time_01_2018_122718.txt", sep="/")

dir.name <- "/Users/wug/git/FINetworkBuild/results/2019"
in.file.name <- paste(dir.name, "HotNet_L_matrix_2019.txt", sep="/")
out.file.name <- paste(dir.name, "HeatKernel_HotNet_time_01_2019_021920.txt", sep="/")

calculate.matrix.exponential(in.file.name, out.file.name)

# # Get parameters from the command
# args <- commandArgs(TRUE);
# 
# if (length(args) < 2) {
# 	stop("Two arguments (matrix.in.file.name, matrix.out.file.name) should be provided!");
# }
# 
# calculate.matrix.exponential(args[1], args[2])
