# This R script is used to calculate matrix exponential
require(Matrix)
require(MASS)

calculate.matrix.exponential <- function(in.file.name,
									     out.file.name) {
	# load matrix into R
	in.matrix <- read.delim(in.file.name, header=FALSE, sep="\t");
	in.matrix <- as.matrix(in.matrix);
	out.matrix <- expm(in.matrix);
	# Output matrix
	write.matrix(out.matrix, file = out.file.name, sep="\t");
}

# Get parameters from the command
args <- commandArgs(TRUE);

if (length(args) < 2) {
	stop("Two arguments (matrix.in.file.name, matrix.out.file.name) should be provided!");
}

calculate.matrix.exponential(args[1], args[2]);
