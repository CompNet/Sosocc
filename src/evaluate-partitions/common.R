

#############################################################################################
# It transforms a matrix into a vector. If the matrix is symmetric, 
#   the parameter 'process.half.mtrx' should be set to TRUE.
# 
# dist.mtrx: a matrix of distance scores
# process.half.mtrx: boolean. If TRUE, only lower triangle of the matrix is used (i.e. symmetric matrix).
#
# returns a vector of distance scores
#############################################################################################
retreive.dist.scores.from.matrix = function(dist.mtrx, process.half.mtrx = FALSE){
	dist.values = c()
	counter = 1
	for(i in seq(1,ncol(dist.mtrx)))
		for(j in seq(1,nrow(dist.mtrx))){
			if(process.half.mtrx && j<i){ #  lower triangle
				dist.values[counter] = dist.mtrx[i,j]
				counter = counter + 1  
			} 
			
			if(!process.half.mtrx) { # whole matrix
				dist.values[counter] = dist.mtrx[i,j]
				counter = counter + 1 
			}
		}
	return(dist.values)
}



#############################################################################################
# It is a wrapper function to transform a matrix of distance scores into vector. 
#   The matrix is loaded from a file.
#
# eval.folder: the folder from which the matrix will be loaded
# result.filename: file name containing the matrix of distance scores
# process.half.mtrx: boolean. If TRUE, only lower triangle of the matrix is used (i.e. symmetric matrix).
#
# returns a vector of distance scores
#############################################################################################
retreive.dist.scores.from.matrix2 = function(eval.folder, result.filename, process.half.mtrx = FALSE){
	mtrx.file = file.path(eval.folder, result.filename)
	dist.mtrx = read.csv(mtrx.file, row.names = 1, header= TRUE, check.names=FALSE) # square matrix

	dist.values = retreive.dist.scores.from.matrix(dist.mtrx, process.half.mtrx=process.half.mtrx)
	return(dist.values)
}

