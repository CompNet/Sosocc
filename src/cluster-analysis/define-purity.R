############################################################################
# Functions used in the processing of more than one measure.
#
# Vincent Labatut 2012-16
# The paper: https://arxiv.org/pdf/1303.5441.pdf
# Github code source: https://github.com/CompNet/TopoMeasures/blob/master/src/PurityMeasure.R
############################################################################
# main lib, used to represent and process graphs
library("igraph")




############################################################################
# Processes the confusion matrix for the specified partitions, i.e. the matrix 
# whose cell (i,j) contains the number of elements belonging to part #i in 
# partition #1 and to part #j in partition #2.
# 
# Note: this function amounts to using process.weighted.confusion.matrix with
# a weight vector of 1s. Except it's faster.
# 
# partition1: the first partition to consider, represented as an integer vector.
#			  Each value in the vector represents the id of a part. The parts 
#			  must be counted starting from one (not zero).
# partition2: the second partition to consider. Same representation than for
#			  the first one. Both are inter-exchangeable (symmetric measure).
# returns: the confusion matrix, a.k.a. contingency table.
############################################################################
process.confusion.matrix <- function(partition1, partition2)
{	p1 <- factor(partition1, levels=1:max(partition1))
	p2 <- factor(partition2, levels=1:max(partition2))
	result <- as.matrix(table(p1,p2))
	return(result)
}




# purity.level has 2 options: purity.level="partition" OR purity.level="cluster"
# when purity.level="cluster", the method outputs the purity value of each cluster in the 'partition' variable
# otherwise, it outputs the general putity associated with the 'partition' variable
compute.purity <- function(partition, ground.truth, purity.level="partition")
{	# process the confusion matrix
	conf.matrix <- process.confusion.matrix(partition,ground.truth)
	
	
	# init
	total <- 0
	nb.total.in.clu <- rep(0,nrow(conf.matrix))
	purity.per.cluster = c()
	assoc.clu.indx = c()
	
	# for each part in partition, identify the corresponding parts in ground.truth
	# i.e. those with the largest intersection (there can be several)
	for(r in 1:nrow(conf.matrix)){
		nb.total.in.clu[r] <- sum(conf.matrix[r,])
		max.val = max(conf.matrix[r,])
		mx.indx = which(conf.matrix[r,] == max.val)
		assoc.clu.indx[r] = paste(mx.indx, collapse=",")
		purity.per.cluster[r] = max.val/nb.total.in.clu[r]
	}
	
	if(purity.level == "partition"){
		# compute the purity of the partition
		n = sum(conf.matrix)
		for(r in 1:nrow(conf.matrix))
			total <- total + (nb.total.in.clu[r]/n)*purity.per.cluster[r]
		return(total)
	}
	
	# prepare output for writing into csv
	m = cbind(1:length(purity.per.cluster), nb.total.in.clu, purity.per.cluster, assoc.clu.indx)
	colnames(m) = c("k", "nb.item.in.clu", "purity", "associated cluster id")
	
	return(m)
}

