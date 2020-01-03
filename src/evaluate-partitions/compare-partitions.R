

#################################################################
# It processes all the specified measures in order to compare the partitions 'partition1' and 'partition2'.
#   Since user may use a similarity or distance measure, we unified this part
#       - by converting a similarity measure into distance one (1-value).
#       - normalizing distance maesures (user still can access to unnormalized version of the distance measures)
#   Therefore, all measures are in [0,1] range,  where 0 (resp. 1) represents identical (resp. completely different) partitions.
#   The only exception among measures is ARI, where its original range is [-1,1]. 
#   To do so, we had to set to zero all negative ARI values in order to get a range of [0,1].
#
# partition1: first partition.
# partition2: second partition.
# measures: vector of measure names: c("vi", "nmi", ..).
#
# returns: vector of values corresponding to each specified measures.
#################################################################
compare.partition.pair.as.distance <- function(partition1, partition2, measures="nmi")
{	# init result vector
	result <- rep(NA,length(measures))
	names(result) <- measures
	
	library(igraph)
	library(clues)
	library(NMF)
	library(e1071) # hamming distance
	library(entropy)
	
	
	if(MI %in% measures){ # a distance value
	   result[MI] <- entropy::mi.plugin(table(partition1,partition2))
	}
	if(VI %in% measures){ # a distance value
		result[VI] <- igraph::compare(partition1, partition2, method="vi")
	}
	if(NVI %in% measures){ # a distance value
	    # this normalization is proposed by the following paper:
	    #   Vinh et al. "Information Theoretic Measures for Clusterings Comparison: Variants, Properties, Normalization and Correction for Chance". 2010
	    # The author of the measure VI is also discussed this normalization in her book: 
	    #   Meila. "Handbook of Cluster Analysis, Chapter 27: Criteria for Comparing Clusterings". 2015
	    # Note that there are also other propositions for normalization, such as log(n), which does not allow to compare datasets of different size.
		result[NVI] <- igraph::compare(partition1, partition2, method="vi")
		# then, normalize the result to obtain a value between 0 and 1
		joint.entropy = entropy::entropy.plugin(table(partition1,partition2)) # based on Vinh et al. 2010
		result[NVI] = result[NVI]/joint.entropy
	}
	if(EDIT %in% measures){
		rel.partition2 = create.relative.plot.membership(partition1, partition2)
		d = hamming.distance(partition1,rel.partition2)
		result[EDIT] = d
	}
	if(EDIT.NORM %in% measures){
		rel.partition2 = create.relative.plot.membership(partition1, partition2)
		d = hamming.distance(partition1,rel.partition2)/length(partition1)
		result[EDIT.NORM] = d
	}
	
	# -----------------------------------------------------------
	
	if(RI %in% measures){ # a similarity value
		res = igraph::compare(partition1, partition2, method="rand")
		result[RI] <- 1-res # dissimilarity ~ distance
	}
	if(HA.ARI %in% measures){ # a similarity value
		res = clues::adjustedRand(partition1, partition2)["HA"]
		if(res<0)res=0 # in case ARI yields negativ value => it is really rare, especially in our application
		result[HA.ARI] = 1-res # dissimilarity ~ distance
	}
#	if(MA.ARI %in% measures){ # a similarity value
#		res = clues::adjustedRand(partition1, partition2)["MA"]
#		if(res<0)res=0 # in case ARI yields negativ value => it is really rare, especially in our application
#		result[MA.ARI] = 1-res # dissimilarity ~ distance
#	}
	if(FM %in% measures){ # a similarity value
		res = clues::adjustedRand(partition1, partition2)["FM"]
		result[FM] = 1-res # dissimilarity ~ distance
	}
	if(JACCARD %in% measures){ # a similarity value
		res = clues::adjustedRand(partition1, partition2)["Jaccard"]
		result[JACCARD] = 1-res # dissimilarity ~ distance
	}
	if(HARMONIC.MEAN.PURITY %in% measures){ # a similarity value
		res1 = NMF::purity(as.factor(partition1), as.factor(partition2))
		res2 = NMF::purity(as.factor(partition2), as.factor(partition1))
		res = (2*res1*res2)/(res1+res2)
		result[HARMONIC.MEAN.PURITY] = 1-res # dissimilarity ~ distance
	} 
	if(NMI %in% measures){
		res = igraph::compare(partition1, partition2, method="nmi")
		result[NMI] <- 1-res # dissimilarity ~ distance
	}
#	if(AMI %in% measures){
#		# we need to use octave tool in order to get the result
#		# we provide octave with a string which contains 2 vector def + function call
#		mem1.str = paste0("[",paste(partition1, collapse=" "),"]")
#		mem2.str = paste0("[",paste(partition2, collapse=" "),"]")
#		# addpath('../') in order to recgnize the location of AMI.m file
#		# do not put ';' at the end so that we can get the result of the function
#		eval.str = paste0("addpath(\"src/evaluate-partitions\");","mem1=",mem1.str,"; mem2=",mem2.str,";AMI(mem1,mem2)")
#		# parse the string result
#		cmd = paste0("octave --no-gui --silent --eval '", eval.str,"'")
#		octave.answer = system(cmd, wait=TRUE, intern = TRUE)
#		# octave answer is something like that: "ans =   0.5"
#		tmp = gsub("ans =", "", octave.answer)
#		# remove whitespaces
#		str.res = gsub(" ", "", tmp)
#		# TODO: as.numeric() handles output, for instance as.numeric("Inf") outputs Inf
#		res = as.numeric(str.res)
#		result[AMI] <-  1-res # dissimilarity ~ distance
#	}
	
	# TODO one can add the processing of other measures here if needed
	
	return(result)
}






#################################################################
# It is a helper function to compute dissimilarity/distance scores between a partition and a set of partitions. There are two use cases:
#   - computing dissimilarity for each pair of partitions obtained by the same partitioning method.
#       In this case, the final dissimilarity matrix that we will obtain is symmetric, so 'is.mtrx.symmetric' should be set to TRUE.
#       The parameter 'nb.part2' represents the number of partitions.
#       The parameter 'part.folder2' have no sense here. So, 'part.folder2' and 'part.folder1' are pointed to the same folder
#   - computing dissimilarity between a partition obtained by a partitioning method X and all partitions obtained by another partitioning method Y
#       So, 'is.mtrx.symmetric' should be set to FALSE
#       The parameters 'nb.part2' and 'part.folder2' are related to the method Y
# This function uses the package 'bigmemory' to handle matrix reading and writing. And, it writes the computed dissimilarity values into the matrix
#   accessed by 'dget(BIG.MEMORY.MATRIX.DESCRIPTOR.FILENAME)'. So, the matrix is not given as input the parameter.
# In the end, it performs 'flush', deletes used variables to reduce used RAM.
# A final note about simultaneous writing during parallel mode => source: http://www.stat.yale.edu/~mjk56/temp/bigmemory-vignette.pdf
#   all writing oprerations are user-friendly, i.e locking meachanism is handled implicitely
#
# measure: an evaluation measure (distance or dissimilarity) based on which the dissimilarity between two partitions is calculated
# i: row 'i' of the final matrix, i.e. i.th partition obtained by the method X
# nb.part2: the number of partitions with method Y if is.mtrx.symmetric=TRUE. Otherwise, it represents the number of all partitions
# is.mtrx.symmetric: boolean. TRUE if the first use case is used
# part.folder1: the folder from which the partitions are accessed
# part.folder2: the same value as 'part.folder1' should be provided if  is.mtrx.symmetric=TRUE
#
# returns nothing
################################################################
compute.distance.scores.by.row = function(measure, i, nb.part2, is.mtrx.symmetric, part.folder1, part.folder2){
	
	table.file = file.path(part.folder1, paste0(MBRSHP.FILE.PREFIX,i,".txt"))
	part1 <- as.numeric(as.matrix(read.table(file=table.file, header=FALSE)))
	
	require(bigmemory)
	require(itertools)
	mdesc <- dget(BIG.MEMORY.MATRIX.DESCRIPTOR.FILENAME)
	big.dist.mtrx = NA
	
	it <- ihasNext(iter(seq(0,nb.part2-1)))
	while(hasNext(it)){
		j = nextElem(it)
		# when it is rectangular matrix, process each case of the current row. Otherwise, process upper triangle part of the matrix
		if(!is.mtrx.symmetric || i < j){ # if the matrix is rectangular (i.e. non-symmetric), it enters into this code block in any case
			table.file = file.path(part.folder2, paste0(MBRSHP.FILE.PREFIX,j,".txt"))
			part2 <- as.numeric(as.matrix(read.table(file=table.file, header=FALSE)))
			result = round(compare.partition.pair.as.distance(part1, part2, measure), digits=4)
					
			big.dist.mtrx <- attach.big.matrix(mdesc) # access matrix
			big.dist.mtrx[i+1,j+1] = as.numeric(result[measure])
		}
	}
	
	if(is.big.matrix(big.dist.mtrx))
		flush(big.dist.mtrx) # For a file-backed big.matrix object, flush() forces any modified information to be written to the file-backing.
	
	rm(j)
	rm(it)
	rm(part1)
	rm(table.file)
	rm(mdesc)
	rm(big.dist.mtrx)
	gc() # garbage collector
	return() # to reduce memory
}


#################################################################
# It computes the dissimilarity score for each pair of partitions. It covers two use cases:
#   - computing dissimilarity for each pair of partitions obtained by the same partitioning method, i.e. partitions in the same partition space
#       In this case, the final dissimilarity matrix that we will obtain is symmetric, so 'is.mtrx.symmetric' should be set to TRUE.
#       The parameters 'nb.part2' and 'part.folder2' have no sense here, so no need to provide them.
#   - computing dissimilarity between a partition obtained by a partitioning method X and all partitions obtained by another partitioning method Y
#       So, 'is.mtrx.symmetric' should be set to FALSE
#       The parameters 'nb.part2' and 'part.folder2' are related to the method Y
# This function requires a lot of RAM memory when the size of the partitions is large (i.e. > 1000)
# We suppose that the distance measure is symetric => dist(part1,part2) = dist(part2, part1)  
# We make use of an external library, called 'bigmemory', to store large matrix in a shared memory.
# We do not need to use lock routines when accessing the matrix. This is handled via the package 'bigmemory'.
#
# In parallel mode: chunk size: 250, nb core for parallelism: 8, nb.part1: 10.000, nb.part2: 10.000 ==> requires 33 GB of RAM
#	 When you reduce the number of core, the required RAM memory is reduced as well
#
# measure: an evaluation measure (distance or dissimilarity) based on which the dissimilarity between two partitions is calculated
# part.folder1: The folder from which the partitions are accessed. It represents the rows of the final dissimilarity matrix.
# part.folder2: For the 2nd use case. It represents the columns of the final dissimilarity matrix.
# nb.part1: the number of all partitions
# nb.part2: For the 2nd use case. It represents the number of all partitions obtained by the method Y. It will be omitted if is.mtrx.symmetric=TRUE
# algo.name1: the folder from which the partitions are accessed
# algo.name2: For the 2nd use case. It will be omitted if is.mtrx.symmetric=TRUE
# is.mtrx.symmetric: it should be set to TRUE for the 1st use case. when TRUE, 'part.folder2' and 'nb.part2' will not be used
# parallel.mode: Boolean. Whether to create matrix in parallel mode or not
# nb.core: Used only in parallel mode
# chunk.size: Used only in parallel mode. It should be less than 'nb.part1' (since parallelizing by row)
# add.colnames: Boolean. Adding column names into the final dissimilarity matrix
# add.rownames: Boolean. Adding row names into the final dissimilarity matrix
#
# returns the dissimilarity matrix 
#################################################################
create.distance.matrix = function(measure, part.folder1, part.folder2=NA, nb.part1, nb.part2=NA, algo.name1, algo.name2=NA, 
		is.mtrx.symmetric=TRUE, parallel.mode=FALSE, nb.core=NA, chunk.size=NA, add.rownames=FALSE, add.colnames=FALSE){
	
# You can test 'foreach dopar'
# cl <- makeCluster(2, outfile="") # manually set the number of processors to use
# registerDoParallel(cl)
# iters = iter(seq(1,10))
# n=10
# dist.mtrx2 = foreach(i1=iters, .combine='rbind', .packages=c("itertools2"), .export=c("n")) %dopar% {
#	curr.row = sapply(iter_deepcopy(iters), function(i2) paste(i1,i2,sep=",") )
#	return(c(rep(NA, (n - length(curr.row))), curr.row))
# }
# stopCluster(cl)
	
	library(parallel)
	library(foreach)
	library(doParallel)

	
	if(is.mtrx.symmetric){
		part.folder2 = part.folder1
		nb.part2 = nb.part1
		algo.name2 = algo.name1
		tlog(40, "[symmetric distance matrix] nb partition in folder '", part.folder1,"': ", nb.part1)	
	} else{
		tlog("[unsymmetric distance matrix] nb partition in folder '", part.folder1,"': ", nb.part1, ",\tnb partition in folder '", part.folder2,"': ", nb.part2)	
	}
	
	if(nb.part1 < 1 && nb.part2 < 1){
		tlog(40, "Error in create.distance.matrix(). One of the partiton folder is empty")
		return(0)
	}
	
	if(parallel.mode){
		if(is.na(nb.core))
			nb.core = 6 # detectCores()-1
		if(is.na(chunk.size))
			chunk.size=100
	}
	
	library(iterators)
	library(bigmemory)
	big.dist.mtrx = big.matrix(nrow=nb.part1, ncol=nb.part2, type="double", init=NA,
			descriptorfile=BIG.MEMORY.MATRIX.DESCRIPTOR.FILENAME, backingfile=BIG.MEMORY.MATRIX.BACKING.FILENAME)
	
	# -----------------------------------------------------------------------------------------
	
	if(parallel.mode){ # parallel mode with chunk strategy (by row)	
		last.indx = 0
		while(last.indx < nb.part1){
			first.indx = last.indx + 1
			remaining = nb.part1-last.indx
			last.indx = last.indx + chunk.size
			if(remaining < chunk.size)
				last.indx = nb.part1
						
			tlog("chunk: [",first.indx,",",last.indx,"] (total=",nb.part1,")")
			# --------------------------------------------------------
			
			# source: https://stackoverflow.com/questions/24327137/error-in-unserializesocklistn-error-reading-from-connection-on-unix
			cl <- makeCluster(nb.core, outfile="") # manually set the number of processors to use
			registerDoParallel(cl)
#			clusterExport(cl=cl,
#					varlist=c(
#							"compare.partition.pair.as.distance", "compute.distance.scores.by.row" 
#							#"VI", "Nmi", "RI", "HA's ARI", "MA's ARI", "Jaccard", "FM", "HarmonicMeanPurity", "MBRSHP.FILE.PREFIX",
#							#"BIG.MEMORY.MATRIX.DESCRIPTOR.FILENAME"
#							),
#					envir=globalenv()
#			)
			
			foreach(i=iter(seq(first.indx-1,last.indx-1))) %dopar% { # since partition files ids start from 0
			# .packages=c("iterators", "itertools")
				source("src/define-constants.R")
				source("src/evaluate-partitions/compare-partitions.R")
				source("src/evaluate-partitions/create-relative-plot-membership.R")
				source("src/cluster-analysis/define-purity.R")
				compute.distance.scores.by.row(measure, i, nb.part2, is.mtrx.symmetric, part.folder1, part.folder2)			
			}
			
			stopCluster(cl)
		}
	} else { # sequantial mode
		foreach(i=iter(seq(0, nb.part1 - 1))) %do% { # since partition files ids start from 0
			compute.distance.scores.by.row(measure, i, nb.part2, is.mtrx.symmetric, part.folder1, part.folder2)	
		}
	}
				
	# -----------------------------------------------------------------------------------------
	
	dist.mtrx = as.matrix(big.dist.mtrx)
	if(!is.matrix(dist.mtrx)) # it might be either NA, a numeric value or a vector
		dist.mtrx = matrix(dist.mtrx, nb.part1, nb.part2)
	
	
	if(is.mtrx.symmetric){ # when it is not a rectangular matrix
		# We want to fill in the lower triangle from upper part.
		# symmetry
		for(i in 1:nb.part1){ # row
			dist.mtrx[i,i] = 0 # 0 as distance for identical partitions
			for(j in 1:nb.part2){ # column
				if(j<i)
					dist.mtrx[i,j] = dist.mtrx[j,i]
			}
		}
	}
	
	

	
	if(!add.rownames)
		rownames(dist.mtrx) = paste(algo.name1, "sol", seq(0, nb.part1 - 1)) # "algoname sol 0", "algoname sol 1", ..
	if(!add.colnames)
		colnames(dist.mtrx) = paste(algo.name2, "sol", seq(0, nb.part2 - 1))
	
	return(dist.mtrx)
}

