


#################################################################
# It computes the average distance (and corresponding standard deviation) between one solution and other solutions.
#   The distance matrix contains already the distance values between each pair of solutions.
#   The average distance can be computed by processing either by column or by row. The only thing to do is to remove
#       distance from the solution to itself, which is of course zero, otherwise, it will reduce erroneously the average distance 
#   This is a generic method, it may handle different situations. For instance, if we run two different partitioning methods,
#       we can compute the average distance between the solutions of the first method and the ones of the second one ==> not symmetric
#
# dist.mtrx: distance matrix between solutions. It is symmetric if we consider only the solutions of a single partitioning method.
#               Otherwise, it is not symmetric if we consider two different partitioning methods
# is.mtrx.symmetric: whether the distance matrix is symmetric. In this case, it should be sylletric
# m: Number of partition, i.e. solutions
# by.row: calculating average by taking distance values by row from the distance matrix
# by.col: calculating average by taking distance values by column from the distance matrix
#
#################################################################
compute.mean.and.sd.dist.scores = function(dist.mtrx, is.mtrx.symmetric, m, by.row=TRUE, by.col=FALSE){
	
	mean.sd.result = c()
	for(i in seq(1,m)){
		
		dist.values = NA
		if(by.col)
			dist.values = dist.mtrx[,i]
		else
			dist.values = dist.mtrx[i,]
		
		# if(is.mtrx.symmetric && length(dist.values)>1)
		if(is.mtrx.symmetric) # dist.values[i]=1 when symmetric, exclude it
			dist.values = dist.values[-i] # if there is only 1 element, this causes to have 'Nan' after mean()
		
		mn = mean(dist.values)
		s = sd(dist.values)
		

		# ------------------
		# if there is only 1 solution, 'dist.values' will be empty since we remove the solution itself. And mean of something empty will cause 'Nan'.
		# To stop cascade effect of Nan for further steps, convert Nan into NA. That way, we can easily remove NA values from mean() with its 'na.rm' param
		if(is.nan(mn))
			mn = NA
		if(is.nan(s))
			mn = NA
		# ------------------
		mean.sd.result = rbind(mean.sd.result, c(sprintf("%.4f",mn),sprintf("%.4f",s)))
	}
	colnames(mean.sd.result) = c(MEAN.COL.NAME, SD.COL.NAME)

	return(mean.sd.result)
}





#################################################################
# This method produces from all optimal solutions the following result files:
# - EVAL.NB.SOL.FILENAME: writing the total number of solution
# - EVAL.EXEC.TIME.FILENAME: writing the execution time spent in ExCC
# - EVAL.CLU.INFO.TABLE.FILENAME: writing cluster numbers obtained in the partitions
# - EVAL.IMB.INFO.TABLE.FILENAME: writing imbalance information: count and percentage
# - EVAL.DIST.MATRIX.FILE.PREFIX,"-",measure: computing the distance between each pair of solutions for a given distance measure
# - EVAL.MEAN.SD.DIST.SCORES.FILE.PREFIX,"-",measure: once having a distance matrix, 
#                                                   computing the mean of distance scores for each solution with the rest
#
# Note that in the method 'create.distance.matrix()', there is a threshold value to run the method in parallel mode.
#   You can change this threshold value in the 'define-constants.R' file
#
#
# eval.folder: evaluation folder in which the result fiels are stored
# part.folder: partition folder in which the partitions was stored
# algo.name: the name of correlation clustering algorithm that has been executed for optimal solutions
# net.folder: the network id in  which the evaluation process is performed
# graph.desc.name: graph description, i.e. SIGNED.UNWEIGHTED.FILE, SIGNED.WEIGHTED.FILE, etc.
# comp.measures: distance measure to compare the optimal solutions
# force: whether or not the existing files are overwritten by a fresh call of all corresponding methods (e.g partitioning method)
#
#################################################################
process.algo.evaluation = function(eval.folder, part.folder, algo.name, net.folder, graph.desc.name, 
                                    comp.measures, force){

    graph.name = paste0(graph.desc.name,".G")
	
	is.mtrx.symmetric = TRUE # this comes by definition of the method
	mbrshps = load.membership.files(part.folder)
	m = length(mbrshps) # nb partition
	print(m)
	if(m>0){ # if there is at least 1 optimal solution
		network.path = file.path(net.folder,graph.name)
		g = read.graph.ils(network.path)
		n = vcount(g) # nb node
	
		
		# ---------------------------------------------------------------------
		# nb solution
		tlog(24, "process algo evaluation: nb solution")
		table.file = file.path(eval.folder, paste0(EVAL.NB.SOL.FILENAME,".csv"))
		mtrx = matrix(m, 1, 1) # put the value into matrix, that way we will process all results in the same way (i.e. via 'csv')
		colnames(mtrx) = NB.SOL.COL.NAME
		rownames(mtrx) = algo.name
		write.csv(x=mtrx, file=table.file, row.names=TRUE)
		# ---------------------------------------------------------------------

		
		# ---------------------------------------------------------------------
		# exec time
		table.file = file.path(eval.folder, paste0(EVAL.EXEC.TIME.FILENAME,".csv"))
		if(!file.exists(table.file) || force){	
		    # read txt file, and write it into a csv file
			result.f.path = file.path(part.folder, EXEC.TIME.FILENAME)
			# there is a single 'exec time' value
			exec.time = read.table(result.f.path)$V1
			mtrx = matrix(exec.time, 1, 1) # matrix with 1 element
			colnames(mtrx) = c(EXEC.TIME.COL.NAME)
			rownames(mtrx) = algo.name
			write.csv(x=mtrx, file=table.file, row.names=TRUE)
		}
		# ---------------------------------------------------------------------
		
		
		# ---------------------------------------------------------------------
		# TODO: used RAM
		#
		# ---------------------------------------------------------------------
		
		
		# ---------------------------------------------------------------------
		# clusters
		# we create matrix of size m x n where is m: nb solutions and n: nb nodes
		# since a matrix should be fixed size, and the number of clsuters can vary across solutions, we put NA into (n-nb.clu) case
		tlog(24, "process algo evaluation: clusters")
		table.file = file.path(eval.folder, paste0(EVAL.CLU.INFO.TABLE.FILENAME,".csv"))
		if(!file.exists(table.file) || force){
			clu.result = c()
			for(mbrshp in mbrshps){
				clu.freq = table(mbrshp) # clu freq in terms of nb cluster
				nb.clu = length(clu.freq)
				remaining.freq = rep(NA, n-nb.clu) # at worst, there will be 'n' clusters
				clu.result = rbind(clu.result, c(nb.clu, clu.freq, remaining.freq))
			}
			colnames(clu.result) = c(NB.CLU.COL.NAME, paste(CLUSTER.COL.NAME.PREFIX, seq(1,n), sep=" "))# e.g. "nb cluster", "clu1", "clu2", etc.
			rownames(clu.result) = paste(algo.name, "sol", seq(0, m-1)) # e.g. "ExCC-All sol0", "ExCC-All sol1", etc.
			write.csv(x=clu.result, file=table.file, row.names=TRUE)
		}
		# ---------------------------------------------------------------------
		
		
		# ---------------------------------------------------------------------
		# imbalance value and proportion
		tlog(24, "process algo evaluation: imbalance")
		table.file = file.path(eval.folder, paste0(EVAL.IMB.INFO.TABLE.FILENAME,".csv"))
		if(!file.exists(table.file) || force){
			imb.result = c()
			my.mbrshps = mbrshps
			# if it is the result of an exact algo, take just the first partition, which is enough for imbalance
			# Later, we need to collect these values and take their average for 'summary' folder, average should be correct with only 1 value here
			mbrshp = mbrshps[[1]]

			imb.count = as.numeric(compute.imbalance.from.membership(g, mbrshp, output.type = "count"))
			imb.perc = as.numeric(compute.imbalance.from.membership(g, mbrshp, output.type = "percentage"))

			imb.perc = sprintf("%.4f",imb.perc)
			

			imb.result = rbind(imb.result, c(imb.count, imb.perc))
			colnames(imb.result) = c(IMB.COUNT.COL.NAME, IMB.PERC.COL.NAME)
			rownames(imb.result) = algo.name

			write.csv(x=imb.result, file=table.file, row.names=TRUE)
		}
		# ---------------------------------------------------------------------

		
		# ---------------------------------------------------------------------
		# distance measures
		tlog(24, "process algo evaluation: distance")
		for(measure in comp.measures){
			mtrx.file = file.path(eval.folder, paste0(EVAL.DIST.MATRIX.FILE.PREFIX,"-",measure,".csv"))
			tlog(28, "process algo evaluation: distance matrix for the measure: ", measure)
	
			print(mtrx.file)
			if(!file.exists(mtrx.file) || force){
				dist.mtrx = NA
				
				# ===========================================================================
				# This is the fast way to compute the normalization of VI and EDIT distances.
				# We do not want to spend much time to compute their normalizations, if the raw distance are already computed
				done = FALSE
				if(measure == EDIT.NORM){
					edit.mtrx.file = file.path(eval.folder, paste0(EVAL.DIST.MATRIX.FILE.PREFIX,"-",EDIT,".csv"))
					if(file.exists(edit.mtrx.file)){
						edit.dist.mtrx = as.matrix(read.csv(edit.mtrx.file, row.names = 1, header= TRUE, check.names=FALSE))
						dist.mtrx = matrix(as.numeric(edit.dist.mtrx),nrow=nrow(edit.dist.mtrx), ncol=ncol(edit.dist.mtrx))/n
						done = TRUE
					}
				}
				# else if(measure == NVI){ # OLD CODE => WE NORMALIZE VI BY JOINT ENTROPY
				#                           NORMALIZATION WITH LOG(N) DOES NOT ALLOW TO COMPARE DATASETS OF DIFFERENT SIZES
				# 	vi.mtrx.file = file.path(eval.folder, paste0(EVAL.DIST.MATRIX.FILE.PREFIX,"-",VI,".csv"))
				# 	if(file.exists(vi.mtrx.file)){
				# 	    vi.dist.mtrx = as.matrix(read.csv(vi.mtrx.file, row.names = 1, header= TRUE, check.names=FALSE))
				# 		dist.mtrx = matrix(as.numeric(vi.dist.mtrx),nrow=nrow(vi.dist.mtrx), ncol=ncol(vi.dist.mtrx))/log(n)
				# 		done = TRUE
				# 	}
				# }
				# ===========================================================================
				
				if(!done){
					par.mode = TRUE
					nb.core = PAR.MODE.NB.CORE.DEFAULT
					chunk.size = PAR.MODE.CHUNK.SIZE.DEFAULT
					if(m < PAR.MODE.THRESH.NB.MEM.FILE){ # if there is just a small amount of files, do it sequantially
						par.mode = FALSE
						nb.core = NA
						chunk.size = NA
					}
					dist.mtrx = create.distance.matrix(measure, part.folder1=part.folder, nb.part1=m, algo.name1=algo.name, 
							is.mtrx.symmetric, parallel.mode=par.mode, nb.core=nb.core, chunk.size=chunk.size)
				}
				
				write.csv(x=dist.mtrx, file=mtrx.file)
			}
			
			
			
			tlog(28, "process algo evaluation: mean & sd values for the measure: ", measure)
			# mean & sd by algo.name: compute the mean of distance scores for each solution with the rest
			table.file = file.path(eval.folder, paste0(EVAL.MEAN.SD.DIST.SCORES.FILE.PREFIX,"-",measure,".csv"))
			if(!file.exists(table.file) || force){
				mtrx.file = file.path(eval.folder, paste0(EVAL.DIST.MATRIX.FILE.PREFIX,"-",measure,".csv"))
				# --------------

				dist.mtrx = read.csv(mtrx.file, row.names = 1, header= TRUE, check.names=FALSE) # sqaure matrix
				dist.mtrx = matrix(as.numeric(as.matrix(dist.mtrx)), nrow(dist.mtrx),ncol(dist.mtrx)) # make it integer matrix

				# --------------
				process.half.mtrx = startsWith(graph.name,SIGNED.UNWEIGHTED.FILE) || startsWith(graph.name,SIGNED.WEIGHTED.FILE)
				curr.values = retreive.dist.scores.from.matrix(dist.mtrx, process.half.mtrx)
				if(length(curr.values) == 0) # if there is only 1 sol, this will be empty or not based on 'process.half.mtrx'
					curr.values = NA
				curr.df = as.data.frame(c(sprintf("%.4f",curr.values)))
				rownames(curr.df) = seq(0,length(curr.values)-1)
				colnames(curr.df) = paste0(DIST.SCORE,"-",measure)
				table.file2 = file.path(eval.folder, paste0(EVAL.ALL.DIST.SCORES.FILE.PREFIX,"-",measure,".csv"))
				write.csv(x=curr.df, file=table.file2, row.names=TRUE)
				# --------------
				result = compute.mean.and.sd.dist.scores(dist.mtrx, is.mtrx.symmetric, m, by.row=TRUE, by.col=FALSE)
				rownames(result) = paste(algo.name, "sol", seq(0, m-1))
				write.csv(x=result, file=table.file, row.names=TRUE)
			}
		}
		# ---------------------------------------------------------------------

	}
}





#################################################################
# It evaluates the solutions of a given network based on the considered algorithm name and graph type (weighted or not, etc.).
#
# n: graph size
# l0: number of clusters
# d: density
# prop.mispl: proportion of misplaced links
# prop.neg: proportion of negative links
# network.no: network id (the identifiers start from 1)
# cor.clu.exact.algo: the name of correlation clustering algorithm that has been executed for optimal solutions
# force: whether or not the existing files are overwritten by a fresh call of all corresponding methods (e.g partitioning method)
#
##################################################################
evaluate.partition = function(n, l0, d, prop.mispl, prop.neg, network.no,
		cor.clu.exact.algo, comp.measures, force)
{
	net.folder = get.input.network.folder.path(n, l0, d, prop.mispl, prop.neg, network.no)
	
	tlog(16, "start to evaluate exact solutions")
	for(graph.desc.name in c(SIGNED.UNWEIGHTED.FILE)){ # SIGNED.WEIGHTED.FILE
		tlog(20, "evaluating partitions => graph.desc.name: ", graph.desc.name)
		tlog(20, "evaluating partitions => algo.name: ", cor.clu.exact.algo)
		
		e.part.folder = get.part.folder.path(n, l0, d, prop.mispl, prop.neg, network.no, cor.clu.exact.algo, graph.desc.name)
		e.eval.folder = get.eval.folder.path(n, l0, d, prop.mispl, detected.imb=NA, prop.neg, k=ALL, network.no, cor.clu.exact.algo, graph.desc.name)
		if(!dir.exists(e.eval.folder))
			dir.create(path=e.eval.folder, showWarnings=FALSE, recursive=TRUE)
		
		tlog(20, "proceed only exact solutions in folder: ", e.eval.folder)
		process.algo.evaluation(e.eval.folder, e.part.folder, cor.clu.exact.algo, net.folder, graph.desc.name, comp.measures, force)
	}
}




#################################################################
# It is the starting method in the aim of evaluating the solutions of the considered networks. 
#   It handles all networks by graph.sizes,  prop.mispls, my.prop.negs and in.rand.net.folders
#
# graph.sizes: a vector of values regarding graph sizes to be considered
# d: density (it is a single value)
# l0: number of clusters to be considered (it is a single value)
# prop.mispls: a vector of values regarding proportion of misplaced links
# prop.negs: a vector of values regarding proportion of negative links (for now, it is not operational)
# in.rand.net.folders: a vector of values regarding input random graph folders. Sequantial integers (1, .., 10)
# cor.clu.exact.algo: the name of correlation clustering algorithm that has been executed for optimal solutions
# force: whether or not the existing files are overwritten by a fresh call of all corresponding methods (e.g partitioning method)
#
##################################################################
evaluate.partitions = function(graph.sizes, d, l0, prop.mispls, prop.negs, in.rand.net.folders,
		cor.clu.exact.algo, comp.measures, force)
{
    # check the necessary condition
    if(cor.clu.exact.algo == get.ExCC.code(enum.all=FALSE)){
        tlog(4, "The partitioning method should be the one which enumerating all optimal solutions, i.e. 'ExCC-All' !")
        return(0)
    }
    
    
	tlog("starts evaluating partitions")
	for(n in graph.sizes){
		tlog(4, "evaluating partitions => n: ", n)
		
		for(prop.mispl in prop.mispls){
			tlog(8, "evaluating partitions => prop.mispl: ", prop.mispl)
			
		    if(is.na(prop.negs) && d == 1){
		        prop.negs = compute.prop.neg(n, d, l0, prop.mispl)
		    }		
			
			for(prop.neg in prop.negs){
				tlog(12, "evaluating partitions => prop.neg: ", prop.neg)
				
								
				net.folder = get.input.network.folder.path(n, l0, d, prop.mispl, prop.neg, network.no=NA)
				if(dir.exists(net.folder)){
				
					for(network.no in in.rand.net.folders){
						tlog(16, "evaluating partitions => network.no: ", network.no)
						
						evaluate.partition(n, l0, d, prop.mispl, prop.neg, network.no,
								cor.clu.exact.algo, comp.measures, force)
						
					}
				}
				
			}
			
		}
		
	}
	
}
