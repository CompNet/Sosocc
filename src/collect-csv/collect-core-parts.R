# TODO: Add comment
# 
# Author: nejat
###############################################################################





#################################################################
# It constructs a single row of the csv file in output. It retrieves the statistics 
#   related to core parts for the given values of the input parameters. 
#
# n: graph order
# l0: number of modules
# d: density
# prop.mispl: proportion of misplaced links
# prop.neg: proportion of negative links
# network.no: network id (the identifiers start from 1
# cor.clu.exact.algo: the name of correlation clustering algorithm to run
# comp.measures: comparison/evaluation measures used to obtain core parts.
# force: whether or not the existing files are overwritten by a fresh call of all corresponding methods (e.g partitioning method)
#
##################################################################
retreive.core.parts = function(n, l0, d, prop.mispl, prop.neg, network.no,
		cor.clu.exact.algo, comp.measures, force)
{
	
	data = c()
	
	#for(algo.name in cor.clu.exact.algos){
	tlog(20, "plot transition networks => algo.name: ", cor.clu.exact.algo)
	
		for(graph.desc.name in c(SIGNED.UNWEIGHTED.FILE)){ # OPPOSITE.SIGNED.UNWEIGHTED.FILE
			
			net.folder = get.input.network.folder.path(n, l0, d, prop.mispl, prop.neg, network.no)
			part.folder = get.part.folder.path(n, l0, d, prop.mispl, prop.neg, network.no, cor.clu.exact.algo, graph.desc.name)
			eval.algo.folder = get.eval.folder.path(n, l0, d, prop.mispl, NA, prop.neg, k=NA, network.no, cor.clu.exact.algo, graph.desc.name)
			clu.analysis.folder = get.clu.analysis.folder.path(n, l0, d, prop.mispl, prop.neg, network.no, cor.clu.exact.algo, graph.desc.name)
			clu.characterization.folder = get.clu.characterization.folder.path(n, l0, d, prop.mispl, prop.neg, network.no, cor.clu.exact.algo, graph.desc.name)
			
			if(dir.exists(net.folder)){
				
				partitions = load.membership.files(part.folder)
				m = length(partitions) # nb partition
				
				for(measure in comp.measures){
					tlog(20, "performing cluster characterization => measure: ", measure)
			
					#matrix.file = file.path(eval.algo.folder,paste0(EVAL.DIST.MATRIX.FILE.PREFIX,"-",NVI,".csv")) # ehre 
					#dist.mtrx = as.matrix(read.csv(matrix.file, row.names = 1, header= TRUE, check.names=FALSE)) # square matrix
					#dist.mtrx = matrix(as.numeric(sprintf("%.4f", dist.mtrx)), m, m)
						
					clu.analysis.measure.folder = file.path(clu.analysis.folder,measure)
					clu.characterization.measure.folder = file.path(clu.characterization.folder,measure)				
					
					silhouette.scores.file = file.path(clu.analysis.measure.folder, paste0(SILH.SCORE.FILENAME,".csv"))
					clu.result = read.csv(file=silhouette.scores.file, header=TRUE)
					avg.sil.scores = clu.result[["avg.sil.score"]]
					m = length(avg.sil.scores) # nb partitions
					k.values = clu.result[["k"]]
					silh.folder.names = prepare.cluster.analysis.result.folders(avg.sil.scores, k.values, "silh")
			
					# ====================================================================================
					max.indx = 1
					k = 1
					if(length(avg.sil.scores) == 2){
					    d = compare.partition.pair.as.distance(partitions[[1]], partitions[[2]], measures=NVI)
					    if(d >= 0.50){
    						max.indx = 2
    						k = 2
					    }
					}
					if(length(avg.sil.scores) > 2 &&  any( (as.numeric(avg.sil.scores) > 0.50)[-c(1,m)] )){
						max.indx = which.max(as.numeric(avg.sil.scores))
						k = k.values[max.indx]
					}	
					# ====================================================================================
					
					
					print(avg.sil.scores)
					print(k)
					if(length(avg.sil.scores) != 1){ # if there are more than 1 solution
						silh.folder.name = silh.folder.names[max.indx]
						
						analysis.silh.folder = file.path(clu.analysis.measure.folder,silh.folder.name)
						clu.analysis.mbrshp = read.table(file=file.path(analysis.silh.folder,paste0(MBRSHP.FILE.PREFIX,".txt")))$V1
						characterization.silh.folder = file.path(clu.characterization.measure.folder,silh.folder.name)
						
						curr.line.all = c()
						for(core.part.threshold in CORE.PART.THRESHOLDS){
							cat("core parth threshold: ", core.part.threshold, "\n")
							clu.data = c() # put membership vecgtor as column into the matrix
							curr.line = c()
							for(clu.no in seq(1,k)){ # for each cluster of similar solutions
								cat("clu no: ", clu.no, "\n")
								characterization.core.folder = file.path(characterization.silh.folder, paste0(CORE.PART.FOLDER.NAME,"=",sprintf("%.2f",core.part.threshold)))
								fpath = file.path(characterization.core.folder, paste0(MBRSHP.FILE.PREFIX,"-",CLUSTER.ID.FOLDER.PREFIX,clu.no,".txt"))
								mem = read.table(file=file.path(characterization.core.folder, paste0(MBRSHP.FILE.PREFIX,"-",CLUSTER.ID.FOLDER.PREFIX,clu.no,".txt")))$V1
								clu.data = cbind(clu.data,mem)	
								
								curr.line = c( curr.line, paste0("clu",clu.no,":",length(which(!is.na(mem))) ))		
							}
							# make it a string, then put into a vector
							curr.line = c(paste(curr.line, collapse=";"))
							
							# get nb common unstable nodes
							if(k!=1){
								vec = sapply(1:nrow(clu.data), function(i) all(is.na(clu.data[i,])) )
								curr.line = c( curr.line, paste0("all-unstable",":",length(which(vec == TRUE)) ))
							} else {
								curr.line = c( curr.line, NA)
							}
						
							# get nb stable node in general core part (so, intersection of all solution classes)
							fpath = file.path(characterization.core.folder, paste0(MBRSHP.FILE.PREFIX,"-","clu",ALL,".txt"))
							if(file.exists(fpath)){
								general.mem = read.table(file=fpath)$V1
								curr.line = c( curr.line, paste0("generalCore",":",length(which(!is.na(general.mem))) ))	
							} else {
								curr.line = c( curr.line, NA )
							}

							curr.line = paste(curr.line, collapse="---")
							curr.line.all = c(curr.line.all, curr.line)
						}
						
						line.matrix = matrix(curr.line.all, nrow=1)
						rownames(line.matrix) = paste0("n=", n, "_l0=", l0, "_d=", d, "_prop.mispl=", prop.mispl, "_prop.neg=", prop.neg, "_network.no=", network.no)
						data = rbind(data, line.matrix)
					
					}
					
				}
			
			}
		}
		
	#}
	
	
	return(data)
}





#################################################################
# It is the starting method in the aim of collecting the statistics related to core parts. 
#   It handles all networks by graph.sizes,  prop.mispls, my.prop.negs and in.rand.net.folders
#   For the given value of the parameter l0, it generates a csv file.
#
# graph.sizes: a vector of values regarding graph orders to be considered
# d: density (it is a single value)
# l0: number of modules to be considered (it is a single value)
# prop.mispls: a vector of values regarding proportion of misplaced links
# prop.negs: a vector of values regarding proportion of negative links (for now, it is not operational)
# in.rand.net.folders: a vector of values regarding input random graph folders. Sequantial integers (1, .., 10)
# cor.clu.exact.algo: the name of correlation clustering algorithm to run
# comp.measures: comparison/evaluation measures used to obtain core parts.
# force: whether or not the existing files are overwritten by a fresh call of all corresponding methods (e.g partitioning method)
#
##################################################################
collect.core.parts = function(graph.sizes, d, l0, prop.mispls, prop.negs, in.rand.net.folders,
		cor.clu.exact.algo, comp.measures, force)
{
	all.data = c()
	
#	fin =c()
	tlog("starts collecting core parts")
	for(n in graph.sizes){
		tlog(8, "partitioning networks => n: ", n)
		
		for(prop.mispl in prop.mispls){
			tlog(8, "collecting core parts => prop.mispl: ", prop.mispl)
			
		    if(is.na(prop.negs) && d == 1){
		        prop.negs = compute.prop.neg(n, d, l0, prop.mispl)
		    }
			
			for(prop.neg in prop.negs){
				tlog(12, "collecting core parts => prop.neg: ", prop.neg)
				
#				temp = c()
				
				for(network.no in in.rand.net.folders){
					tlog(16, "collecting core parts => network.no: ", network.no)
					
					data = retreive.core.parts(n, l0, d, prop.mispl, prop.neg, network.no,
							cor.clu.exact.algo, comp.measures, force)
					if(length(data)>0)
						all.data = rbind(all.data, data)
#					temp = c(temp, data)
				}
				
#				fin = c(fin, paste0("n=",n,", prop.mispl=",prop.mispl,", network.no=",network.no,", max=", max(temp)))
				
			}
			
		}
		
	}
#	print(fin)
	
	if(length(all.data)>0){
		colnames(all.data) = paste("core part threshold: ", sprintf("%.2f",CORE.PART.THRESHOLDS))
		#rownames(all.conn.comp.data) = NULL
		if(!dir.exists(OUTPUT.CSV.FOLDER))
		    dir.create(OUTPUT.CSV.FOLDER, recursive=FALSE, showWarnings=FALSE)
		
		write.csv(file=file.path(OUTPUT.CSV.FOLDER,paste0("core-part-l0=",l0,".csv")), x=all.data)
	}
}
