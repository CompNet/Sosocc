


##################################################################
# For each considered network it reads the matrix values from the file 'out.filename' and stores them in a vector. If matrix is symmetric, only lower triangle will be considered.
#  In the end, it combines the vector values of each network by row.
#  So, in the end, we have one huge vector of column. The final vector has row names so that we can distinguish which network a value belongs to. 
#  The final vector is written into file 'out.filename'.
# Below, we see an example of output:
#
#               Distance-score-VI
# network1 0	0.1251
# network1 1	0.1251
# network1 2	0.2502
# network2 0	0.6931
# network2 1	0.5022
# network2 2	0.191
# network3 0	0.0955
# network3 1	0.0955

##################################################################
collect.and.record.for.matrix.values.per.network = function(n, l0, d, prop.mispl, detected.imb.interval, prop.neg, folder.at.k.level, in.rand.net.folders,
		 algo.name, graph.desc.name, in.filename, out.filename, desc.name, force){

    result.folder = get.eval.folder.path(n, l0, d, prop.mispl, detected.imb.interval, prop.neg, folder.at.k.level, SUMMARY.FOLDER.NAME, algo.name, graph.desc.name)
    print("---")
    print(result.folder)
    
    if(!dir.exists(result.folder))
        dir.create(path=result.folder, showWarnings=FALSE, recursive=TRUE)
    tlog(24, paste0("writing values for ",desc.name," into '"), result.folder, "' folder")
    
    out.file.path = file.path(result.folder, out.filename)			
    if(!file.exists(out.file.path) || force){
        
        df = c()
    	tlog(24, paste0("taking average for ",desc.name))
    	for(i in seq(1,length(in.rand.net.folders))){
    		network.no = in.rand.net.folders[i]
    		tlog(28, paste0("taking average for ",desc.name," => network.no: "), network.no)
    		eval.folder = get.eval.folder.path(n, l0, d, prop.mispl, detected.imb.interval, prop.neg, folder.at.k.level, network.no, algo.name, graph.desc.name)
    		
    		process.half.mtrx = (graph.desc.name == SIGNED.UNWEIGHTED.FILE)
            curr.values = retreive.dist.scores.from.matrix2(eval.folder, in.filename, process.half.mtrx)
    
    		
    		if(length(curr.values) == 0) # if there is only 1 sol, this will be empty or not based on 'process.half.mtrx'
    			curr.values = NA
    		curr.df = as.data.frame(c(sprintf("%.4f",curr.values)))
    		rownames(curr.df) = paste0("network",network.no," ",seq(0,length(curr.values)-1))
    		df = rbind(df, curr.df)
    	}
    	
    	colnames(df) = desc.name
    
    	write.csv(x=df, file=out.file.path, row.names=TRUE)
	
    }
    
}		 


##################################################################
# For each considered network it reads the vector values from the column(s) 'colnames' of the file 'out.filename'.
#  In the end, it combines the vector values of each network by row.
# So, in the end, we have one huge vector of column. The final vector has row names so that we can distinguish which network a value belongs to. 
#  The final vector is written into file 'out.filename'.
# Below, we see an example of output:
#
#            nb cluster
# network1 sol0	7
# network1 sol1	6
# network1 sol2	6
# network2 sol0	6
# network2 sol1	5
# network2 sol2	6
# network3 sol0	5
# network3 sol1	6
# network4 sol0	6
# network4 sol1	5
# network4 sol2	6
# network4 sol3	6
# network4 sol4	4

##################################################################
collect.and.record.for.vector.values.per.network = function(n, l0, d, prop.mispl, detected.imb.interval, prop.neg, folder.at.k.level, in.rand.net.folders,
		 algo.name, graph.desc.name, in.filename, out.filename, colnames, force){

    result.folder = get.eval.folder.path(n, l0, d, prop.mispl, detected.imb.interval, prop.neg, folder.at.k.level, SUMMARY.FOLDER.NAME, algo.name, graph.desc.name)
    print("---")
    print(result.folder)
    
    if(!dir.exists(result.folder))
        dir.create(path=result.folder, showWarnings=FALSE, recursive=TRUE)
    tlog(24, paste("writing values for ",colnames," into '"), result.folder, "' folder")
    
    out.file.path = file.path(result.folder, out.filename)	
    if(!file.exists(out.file.path) || force){
    
    	df = c()
    	#tlog(24, paste("taking average for ",colnames))
    	for(i in seq(1,length(in.rand.net.folders))){
    		network.no = in.rand.net.folders[i]
    		tlog(28, paste("taking average for ",colnames," => network.no: "), network.no)
    		eval.folder = get.eval.folder.path(n, l0, d, prop.mispl, detected.imb.interval, prop.neg, folder.at.k.level, network.no, algo.name, graph.desc.name)
    		
    		# ------------------------------------------------------------------
    					
    		table.file = file.path(eval.folder, in.filename)
    		#print(table.file)
    		print("/*/*/!!!!!!!!!!!!!!!!!!")
    		print(table.file)
    	
    		if(file.exists(table.file))
    			result.mtrx = as.matrix(read.csv(table.file, row.names = 1, header= TRUE, check.names=FALSE))
    		else{
    			result.mtrx = matrix(NA,nrow=1, ncol=length(colnames))
    			colnames(result.mtrx) = colnames
    		}
    		
    		# curr.df = result.df
    		curr.mtrx = result.mtrx[,colnames]
    		
    		curr.mtrx = matrix(result.mtrx[,colnames], ncol=length(colnames)) # in case of 1 colname, it becomes vector => use data.frame() to make it df
    		colnames(curr.mtrx) = colnames
    		
    		
    		rownames(curr.mtrx) = paste0("network",network.no," sol",seq(0,nrow(curr.mtrx)-1))
    		if(PROP.CLASS.CORE.PART.SIZE.COL.NAME %in% colnames)
    			rownames(curr.mtrx) = paste0("network",network.no," k=",seq(1,nrow(curr.mtrx)))
    		df = rbind(df, curr.mtrx)
    		# ------------------------------------------------------------------
    	}
    	colnames(df) = colnames

    	write.csv(x=df, file=out.file.path, row.names=TRUE)
	}
	
}




#################################################################
# It is the wrapper function to get summary of the considered networks. The goal is to combine all results of the same kind into the same file.
#   It iterates over exact algorithms used to obtain partitions, graph types and dissimilarity measures.
# The list of generated files are in the following an they are consistent with 'process.algo.evaluation()' in the file 'evaluate-partitions.R': 
# - EVAL.NB.SOL.FILENAME: writing the total number of solution
# - EVAL.EXEC.TIME.FILENAME: writing the execution time spent in ExCC
# - EVAL.CLU.INFO.TABLE.FILENAME: writing cluster numbers obtained in the partitions
# - EVAL.IMB.INFO.TABLE.FILENAME: writing imbalance information: count and percentage
# - EVAL.BEST.K.FOR.KMEDOIDS: writing the best number of clusters in kmedoids
# - EVAL.CLASS.CORE.PART.SIZE.FILE.NAME: writing the core part sizes of each clsuter of kmedoids
# - EVAL.MEAN.SD.DIST.SCORES.FILE.PREFIX,"-",measure: once having a distance matrix, 
#                                                   computing the mean of distance scores for each solution with the rest
# - EVAL.DIST.MATRIX.FILE.PREFIX,"-",measure: computing the distance between each pair of solutions for a given distance measure
#
# n: graph size
# l0: number of cluster
# d: density
# prop.mispl: proportion of misplaced links
# prop.neg: proportion of negative links
# folder.at.k.level
# network.no: network id (the identifiers start from 1)
# cor.clu.exact.algo: the name of correlation clustering algorithm to run
# comp.measures: dissimilarity measure to compare the optimal solutions
# force: whether or not the existing files are overwritten by a fresh call of all corresponding methods (e.g partitioning method)
#
##################################################################
take.average = function(n, l0, d, prop.mispl, detected.imb.interval, prop.neg, folder.at.k.level, in.rand.net.folders,
	 cor.clu.exact.algo, comp.measures, force)
{
	tlog(16, "start to collect values for networks with exact algorithms")
	for(graph.desc.name in c(SIGNED.UNWEIGHTED.FILE)){ # OPPOSITE.SIGNED.UNWEIGHTED.FILE
		tlog(20, "collecting values for networks => graph.desc.name: ", graph.desc.name)
		
		
		#for(e.algo.name in cor.clu.exact.algos){
			tlog(24, "taking average over networks => algo.name: ", cor.clu.exact.algo)
			
#			# TODO temporary
#			#result.folder = get.eval.folder.path(n, l0, d, prop.mispl, prop.neg, SUMMARY.FOLDER.NAME, cor.clu.exact.algo, graph.desc.name)
#			#unlink(result.folder, recursive=TRUE)

			
			# -------------------------------------------------------------		
			# only exact solutions
			tlog(24, "collecting values for networks => only exact solutions")
			
			in.filename = paste0(EVAL.NB.SOL.FILENAME,".csv")
			out.filename = paste0(EVAL.NB.SOL.FILENAME,".csv")
			collect.and.record.for.vector.values.per.network(n, l0, d, prop.mispl, detected.imb.interval, prop.neg, folder.at.k.level,
					 in.rand.net.folders, cor.clu.exact.algo, graph.desc.name, in.filename, out.filename, NB.SOL.COL.NAME, force)
					 
			in.filename = paste0(EVAL.EXEC.TIME.FILENAME,".csv")
			out.filename = paste0(EVAL.EXEC.TIME.FILENAME,".csv")
			collect.and.record.for.vector.values.per.network(n, l0, d, prop.mispl, detected.imb.interval, prop.neg, folder.at.k.level,
					 in.rand.net.folders, cor.clu.exact.algo, graph.desc.name, in.filename, out.filename, EXEC.TIME.COL.NAME, force)
			
			in.filename = paste0(EVAL.CLU.INFO.TABLE.FILENAME,".csv")
			out.filename = paste0(EVAL.CLU.INFO.TABLE.FILENAME,".csv")
			collect.and.record.for.vector.values.per.network(n, l0, d, prop.mispl, detected.imb.interval, prop.neg, folder.at.k.level,
					 in.rand.net.folders, cor.clu.exact.algo, graph.desc.name, in.filename, out.filename, NB.CLU.COL.NAME, force)
			
			in.filename = paste0(EVAL.IMB.INFO.TABLE.FILENAME,".csv")	 
			out.filename = paste0(EVAL.IMB.INFO.TABLE.FILENAME,".csv")
			collect.and.record.for.vector.values.per.network(n, l0, d, prop.mispl, detected.imb.interval, prop.neg, folder.at.k.level, in.rand.net.folders, 
			                                                 cor.clu.exact.algo, graph.desc.name, in.filename, out.filename, c(IMB.COUNT.COL.NAME,IMB.PERC.COL.NAME), force)
			

			# # ----------
			# result.folder = get.eval.folder.path(n, l0, d, prop.mispl, detected.imb.interval, prop.neg, SUMMARY.FOLDER.NAME, cor.clu.exact.algo, graph.desc.name)
			# edit.dist.comps.table.file = file.path(result.folder, paste0(EVAL.EDIT.DIST.COMPS,".csv"))
			# edit.dist.nb.comp.table.file = file.path(result.folder, paste0(EVAL.EDIT.DIST.NB.COMP,".csv"))
			# if(file.exists(file.path(result.folder,edit.dist.comps.table.file))
			# 	&& file.exists(file.path(result.folder,edit.dist.nb.comp.table.file)) )
			# {
			# 	print("GIRDI !!!!")
			# 	# new: nb component in transition graph
			# 	in.filename = paste0(EVAL.EDIT.DIST.NB.COMP,".csv")	 
			# 	out.filename = paste0(EVAL.EDIT.DIST.NB.COMP,".csv")
			# 	collect.and.record.for.vector.values.per.network(n, l0, d, prop.mispl, prop.neg, in.rand.net.folders,
			# 		 cor.clu.exact.algo, graph.desc.name, in.filename, out.filename, 
			# 		 c(EDIT1.COMP.MEM.COL.NAME, EDIT12.COMP.MEM.COL.NAME, EDIT123.COMP.MEM.COL.NAME, EDIT1234.COMP.MEM.COL.NAME), force)
			# 	 
			# 	# based on the previous result file, create another result file
			# 	result.folder = get.eval.folder.path(n, l0, d, prop.mispl, detected.imb.interval, prop.neg, SUMMARY.FOLDER.NAME, cor.clu.exact.algo, graph.desc.name)
			# 	table.file = file.path(result.folder, paste0(EVAL.EDIT.DIST.NB.COMP,".csv"))
			# 	result.df = read.csv(table.file, row.names = 1, header= TRUE, check.names=FALSE)
			# 	prop.vals = sapply(1:ncol(result.df),
			# 			function(i){
			# 				nb.comp1 = length(which(as.numeric(result.df[,i]) == 1))
			# 				nb.all = length(which(!is.na(result.df[,i]))) # exclude NAs, that means only 1 opt sol
			# 				return(nb.comp1/nb.all)
			# 			})
			# 	res = matrix(prop.vals, nrow=1)
			# 	colnames(res) = c(EDIT1.COMP.MEM.COL.NAME, EDIT12.COMP.MEM.COL.NAME, EDIT123.COMP.MEM.COL.NAME, EDIT1234.COMP.MEM.COL.NAME)
			# 	out.file.path = file.path(result.folder, paste0(EVAL.EDIT1.DIST.NB.COMP.PROP,".csv"))	
			# 	write.csv(x=res, file=out.file.path, row.names=TRUE)
			# }
			# # ----------


			# ----------
			# new: best k for kmedoids
			for(measure in comp.measures){
			    tlog(28, "collecting values ==> measure: ", measure)
			    
    			in.filename = paste0(EVAL.BEST.K.FOR.KMEDOIDS,"-",measure,".csv")
    			out.filename = paste0(EVAL.BEST.K.FOR.KMEDOIDS,"-",measure,".csv")
    			collect.and.record.for.vector.values.per.network(n, l0, d, prop.mispl, detected.imb.interval, prop.neg, folder.at.k.level, in.rand.net.folders, 
    			              cor.clu.exact.algo, graph.desc.name, in.filename, out.filename, c(BEST.K.FOR.SILH.COL.NAME), force)
    			 
    			# based on the previous result file, create another result file
    			result.folder = get.eval.folder.path(n, l0, d, prop.mispl, detected.imb.interval, prop.neg, folder.at.k.level, SUMMARY.FOLDER.NAME, cor.clu.exact.algo, graph.desc.name)
    			table.file = file.path(result.folder, paste0(EVAL.BEST.K.FOR.KMEDOIDS,"-",measure,".csv"))
    			result.df = read.csv(table.file, row.names = 1, header= TRUE, check.names=FALSE)
    			prop.vals = sapply(1:ncol(result.df),
    					function(i){
    						nb.single.clu = length(which(as.numeric(result.df[,i]) == 1))
    						nb.all = length(which(!is.na(result.df[,i]))) # exclude NAs, that means only 1 opt sol
    						return(nb.single.clu/nb.all)
    					})
    			res = matrix(prop.vals, nrow=1)
    			colnames(res) = c(PROP.SINGLE.CLU.FOR.SILH.COL.NAME)
    			out.file.path = file.path(result.folder, paste0(EVAL.SINGLE.CLU.PROP.FOR.KMEDOIDS,".csv"))	
    			write.csv(x=res, file=out.file.path, row.names=TRUE)
			}
			# ----------


			# ----------
			# new: intersection size in core part analysis
			for(core.part.threshold in CORE.PART.THRESHOLDS){
				in.filename = paste0(EVAL.CLASS.CORE.PART.SIZE.FILE.NAME,"-tresh=",sprintf("%.2f",core.part.threshold),".csv")
				out.filename = paste0(EVAL.CLASS.CORE.PART.SIZE.FILE.NAME,"-tresh=",sprintf("%.2f",core.part.threshold),".csv")
				collect.and.record.for.vector.values.per.network(n, l0, d, prop.mispl, detected.imb.interval, prop.neg, folder.at.k.level, in.rand.net.folders,
				             cor.clu.exact.algo, graph.desc.name, in.filename, out.filename, c(PROP.CLASS.CORE.PART.SIZE.COL.NAME), force)
			 }
			# ----------	 
			 
 			 
			 
			 
			for(measure in comp.measures){
				tlog(28, "collecting values ==> measure: ", measure)
		
#				colname = MEAN.COL.NAME
				desc = MEAN.COL.NAME
				in.filename = paste0(EVAL.DIST.MATRIX.FILE.PREFIX,"-",measure,".csv")	 
				#in.filename = paste0(EVAL.MEAN.SD.DIST.SCORES.FILE.PREFIX,"-",measure,".csv")	 
				out.filename = paste0(EVAL.MEAN.SD.DIST.SCORES.FILE.PREFIX,"-",measure,".csv")
#				collect.and.record.for.vector.values.per.network(n, l0, d, prop.mispl, prop.neg,
#					 in.rand.net.folders, cor.clu.exact.algo, graph.desc.name, in.filename, out.filename, colname)
				collect.and.record.for.matrix.values.per.network(n, l0, d, prop.mispl, detected.imb.interval, prop.neg, folder.at.k.level,
					 in.rand.net.folders, cor.clu.exact.algo, graph.desc.name, in.filename, out.filename, desc, force)
					 			
				desc = paste0(DIST.SCORE,"-",measure) # provide it s conme, butuse as description
				in.filename = paste0(EVAL.DIST.MATRIX.FILE.PREFIX,"-",measure,".csv")	 
				out.filename = paste0(EVAL.ALL.DIST.SCORES.FILE.PREFIX,"-",measure,".csv")
				collect.and.record.for.matrix.values.per.network(n, l0, d, prop.mispl, detected.imb.interval, prop.neg, folder.at.k.level,
					 in.rand.net.folders, cor.clu.exact.algo, graph.desc.name, in.filename, out.filename, desc, force)
	 
			}
			# ------------------------------------------------------------

			# -------------------------------------------------------------		
		#}	
	}
	
}







#################################################################
# It is the wrapper function to get summary of the considered networks. The goal is to combine all results of the same kind into the same file.
#  It iterates over graph sizes, proportion of misplaced links, proportion of negative links and generated different network sizes.
#
# graph.sizes: a vector of values regarding graph sizes to be considered
# d: density (it is a single value)
# l0: number of clusters to be considered (it is a single value)
# prop.mispls: a vector of values regarding proportion of misplaced links
# prop.negs: a vector of values regarding proportion of negative links (for now, it is not operational)
# in.rand.net.folders: a vector of values regarding input random graph folders. Sequantial integers (1, .., 10)
# cor.clu.exact.algo: the name of correlation clustering algorithm to run
# comp.measures: dissimilarity measure to compare the optimal solutions
# force: whether or not the existing files are overwritten by a fresh call of all corresponding methods (e.g partitioning method)
#
##################################################################
take.average.over.networks = function(graph.sizes, d, l0, prop.mispls, detected.imb.intervals, prop.negs, in.rand.net.folders,
		cor.clu.exact.algo, comp.measures, force)
{
	tlog("starts taking average over networks")
	for(n in graph.sizes){
		tlog(4, "taking average over networks => n: ", n)
		
	    # ===========================================================

		for(prop.mispl in prop.mispls){
			tlog(8, "taking average over networks => prop.mispl: ", prop.mispl)
			
			# my.prop.negs = NA
			# if(is.na(prop.negs)){
			# 	prop.neg = compute.prop.neg(n, d, l0, prop.mispl)
			# 	my.prop.negs = prop.neg # we show explicitely that we get only 1 value since density=1.0
			# }		
		    eval.upper.folder = get.eval.folder.path(n, l0, d, prop.mispl, NA, prop.neg=NA, network.no=NA)
		    if(dir.exists(eval.upper.folder)){
    		    existing.folders = list.dirs(path = eval.upper.folder, full.names = FALSE, recursive = FALSE)
    		    my.prop.negs = as.numeric(gsub("propNeg=","",existing.folders))
    			
    			for(prop.neg in my.prop.negs){
    				tlog(12, "taking average over networks => prop.neg: ", prop.neg)
    				
    				net.folder = get.input.network.folder.path(n, l0, d, prop.mispl, prop.neg, network.no=NA)
    				if(dir.exists(net.folder)){
    				    
    				    upper.eval.folder = get.eval.folder.path(n, l0, d, prop.mispl, NA, prop.neg, NA, NA, NA, NA)
    				    existing.dirs = list.dirs(path = upper.eval.folder, full.names = FALSE, recursive = FALSE)
    	                # folder names between network numbers and prop.neg. Ex: c("All","k=1","k=2")
    				    #               in "/home/nejat/eclipse/workspace-neon/CC-Opti-Heur-Analysis/out/evaluate-partitions/n=20_l0=4_dens=1.0000/propMispl=0.1000/propNeg=0.7895/k=1"
    				    for(folder.at.k.level in existing.dirs){
    				        
    				        # retreive network numbers
    				        files = list.files(file.path(upper.eval.folder,folder.at.k.level))
    				        indx = which(startsWith(files,"network"))
    				        net.folders = as.integer(gsub("network=","",files[indx])) # a subset of 'in.rand.net.folders' if 'folder.at.k.level' starts with "k=". Otherwise, the whole list
        					take.average(n, l0, d, prop.mispl, NA, prop.neg, folder.at.k.level, net.folders,
        							cor.clu.exact.algo, comp.measures, force)
    				    }
    				}
    			}
		    }
			
		}
	    
	    # === end
	        
        # # retreive the detected imbalance intervals used for n, l0 and d
        # eval.upper.folder = get.eval.folder.path(n, l0, d, prop.mispl=NA, detected.imb=NA, prop.neg=NA, network.no=NA)
        # existing.folders = list.dirs(path = eval.upper.folder, full.names = FALSE, recursive = FALSE)
        # indx = which(startsWith(existing.folders, DETECTED.IMB.PROP.PARAM.NAME))
        # detected.imb.intervals = gsub(paste0(DETECTED.IMB.PROP.PARAM.NAME,"="),"",existing.folders[indx])
        
        for(detected.imb.interval in detected.imb.intervals){
            tlog(8, "evaluating partitions => detected.imb.interval: ", detected.imb.interval)
            
            eval.upper.folder = get.eval.folder.path(n, l0, d, NA, detected.imb.interval, prop.neg=NA, network.no=NA)
            
            if(dir.exists(eval.upper.folder)){
	            existing.folders = list.dirs(path = eval.upper.folder, full.names = FALSE, recursive = FALSE)
	            my.prop.negs = as.numeric(gsub("propNeg=","",existing.folders))
	            
	            for(prop.neg in my.prop.negs){
	                tlog(12, "evaluating partitions => prop.neg: ", prop.neg)
	                
	                upper.eval.folder = get.eval.folder.path(n, l0, d, NA, detected.imb.interval, prop.neg, NA, NA, NA, NA)
	                existing.dirs = list.dirs(path = upper.eval.folder, full.names = FALSE, recursive = FALSE)
	                # folder names between network numbers and prop.neg. Ex: c("All","k=1-VI","k=2-VI")
	                #               in "/home/nejat/eclipse/workspace-neon/CC-Opti-Heur-Analysis/out/evaluate-partitions/n=20_l0=4_dens=1.0000/propMispl=0.1000/propNeg=0.7895/k=1"
	                for(folder.at.k.level in existing.dirs){
	                    
	                    # retreive network numbers
	                    files = list.files(file.path(upper.eval.folder,folder.at.k.level))
	                    indx = which(startsWith(files,"network"))
	                    net.folders = as.integer(gsub("network=","",files[indx])) # a subset of 'in.rand.net.folders' if 'folder.at.k.level' starts with "k=". Otherwise, the whole list
	                    take.average(n, l0, d, NA, detected.imb.interval, prop.neg, folder.at.k.level, net.folders,
	                                 cor.clu.exact.algo, comp.measures, force)
	                }
	                
	            }
            }
        }
	        
	    
	    
		
	}
	
}