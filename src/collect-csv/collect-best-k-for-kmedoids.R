# TODO: Add comment
# 
# Author: nejat
###############################################################################






#################################################################
# It constructs a single row of the csv file in output. It retrieves the statistics 
#   related to the number of clusters obtained in the kmedoids method for the given values of the input parameters. 
#
# n: graph order
# l0: number of cluster
# d: density
# prop.mispl: proportion of misplaced links
# prop.neg: proportion of negative links
# network.no: network id (the identifiers start from 1
# cor.clu.exact.algo: the name of correlation clustering algorithm to run
# measure: comparison/evaluation measure used to obtain a distance matrix (before applying the kmedoids method).
# force: whether or not the existing files are overwritten by a fresh call of all corresponding methods (e.g partitioning method)
#
##################################################################
retreive.best.k.for.kmedoids = function(n, l0, d, prop.mispl, prop.neg, network.no,
		cor.clu.exact.algo, measure, force)
{
	
	curr.data =c()
		
	tlog(20, "collecting => algo.name: ", cor.clu.exact.algo)
	
	for(graph.desc.name in c(SIGNED.UNWEIGHTED.FILE)){ # OPPOSITE.SIGNED.UNWEIGHTED.FILE
		tlog(24, "plot transition networks => graph.desc.name: ", graph.desc.name)
		part.folder = part.folder = get.part.folder.path(n, l0, d, prop.mispl, prop.neg, network.no, cor.clu.exact.algo, graph.desc.name)
		eval.folder = get.eval.folder.path(n, l0, d, prop.mispl, NA, prop.neg, k=ALL, network.no, cor.clu.exact.algo, graph.desc.name)

		mems = load.membership.files(part.folder)
		nb.sol = length(mems)
		#print(nb.sol)		
		if(nb.sol>1){
			edit.dist.nb.comp.table.file = file.path(eval.folder, paste0(EVAL.BEST.K.FOR.KMEDOIDS,"-",measure,".csv"))
			df = read.csv(edit.dist.nb.comp.table.file, row.names = 1, header= TRUE, check.names=FALSE, stringsAsFactors=FALSE)
			curr.line = matrix(df[,c(BEST.K.FOR.SILH.COL.NAME)])
			rownames(curr.line) = paste0("n=", n, "l0=", l0, "d=", d, "prop.mispl=", prop.mispl, "prop.neg=", prop.neg, "network.no=", network.no)
			curr.data = rbind(curr.data, curr.line)
		} else {
		    curr.line = matrix(NA)
		    rownames(curr.line) = paste0("n=", n, "l0=", l0, "d=", d, "prop.mispl=", prop.mispl, "prop.neg=", prop.neg, "network.no=", network.no)
		    curr.data = rbind(curr.data, curr.line)
		}
	}
	
	
	return(curr.data)
}









#################################################################
#
# It is the starting method in the aim of collecting the statistics related to the number of clusters obtained in the kmedoids method. 
#   It handles all networks by graph.sizes, prop.mispls, my.prop.negs and in.rand.net.folders
#   For the given value of the parameter l0, it generates a csv file.
#
# graph.sizes: a vector of values regarding graph orders to be considered
# d: density (it is a single value)
# l0: number of clusters to be considered (it is a single value)
# prop.mispls: a vector of values regarding proportion of misplaced links
# prop.negs: a vector of values regarding proportion of negative links (for now, it is not operational)
# in.rand.net.folders: a vector of values regarding input random graph folders. Sequantial integers (1, .., 10)
# cor.clu.exact.algo: the name of correlation clustering algorithm to run
# comp.measures: comparison/evaluation measures used to obtain a distance matrix (before applying the kmedoids method).
# force: whether or not the existing files are overwritten by a fresh call of all corresponding methods (e.g partitioning method)
#
##################################################################
collect.all.best.k.for.kmedoids = function(graph.sizes, d, l0, prop.mispls, prop.negs, in.rand.net.folders,
		cor.clu.exact.algo, comp.measures, force)
{
    
    for(measure in comp.measures){
        tlog(20, "collecting => measure: ", measure)
        
    	all.data = c()
    	
    	tlog("starts collecting")
    	for(n in graph.sizes){
    		tlog(8, "partitioning networks => n: ", n)
    		
    		for(prop.mispl in prop.mispls){
    			tlog(8, "collecting => prop.mispl: ", prop.mispl)
    			
    		    if(is.na(prop.negs) && d == 1){
    		        prop.negs = compute.prop.neg(n, d, l0, prop.mispl)
    		    }	
    			
    			for(prop.neg in prop.negs){
    				tlog(12, "collecting => prop.neg: ", prop.neg)
    				
    				for(network.no in in.rand.net.folders){
    					tlog(16, "collecting => network.no: ", network.no)
    					
    				    curr.data = retreive.best.k.for.kmedoids(n, l0, d, prop.mispl, prop.neg, network.no,
    							cor.clu.exact.algo, measure, force)
    					all.data = rbind(all.data, curr.data)
    				}
    				
    			}
    			
    		}
    		
    	}

    	if(length(all.data)>0){
    		colnames(all.data) = c(BEST.K.FOR.SILH.COL.NAME)
    		#rownames(all.data) = NULL
    		if(!dir.exists(OUTPUT.CSV.FOLDER))
    		    dir.create(OUTPUT.CSV.FOLDER, recursive=FALSE, showWarnings=FALSE)
    		write.csv(file=file.path(OUTPUT.CSV.FOLDER,paste0("cluster-analysis-",measure,"-l0=",l0,".csv")), x=all.data)
    	}
    }
}
