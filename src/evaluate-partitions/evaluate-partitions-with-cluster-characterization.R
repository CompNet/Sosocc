# TODO: Add comment
# 
# Author: nejat
###############################################################################



#################################################################
# It prepares the statistics about the core part size of the detected clusters of kmedoids
#       for each cluster of kmedoids and for each value of core part thresholds 'CORE.PART.THRESHOLDS'. 
#   This method requires the presence of the file containing the membership information about core part for each cluster of kmedoids. 
#   In other words, the core part information in each cluster of kmedoids is already known.
#
# clu.analysis.measure.folder: the folder related to cluster analysis (i.e. kmedoids) based on a known dissimilarity measure
# clu.characterization.measure.folder: the folder related to cluster characterization analysis based on a known dissimilarity measure (the folder of core part and representatives)
# eval.algo.folder: the folder containing evaluation process
# measure
# n: the number of nodes in the considered graph
#
#################################################################
process.core.part.size.in.solution.classes = function(clu.analysis.measure.folder, clu.characterization.measure.folder, eval.algo.folder, measure, n){
    
    table.file = file.path(eval.algo.folder, paste0(EVAL.BEST.K.FOR.KMEDOIDS,"-",measure,".csv"))
    if(file.exists(table.file)){
        df = read.csv(table.file, check.names=FALSE, stringsAsFactors=FALSE)
        k = as.numeric(df[,BEST.K.FOR.SILH.COL.NAME])
        best.silh = df[,BEST.SILH.COL.NAME]
        silh.folder.name = paste0("k=",k,"_silh=", sprintf("%.4f",best.silh))
        clu.analysis.mbrshp = as.numeric(unlist(strsplit(df[,BEST.MEM.FOR.SILH.COL.NAME], ",")))
        analysis.silh.folder = file.path(clu.analysis.measure.folder,silh.folder.name)
        characterization.silh.folder = file.path(clu.characterization.measure.folder,silh.folder.name)
        
        for(core.part.threshold in CORE.PART.THRESHOLDS){
            characterization.thresh.folder = file.path(characterization.silh.folder,
                                                       paste0("core-part=",sprintf("%.2f",core.part.threshold)))
            print(characterization.thresh.folder)
            vals = c()
            for(clu.no in 1:k){
                fpath = file.path(characterization.thresh.folder,paste0(MBRSHP.FILE.PREFIX,"-",CLUSTER.ID.FOLDER.PREFIX,clu.no,".txt"))
                print(fpath)
                if(file.exists(fpath)){
                    mem = read.table(file=fpath)$V1
                    core.part.size = n - length(which(is.na(mem)))
                    vals = c(vals,core.part.size)
                }
            }
            
            print(vals)
            norm.vals = vals/n
            class.core.part.size.prop.file = file.path(eval.algo.folder,
                                                       paste0(EVAL.CLASS.CORE.PART.SIZE.FILE.NAME,"-tresh=",sprintf("%.2f",core.part.threshold),".csv"))
            mtrx = matrix(sprintf("%.4f",norm.vals), ncol=1, nrow=length(vals)) # normalized value
            colnames(mtrx) = PROP.CLASS.CORE.PART.SIZE.COL.NAME
            rownames(mtrx) = paste0("k=",1:k)
            print(class.core.part.size.prop.file)
            write.csv(x=mtrx, file=class.core.part.size.prop.file, row.names=TRUE)
        }
        
    }
}




#################################################################
# It is the wrapper function to prepare the statistics about the core part size of each cluster of kmedoids results.
#   It iterates over exact algorithms used to obtain partitions, graph types and dissimilarity measures.
#
#
# n: graph size
# l0: number of cluster
# d: density
# prop.mispl: proportion of misplaced links
# prop.neg: proportion of negative links
# network.no: network id (the identifiers start from 1)
# cor.clu.exact.algo: the name of correlation clustering algorithm to run
# comp.measures: dissimilarity measure to compare the optimal solutions
# force: whether or not the existing files are overwritten by a fresh call of all corresponding methods (e.g partitioning method)
#
##################################################################
evaluate.core.part.size.in.solution.classes = function(n, l0, d, prop.mispl, prop.neg, network.no,
		cor.clu.exact.algo, comp.measures, force)
{
			
	# only exact solutions
	tlog(16, "start to evaluate exact solutions")
	for(graph.desc.name in c(SIGNED.UNWEIGHTED.FILE)){ # OPPOSITE.SIGNED.UNWEIGHTED.FILE
		tlog(24, "evaluating partitions => graph.desc.name: ", graph.desc.name)
			
		# for(algo.name in cor.clu.exact.algos){
			tlog(20, "evaluating partitions => algo.name: ", cor.clu.exact.algo)
			
			net.folder = get.input.network.folder.path(n, l0, d, prop.mispl, prop.neg, network.no)
			part.folder = get.part.folder.path(n, l0, d, prop.mispl, prop.neg, network.no, cor.clu.exact.algo, graph.desc.name)
			eval.algo.folder = get.eval.folder.path(n, l0, d, prop.mispl, detected.imb=NA, prop.neg, k=ALL, network.no, cor.clu.exact.algo, graph.desc.name)
			clu.analysis.folder = get.clu.analysis.folder.path(n, l0, d, prop.mispl, prop.neg, network.no, cor.clu.exact.algo, graph.desc.name)
			clu.characterization.folder = get.clu.characterization.folder.path(n, l0, d, prop.mispl, prop.neg, network.no, cor.clu.exact.algo, graph.desc.name)
			
			if(dir.exists(net.folder)){
				
#				partitions = load.membership.files(part.folder)
#				m = length(partitions) # nb partition
#				n = length(partitions[[1]]) # node size
				
				for(measure in comp.measures){
    				tlog(20, "performing cluster characterization => measure: ", measure)
    				clu.analysis.measure.folder = file.path(clu.analysis.folder,measure)
    				clu.characterization.measure.folder = file.path(clu.characterization.folder,measure)
    				
    #				matrix.file = file.path(eval.algo.folder,paste0(EVAL.DIST.MATRIX.FILE.PREFIX,"-",measure,".csv"))
    #				dist.mtrx = as.matrix(read.csv(matrix.file, row.names = 1, header= TRUE, check.names=FALSE)) # square matrix
    #				dist.mtrx = matrix(as.numeric(sprintf("%.4f", dist.mtrx)), m, m)
    
    				process.core.part.size.in.solution.classes(clu.analysis.measure.folder, clu.characterization.measure.folder, eval.algo.folder, measure, n)

				
				}
				
	
			}
			
		# }
			
	}
	
}






#################################################################
# It is the wrapper function to prepare the statistics about the core part size of each cluster of kmedoids results.
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
evaluate.partitions.with.cluster.characterization = function(graph.sizes, d, l0, prop.mispls, prop.negs, in.rand.net.folders,
		 cor.clu.exact.algo, comp.measures, force)
{
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
						
						evaluate.core.part.size.in.solution.classes(n, l0, d, prop.mispl, prop.neg, network.no,
								cor.clu.exact.algo, comp.measures, force)
							
						
					}
				}
				
			}
			
		}
		
	}
	
}


