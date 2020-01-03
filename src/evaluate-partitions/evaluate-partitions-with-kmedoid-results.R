# TODO: Add comment
# 
# Author: nejat
###############################################################################


#################################################################
# It determines the best number of kmedoids clustering analysis. There are two cases:
#   - only two partitions:
#       We do not have any calculated silhouette scores in this case due to the definition of silhouette.
#       So, we compute a dissimilarity score between these two partitions. 
#       If the score is less than 0.5, we conclude a single cluster, otherwise 2 single cluster.
#   - more than two partitions: 
#       in this case, at least one silhouette score is guaranteed. We retreive the max silhouette score amon the existing values.
#       Then, we check if the max vcalue is greater than 'SILH.THRESH'. If so, we can be confident about the silhouette score.
#       Otherwise, there is a single cluster. Note that in this case, the first and last values (i.e. k=1 or k=n) can not be computed
#       due to the definition of silhouete. User may want to get those values for k=3 or k=4. Unforetunately, we omit this point for now.
# 
# 
# mems: a vector containing all partitions
# nb.sol: the number of partitions
# avg.sil.scores: a vector of silhouette scores for each cluster value (e.g. 2, 3, 4, etc..)
# clu.analysis.measure.folder: the folder related to cluster analysis (i.e. kmedoids) based on a known dissimilarity measure
# measure:
#
#################################################################
process.best.k.in.kmedoids = function(mems, nb.sol, avg.sil.scores, eval.folder, clu.analysis.measure.folder, measure){
    
    # ----------------------------------
    # silhouette
    best.k.silh = NA
    if(nb.sol == 2){ # i.e. nb.sol=2
        d = compare.partition.pair.as.distance(mems[[1]], mems[[2]], measures=NVI) # we needed a normalized measure, NVI is arbitrary
        print("dis")
        print(d)
        if(d < 0.5) # ===> this 0.5 value is conceptually different than SILH.THRESH. Here, we do not evaluate any silhouette value
            best.k.silh = 1
        else
            best.k.silh = 2
    }
    else{
        best.k.silh = which.max(avg.sil.scores)
        best.silh.score = avg.sil.scores[best.k.silh]
        if(best.silh.score < SILH.THRESH)
            best.k.silh = 1
    }
    

    k.desc = paste0("k=", best.k.silh)
    best.silh = sprintf("%.4f",as.numeric(avg.sil.scores[best.k.silh]))
    desc = paste0("silh","=", best.silh)
    folder.name = paste(k.desc, desc, sep="_")
    analysis.silh.folder = file.path(clu.analysis.measure.folder,folder.name)
    best.silh.mbrshp = read.table(file=file.path(analysis.silh.folder,paste0(MBRSHP.FILE.PREFIX,".txt")))$V1
    # ----------------------------------
    
    table.file = file.path(eval.folder, paste0(EVAL.BEST.K.FOR.KMEDOIDS,"-",measure,".csv"))
    res = matrix(c(best.k.silh, best.silh, paste0(best.silh.mbrshp,collapse=",")), nrow=1)
    rownames(res) = "kmedoids"
    colnames(res) = c(BEST.K.FOR.SILH.COL.NAME, BEST.SILH.COL.NAME, BEST.MEM.FOR.SILH.COL.NAME)
    write.csv(x=res, file=table.file, row.names=TRUE)
    
}


#################################################################
# It is the wrapper function to determine the best number of cluster number in kmedoids results.
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
evaluate.best.k.in.kmedoids = function(n, l0, d, prop.mispl, prop.neg, network.no,
		cor.clu.exact.algo, comp.measures, force)
{
	
	net.folder = get.input.network.folder.path(n, l0, d, prop.mispl, prop.neg, network.no)
	
	# only exact solutions
	tlog(16, "start to evaluate exact solutions")
	for(graph.desc.name in c(SIGNED.UNWEIGHTED.FILE)){ # OPPOSITE.SIGNED.UNWEIGHTED.FILE
		tlog(24, "evaluating partitions => graph.desc.name: ", graph.desc.name)
			
		tlog(20, "evaluating partitions => algo.name: ", cor.clu.exact.algo)
		
#			graph.name =paste0(graph.desc.name,".G")
		r = part.folder = get.part.folder.path(n, l0, d, prop.mispl, prop.neg, network.no, cor.clu.exact.algo, graph.desc.name)
		plot.folder = get.transition.graph.plot.folder.path(n, l0, d, prop.mispl, prop.neg, network.no, cor.clu.exact.algo, graph.desc.name)
		clu.analysis.folder = get.clu.analysis.folder.path(n, l0, d, prop.mispl, prop.neg, network.no, cor.clu.exact.algo, graph.desc.name)
		eval.folder = get.eval.folder.path(n, l0, d, prop.mispl, detected.imb=NA, prop.neg, k=ALL, network.no, cor.clu.exact.algo, graph.desc.name)
		
		mems = load.membership.files(part.folder)
		nb.sol = length(mems)
		print(nb.sol)		
		
		if(nb.sol > 1){
			for(measure in comp.measures){
				tlog(20, "performing cluster characterization => measure: ", measure)
		
				clu.analysis.measure.folder = file.path(clu.analysis.folder,measure)
				silhouette.scores.file = file.path(clu.analysis.measure.folder, paste0(SILH.SCORE.FILENAME,".csv"))
				clu.result.silh = read.csv(file=silhouette.scores.file, header=TRUE)
				avg.sil.scores = clu.result.silh[["avg.sil.score"]]
	

				process.best.k.in.kmedoids(mems, nb.sol, avg.sil.scores, eval.folder, clu.analysis.measure.folder, measure)
					
			}
			
		}
					

	}
	
}





#################################################################
# It is the wrapper function to determine the best number of cluster number in kmedoids results.
#  It iterates over graph sizes, proportion of misplaced links, proportion of negative links and generated different network sizes.
#
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
evaluate.partitions.with.kmedoid.results = function(graph.sizes, d, l0, prop.mispls, prop.negs, in.rand.net.folders,
		 cor.clu.exact.algo, comp.measures, force)
{
	tlog("starts evaluating partitions")
	for(n in graph.sizes){
		tlog(4, "evaluating partitions => n: ", n)
		
		for(prop.mispl in prop.mispls){
			tlog(8, "evaluating partitions => prop.mispl: ", prop.mispl)
			
		    if(is.na(prop.negs) && d == 1){
		        prop.negs = compute.prop.neg(n, d, k, prop.mispl)
		    }		
			
			for(prop.neg in prop.negs){
				tlog(12, "evaluating partitions => prop.neg: ", prop.neg)
				
								
				net.folder = get.input.network.folder.path(n, l0, d, prop.mispl, prop.neg, network.no=NA)
				if(dir.exists(net.folder)){
				
					for(network.no in in.rand.net.folders){
						tlog(16, "evaluating partitions => network.no: ", network.no)
						
						evaluate.best.k.in.kmedoids(n, l0, d, prop.mispl, prop.neg, network.no,
								cor.clu.exact.algo, comp.measures, force)
								
					}
				}
				
			}
			
		}
		
	}
	
}


