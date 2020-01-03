# TODO: Add comment
# 
# Author: nejat
###############################################################################








#################################################################
#
# d: density
# k: nb cluster
#
##################################################################
retreive.nb.opt.solution = function(n, l0, d, prop.mispl, prop.neg, network.no,
		cor.clu.exact.algo, comp.measures, force, plot.formats)
{
	
	data =c()
		
	#tlog(16, "start to plot transition with exact algorithms")	
	tlog(20, "plot transition networks => algo.name: ", cor.clu.exact.algo)
	
	for(graph.desc.name in c(SIGNED.UNWEIGHTED.FILE)){ # OPPOSITE.SIGNED.UNWEIGHTED.FILE
		tlog(24, "plot transition networks => graph.desc.name: ", graph.desc.name)
		part.folder = part.folder = get.part.folder.path(n, l0, d, prop.mispl, prop.neg, network.no, cor.clu.exact.algo, graph.desc.name)
		#eval.folder = get.eval.folder.path(n, l0, d, prop.mispl, prop.neg, network.no, cor.clu.exact.algo, graph.desc.name)

		curr.line = matrix(NA,1,1)
		if(dir.exists(part.folder)){
			mems = load.membership.files(part.folder)
			curr.line = matrix(length(mems),1,1)
		}
		rownames(curr.line) = paste0("n=", n, "l0=", l0, "d=", d, "prop.mispl=", prop.mispl, "prop.neg=", prop.neg, "network.no=", network.no)
		data = rbind(data, curr.line)
	}
	
	
	return(data)
}









#################################################################
#
# graph.sizes
# d: density
# k: nb cluster
# prop.mispls
# prop.negs
# cor.clu.heur.algos
# heur.reps: heuristic repetitions. Sequantial integers (1, .., 100)
# cor.clu.exact.algos
# in.rand.g.folders: input random graph folders. Sequantial integers (1, .., 10)
# force
# plot.formats
#
##################################################################
collect.all.nb.opt.solution = function(graph.sizes, d, l0, prop.mispls, prop.negs, in.rand.net.folders,
		cor.clu.exact.algo, comp.measures, force, plot.formats)
{
	all.data = c()
	
	tlog("starts collecting networks")
	for(n in graph.sizes){
		tlog(8, "collecting networks => n: ", n)
		
		for(prop.mispl in prop.mispls){
			tlog(8, "collecting networks => prop.mispl: ", prop.mispl)
			
		    if(is.na(prop.negs) && d == 1){
		        prop.negs = compute.prop.neg(n, d, l0, prop.mispl)
		    }
			
			for(prop.neg in prop.negs){
				tlog(12, "collecting networks => prop.neg: ", prop.neg)
				
				for(network.no in in.rand.net.folders){
					tlog(16, "partitioning networks => network.no: ", network.no)
					
					data = retreive.nb.opt.solution(n, l0, d, prop.mispl, prop.neg, network.no,
							cor.clu.exact.algo, comp.measures, force, plot.formats)
					if(!is.na(data))
						all.data = rbind(all.data, data)
				}
			}
			
		}
		
	}
	
	if(length(all.data)>0){
		colnames(all.data) = c(NB.SOL.COL.NAME)
		#rownames(all.conn.comp.data) = NULL
		if(!dir.exists(OUTPUT.CSV.FOLDER))
		    dir.create(OUTPUT.CSV.FOLDER, recursive=FALSE, showWarnings=FALSE)
		
		write.csv(file=file.path(OUTPUT.CSV.FOLDER,paste0("nb-sol-l0=",l0,".csv")), x=all.data)
	}
}
