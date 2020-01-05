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
retreive.nb.comp.in.transition.graphs = function(n, l0, d, prop.mispl, prop.neg, network.no,
		cor.clu.exact.algo, comp.measures, force)
{
	
	conn.comp.data =c()
		
	tlog(16, "start to plot transition with exact algorithms")	
	#for(algo.name in cor.clu.exact.algo){
		tlog(20, "plot transition networks => algo.name: ", cor.clu.exact.algo)
		
		for(graph.desc.name in c(SIGNED.UNWEIGHTED.FILE)){ # OPPOSITE.SIGNED.UNWEIGHTED.FILE
			tlog(24, "plot transition networks => graph.desc.name: ", graph.desc.name)
			part.folder = part.folder = get.part.folder.path(n, l0, d, prop.mispl, prop.neg, network.no, cor.clu.exact.algo, graph.desc.name)
			eval.folder = get.eval.folder.path(n, l0, d, prop.mispl, NA, prop.neg, k=ALL, network.no, cor.clu.exact.algo, graph.desc.name)
			
			mems = load.membership.files(part.folder)
			nb.sol = length(mems)
			print(nb.sol)		
			if(nb.sol>1){

			
				edit.dist.nb.comp.table.file = file.path(eval.folder, paste0(EVAL.EDIT.DIST.NB.COMP,".csv"))
				df = read.csv(edit.dist.nb.comp.table.file, row.names = 1, header= TRUE, check.names=FALSE, stringsAsFactors=FALSE)
				curr.line = df[,c(EDIT1.COMP.MEM.COL.NAME, EDIT1234.COMP.MEM.COL.NAME)]
				rownames(curr.line) = paste0("n=", n, "l0=", l0, "d=", d, "prop.mispl=", prop.mispl, "prop.neg=", prop.neg, "network.no=", network.no)
				conn.comp.data = rbind(conn.comp.data, curr.line)

			}
		}
	#}
	
	return(conn.comp.data)
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
collect.nb.comp.in.all.transition.graphs = function(graph.sizes, d, l0, prop.mispls, prop.negs, in.rand.net.folders,
		cor.clu.exact.algo, comp.measures, force)
{
	all.conn.comp.data = c()
	
	tlog("starts partitioning networks")
	for(n in graph.sizes){
		tlog(8, "partitioning networks => n: ", n)
		
		for(prop.mispl in prop.mispls){
			tlog(8, "partitioning networks => prop.mispl: ", prop.mispl)
			
		    if(is.na(prop.negs) && d == 1){
		        prop.negs = compute.prop.neg(n, d, l0, prop.mispl)
		    }
			
			for(prop.neg in prop.negs){
				tlog(12, "partitioning networks => prop.neg: ", prop.neg)
				
				for(network.no in in.rand.net.folders){
					tlog(16, "partitioning networks => network.no: ", network.no)
					
					conn.comp.data = retreive.nb.comp.in.transition.graphs(n, l0, d, prop.mispl, prop.neg, network.no,
							cor.clu.exact.algo, comp.measures, force)
					all.conn.comp.data = rbind(all.conn.comp.data, conn.comp.data)
				}
				
			}
			
		}
		
	}
	
	if(length(all.conn.comp.data)>0){
		colnames(all.conn.comp.data) = c("nbComp-Edit1", "nbComp-Edit1&Edit2&Edit3&Edit4")
		#rownames(all.conn.comp.data) = NULL
		if(!dir.exists(OUTPUT.CSV.FOLDER))
		    dir.create(OUTPUT.CSV.FOLDER, recursive=FALSE, showWarnings=FALSE)
		
		write.csv(file=file.path(OUTPUT.CSV.FOLDER,paste0("transition-graph-conn-compenents-stats-l0=",l0,".csv")), x=all.conn.comp.data)
	}
}
