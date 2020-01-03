
library(igraph)





#################################################################
#
# d: density
# k: nb cluster
#
##################################################################
add.layouts = function(n, k, d, prop.mispl, prop.neg, network.no, unsigned.layout.methods, signed.layout.methods)
{
	net.folder = get.input.network.folder.path(n, k, d, prop.mispl, prop.neg, network.no)
	
	graph.name = paste0(SIGNED.UNWEIGHTED.FILE, ".graphml") # our reference graph to get optimal and worst imbalance count value
	layout.graph.name = paste0(GRAPH.FILE.LAYOUT.PREFIX,"-",graph.name)
	network.path = file.path(net.folder, graph.name)
	g = suppressWarnings(read.graph(file=network.path, format="graphml"))
	gpos <- delete.edges(graph=g,edges=which(E(g)$weight<0))
	
	tlog(20, "adding unsigned graph based layouts for all networks")
	for(layout.method in unsigned.layout.methods){
		tlog(24,"layout method: ", layout.method)
		lyt = NA
		if(layout.method == FRUCHTERMAN.REINGOLD){
			lyt <- layout.fruchterman.reingold(graph=gpos)	
		}
		if(layout.method == KAMADA.KAWAI){
			lyt <- layout.kamada.kawai(graph=gpos)
		}
		
		g <- set.vertex.attribute(graph=g, name=paste0("x_",layout.method),value=lyt[,1])
		g <- set.vertex.attribute(graph=g, name=paste0("y_",layout.method),value=lyt[,2])
	}
	
	tlog(20, "adding signed graph based layouts for all networks")
	for(layout.method in signed.layout.methods){
		tlog(24,"layout method: ", layout.method)
		lyt = layout.signed.laplacian(g, method=layout.method)
		g <- set.vertex.attribute(graph=g, name=paste0("x_",layout.method),value=lyt[,1])
		g <- set.vertex.attribute(graph=g, name=paste0("y_",layout.method),value=lyt[,2])
	}
	
	write.graph(graph=g,file=file.path(net.folder,layout.graph.name),format="graphml")
		
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
add.layouts.for.all.networks = function(graph.sizes, d, k, prop.mispls, prop.negs, in.rand.net.folders,
		 unsigned.layout.methods, signed.layout.methods)
{
	tlog("starts add layouts for all networks")
	for(n in graph.sizes){
		tlog(4, "adding layouts for all networks => n: ", n)
		
		for(prop.mispl in prop.mispls){
			tlog(8, "adding layouts for all networks => prop.mispl: ", prop.mispl)
			
			if(is.na(prop.negs) && d == 1){
			    prop.negs = compute.prop.neg(n, d, k, prop.mispl)
			}
			
			for(prop.neg in prop.negs){
				tlog(12, "adding layouts for all networks => prop.neg: ", prop.neg)
				
				net.folder = get.input.network.folder.path(n, k, d, prop.mispl, prop.neg, network.no=NA)
				if(dir.exists(net.folder)){
	
					for(network.no in in.rand.net.folders){
						tlog(16, "adding layouts for all networks => network.no: ", network.no)
						
						add.layouts(n, k, d, prop.mispl, prop.neg, network.no, unsigned.layout.methods, signed.layout.methods)
						
					}
				
				}
				
			}
			
		}
		
	}
	
}
