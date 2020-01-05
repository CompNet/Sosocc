

##################################################################
##
##
##################################################################
#compute.imbalance.from.membership = function(g, membership, output.type = "count"){	
#	membership = as.integer(membership)
#	edge.mat <- get.edgelist(g)
#	clus.mat <- cbind(membership[as.integer(edge.mat[,1])+1], membership[as.integer(edge.mat[,2])+1])  # +1 since node ids start from 0
#	
#	same.clus <- clus.mat[,1]==clus.mat[,2]
#	
#	#compare link signs and positions 
#	neg.links <- E(g)$weight<0
#	pos.links <- E(g)$weight>0
#	neg.misplaced <- same.clus & neg.links
#	pos.misplaced <- !same.clus & pos.links
#	all.misplaced <- neg.misplaced | pos.misplaced
#	
#	imb.value = sum(abs(E(g)$weight[all.misplaced]))
#	
#	
#	if(output.type == "count"){
#		return(format(round(imb.value, 5), nsmall = 5)) # 5 decimal floating
#	}
#	else{ # if(output.type == "percentage")
#		perc = (imb.value/ sum(abs(E(g)$weight)))*100
#		return(format(round(perc, 5), nsmall = 5))
#	}
#}



#################################################################
#
# d: density
# k: nb cluster
#
##################################################################
create.stats.for.instance = function(n, k, d, prop.mispl, prop.neg, network.no, cor.clu.exact.algos, comp.measures, force)
{
	# NETWORKS.STAT.FOLDER

	
	net.folder = get.input.network.folder.path(n, k, d, prop.mispl, prop.neg, network.no)
	tlog(16, "start to perform networks stats with exact algorithms")	
	for(algo.name in cor.clu.exact.algos){
		tlog(20, "performing networks stats => algo.name: ", algo.name)
		
		for(graph.desc.name in c(SIGNED.UNWEIGHTED.FILE)){ # OPPOSITE.SIGNED.UNWEIGHTED.FILE
			tlog(24, "performing networks stats => graph.desc.name: ", graph.desc.name)
			
			part.folder = get.part.folder.path(n, k, d, prop.mispl, prop.neg, network.no, algo.name, graph.desc.name)
			if(!dir.exists(part.folder))
				dir.create(path=part.folder, showWarnings=FALSE, recursive=TRUE)
			
			tlog(24, "performing networks stats in ", part.folder)		
			g.name = paste0(graph.desc.name,".G")
			network.path = file.path(net.folder,g.name)
			g = read.graph.ils(network.path)
			# BE CAREFUL THAT node ids start from 0
			
			k=K # this is the paramter used in main
			init.mem = rep(c(1,k),each=n/2) # membership info when graph generation is performed
			init.in.clu.pos.degrees = rep((n/k)-1, times=n) # each cluster there are n/k nodes, so (n/k - 1) edges for each node in the cluster
			init.in.clu.neg.degrees = rep(0, times=n)
			init.betw.clu.pos.degrees = rep(0, times=n) # each cluster there are n/k nodes, so (n/k - 1) edges for each node in the cluster
			init.betw.clu.neg.degrees = rep((n/k), times=n)
			g.pos <- delete.edges(graph=g,edges=which(E(g)$weight<0))
			g.pos.degrees = degree(g.pos)
			
			final.in.clu.pos.degrees = c()
			final.in.clu.neg.degrees = c()
			final.betw.clu.pos.degrees = c()
			final.betw.clu.neg.degrees = c()
			for(clu.no in 1:k){
				indx = which(init.mem == clu.no)
				clu.g = induced.subgraph(graph=g,vids=indx)
				clu.pos.g <- delete.edges(graph=clu.g,edges=which(E(clu.g)$weight<0))
				in.clu.pos.degrees = degree(clu.pos.g)
				in.clu.neg.degrees = init.in.clu.pos.degrees[indx] - in.clu.pos.degrees
				betw.clu.pos.degrees = g.pos.degrees[indx] - in.clu.pos.degrees
				betw.clu.neg.degrees = init.betw.clu.neg.degrees[indx] - betw.clu.pos.degrees
				
				final.in.clu.pos.degrees = c(final.in.clu.pos.degrees,in.clu.pos.degrees)
				final.in.clu.neg.degrees = c(final.in.clu.neg.degrees,in.clu.neg.degrees)
				final.betw.clu.pos.degrees = c(final.betw.clu.pos.degrees,betw.clu.pos.degrees)
				final.betw.clu.neg.degrees = c(final.betw.clu.neg.degrees,betw.clu.neg.degrees)
			}
			
			print("final.in.clu.pos.degrees")
			print(final.in.clu.pos.degrees)
			print("final.in.clu.neg.degrees")
			print(final.in.clu.neg.degrees)
			print("final.betw.clu.pos.degrees")
			print(final.betw.clu.pos.degrees)
			print("final.betw.clu.neg.degrees")
			print(final.betw.clu.neg.degrees)
			
			print("------------------")
			
			print("init imb")
			print(compute.imbalance.from.membership(g, init.mem, output.type = "count"))
			print("actual imb")
			partitions = load.membership.files(part.folder)
			any.solution = partitions[[1]]
			print(any.solution)
			print(compute.imbalance.from.membership(g, any.solution, output.type = "count"))
		}
	}
	
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
#
##################################################################
perform.networks.stats = function(graph.sizes, d, k, prop.mispls, prop.negs, in.rand.net.folders, cor.clu.exact.algos, comp.measures, force)
{
	tlog("starts performing networks stats")
	for(n in graph.sizes){
		tlog(8, "performing networks stats => n: ", n)
		
		for(prop.mispl in prop.mispls){
			tlog(8, "performing networks stats => prop.mispl: ", prop.mispl)
			
			my.prop.negs = NA
			if(is.na(prop.negs)){
				prop.neg = compute.prop.neg(n, d, k, prop.mispl)
				my.prop.negs = prop.neg # we show explicitely that we get only 1 value since density=1.0
			}
			
			for(prop.neg in my.prop.negs){
				tlog(12, "performing networks stats => prop.neg: ", prop.neg)
				
				for(network.no in in.rand.net.folders){
					tlog(16, "performing networks stats => network.no: ", network.no)
					
					create.stats.for.instance(n, k, d, prop.mispl, prop.neg, network.no, cor.clu.exact.algos, comp.measures, force)
				}
				
			}
			
		}
		
	}
	
}