# TODO: Add comment
# 
# Author: nejat
###############################################################################


###################################################################
# It builds the adjacency matrix from edit distance matrix for a given edit distance value.
#  For instance, if nb.edit = 2, the adjacency matrix will have the values of 1 for each pair of partitions that are close each other with 2-Edits.
#
# returns adjacency matrix for a given edit distance value.
###################################################################
build.adj.matrix.from.edit.dist.matrix = function(edit.dist.mtrx, nb.edit){
    nb.sol = nrow(edit.dist.mtrx)
    adj.mtrx = as.numeric(as.matrix(edit.dist.mtrx)) # as.numeric() converts the matrix into vector
    adj.mtrx[which(adj.mtrx != nb.edit)] = 0
    adj.mtrx = matrix(adj.mtrx, nb.sol,nb.sol)
    # since the matrix contains only nb edit dist = nb.edit, we need to make it 1s matrix
    adj.mtrx = adj.mtrx/nb.edit
    return(adj.mtrx)
}


#################################################################
# This function requires the existence of a Edit-distance matrix, that is,
#   a matrix A where A_ij indicates the number of edits to obtain the partition j from partition i.
# So, this function first creates the adjacency matrix for a given edit-distance value, then creates the corresponding graph.
#   Finally, it obtains the statistics of the number and size of the connected components, and writes them into file .
#   The following Edit-distance values are used:
#   - 1-Edit, 2-Edit, 3-Edit and 4-Edit distance matrices
#   - 1&2-Edit, 1&2&3-Edit and 1&2&3&4-Edit distance matrices 
#
#
# eval.folder: the folder where the Edit-distance matrix is stored
# part.folder: the folder containing all the partitions. This is used to know the number of partitions
# force: Boolean. if TRUE, the statistics about connected components will be re-calculated
#
##################################################################
process.edit.dist.connected.comps = function(eval.folder, part.folder, force){
	
	mbrshps = load.membership.files(part.folder)
	m = length(mbrshps) # nb partition
#	network.path = file.path(net.folder,graph.name)
#	g = read.graph.ils(network.path)
#	n = vcount(g) # nb node
	
	tlog(24, "process edit dist connected comps")
	
	edit.dist.comps.table.file = file.path(eval.folder, paste0(EVAL.EDIT.DIST.COMPS,".csv"))
	edit.dist.nb.comp.table.file = file.path(eval.folder, paste0(EVAL.EDIT.DIST.NB.COMP,".csv"))
	unlink(edit.dist.comps.table.file)
	unlink(edit.dist.nb.comp.table.file)
	
	edit.mtrx.file = file.path(eval.folder, paste0(EVAL.DIST.MATRIX.FILE.PREFIX,"-",EDIT,".csv"))
	
	if(force || !(file.exists(edit.dist.comps.table.file) && file.exists(edit.dist.comps.table.file))){
	    # if at least 2 solutions
    	if(m>1 && file.exists(edit.mtrx.file)){
    		edit.dist.matrix = as.matrix(read.csv(edit.mtrx.file, row.names = 1, header= TRUE, check.names=FALSE))
    					
    		res.comps = c()
    		res.nb.comps = c()
    		adj.matrix1 = build.adj.matrix.from.edit.dist.matrix(edit.dist.matrix, nb.edit=1)
    		g1 = graph_from_adjacency_matrix(adj.matrix1, mode="undirected", weighted=NULL, diag=FALSE)
    		comp.mem1 = components(g1)$membership
    		res.comps = cbind(res.comps, paste0(comp.mem1, collapse=","))
    		res.nb.comps = cbind(res.nb.comps, length(unique(comp.mem1)))
    		
    		adj.matrix2 = build.adj.matrix.from.edit.dist.matrix(edit.dist.matrix, nb.edit=2)	
    		adj.matrix3 = build.adj.matrix.from.edit.dist.matrix(edit.dist.matrix, nb.edit=3)
    		adj.matrix4 = build.adj.matrix.from.edit.dist.matrix(edit.dist.matrix, nb.edit=4)
    		
    		adj.matrix12 = adj.matrix1 + adj.matrix2
    		g12 = graph_from_adjacency_matrix(adj.matrix12, mode="undirected", weighted=NULL, diag=FALSE)
    		comp.mem12 = components(g12)$membership
    		res.comps = cbind(res.comps, paste0(comp.mem12, collapse=","))
    		res.nb.comps = cbind(res.nb.comps, length(unique(comp.mem12)))
    		
    		adj.matrix123 = adj.matrix1 + adj.matrix2 + adj.matrix3
    		g123 = graph_from_adjacency_matrix(adj.matrix123, mode="undirected", weighted=NULL, diag=FALSE)
    		comp.mem123 = components(g123)$membership
    		res.comps = cbind(res.comps, paste0(comp.mem123, collapse=","))
    		res.nb.comps = cbind(res.nb.comps, length(unique(comp.mem123)))
    		
    		adj.matrix1234 = adj.matrix1 + adj.matrix2 + adj.matrix3 + adj.matrix4
    		g1234 = graph_from_adjacency_matrix(adj.matrix1234, mode="undirected", weighted=NULL, diag=FALSE)
    		comp.mem1234 = components(g1234)$membership
    		res.comps = cbind(res.comps, paste0(comp.mem1234, collapse=","))
    		res.nb.comps = cbind(res.nb.comps, length(unique(comp.mem1234)))
    		
    		
    		colnames(res.comps) = c(EDIT1.COMP.MEM.COL.NAME, EDIT12.COMP.MEM.COL.NAME, EDIT123.COMP.MEM.COL.NAME,
    			EDIT1234.COMP.MEM.COL.NAME)
    		rownames(res.comps) = "Component membership"
    		colnames(res.nb.comps) = c(EDIT1.COMP.MEM.COL.NAME, EDIT12.COMP.MEM.COL.NAME, EDIT123.COMP.MEM.COL.NAME,
    			EDIT1234.COMP.MEM.COL.NAME)	
    		rownames(res.nb.comps) = "Nb component"
    		
    		write.csv(x=res.comps, file=edit.dist.comps.table.file, col.names=TRUE, row.names=TRUE)
    		write.csv(x=res.nb.comps, file=edit.dist.nb.comp.table.file, col.names=TRUE, row.names=TRUE)
    	}
    }

}



#################################################################
# The wrapper function to create the statistics about the connected components of the graphs created by Edit-distance matrix.
# It prepares such statistics based on the considered algorithm name and graph type (weighted or not, etc.).
# 
# n: graph size
# k: number of cluster
# d: density
# prop.mispl: proportion of misplaced links
# prop.neg: proportion of negative links
# network.no: network id (the identifiers start from 1)
# cor.clu.exact.algo: the name of correlation clustering algorithm to run
# force: whether or not the existing files are overwritten
#
##################################################################
create.edit.dist.connected.comps = function(n, l0, d, prop.mispl, prop.neg, network.no,
		cor.clu.exact.algo, force)
{
	
	# only exact solutions
	tlog(16, "start to evaluate exact solutions")
	for(graph.desc.name in c(SIGNED.UNWEIGHTED.FILE)){ # OPPOSITE.SIGNED.UNWEIGHTED.FILE
		tlog(24, "evaluating partitions => graph.desc.name: ", graph.desc.name)
			
		#for(e.algo.name in cor.clu.exact.algos){
			tlog(20, "evaluating partitions => algo.name: ", cor.clu.exact.algo)
			
			e.part.folder = get.part.folder.path(n, l0, d, prop.mispl, prop.neg, network.no, cor.clu.exact.algo, graph.desc.name)
			e.eval.folder = get.eval.folder.path(n, l0, d, prop.mispl, NA, prop.neg, k=ALL, network.no, cor.clu.exact.algo, graph.desc.name)
			if(!dir.exists(e.eval.folder))
				dir.create(path=e.eval.folder, showWarnings=FALSE, recursive=TRUE)
			
			tlog(24, "proceed only exact solutions in folder: ", e.eval.folder)
			process.edit.dist.connected.comps(e.eval.folder, e.part.folder, force)
		#}
	}
	
}




#################################################################
# The wrapper function to create the statistics about the connected components of the graphs created by Edit-distance matrix.
#   It handles all networks by graph.sizes, prop.mispls, my.prop.negs and in.rand.net.folders
#
# graph.sizes: a vector of values regarding graph sizes to be considered
# d: density (it is a single value)
# l0: number of clusters to be considered (it is a single value)
# prop.mispls: a vector of values regarding proportion of misplaced links
# prop.negs: a vector of values regarding proportion of negative links (for now, it is not operational)
# cor.clu.exact.algo: the name of correlation clustering algorithm to run
# in.rand.g.folders: input random graph folders. Sequantial integers (1, .., 10)
# force: whether or not the existing files are overwritten.
#
##################################################################
create.all.edit.dist.connected.comps = function(graph.sizes, d, l0, prop.mispls, prop.negs, in.rand.net.folders,
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
						
						create.edit.dist.connected.comps(n, l0, d, prop.mispl, prop.neg, network.no,
								cor.clu.exact.algo, force)
						
					}
				}
				
			}
			
		}
		
	}
	
}

