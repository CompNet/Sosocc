# File-related constants and functions.
# 
# Author: Vincent Labatut 07/2017
###############################################################################
library("igraph")





#################################################################
#
# workaround: since when density=1, prop.neg=NA.
# So, we need to handle this exceptional case here since prop.neg is changed in terms of 'prop.mispl' (and 'n')
#################################################################
compute.prop.neg.if.not.available = function(g.params, by, param.value){
	
	prop.neg = g.params$prop.neg
	if(is.na(prop.neg)){
		if(by == GRAPH.SIZE)
			prop.neg = compute.prop.neg(param.value, g.params$d, g.params$l0, g.params$prop.mispl) # param.value for: n
		else if(by == PROP.MISPL || by == REAL.IMB.PERC)
			prop.neg = compute.prop.neg(g.params$n, g.params$d, g.params$l0, param.value) # param.value for: prop.mispl
		else
			prop.neg = compute.prop.neg(g.params$n, g.params$d, g.params$l0, g.params$prop.mispl)
	}
	return(prop.neg)
}



#################################################################
# by: one of these values GRAPH.SIZE, PROP.MISPL, PROP.NEG, NETWORK.NO, ALGO.NAME, GRAPH.DESC.NAME
#
##################################################################
get.part.folder.path.by = function(by, param.value, g.params)
{	
	
	prop.neg = compute.prop.neg.if.not.available(g.params, by, param.value)
		
	n = g.params$n
	l0 = g.params$l0
	d = g.params$d
	prop.mispl = g.params$prop.mispl
#	prop.neg = g.params$prop.neg
	# k = g.params$k
	network.no = g.params$network.no
	algo.name = g.params$algo.name
	graph.desc.name = g.params$graph.desc.name
	
	
	folder = NA
	if(by == GRAPH.SIZE)
		folder = get.part.folder.path(param.value, l0, d, prop.mispl, prop.neg, network.no, algo.name, graph.desc.name)
	else if(by == PROP.MISPL || by == REAL.IMB.PERC)
		folder = get.part.folder.path(n, l0, d, param.value, prop.neg, network.no, algo.name, graph.desc.name)
	else if(by == PROP.NEG)
		folder = get.part.folder.path(n, l0, d, prop.mispl, param.value, network.no, algo.name, graph.desc.name)
	else if(by == NETWORK.NO)
		folder = get.part.folder.path(n, l0, d, prop.mispl, prop.neg, param.value, algo.name, graph.desc.name)
	else if(by == ALGO.NAME)
		folder = get.part.folder.path(n, l0, d, prop.mispl, prop.neg, network.no, param.value, graph.desc.name)
	else if(by == GRAPH.DESC.NAME)
		folder = get.part.folder.path(n, l0, d, prop.mispl, prop.neg, network.no, algo.name, param.value)
	
	return(folder)
}


#################################################################
# by: one of these values GRAPH.SIZE, PROP.MISPL, PROP.NEG, NETWORK.NO, ALGO.NAME, GRAPH.DESC.NAME
#
##################################################################
get.eval.folder.path.by = function(by, param.value, g.params)
{	
	prop.neg = compute.prop.neg.if.not.available(g.params, by, param.value)
	
	n = g.params$n
	l0 = g.params$l0
	d = g.params$d
	prop.mispl = g.params$prop.mispl
#	prop.neg = g.params$prop.neg
	k = g.params$k
	network.no = g.params$network.no
	algo.name = g.params$algo.name
	graph.desc.name = g.params$graph.desc.name
	
	
	folder = NA
	if(by == GRAPH.SIZE)
		folder = get.eval.folder.path(param.value, l0, d, prop.mispl, prop.neg, k, network.no, algo.name, graph.desc.name)
	else if(by == PROP.MISPL || by == REAL.IMB.PERC)
		folder = get.eval.folder.path(n, l0, d, param.value, prop.neg, k, network.no, algo.name, graph.desc.name)
	else if(by == PROP.NEG)
		folder = get.eval.folder.path(n, l0, d, prop.mispl, param.value, k, network.no, algo.name, graph.desc.name)
	else if(by == NETWORK.NO)
		folder = get.eval.folder.path(n, l0, d, prop.mispl, prop.neg, k, param.value, algo.name, graph.desc.name)
	else if(by == ALGO.NAME)
		folder = get.eval.folder.path(n, l0, d, prop.mispl, prop.neg, k, network.no, param.value, graph.desc.name)
	else if(by == GRAPH.DESC.NAME)
		folder = get.eval.folder.path(n, l0, d, prop.mispl, prop.neg, k, network.no, algo.name, param.value)
	
	return(folder)
}


#################################################################
# by: 
#
##################################################################
get.plot.folder.path.by = function(g.params, by, subplot.by=NA){

	prop.neg = g.params$prop.neg

	# -------------------------------------
	if(is.na(prop.neg)){
		# since prop.neg = NA, we know that prop.neg is actually a single value (i.e. not a range) for a given (graph.size, prop.mispl) pair
		# but, we need to know both graph.size and prop.mispl to compute prop.neg
		if( (!is.na(by) && (by == GRAPH.SIZE || by == PROP.MISPL || by == DETECTED.IMB.PROP.PARAM.NAME))
				|| (!is.na(subplot.by) && (subplot.by == GRAPH.SIZE || subplot.by == PROP.MISPL || subplot.by == DETECTED.IMB.PROP.PARAM.NAME)) )
			prop.neg = NOT.AVAILABLE.FOLDER.NAME
		else
			prop.neg = sprintf("%.4f", compute.prop.neg(g.params$n, g.params$d, g.params$l0, g.params$prop.mispl))
	}
	# -------------------------------------
	
	plot.folder = get.plot.folder.path(g.params$n, g.params$l0, g.params$d, g.params$prop.mispl, g.params$detected.imb.interval, paste0("propNeg=",prop.neg),
	                                   g.params$k, g.params$network.no, g.params$algo.name, g.params$graph.desc.name)
	return(plot.folder)
}




# ============================================================================================
# ============================================================================================



###############################################################################
# Builds the path of a subfolder corresponding to one specific parameter set.
#
# n: number of nodes in the graph.
# l0: number of clusters in the graph.
# d: density of the graph.
# prop.mispls: proportion of misplaced links.
# prop.neg: proportion of negative links in the network.
# network.no
#
# returns: the folder path defined to store the network and the associated result files.
###############################################################################
get.input.network.folder.path <- function(n, l0, d, prop.mispl=NA, prop.neg=NA, network.no=NA)
{	result <- file.path(RANDOM.NETWORKS.FOLDER,paste0(paste0("n=",n),paste0("_l0=",l0),paste0("_dens=",sprintf("%.4f",d))))
	if(!is.na(prop.mispl))
		result <- file.path(result, paste0("propMispl=",sprintf("%.4f",prop.mispl)))
	if(!is.na(prop.neg))
		result <- file.path(result, paste0("propNeg=",sprintf("%.4f",prop.neg)))
	
	if(!is.na(network.no))
		result <- file.path(result, paste0("network=", network.no))
		
	return(result)
}




###############################################################################
# Builds the path of a subfolder corresponding to one specific parameter set.
#
# n: number of nodes in the graph.
# l0: number of clusters in the graph.
# d: density of the graph.
# prop.mispls: proportion of misplaced links.
# prop.neg: proportion of negative links in the network.
# network.no
# algo.name
# graph.desc.name
#
# returns: the folder path defined to store the network and the associated result files.
###############################################################################
get.part.folder.path <- function(n, l0, d, prop.mispl=NA, prop.neg=NA, network.no=NA, algo.name=NA, graph.desc.name=NA)
{	
	result <- file.path(PARTITIONS.FOLDER,paste0(paste0("n=",n),paste0("_l0=",l0),paste0("_dens=",sprintf("%.4f",d))))
	if(!is.na(prop.mispl))
		result <- file.path(result, paste0("propMispl=",sprintf("%.4f",prop.mispl)))
	if(!is.na(prop.neg))
		result <- file.path(result, paste0("propNeg=",sprintf("%.4f",prop.neg)))
	
	if(!is.na(network.no))
		result <- file.path(result, paste0("network=", network.no))		
	if(!is.na(algo.name))
		result <- file.path(result, algo.name)
	if(!is.na(graph.desc.name))
		result <- file.path(result, graph.desc.name)
	
	return(result)
}



###############################################################################
# Builds the path of a subfolder corresponding to one specific parameter set.
#
# n: number of nodes in the graph.
# l0: number of clusters in the graph.
# d: density of the graph.
# prop.mispls: proportion of misplaced links.
# prop.neg: proportion of negative links in the network.
# k: number of classes found in the clustering process of optimal partitions (e.g. k-medoids). 
# network.no
# algo.name
# graph.desc.name
#
# returns: the folder path defined to store the network and the associated result files.
###############################################################################
get.eval.folder.path <- function(n, l0, d, prop.mispl=NA, detected.imb=NA, prop.neg=NA, k=NA, network.no=NA, algo.name=NA, graph.desc.name=NA)
{	
	result <- file.path(EVALUATION.FOLDER,paste0(paste0("n=",n),paste0("_l0=",l0),paste0("_dens=",sprintf("%.4f",d))))
	if(!is.na(prop.mispl))
		result <- file.path(result, paste0("propMispl=",sprintf("%.4f",prop.mispl)))
	if(!is.na(detected.imb))
	    result <- file.path(result, paste0(DETECTED.IMB.PROP.PARAM.NAME,"=",detected.imb)) # ex: [0.10,0.15)]
	if(!is.na(prop.neg))
		result <- file.path(result, paste0("propNeg=",sprintf("%.4f",prop.neg)))
	

	if(is.character(k) && !is.na(k))
	    result <- file.path(result, k) # paste0("k=", k)
	# else if(is.numeric(k) && !is.na(k))
	#     result <- file.path(result, paste0("k=", k))
	
	
	# 'network.no' might be a string (like 'Summary') when we process all those values for plots
	if(is.character(network.no) && !is.na(network.no))
		result <- file.path(result, network.no)
	else if(is.numeric(network.no) && !is.na(network.no))
		result <- file.path(result, paste0("network=", network.no))
		
	if(!is.na(algo.name))
		result <- file.path(result, algo.name)
	if(!is.na(graph.desc.name))
		result <- file.path(result, graph.desc.name)
		
	return(result)
}



###############################################################################
# Builds the path of a subfolder corresponding to one specific parameter set.
#
# n: number of nodes in the graph.
# l0: number of clusters in the graph.
# d: density of the graph.
# prop.mispls: proportion of misplaced links.
# prop.neg: proportion of negative links in the network. The difference between this method and other 'get' methods. It can be 'NotAvailable' 
# k: number of classes found in the clustering process of optimal partitions (e.g. k-medoids). The defautl value is 1, i.E. all networks are together.
# network.no
# algo.name
# graph.desc.name
#
# returns: the folder path defined to store the network and the associated result files.
###############################################################################
get.plot.folder.path <- function(n, l0, d, prop.mispl=NA, detected.imb=NA, prop.neg=NA, k=NA, network.no=NA, algo.name=NA, graph.desc.name=NA)
{		
	# 'n' might be a string (like 'All') when we process all those values for plots
	result <- file.path(PLOT.FOLDER,paste0(paste0("n=",n),paste0("_l0=",l0),paste0("_dens=",sprintf("%.4f",d))))
	

	# 'prop.mispl' might be a string (like 'All') when we process all those values for plots
	if(!is.na(prop.mispl) && is.character(prop.mispl))
		result <- file.path(result, paste0("propMispl=",prop.mispl))
	else if(!is.na(prop.mispl) && is.numeric(prop.mispl))
		result <- file.path(result, paste0("propMispl=",sprintf("%.4f",prop.mispl)))
	
	if(!is.na(detected.imb))
	    result <- file.path(result, paste0(DETECTED.IMB.PROP.PARAM.NAME,"=",detected.imb)) # ex: [0.10,0.15)]
	
	if(!is.na(prop.neg) && is.character(prop.neg))
		result <- file.path(result, prop.neg)
	else if(!is.na(prop.neg) && is.numeric(prop.neg))
		result <- file.path(result, paste0("propNeg=",sprintf("%.4f",prop.neg)))
	
	if(is.character(k) && !is.na(k))
	    result <- file.path(result, k) # paste0("k=", k)
	# else if(is.numeric(k) && !is.na(k))
	#     result <- file.path(result, paste0("k=", k))
		
	# 'network.no' might be a string (like 'Summary') when we process all those values for plots
	if(!is.na(network.no) && is.character(network.no))
		result <- file.path(result, paste0("network=",network.no))
	else if(!is.na(network.no) && is.numeric(network.no))
		result <- file.path(result, paste0("network=", network.no))
		
	if(!is.na(algo.name))
		result <- file.path(result, algo.name)
	if(!is.na(graph.desc.name))
		result <- file.path(result, graph.desc.name)
		
	return(result)
}


###############################################################################
#
#
# n: number of nodes in the graph.
# l0: number of clusters in the graph.
# d: density of the graph.
# prop.mispls: proportion of misplaced links.
# prop.neg: proportion of negative links in the network. The difference between this method and other 'get' methods. It can be 'NotAvailable' 
# k: number of classes found in the clustering process of optimal partitions (e.g. k-medoids). The defautl value is 1, i.E. all networks are together.
# network.no
# algo.name
# graph.desc.name
#
# returns:
###############################################################################
get.transition.graph.plot.folder.path <- function(n, l0, d, prop.mispl=NA, prop.neg=NA, k=NA, network.no=NA, algo.name=NA, graph.desc.name=NA)
{		
	# 'n' might be a string (like 'All') when we process all those values for plots
	result <- file.path(TRANSITION.GRAPH.PLOT.FOLDER,paste0(paste0("n=",n),paste0("_l0=",l0),paste0("_dens=",sprintf("%.4f",d))))
	

	# 'prop.mispl' might be a string (like 'All') when we process all those values for plots
	if(!is.na(prop.mispl) && is.character(prop.mispl))
		result <- file.path(result, prop.mispl)
	else if(!is.na(prop.mispl) && is.numeric(prop.mispl))
		result <- file.path(result, paste0("propMispl=",sprintf("%.4f",prop.mispl)))
	
	if(!is.na(prop.neg) && is.character(prop.neg))
		result <- file.path(result, prop.neg)
	else if(!is.na(prop.neg) && is.numeric(prop.neg))
		result <- file.path(result, paste0("propNeg=",sprintf("%.4f",prop.neg)))
	
	if(!is.na(k))
	    result <- file.path(result, paste0("k=", k))
	else
	    result <- file.path(result, paste0("k=", 1)) # default value => one solution class
		
	# 'network.no' might be a string (like 'Summary') when we process all those values for plots
	if(!is.na(network.no) && is.character(network.no))
		result <- file.path(result, network.no)
	else if(!is.na(network.no) && is.numeric(network.no))
		result <- file.path(result, paste0("network=", network.no))
		
	if(!is.na(algo.name))
		result <- file.path(result, algo.name)
	if(!is.na(graph.desc.name))
		result <- file.path(result, graph.desc.name)
		
	return(result)
}




###############################################################################
# Builds the path of a subfolder corresponding to one specific parameter set.
#
# n: number of nodes in the graph.
# l0: number of clusters in the graph.
# d: density of the graph.
# prop.mispls: proportion of misplaced links.
# prop.neg: proportion of negative links in the network.
# network.no
# algo.name
# graph.desc.name
#
# returns: the folder path defined to store the network and the associated result files.
###############################################################################
get.clu.analysis.folder.path <- function(n, l0, d, prop.mispl=NA, prop.neg=NA, network.no=NA, algo.name=NA, graph.desc.name=NA)
{	
	result <- file.path(CLU.ANALYSIS.FOLDER,paste0(paste0("n=",n),paste0("_l0=",l0),paste0("_dens=",sprintf("%.4f",d))))
	if(!is.na(prop.mispl))
		result <- file.path(result, paste0("propMispl=",sprintf("%.4f",prop.mispl)))
	if(!is.na(prop.neg))
		result <- file.path(result, paste0("propNeg=",sprintf("%.4f",prop.neg)))

	# 'network.no' might be a string (like 'Summary') when we process all those values for plots
	if(is.character(network.no) && !is.na(network.no))
		result <- file.path(result, network.no)
	else if(is.numeric(network.no) && !is.na(network.no))
		result <- file.path(result, paste0("network=", network.no))
		
	if(!is.na(algo.name))
		result <- file.path(result, algo.name)
	if(!is.na(graph.desc.name))
		result <- file.path(result, graph.desc.name)
		
	return(result)
}



###############################################################################
# Builds the path of a subfolder corresponding to one specific parameter set.
#
# n: number of nodes in the graph.
# l0: number of clusters in the graph.
# d: density of the graph.
# prop.mispls: proportion of misplaced links.
# prop.neg: proportion of negative links in the network.
# network.no
# algo.name
# graph.desc.name
#
# returns: the folder path defined to store the network and the associated result files.
###############################################################################
get.clu.characterization.folder.path <- function(n, l0, d, prop.mispl=NA, prop.neg=NA, network.no=NA, algo.name=NA, graph.desc.name=NA)
{	
	result <- file.path(CLUSTER.CHARACTERIZATION.FOLDER,paste0(paste0("n=",n),paste0("_l0=",l0),paste0("_dens=",sprintf("%.4f",d))))
	if(!is.na(prop.mispl))
		result <- file.path(result, paste0("propMispl=",sprintf("%.4f",prop.mispl)))
	if(!is.na(prop.neg))
		result <- file.path(result, paste0("propNeg=",sprintf("%.4f",prop.neg)))
	
	
	# 'network.no' might be a string (like 'Summary') when we process all those values for plots
	if(is.character(network.no) && !is.na(network.no))
		result <- file.path(result, network.no)
	else if(is.numeric(network.no) && !is.na(network.no))
		result <- file.path(result, paste0("network=", network.no))
		
	if(!is.na(algo.name))
		result <- file.path(result, algo.name)
	if(!is.na(graph.desc.name))
		result <- file.path(result, graph.desc.name)
		
	return(result)
}



###############################################################################
# Builds the path of a subfolder corresponding to one specific parameter set.
#
# n: number of nodes in the graph.
# l0: number of clusters in the graph.
# d: density of the graph.
# prop.mispls: proportion of misplaced links.
# prop.neg: proportion of negative links in the network.
# k: number of classes found in the clustering process of optimal partitions (e.g. k-medoids). The defautl value is 1, i.E. all networks are together.
# network.no
# algo.name
# graph.desc.name
#
# returns: the folder path defined to store the network and the associated result files.
###############################################################################
get.clu.analysis.difference.folder.path <- function(n, l0, d, prop.mispl=NA, prop.neg=NA, network.no=NA, algo.name=NA, graph.desc.name=NA)
{	
	result <- file.path(CLU.ANALYSIS.DIFFERENCE.FOLDER,paste0(paste0("n=",n),paste0("_l0=",l0),paste0("_dens=",sprintf("%.4f",d))))
	if(!is.na(prop.mispl))
		result <- file.path(result, paste0("propMispl=",sprintf("%.4f",prop.mispl)))
	if(!is.na(prop.neg))
		result <- file.path(result, paste0("propNeg=",sprintf("%.4f",prop.neg)))
	
	
	# 'network.no' might be a string (like 'Summary') when we process all those values for plots
	if(is.character(network.no) && !is.na(network.no))
		result <- file.path(result, network.no)
	else if(is.numeric(network.no) && !is.na(network.no))
		result <- file.path(result, paste0("network=", network.no))
		
	if(!is.na(algo.name))
		result <- file.path(result, algo.name)
	if(!is.na(graph.desc.name))
		result <- file.path(result, graph.desc.name)
		
	return(result)
}



###############################################################################
# Builds the path of a subfolder corresponding to one specific parameter set.
#
# n: number of nodes in the graph.
# l0: number of clusters in the graph.
# d: density of the graph.
# prop.mispls: proportion of misplaced links.
# prop.neg: proportion of negative links in the network.
# k: number of classes found in the clustering process of optimal partitions (e.g. k-medoids). The defautl value is 1, i.E. all networks are together.
# network.no
# algo.name
# graph.desc.name
#
# returns: the folder path defined to store the network and the associated result files.
###############################################################################
get.networks.stats.folder.path <- function(n, l0, d, prop.mispl=NA, prop.neg=NA, k=NA, network.no=NA, algo.name=NA, graph.desc.name=NA)
{	
	result <- file.path(NETWORKS.STAT.FOLDER,paste0(paste0("n=",n),paste0("_l0=",l0),paste0("_dens=",sprintf("%.4f",d))))
	if(!is.na(prop.mispl))
		result <- file.path(result, paste0("propMispl=",sprintf("%.4f",prop.mispl)))
	if(!is.na(prop.neg))
		result <- file.path(result, paste0("propNeg=",sprintf("%.4f",prop.neg)))
	
	if(!is.na(k))
	    result <- file.path(result, paste0("k=", k))
	else
	    result <- file.path(result, paste0("k=", 1)) # default value => one solution class
	
	
	# 'network.no' might be a string (like 'Summary') when we process all those values for plots
	if(is.character(network.no) && !is.na(network.no))
		result <- file.path(result, network.no)
	else if(is.numeric(network.no) && !is.na(network.no))
		result <- file.path(result, paste0("network=", network.no))
	
	if(!is.na(algo.name))
		result <- file.path(result, algo.name)
	if(!is.na(graph.desc.name))
		result <- file.path(result, graph.desc.name)
	
	return(result)
}
