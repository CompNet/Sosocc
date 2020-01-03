

#################################################################
# Based on individual graph parmeters it creates a new list covering all. 
#   This eases the use and transfter of parameters from one function to another.
# 
# n: graph size
# l0: the number of initial cluster numbers when generating first a graph.
#       It is a beforehand value. The real cluster number can be different than 'l0'.
# d: density
# prop.mispl: propotion of misplaced links. It is beforehand value. The detected value of imbalance can be different.
# detected.imb.interval: detected value of imbalance. It is used when evaluation results are organized by this value (and not prop.mispl)
# prop.neg: proportion of negative links
# k: number of solution classes in k-medoids result
# network.no
# algo.name: ExCC or another algo name
# graph.desc.name: graph type. Ex: weighted, unweighted, etc.
#
##################################################################
build.g.params = function(n, l0, d, prop.mispl, detected.imb.interval, prop.neg, k, network.no, algo.name, graph.desc.name){
	
	g.params = list(n=n, l0=l0, d=d, prop.mispl=prop.mispl, detected.imb.interval=detected.imb.interval, prop.neg=prop.neg, k=k, network.no=network.no,
					algo.name=algo.name, graph.desc.name=graph.desc.name)
	return(g.params)
}


#################################################################
# It assigns a value to a specific key of the list called 'g.params'.
#   The key we search for is found thanks to the variable 'by'.
#
# g.params: a list of graph parameters.
# by: a stringenabling to find the right key in the list. Ex: GRAPH.SIZE
# value: the value which will be used for assignemnt. Ex: 24
#
##################################################################
update.g.params.by = function(g.params, by, value){
	
	if(!is.na(by)){
		if(by == GRAPH.SIZE){
			g.params$n = value
		}
		else if(by == PROP.MISPL){
			g.params$prop.mispl = value
		}
	    else if(by == DETECTED.IMB.PROP.PARAM.NAME){
	        g.params$detected.imb.interval = value
	    }
		else if(by == PROP.NEG){
			g.params$prop.neg = value
		}
		else if(by == NETWORK.NO){
			g.params$network.no = value
		}
		else if(by == ALGO.NAME){
			g.params$algo.name = value
		}
		else if(by == GRAPH.DESC.NAME){
			g.params$graph.desc.name = value
		}

	}
	
	return(g.params)
}



#################################################################
# It creates plot title based on graph parameters.
#
# g.params: a list of graph parameters.
#
##################################################################
create.plot.title = function(g.params){
    # we did not consider g.params$k ==> number of solution classes (useless here)
    # g.params$l0, neither
	plot.title = prepare.plot.title(g.params$n, g.params$prop.mispl, g.params$detected.imb.interval, g.params$prop.neg, g.params$network.no,
			g.params$algo.name, g.params$graph.desc.name)			
	return(plot.title)
}



