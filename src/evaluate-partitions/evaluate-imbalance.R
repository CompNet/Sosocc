



############################################################################
# It calculates the imbalance in terms of count or percentage for Correlation Clustering problem.
#
# g: graph
# membership: a membership, i.e. partition vector
# output.type: either "count" or "percentage" or "node.imbalance"
#
# returns imbalance
############################################################################
compute.imbalance.from.membership = function(g, membership, output.type = "count"){

	membership = as.integer(membership)
	edge.mat <- get.edgelist(g)
	edge.mat <- matrix(as.integer(edge.mat), nrow(edge.mat), ncol(edge.mat))
	if(edge.mat[1,1] == 0) # check if node ids start from 0 or 1. This affects directly 'clus.mat' variable
		edge.mat = edge.mat + 1
	
	clus.mat <- cbind(membership[edge.mat[,1]], membership[edge.mat[,2]])
		
	neg.links <- E(g)$weight<0
	pos.links <- E(g)$weight>=0

	misplaced <- (clus.mat[,1]==clus.mat[,2] & neg.links) | (clus.mat[,1]!=clus.mat[,2] & pos.links)
	imb.val = sum(abs(E(g)$weight[misplaced]))
	
	# --------------------------------
	lpos.imbalance <- E(g)$weight * as.numeric(misplaced & pos.links)
	lneg.imbalance <- abs(E(g)$weight) * as.numeric(misplaced & neg.links)
	npos.imbalance <- sapply(1:vcount(g), function(u) 
			{	idx <- which(edge.mat[,1]==u | edge.mat[,2]==u)
				result <- sum(lpos.imbalance[idx])
				return(result)
			})
	nneg.imbalance <- sapply(1:vcount(g), function(u) 
			{	idx <- which(edge.mat[,1]==u | edge.mat[,2]==u)
				result <- sum(lneg.imbalance[idx])
				return(result)
			})
	
	max.val = max(c(npos.imbalance,nneg.imbalance))
	if(max.val != 0){ # if the situation has some imbalance
		npos.imbalance <- npos.imbalance / max.val # normalized
		nneg.imbalance <- nneg.imbalance / max.val # normalized
	}
	# --------------------------------
	
	# make them explicit
	n.in.clu.imb = nneg.imbalance # negative misplaced links are the misplaced link insde clusters
	n.betw.clu.imb = npos.imbalance # pisitive misplaced links are the misplaced link between clusters
	
	# ========================
	# ========================
	
	if(output.type == "count")
		return(format(round(imb.val, 3), nsmall = 3)) # 3 decimal floating
	else if(output.type == "percentage"){
		perc = (imb.val/ sum(abs(E(g)$weight)))*100
		return(format(round(perc, 3), nsmall = 3))
	} else if(output.type == "node.imbalance") # normalized
		return( list(in.imb=n.in.clu.imb, betw.imb=n.betw.clu.imb) )
	else
		return(NA)
	
}



############################################################################
# It calculates the imbalance in terms of count or percentage for the relaxed version of Correlation Clustering problem.
#
# g: graph
# membership: a membership, i.e. partition vector
# output.type: either "count" or "percentage" or "node.imbalance"
#
# returns imbalance
############################################################################
compute.relaxed.imbalance.from.membership = function(g, membership, output.type = "count"){
	
	edge.mat <- get.edgelist(g)
	edge.mat <- matrix(as.integer(edge.mat), nrow(edge.mat), ncol(edge.mat))
	if(edge.mat[1,1] == 0) # check if node ids start from 0 or 1. This affects directly 'clus.mat' variable
		edge.mat = edge.mat + 1
	
	clus.mat <- cbind(membership[edge.mat[,1]], membership[edge.mat[,2]])
	
	#compare link signs and positions 
	neg.links <- E(g)$weight<0
	pos.links <- E(g)$weight>=0
	
	nb.clu = length(unique(membership))
	
	# compute the imbalance (i.e. cost) of intra-edges
	imb.val=0
	n.in.clu.imb = rep(0, vcount(g))
	for(clu in seq_len(nb.clu)){
		in.edges = (clus.mat[, 1] == clu & clus.mat[, 2] == clu)
		
		pos.misplaced = pos.links & in.edges
		neg.misplaced = neg.links & in.edges
		pos.cost = sum(E(g)$weight[pos.misplaced])
		neg.cost = sum(abs(E(g)$weight[neg.misplaced]))
		
		
		l.imbalance=NA
		if(neg.cost > pos.cost){ # if positive links are dominant, take the negative ones
			l.imbalance <- E(g)$weight * as.numeric(pos.misplaced)
			imb.val = imb.val + pos.cost
		}
		else{ # if neg.cost=pos.cost, the neg.cost will be chosen
			l.imbalance <- abs(E(g)$weight) * as.numeric(neg.misplaced)
			imb.val = imb.val + neg.cost
		}
		
		
		# node imbalance
		n.imbalance <- sapply(1:vcount(g), function(u) 
				{	idx <- which(edge.mat[,1]==u | edge.mat[,2]==u)
					result <- sum(l.imbalance[idx])
					return(result)
				})
		n.in.clu.imb = n.in.clu.imb + n.imbalance
	}
	
	
	# --------------------------------------------------------------
	# --------------------------------------------------------------
	
	
	# compute the imbalance (i.e. cost) of inter-edges if nb.clu > 1
	n.betw.clu.imb = rep(0, vcount(g))
	if(nb.clu > 1){
		pair.list = combn(x=nb.clu, m=2) # x'in m'li combinasyonu
		nb.pair = length(pair.list)/2
		
		for(i in seq_len(nb.pair)){
			clu1 = pair.list[1, i]
			clu2 = pair.list[2, i]
			
			betw.edges = (clus.mat[, 1] == clu1 & clus.mat[, 2] == clu2) | (clus.mat[, 1] == clu2 & clus.mat[, 2] == clu1)
			
			pos.misplaced = pos.links & betw.edges
			neg.misplaced = neg.links & betw.edges
			pos.cost = sum(E(g)$weight[pos.misplaced])
			neg.cost = sum(abs(E(g)$weight[neg.misplaced]))
			
			
			l.imbalance=NA
			if(pos.cost > neg.cost){ # if positive links are dominant, take the negative ones
				l.imbalance <- abs(E(g)$weight) * as.numeric(neg.misplaced)
				imb.val = imb.val + neg.cost
			}
			else{  # if neg.cost=pos.cost, the pos.cost will be chosen
				l.imbalance <- E(g)$weight * as.numeric(pos.misplaced)
				imb.val = imb.val + pos.cost
			}
			
			# node imbalance
			n.imbalance <- sapply(1:vcount(g), function(u) 
					{	idx <- which(edge.mat[,1]==u | edge.mat[,2]==u)
						result <- sum(l.imbalance[idx])
						return(result)
					})
			n.betw.clu.imb = n.betw.clu.imb + n.imbalance
		}
	}
	
	max.val = max(c(n.in.clu.imb,n.betw.clu.imb))
	if(max.val != 0){ # if the situation has some imbalance
		n.in.clu.imb <- n.in.clu.imb / max.val # normalized
		n.betw.clu.imb <- n.betw.clu.imb / max.val # normalized
	}
	
	
	# ========================
	# ========================
	
	if(output.type == "count")
		return(format(round(imb.val, 3), nsmall = 3)) # 3 decimal floating
	else if(output.type == "percentage"){
		perc = (imb.val/ sum(abs(E(g)$weight)))*100
		return(format(round(perc, 3), nsmall = 3))
	} else if(output.type == "node.imbalance") # normalized
		return( list(in.imb=n.in.clu.imb, betw.imb=n.betw.clu.imb) )
	else
		return(NA)
}


