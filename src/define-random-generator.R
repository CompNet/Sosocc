
#################################################################
# It generates a membership vector associated with a partition of (almost) equal-sized modules. 
# This generation respects the given values of the parameters n and l0.
#################################################################
generate.param.membership = function(n, l0){
	nl0 = n%/%l0
	clu.sizes = rep(nl0, l0)
	nb.remaining = n%%l0
	remaining = rep(0, l0)
	if(nb.remaining > 0)
		remaining[1:nb.remaining] = 1
	clu.sizes = clu.sizes + remaining
	membership <- rep(1:l0,clu.sizes)
	return(membership)
}



#################################################################
# It computes the proportion of negative links given four parameters: n, d, l0 and prop.mispl
# This calculation is specific to our random signed graph generator. 
#################################################################
compute.prop.neg = function(n, d, l0, prop.mispl){

	prop.mispl = as.numeric(prop.mispl)
	n = as.numeric(n)
	l0 = as.numeric(l0)
	
	# membership <- rep(1:l0,each=n%/%l0)
	membership <- generate.param.membership(n, l0)
	pext <- sum(apply(t(combn(x=max(membership),m=2,simplify=TRUE)), 1, function(r)
					{	n1 <- length(which(membership==r[1]))
						n2 <- length(which(membership==r[2]))
						n1 * n2
					})) / (n*(n-1)/2)
	
			
	if(d == 1){
		prop.neg = pext
	} else {
		prop.neg = ((1-prop.mispl) - (1-pext))/(1-2*prop.mispl) # proportion of negative links
	}
	
	return(prop.neg)
}
