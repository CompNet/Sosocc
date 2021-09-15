

########################################################################
# It converts a weighted signed graph into an unweighted signed graph.
########################################################################
convert.weight.into.unweight.input.graph = 
		function(network.path, unweighted.network.path)
{
	network = file(network.path, 'r')  # connection for reading 
	first.line = readLines(network, n = 1) 
	t <- read.table(network, header = FALSE, skip=1) # skip first line
	close(network) 
	
	weights = t$V3
	converted.weights = sapply(weights, 
			function(w){
				if(w > 0)
					return(1)
				else if(w < 0)
					return(-1)
			}
	)
	t$V3 = converted.weights
	
	write(first.line, unweighted.network.path)
	write.table(
			t, 
			unweighted.network.path, 
			sep="\t", 
			append=TRUE, 
			row.names=FALSE, 
			col.names=FALSE
	)
	
	return(t)
}

