# TODO: Add comment
# 
# Author: nejat
###############################################################################






#################################################################
#
# It constructs a single row of the csv file in output. It retrieves the statistics 
#   related to the execution times of a partitioning method for the given values of the input parameters. 
#
# n: graph order
# l0: number of modules
# d: density
# prop.mispl: proportion of misplaced links
# prop.neg: proportion of negative links
# network.no: network id (the identifiers start from 1
# cor.clu.exact.algo: the name of correlation clustering algorithm to run
# colnames: column names of the csv file in output
# force: whether or not the existing files are overwritten by a fresh call of all corresponding methods (e.g partitioning method)
#
##################################################################
retreive.exec.time = function(n, l0, d, prop.mispl, prop.neg, network.no,
		cor.clu.exact.algo, colnames, force)
{
	
	curr.data =c()
		
	tlog(20, "collecting => algo.name: ", cor.clu.exact.algo)
	
	for(graph.desc.name in c(SIGNED.UNWEIGHTED.FILE)){ # OPPOSITE.SIGNED.UNWEIGHTED.FILE
		tlog(24, "collecting => graph.desc.name: ", graph.desc.name)
		part.folder = part.folder = get.part.folder.path(n, l0, d, prop.mispl, prop.neg, network.no, cor.clu.exact.algo, graph.desc.name)
		eval.folder = get.eval.folder.path(n, l0, d, prop.mispl, NA, prop.neg, k=ALL, network.no, cor.clu.exact.algo, graph.desc.name)

		#mems = load.membership.files(part.folder)
		#nb.sol = length(mems)
        print("----")
		#print(nb.sol)		
		#if(nb.sol>1){
            table.file = file.path(eval.folder, paste0(EVAL.EXEC.TIME.FILENAME,".csv"))
            print(table.file)

            #curr.line = matrix("NA",1,1)
            #colnames(curr.line) = EXEC.TIME.COL.NAME
            if(dir.exists(part.folder) && file.exists(table.file))
		        curr.data = read.csv(table.file, row.names = 1, header= TRUE, check.names=FALSE, stringsAsFactors=FALSE)
            #print(curr.data)

		#} else {
		#    curr.line = matrix(NA,nrow=1,ncol=length(colnames))
        #    colnames(curr.line) = colnames
		#    rownames(curr.line) = paste0("n=", n, "l0=", l0, "d=", d, "prop.mispl=", prop.mispl, "prop.neg=", prop.neg, "network.no=", network.no)
		#    curr.data = rbind(curr.data, curr.line)
		#}
	}
	
	
	return(curr.data)
}









#################################################################
#
# It is the starting method in the aim of collecting the statistics related to execution times of an partitioning method. 
#   It handles all networks by graph.sizes,  prop.mispls, my.prop.negs and in.rand.net.folders
#   For the given value of the parameter l0, it generates a csv file.
#
# graph.sizes: a vector of values regarding graph orders to be considered
# d: density (it is a single value)
# l0: number of modules to be considered (it is a single value)
# prop.mispls: a vector of values regarding proportion of misplaced links
# prop.negs: a vector of values regarding proportion of negative links (for now, it is not operational)
# in.rand.net.folders: a vector of values regarding input random graph folders. Sequantial integers (1, .., 10)
# cor.clu.exact.algos: the names of correlation clustering algorithms to run
# force: whether or not the existing files are overwritten by a fresh call of all corresponding methods (e.g partitioning method)
#
##################################################################
collect.all.exec.times = function(graph.sizes, d, l0, prop.mispls, prop.negs, in.rand.net.folders,
		cor.clu.exact.algos, force)
{
    

    	all.data = c()
    	
    	tlog("starts collecting")
    	for(n in graph.sizes){
    		tlog(8, "collecting networks => n: ", n)
    		
    		for(prop.mispl in prop.mispls){
    			tlog(8, "collecting => prop.mispl: ", prop.mispl)
    			
                my.prop.negs = prop.negs # if we do not do that, for each n value, prop negs will not be the initial value(s)
                if(is.na(my.prop.negs) && d == 1){
                    my.prop.negs = compute.prop.neg(n, d, l0, prop.mispl)
                }
                
                for (prop.neg in my.prop.negs) {
    				tlog(12, "collecting => prop.neg: ", prop.neg)
    				
    				for(network.no in in.rand.net.folders){
    					tlog(16, "collecting => network.no: ", network.no)
    					
						curr.data = c()
                        for(cor.clu.exact.algo in cor.clu.exact.algos){
                            data = retreive.exec.time(n, l0, d, prop.mispl, prop.neg, network.no,
        							cor.clu.exact.algo, force)

						    if(length(data)>0)
							    curr.data = c(curr.data, data)
                            else
                                curr.data = c(curr.data, NA)

                            
                        }
                        curr.data = matrix(curr.data, nrow=1)
                        rownames(curr.data) = paste0("n=", n, ",l0=", l0, ",d=", d, ",prop.mispl=", prop.mispl, ",prop.neg=", prop.neg, ",network.no=", network.no)
    					all.data = rbind(all.data, curr.data)
    				}
    				
    			}
    			
    		}
    		
    	}

    	if(length(all.data)>0){
    		colnames(all.data) = cor.clu.exact.algos
    		#rownames(all.data) = NULL
    		if(!dir.exists(OUTPUT.CSV.FOLDER))
    		    dir.create(OUTPUT.CSV.FOLDER, recursive=FALSE, showWarnings=FALSE)
    		write.csv(file=file.path(OUTPUT.CSV.FOLDER,paste0("exec-time-l0=",l0,"_n=",paste0(graph.sizes, collapse=","),"_d=",d,".csv")), x=all.data)
    	}

}
