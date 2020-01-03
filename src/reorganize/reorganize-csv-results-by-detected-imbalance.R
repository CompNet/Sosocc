

#################################################################
#
#################################################################
reorganize.network.folder.by.detected.imbalance = function(counter.per.detected.imb.interval, detected.imb.val, n, l0, d, prop.mispl,
                                                           prop.neg, k, network.no, cor.clu.exact.algo, graph.desc.name)
{
    detected.imb.interval.seq.vals = seq(DETECTED.IMB.LOWER.BOUND, DETECTED.IMB.UPPER.BOUND, 0.05)
    
	for(i in 1:(length(detected.imb.interval.seq.vals)-1)){
	    lower.bound = detected.imb.interval.seq.vals[i]
	    upper.bound = detected.imb.interval.seq.vals[i+1]
	    desc = paste0("[",lower.bound,",",upper.bound,")")
		if(detected.imb.val>=lower.bound && detected.imb.val<upper.bound){
		    if(desc %in% names(counter.per.detected.imb.interval))
		        counter.per.detected.imb.interval[[desc]] = counter.per.detected.imb.interval[[desc]] + 1
		    else
		        counter.per.detected.imb.interval[[desc]] = 1
		    
		    new.network.no = counter.per.detected.imb.interval[[desc]]
		    new.eval.folder = get.eval.folder.path(n, l0, d, NA, desc, prop.neg, k=ALL, new.network.no, cor.clu.exact.algo, NA)
            dir.create(new.eval.folder, showWarnings=FALSE, recursive=TRUE)

            source.eval.folder = get.eval.folder.path(n, l0, d, prop.mispl, NA, prop.neg, k=ALL, network.no, cor.clu.exact.algo, graph.desc.name)
            file.copy(source.eval.folder, new.eval.folder, recursive=TRUE)
		    
            # at network level, put the source information, i.e. association of initial network ids with new networks ids 
            new.network.eval.folder = get.eval.folder.path(n, l0, d, NA, desc, prop.neg, k=ALL, new.network.no, NA, NA)
            write.table(x=source.eval.folder, file=file.path(new.network.eval.folder,"source.txt"), row.names=F, col.names=F)
		}
	}
    return(counter.per.detected.imb.interval)
}


#################################################################
# 
#
# n: graph size
# l0: number of clusters
# d: density
# prop.mispl: proportion of misplaced links
# prop.neg: proportion of negative links
# network.no: network id (the identifiers start from 1)
# cor.clu.exact.algo: the name of correlation clustering algorithm that has been executed for optimal solutions
# force: whether or not the existing files are overwritten by a fresh call of all corresponding methods (e.g partitioning method)
#
##################################################################
reorganize.csv.results.by.detected.imbalance = function(n, l0, d, prop.mispls, prop.negs, in.rand.net.folders,
		cor.clu.exact.algo, force)
{
    k=ALL
    
    notExist = TRUE
    upper.folder = get.eval.folder.path(n, l0, d, NA, NA, NA, NA, NA, NA)
    existing.dirs = list.dirs(path = upper.folder, full.names = FALSE, recursive = FALSE)
    if(any(startsWith(existing.dirs, DETECTED.IMB.PROP.PARAM.NAME)))
        notExist = FALSE
    
    if(notExist || force){
        if(!notExist){ # force=TRUE and dirs do exist
            existing.dirs.path = list.dirs(path = upper.folder, full.names = TRUE, recursive = FALSE)
            indx = which(startsWith(existing.dirs, DETECTED.IMB.PROP.PARAM.NAME))
            print(existing.dirs.path[indx])
            unlink(existing.dirs.path[indx], recursive=TRUE)
        }
    
    	tlog(16, "start to evaluate exact solutions")
    	for(graph.desc.name in c(SIGNED.UNWEIGHTED.FILE)){ # SIGNED.WEIGHTED.FILE
    		tlog(20, "evaluating partitions => graph.desc.name: ", graph.desc.name)
    		tlog(20, "evaluating partitions => algo.name: ", cor.clu.exact.algo)

    		
    		
    		for(prop.mispl in prop.mispls){
    		    tlog(8, "evaluating partitions => prop.mispl: ", prop.mispl)
    		    
    		    if(is.na(prop.negs) && d==1){
    		        prop.negs = compute.prop.neg(n, d, l0, prop.mispl)
    		    }
    		    
    		    for(prop.neg in prop.negs){
    		        tlog(12, "evaluating partitions => prop.neg: ", prop.neg)
    		        
    		        counter.per.detected.imb.interval = list()

    		        net.folder = get.input.network.folder.path(n, l0, d, prop.mispl, prop.neg, network.no=NA)
    		        if(dir.exists(net.folder)){
    		            
                	    for(network.no in in.rand.net.folders){
                	        tlog(16, "evaluating partitions => network.no: ", network.no)
            
                	        #net.folder = get.input.network.folder.path(n, l0, d, prop.mispl, prop.neg, network.no)
                	        #e.part.folder = get.part.folder.path(n, l0, d, prop.mispl, prop.neg, network.no, cor.clu.exact.algo, graph.desc.name)
                	        e.eval.folder = get.eval.folder.path(n, l0, d, prop.mispl, NA, prop.neg, k=k, network.no, cor.clu.exact.algo, graph.desc.name)
                	        if(!dir.exists(e.eval.folder))
                	            dir.create(path=e.eval.folder, showWarnings=FALSE, recursive=TRUE)
            
                	        table.file = file.path(e.eval.folder, paste0(EVAL.IMB.INFO.TABLE.FILENAME,".csv"))
                	        print(table.file)
                	        if(file.exists(table.file)){
                	            print("exists")
                	            df = as.matrix(read.csv(table.file, row.names = 1, header= TRUE, check.names=FALSE, stringsAsFactors=FALSE))
                	            detected.imb.val = as.numeric(df[,IMB.PERC.COL.NAME])
                	            detected.imb.val = as.numeric(detected.imb.val)/100 # convert to proportion from percentage
                	            
                	            counter.per.detected.imb.interval =
                	                reorganize.network.folder.by.detected.imbalance(counter.per.detected.imb.interval, detected.imb.val, n, l0, d, prop.mispl,
                	                                                           prop.neg, k, network.no, cor.clu.exact.algo, graph.desc.name)
                	        } else {
                	            # if file does not exist, then there supposed to be an error !!!
                	            print("ERROR !!!")
                	        }
            
                	    } # end of updating 'detected.imb.data.str'
            
                		# print(counter.per.detected.imb.interval)
    		        }
    
    		    }
    		}
    	}
        
    }
}




#################################################################
#
#
# graph.sizes: a vector of values regarding graph sizes to be considered
# d: density (it is a single value)
# l0: number of clusters to be considered (it is a single value)
# prop.mispls: a vector of values regarding proportion of misplaced links
# prop.negs: a vector of values regarding proportion of negative links (for now, it is not operational)
# in.rand.net.folders: a vector of values regarding input random graph folders. Sequantial integers (1, .., 10)
# cor.clu.exact.algo: the name of correlation clustering algorithm that has been executed for optimal solutions
# force: whether or not the existing files are overwritten by a fresh call of all corresponding methods (e.g partitioning method)
#
##################################################################
reorganize.all.csv.results.by.detected.imbalance = function(graph.sizes, d, l0, prop.mispls, prop.negs, in.rand.net.folders,
		cor.clu.exact.algo, force)
{

	tlog("starts evaluating partitions")
	for(n in graph.sizes){
		tlog(4, "evaluating partitions => n: ", n)

	    reorganize.csv.results.by.detected.imbalance(n, l0, d, prop.mispls, prop.negs, in.rand.net.folders,
				cor.clu.exact.algo, force)
	        
	}
	
}
