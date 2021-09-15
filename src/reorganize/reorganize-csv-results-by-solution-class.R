





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
reorganize.csv.results.by.sol.class = function(n, l0, d, prop.mispl, detected.imb.interval, prop.neg, in.rand.net.folders,
		cor.clu.exact.algo, comp.measures, force)
{
    notExist = TRUE
    upper.folder = get.eval.folder.path(n, l0, d, prop.mispl, detected.imb.interval, prop.neg, NA, NA, NA, NA)
    existing.dirs = list.dirs(path = upper.folder, full.names = FALSE, recursive = FALSE)
    if(any(startsWith(existing.dirs, "k=")))
        notExist = FALSE
    
    if(notExist || force){
        if(!notExist){ # force=TRUE and dirs do exist
            existing.dirs.path = list.dirs(path = upper.folder, full.names = TRUE, recursive = FALSE)
            indx = which(startsWith(existing.dirs, "k="))
            print(existing.dirs.path[indx])
            unlink(existing.dirs.path[indx], recursive=TRUE)
        }
    
    	tlog(16, "start to evaluate exact solutions")
    	for(graph.desc.name in c(SIGNED.UNWEIGHTED.FILE)){ # SIGNED.WEIGHTED.FILE
    		tlog(20, "reorganizing => graph.desc.name: ", graph.desc.name)
    		tlog(20, "reorganizing => algo.name: ", cor.clu.exact.algo)
    		
    		for(measure in comp.measures){
    		    tlog(20, "reorganizing => measure: ", measure)

        		networks.clustering.list = list()
    
        	    for(network.no in in.rand.net.folders){
        	        tlog(16, "reorganizing => network.no: ", network.no)
    
        	        e.eval.folder = get.eval.folder.path(n, l0, d, prop.mispl, detected.imb.interval, prop.neg, k=ALL, network.no, cor.clu.exact.algo, graph.desc.name)
        	        if(!dir.exists(e.eval.folder))
        	            dir.create(path=e.eval.folder, showWarnings=FALSE, recursive=TRUE)
    
        	        table.file = file.path(e.eval.folder, paste0(EVAL.BEST.K.FOR.KMEDOIDS,"-",measure,".csv"))
        	        #print(table.file)
        	        if(file.exists(table.file)){
        	            #print("exists")
        	            df = as.matrix(read.csv(table.file, row.names = 1, header= TRUE, check.names=FALSE, stringsAsFactors=FALSE))
        	            k = as.numeric(df[,BEST.K.FOR.SILH.COL.NAME])
    
        	            # think networks.clustering.list[[as.character(k)]] as an empty vector, so we can append
        	            networks.clustering.list[[as.character(k)]] = c(networks.clustering.list[[as.character(k)]], network.no)
        	        } else {
        	            # if file does not exist, then there supposed to be single optimal partition) ==> it is different than k=1
        	            networks.clustering.list[[NOT.AVAILABLE.FOLDER.NAME]] = c(networks.clustering.list[[NOT.AVAILABLE.FOLDER.NAME]], network.no)
        	        }
    
        	    } # end of updating 'networks.clustering.list'
    
        		# EXAMPLE OF 'networks.clustering.list': keys represent k values, and values represent network ids.
        		# $`1`
        		# [1] 5 4
        		# $`2`
        		# [1] 2
    
        		keys = names(networks.clustering.list)
        		for(k in keys){
        	        network.no.vec = networks.clustering.list[[k]]
        	        for(network.no in network.no.vec){
        	            # k.val = as.integer(k)
        	            # if(k == NOT.AVAILABLE.FOLDER.NAME)
        	            #     k.val = as.character(paste0("k=",NOT.AVAILABLE.FOLDER.NAME))
        	            k.val = paste0("k=",k,"-",measure)
        	            new.eval.folder = get.eval.folder.path(n, l0, d, prop.mispl, detected.imb.interval, prop.neg, k=k.val, network.no, cor.clu.exact.algo, NA)
        	            dir.create(new.eval.folder, showWarnings=FALSE, recursive=TRUE)
    
        	            source.eval.folder = get.eval.folder.path(n, l0, d, prop.mispl, detected.imb.interval, prop.neg, k=ALL, network.no, cor.clu.exact.algo, graph.desc.name)
        	            file.copy(source.eval.folder, new.eval.folder, recursive=TRUE)
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
reorganize.all.csv.results.by.sol.class = function(graph.sizes, d, l0, prop.mispls, detected.imb.intervals, prop.negs, in.rand.net.folders,
		cor.clu.exact.algo, comp.measures, force)
{

	tlog("starts reorganizing")
	for(n in graph.sizes){
		tlog(4, "reorganizing => n: ", n)

	    # ===========================================================
	   
		for(prop.mispl in prop.mispls){
			tlog(8, "reorganizing => prop.mispl: ", prop.mispl)
			# my.prop.negs = NA
			# if(is.na(prop.negs)){
			# 	prop.neg = compute.prop.neg(n, d, l0, prop.mispl)
			# 	my.prop.negs = prop.neg # we show explicitely that we get only 1 value since density=1.0
			# }
		    
		    eval.upper.folder = get.eval.folder.path(n, l0, d, prop.mispl, NA, prop.neg=NA, network.no=NA)
		    if(dir.exists(eval.upper.folder)){
    		    existing.folders = list.dirs(path = eval.upper.folder, full.names = FALSE, recursive = FALSE)
    		    my.prop.negs = as.numeric(gsub("propNeg=","",existing.folders))
    
    			for(prop.neg in my.prop.negs){
    				tlog(12, "reorganizing => prop.neg: ", prop.neg)
    
    			    net.folder = get.input.network.folder.path(n, l0, d, prop.mispl, prop.neg, network.no=NA)
    			    if(dir.exists(net.folder)){
    			        
        			    reorganize.csv.results.by.sol.class(n, l0, d, prop.mispl, NA, prop.neg, in.rand.net.folders,
        						cor.clu.exact.algo, comp.measures, force)
    			    }
    			}
		    }
		}
	    
	    # ======== end
	    # ===========================================================	    
	  
       # # retreive the detected imbalance intervals used for n, l0 and d
       # eval.upper.folder = get.eval.folder.path(n, l0, d, prop.mispl=NA, detected.imb=NA, prop.neg=NA, network.no=NA)
       # existing.folders = list.dirs(path = eval.upper.folder, full.names = FALSE, recursive = FALSE)
       # indx = which(startsWith(existing.folders, DETECTED.IMB.PROP.PARAM.NAME))
       # detected.imb.intervals = gsub(paste0(DETECTED.IMB.PROP.PARAM.NAME,"="),"",existing.folders[indx])

       for(detected.imb.interval in detected.imb.intervals){
           tlog(8, "reorganizing => detected.imb.interval: ", detected.imb.interval)
           
           eval.upper.folder = get.eval.folder.path(n, l0, d, NA, detected.imb.interval, prop.neg=NA, network.no=NA)
           if(dir.exists(eval.upper.folder)){
               existing.folders = list.dirs(path = eval.upper.folder, full.names = FALSE, recursive = FALSE)
               my.prop.negs = as.numeric(gsub("propNeg=","",existing.folders))
               
               for(prop.neg in my.prop.negs){
                   tlog(12, "reorganizing => prop.neg: ", prop.neg)
                   
                   eval.upper.folder = get.eval.folder.path(n, l0, d, NA, detected.imb.interval, prop.neg=prop.neg, k=ALL, network.no=NA)
                   existing.folders = list.dirs(path = eval.upper.folder, full.names = FALSE, recursive = FALSE)
                   print(existing.folders)
                   my.in.rand.net.folders = as.integer(gsub("network=","",existing.folders)) # different than 'in.rand.net.folders'
    
                   reorganize.csv.results.by.sol.class(n, l0, d, NA, detected.imb.interval, prop.neg, my.in.rand.net.folders,
                                                       cor.clu.exact.algo, comp.measures, force)
               }
           }
       }
	   
	   # ======== end
	    
	    

	}
	
}
