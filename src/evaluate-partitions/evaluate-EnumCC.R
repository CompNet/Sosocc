


#################################################################
# It evaluates the solutions of a given network based on the considered algorithm name and graph type (weighted or not, etc.).
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
process.evaluation.EnumCC = function(n, l0, d, prop.mispl, prop.neg, network.no,
		exact.enum.algo, force)
{
	net.folder = get.input.network.folder.path(n, l0, d, prop.mispl, prop.neg, network.no)
	
	tlog(16, "start to evaluate exact solutions")
	for(graph.desc.name in c(SIGNED.UNWEIGHTED.FILE)){ # SIGNED.WEIGHTED.FILE
		tlog(20, "evaluating partitions => graph.desc.name: ", graph.desc.name)
		tlog(20, "evaluating partitions => algo.name: ", exact.enum.algo)
		
		part.folder = get.part.folder.path(n, l0, d, prop.mispl, prop.neg, network.no, exact.enum.algo, graph.desc.name)
		eval.folder = get.eval.folder.path(n, l0, d, prop.mispl, detected.imb=NA, prop.neg, k=ALL, network.no, exact.enum.algo, graph.desc.name)
		if(!dir.exists(eval.folder))
			dir.create(path=eval.folder, showWarnings=FALSE, recursive=TRUE)
		
		tlog(20, "proceed only exact solutions in folder: ", eval.folder)

        m = read.csv(file.path(eval.folder, paste0(EVAL.NB.SOL.FILENAME,".csv")), check.names=F)[,NB.SOL.COL.NAME] # nb solutions

        if(file.exists(file.path(part.folder, EXEC.TIME.FILENAME))){

		    # ---------------------------------------------------------------------
		    # nb enum
		    tlog(24, "process algo evaluation: nb enum")
            # there are as many enumerations as number of directories (e.g. 3, when there are "1", "2", "3" directories) 
            nb.enum = length(list.dirs(path=part.folder, full.names=FALSE))-1 # -1 for excluding "." from the results
		    table.file = file.path(eval.folder, paste0(EVAL.NB.ENUM.FILENAME,".csv"))
		    mtrx = matrix(nb.enum, 1, 1) # put the value into matrix, that way we will process all results in the same way (i.e. via 'csv')
		    colnames(mtrx) = NB.ENUM.COL.NAME
		    rownames(mtrx) = exact.enum.algo
		    write.csv(x=mtrx, file=table.file, row.names=TRUE)
		    # ---------------------------------------------------------------------


		    # ---------------------------------------------------------------------
		    # exec time
		    table.file = file.path(eval.folder, paste0(EVAL.EXEC.TIME.FILENAME,".csv"))
		    if(!file.exists(table.file) || force){
        
                exec.time = 0
                if(!file.exists(file.path(part.folder,"jump-exec-time1.txt")) || m>=ENUM.POPULATE.CC.MAX.NB.SOLS){
                    ## read txt file, and write it into a csv file
			        result.f.path = file.path(part.folder, EXEC.TIME.FILENAME)
			        ## there is a single 'exec time' value
			        exec.time = read.table(result.f.path)$V1 # this is total exec time
                } else {
                    for(i in 1:nb.enum){
                        result.f.path = file.path(part.folder, i, "execTime.txt")
                        exec.time = exec.time + read.table(result.f.path)$V1 # this is total exec time    
                        exec.time = exec.time + read.table(file.path(part.folder,paste0("jump-exec-time",i,".txt")))$V1
                    }
                }


                ExCC.part.folder = file.path(part.folder,"..","..",get.ExCC.code(enum.all=FALSE),"signed-unweighted")
                ExCC.table.file = file.path(ExCC.part.folder, paste0(EVAL.EXEC.TIME.FILENAME,".txt"))
                ExCC.exec.time = read.table(ExCC.table.file)$V1 # this is total exec time

				tot.exec.time = exec.time+ExCC.exec.time
				max.exec.time = ExCC.exec.time+ENUM.POPULATE.CC.MAX.TIME.LIMIT
				if(tot.exec.time>max.exec.time)
					tot.exec.time = max.exec.time
                mtrx = matrix(tot.exec.time, 1, 1) # matrix with 1 element
			    colnames(mtrx) = c(EXEC.TIME.COL.NAME)
			    rownames(mtrx) = exact.enum.algo
			    write.csv(x=mtrx, file=table.file, row.names=TRUE)
            }
            # ---------------------------------------------------------------------


            # ---------------------------------------------------------------------
		    # delay exec time
            table.file = file.path(eval.folder, paste0(EVAL.DELAY.EXEC.TIME.FILENAME,"-",get.ExCC.code(enum.all=TRUE),".csv"))
		    if(!file.exists(table.file) || force){
                # EnumCC
                table.file.exec.time = file.path(eval.folder, paste0(EVAL.EXEC.TIME.FILENAME,".csv"))
                exec.time = read.csv(table.file.exec.time, header=T, check.names=F)[,EXEC.TIME.COL.NAME] # this is total exec time of EnumCC

                # ExCC-All
                ExCC.all.eval.folder = file.path(eval.folder,"..","..",get.ExCC.code(enum.all=TRUE),"signed-unweighted")
                ExCC.all.table.file = file.path(ExCC.all.eval.folder, paste0(EVAL.EXEC.TIME.FILENAME,".csv"))
                ExCC.all.exec.time = read.csv(ExCC.all.table.file, header=T, check.names=F)[,EXEC.TIME.COL.NAME]

                ## ExCC.all: we extract only thime spent during the second phase
                #ExCC.all.part.folder = file.path(part.folder,"..","..",get.ExCC.code(enum.all=TRUE),"signed-unweighted")
                #ExCC.all.table.file = file.path(ExCC.all.part.folder, "logcplex.txt")
                ##print(table.file2)
                #ExCC.all.second.phase.exec.time = retreive.ExCC.All.second.phase.exec.time.from.log.file(ExCC.all.table.file)
		        #ExCC.all.exec.time = ExCC.exec.time + ExCC.all.second.phase.exec.time

                delay.exec.time = exec.time - ExCC.all.exec.time
                mtrx = matrix(delay.exec.time, 1, 1) # matrix with 1 element
		        colnames(mtrx) = c(DELAY.EXEC.TIME.COL.NAME)
		        rownames(mtrx) = exact.enum.algo
                table.file3 = file.path(eval.folder, paste0(EVAL.DELAY.EXEC.TIME.FILENAME,".csv"))
		        write.csv(x=mtrx, file=table.file, row.names=TRUE)

                # do the same process in the folder of ExCC.all (with symmetry)
                ExCC.all.eval.folder = file.path(eval.folder,"..","..",get.ExCC.code(enum.all=TRUE),"signed-unweighted")
                delay.exec.time = ExCC.all.exec.time - exec.time 
                mtrx = matrix(delay.exec.time, 1, 1) # matrix with 1 element
		        colnames(mtrx) = c(DELAY.EXEC.TIME.COL.NAME)
		        rownames(mtrx) = get.ExCC.code(enum.all=TRUE)
                table.file = file.path(ExCC.all.eval.folder, paste0(EVAL.DELAY.EXEC.TIME.FILENAME,"-",exact.enum.algo,".csv"))
		        write.csv(x=mtrx, file=table.file, row.names=TRUE)
            }
		    # ---------------------------------------------------------------------

        }

	}
}






#################################################################
# It is the starting method in the aim of evaluating the solutions of the considered networks. 
#   It handles all networks by graph.sizes,  prop.mispls, my.prop.negs and in.rand.net.folders
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
evaluate.EnumCC = function(graph.sizes, d, l0, prop.mispls, prop.negs, in.rand.net.folders,
		exact.enum.algo, force)
{
    ## check the necessary condition
    #if(cor.clu.exact.algo == get.ExCC.code(enum.all=FALSE)){
    #    tlog(4, "The partitioning method should be the one which enumerating all optimal solutions, i.e. 'ExCC-All' !")
    #    return(0)
    #}
    
    
	tlog("starts evaluating partitions")
	for(n in graph.sizes){
		tlog(4, "evaluating partitions => n: ", n)
		
		for(prop.mispl in prop.mispls){
			tlog(8, "evaluating partitions => prop.mispl: ", prop.mispl)
			
            my.prop.negs = prop.negs # if we do not do that, for each n value, prop negs will not be the initial value(s)
            if(is.na(my.prop.negs) && d == 1){
                my.prop.negs = compute.prop.neg(n, d, l0, prop.mispl)
            }
            
            for (prop.neg in my.prop.negs) {
				tlog(12, "evaluating partitions => prop.neg: ", prop.neg)
				
								
				net.folder = get.input.network.folder.path(n, l0, d, prop.mispl, prop.neg, network.no=NA)
				if(dir.exists(net.folder)){
				
					for(network.no in in.rand.net.folders){
						tlog(16, "evaluating partitions => network.no: ", network.no)
						
						process.evaluation.EnumCC(n, l0, d, prop.mispl, prop.neg, network.no,
								exact.enum.algo, force)
						
					}
				}
				
			}
			
		}
		
	}
	
}
