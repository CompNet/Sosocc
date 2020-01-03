

#################################################################
# It checks if node 'vertex' can be placed into cluster 'cluster'.
#  The function verifies if all the nodes of 'cluster' and 'vertex' have the value of TRUE in comembership
#
# el: edge list
# comembership: binary vector values indicating if pairs of nodes are together or not in the partition of core part 
# cluster: ids of nodes located at the cluster for whcih merge is desired (so, not cluster id)
# vertex: a node id for which merge operation is investigated into 'cluster'
#
#################################################################
check.cluster.consistency = function(el, comembership, cluster, vertex){
    result = TRUE
    for(v in cluster){
        indx1 = which(el[,1] == v & el[,2] == vertex)
        indx2 = which(el[,1] == vertex & el[,2] == v)
        if( !(length(indx1)>0 && comembership[indx1]) && !(length(indx2)>0 && comembership[indx2]) ){
            result = FALSE
            return(result)
        }
    }
    return(result)
}

#################################################################
# It checks if 'partition2' is already found in 'partitions'.
# Note that 'partitions' is a matrix whose columns correspond to partitions and rows correspond to nodes.
# If there is any NA value, the function removes it for processing. If the number of NA values in 'partition2' and any partition in 'partitions' 
#   is not the same, then, it stops (because it is not normal)
#
# partitions: matrix whose columns correspond to partitions and rows correspond to nodes
# partition2: a partition candidate that will be inserted into 'partitions2' if an identical partitions is not already there
#
#################################################################
are.partitions.identical = function(partitions, partition2){
    # print("----- BEGIN ----")
    # print(partitions)
    # print(partition2)
    na.indx = which(is.na(partition2))
    if(length(na.indx>0))
        partition2 = partition2[-na.indx]
    
    for(c in 1:ncol(partitions)){
        partition1 = partitions[,c]
        if(length(na.indx>0))
            partition1 = partition1[-na.indx]
        na.indx1 = which(is.na(partition1))
        if(length(na.indx1) == 0){
            partition = create.relative.plot.membership(partition1, partition2)
            if(all(partition == partition1))
                return(TRUE)
        }
    }
    return(FALSE)
}


#################################################################
# This function is run when CORE.PARTH.THRESHOLD < 1. This is because, the edgelist 'comembership' does not allow to create a non-overlapping partition.
#   So, we end up having fuzzy partition. This means that there might be several non-overlapping partitions (after exploring the corresponding fuzzy partition)
#   So, the output of this function can be visualized as a fuzzy partition.
#   Of course, this function can be used for CORE.PARTH.THRESHOLD = 1, but we think that reading the function 'determine.core.part.strict.agreement.partition()'
#       is much more easier!!
#   Before explaining the steps, lets give an example (we suppose that the vector 'comembership' is prepared with a specific value of CORE.PARTH.THRESHOLD):
#       edge.list     comembership
#       [,1] [,2]
# [1,]    1    2          TRUE
# [2,]    1    3          FALSE
# [3,]    1    4          TRUE
# [4,]    1    5          FALSE
# [5,]    2    3          TRUE
# [6,]    2    4          FALSE
# [7,]    2    5          TRUE
# [8,]    3    4          FALSE
# [9,]    3    5          TRUE
# [10,]   4    5          FALSE
#       In this example, wee see that nodes 1 and 2 are together, and nodes1 and 4 are together. But, nodes 2 and 4 are NOT together ==> fuzzy nature
#       So, the results that we expect are the followings:
#       partition1 = c(1, 1, 2, 3, 2)
#       partition2 = c(1, 2, 2, 1, 2)
#   We see that we may obtain several partitions for core part. This is only possible when CORE.PARTH.THRESHOLD < 1.
#
#
#  THE STEPS:
#   we iterate over the edgelist and comemberhsip at the same time. 
#   => for each pair of nodes i and j that will be together in the partition of the core part:
#       1) if we check if nodes i and j are not already assigned to any cluster, all is ok, so put them together into a new cluster
#       2) if one of them is assigned to a cluster, let us say node i is assigned to cluster.i, and node j is not assigned yet to any cluster
#           2.1) Can we place node j into the cluster that node i belongs to, i.e. cluster.i ?
#                2.1.1) If yes (according to cluster consistency), we check 
#                       if this updated partition (i.e. after adding j into cluster.i) is already in the variable 'partitions'
#                2.1.2) If no, this causes an inconsistency. We will do two operations
#                   2.1.2.1) keep node i in its cluster, which is cluster.i. And, search for another cluster for node j
#                       2.1.2.1.1)  If node j can be placed into another existing cluster, put it into that cluster   
#                       2.1.2.1.2)  If node j can NOT be placed into another existing cluster, put it into a new cluster, i.e. singleton cluster
#                   2.1.2.2) keep nodes i and j in the same cluster and that cluster should be different than cluster.i. 
#                       2.1.2.2.1) Since we know the cluster that node j belongs to in 2.1.2.1), we put node i into the cluster of node j, i.e. cluster.j
#                                   Note that, after this step the resting nodes in cluster.i may be merged with another existing cluster.
#                                   We deal with it in step 4). 
#                                       => An example can be seen in our toy example in the iteration of nodes i=2 and j=3
#                                           Therein, nodes 1 and 4 become singletons, nonetheless, they can be merged
#       3) if both nodes are already assigned to different clusters, but if they can be together according to comembership
#           3.1) if possible (according to cluster consistency), place node i into the cluster that node j belongs to, i.e. cluster.j
#           3.2) if possible (according to cluster consistency), place node j into the cluster that node i belongs to, i.e. cluster.i
#       4) In the end of the loop => for each partition,
#           try to merge clusters that are consistent (sometimes, we end up having 2 diff clusters that can be merged)
#
# To test:
# n = 5 # number of nodes in the graph
# el <- t(combn(x=n, m=2, simplify=TRUE)) # edgelist
# comembership = (round(runif(nrow(el))) == TRUE)
# determine.core.part.relax.agreement.partition(comembership, n)
# print(cbind(el,comembership))
#
#
# OUTPUT: a matrix whose colums correspond to partitions and rows correspond to nodes
#  Example:
#         partition   partition
# [1,]         1         1
# [2,]         1         2
# [3,]         2         2
# [4,]         3         1
# [5,]         1         2
#################################################################
determine.core.part.relax.agreement.partition = function(comembership, n)
{
    get.max.cluster.no = function(partition){
        indx = which(!is.na(partition))
        if(length(indx) == 0)
            return(0)
        return(max(partition[indx]))
    }
    
    el <- t(combn(x=n, m=2, simplify=TRUE)) # edgelist
    
    partitions = c()
    partition = rep(NA,n) # init
    partitions = cbind(partitions, partition)
    same.cluster.indx = which(comembership == TRUE) # i.e. core part index
    #print(same.cluster.indx)	
    
    # handle if there is not any common 'same cluster' info ==> i.e. no core part detected
    if(length(same.cluster.indx) != 0){
        el <- t(combn(x=n, m=2, simplify=TRUE)) # edgelist
        same.cluster.el = el[same.cluster.indx,] # edgelist in which edge pair in the same cluster
        #print(same.cluster.el)
        if(length(same.cluster.indx)==1 && !is.matrix(same.cluster.el)) # handle when length(same.cluster.indx) = 1
            same.cluster.el = as.matrix(t(same.cluster.el))
        
        for(r in 1:nrow(same.cluster.el)){ # for each pair of nodes that will be together in the membership of the core part
            i = same.cluster.el[r,1]
            j = same.cluster.el[r,2]
            #cat("i:",i,", j:",j,"\n")
            
            
            for(c in 1:ncol(partitions)){
                partition = partitions[,c]
                    if(!is.na(partition[i]) && !is.na(partition[j]) && partition[i] != partition[j]){
                        # print("none of i and j is NA => inconsistency")

                        # there is an inconsistency: create two different partitions to handle both cases
                        
                        #  >> 1) if possible (according to cluster consistency), place node i into the cluster that node j belongs to
                        partition2 = partition
                        cluster.j = which(partition2 == partition2[j])
                        if(check.cluster.consistency(el, comembership, cluster.j, i)){
                            partition2[i] = partition2[j]
                            if(!are.partitions.identical(partitions, partition2))
                                partitions = cbind(partitions, partition2)
                        }
                            
                        #  >> 2) if possible (according to cluster consistency), place node j into the cluster that node i belongs to
                        partition2 = partition
                        cluster.i = which(partition2 == partition2[i])
                        if(check.cluster.consistency(el, comembership, cluster.i, j)){
                            partition2[j] = partition2[i]
                            if(!are.partitions.identical(partitions, partition2))
                                partitions = cbind(partitions, partition2)
                        }

                    }
                    else if(is.na(partition[i]) && is.na(partition[j])){ # no inconsistency
                        #print("i  and j are NA")
                        partition[c(i,j)] = get.max.cluster.no(partition) + 1
                        partitions[,c] = partition
                    }
                    else if( (!is.na(partition[i]) && is.na(partition[j])) || (is.na(partition[i]) && !is.na(partition[j])) ) {
                        # print("one of i and j is NA")
                        # ======================================================
                        # the code block works when !is.na(partition[i]) && is.na(partition[j]) = TRUE
                        # so, if is.na(partition[i])=TRUE, then swap i and j
                       if(is.na(partition[i])){
                           tmp = i
                           i = j
                           j = tmp
                       }                        
                       # ======================================================
                        
                       # 1) Can we place node j into the cluster that node i belongs to ?
                        cluster.i = which(partition == partition[i])
                        if(check.cluster.consistency(el, comembership, cluster.i, j)){ # if TRUE, the answer is yes
                            partition[j] = partition[i]
                            # then, we check if this updated partition is already in the variable 'partitions'
                            if(!are.partitions.identical(partitions, partition)) # it is not there, then put it into 'partitions'
                                partitions[,c] = partition # update
                            else # otherwise, remove the current investigated partition, since there is an existing one
                                partitions = partitions[,-c] # update
                        }
                        else { # Otherwise => the answer is no. So, this causes an inconsistency. We will do two operations
    
                            # 1) keep node i in its cluster, which is partition[i]. And, search for another cluster for node j
                            #      1.1)  If node j can be placed into another cluster, put it into that cluster
                            tmp = partition[-cluster.i]
                            other.cluster.labels = unique(tmp[which(!is.na(tmp))])
                            ok=FALSE
                            for(clu in other.cluster.labels){
                                cluster.i = which(partition == clu) # another existing cluster label
                                if(check.cluster.consistency(el, comembership, cluster.i, j)){ # if j can be together with all nodes of the cluster that i belongs to
                                    partition[j] = clu
                                    partitions[,c] = partition # update
                                    ok = TRUE
                                    break
                                }
                            }
                            if(!ok){ # 1.2)  If node j can NOT be placed into another cluster, put it into a new cluster, i.e. singleton cluster
                                partition[j] = get.max.cluster.no(partition) + 1
                                partitions[,c] = partition # update
                            }
       
                            # 2) keep nodes i and j in the same cluster and that cluster should be different than partition[i]. 
                            #       Since we know the cluster that node j belongs to in the 1st case, we put node i into partition[j]
                            partition[i] = partition[j]
                            if(!are.partitions.identical(partitions, partition))
                                partitions = cbind(partitions, partition) # update
                        }
                    }
                    
            }
            
            # print(partitions)
            
            # ==============================
            # for each partition, try to merge clusters that are consistent (sometimes, we end up having 2 diff clusters that can be merged)
            for(c in 1:ncol(partitions)){
                #cat("c:",c,"\n")
                partition = partitions[,c]
                #print(partition)

                any.change = FALSE
                remove.indx.for.partitions = c()
                change = TRUE
                while(change){ # contineu until no merge operation is performed
                    change = FALSE
                    cluster.labels = unique(partition[which(!is.na(partition))])
                    cluster.labels = cluster.labels[order(cluster.labels)]
                    if(length(cluster.labels)>1){
                        cluster.pairs = t(combn(x=cluster.labels, m=2, simplify=TRUE)) # edgelist
                        for(r in 1:nrow(cluster.pairs)){
                            clu1 = cluster.pairs[r,1]
                            clu2 = cluster.pairs[r,2]
                            #cat("clu1:",clu1,", clu2:",clu2,"\n")
                            cluster1 = which(partition == clu1) # cluster label
                            cluster2 = which(partition == clu2) # cluster label

                            
                            ok = TRUE
                            for(v in cluster1){
                                if(!check.cluster.consistency(el, comembership, cluster2, v)){
                                    ok = FALSE
                                    break
                                }
                            }
                            if(ok){ # merge cluster 1 and cluster 2
                                #print('OK MERGE !!!!!')
                                change=TRUE
                                any.change = TRUE
                                # add cluster1 into cluster2
                                partition[cluster1] = clu2 
                                # renumber partition vector, i.e. from c(2 2 3 4 2 3) to c(1 1 2 3 1 2)
                                partition = as.integer(as.factor(partition))
                                break
                            }
                        }
                    }
                } # end while
                
                
                if(any.change){
                    if(are.partitions.identical(partitions, partition)) # it is not there, then put it into 'partitions'
                        # we can not update 'partitions' here, as we have not finished yet the iteration over partitions
                        remove.indx.for.partitions = c(remove.indx.for.partitions, c)
                }
                
            }
            
            if(length(remove.indx.for.partitions)>0)
                partitions = partitions[,-remove.indx.for.partitions] # update
            # ==============================
        }
        
      
        # ================================================ 
        
        
        # I put this info as plot title
        na.indx = which(is.na(partition))
        partition[na.indx] = NA # I prefer not to put different clolor for each single cluster node => white
        
    } else {
        tlog("we can not find any core part")
        # I put this info as plot title
        partition = rep(NA,n)
        
    }
    
    return(partitions)
}


#################################################################
# This function is executed only for CORE.PART.THRESHOLD = 1. The reason is that it is easy to read the code (and also faster).
# The idea is to build a partition representing a core part structure, i.e. nodes that are always together over all partitions.
# This function does not treat the case where nodes that are never together over all partitions. They are treated in function "determine.core.part.disagreement.partition()"
#
# comembership: binary vector values indicating if pairs of nodes are together or not in the partition of core part 
# n: number of nodes in the graph
#
#################################################################
determine.core.part.strict.agreement.partition = function(comembership, n) # when threshold = 1
{
	partition = rep(NA,n) # init
	same.cluster.indx = which(comembership == TRUE) # i.e. core part index
	#print(same.cluster.indx)	

	# handle if there is not any common 'same cluster' info ==> i.e. no core part detected
	if(length(same.cluster.indx) != 0){
		el <- t(combn(x=n, m=2, simplify=TRUE)) # edgelist
		same.cluster.el = el[same.cluster.indx,] # edgelist in which edge pair in the same cluster
		#print(same.cluster.el)
		if(length(same.cluster.indx)==1 && !is.matrix(same.cluster.el)) # handle when length(same.cluster.indx) = 1
			same.cluster.el = as.matrix(t(same.cluster.el))
	
		cluster.no.counter = 1
		for(r in 1:nrow(same.cluster.el)){
			i = same.cluster.el[r,1]
			j = same.cluster.el[r,2]
			
			#print("-----")
			if(is.na(partition[i]) && is.na(partition[j])){
				partition[c(i,j)] = cluster.no.counter
				cluster.no.counter = cluster.no.counter + 1
			}
			else if(!is.na(partition[i]) && is.na(partition[j]))
				partition[j] = partition[i]
			else if(is.na(partition[i]) && !is.na(partition[j]))
				partition[i] = partition[j]
		#print(partition)
		}
	
	} else {
		tlog("we can not find any core part")
		# I put this info as plot title
		partition = rep(NA,n)
		
	}
	
	return(partition)
}

##################################################################
# It is wrapper function in order to choose the right function about processing of core part partition.
#
# core.part.threshold: threshold value for core part agreemnt. 1 means that all pairs of nodes is always together over all partitions.
# comembership: binary vector values indicating if pairs of nodes are together or not in the partition of core part 
# n: number of nodes in the graph
#
##################################################################
determine.core.part.agreement.partition = function(core.part.threshold, comembership, n){
    partition = NA
    if(core.part.threshold == 1)
        partition = determine.core.part.strict.agreement.partition(comembership, n)
    else {
        partitions = determine.core.part.strict.agreement.partition(comembership, n)  # a matrix
        # partition = partitions[,1] # take the first one, since it might be many
        partition = NA
    }
    return(partition)
}



##################################################################
# This function is a complement of the function 'determine.core.part.agreement.partition()'. 
#   The latter finds the nodes which are (almost) always together. The former has the aim of finding the nodes which are (almost) never together. 
# These nodes will appear as singleton node cluster in core part membership
#
# comembership: binary vector values indicating if pairs of nodes are together or not in the partition of core part 
# partition: partition obtained in function "determine.core.part.agreement.partition()"
#
##################################################################
determine.core.part.disagreement.partition = function(comembership, partition)
{
	n = length(partition)
	
	el <- t(combn(x=n, m=2, simplify=TRUE)) # edgelist
	clu.counter = length(unique( partition[which(!is.na(partition))] ))
	resting.nodes = which(is.na(partition)) # nodes except those which are always together over all partitions
	for(i in resting.nodes){
		#print(i)
	    
	    # [1,]    1    2
	    # [2,]    1    3
	    # [3,]    1    4
	    # [4,]    1    5
	    # [5,]    2    3
	    # [6,]    2    4
	    # [7,]    2    5
	    # [8,]    3    4
	    # [9,]    3    5
	    # [10,]   4    5
	    # indx should give = c(10, 3, 6, 8) for i=4
		indx = c(which(el[,1] == i), which(el[,2] == i))
		
		v.vec = comembership[indx]
		if(all(v.vec == TRUE)){ # if all values = TRUE, this means that node i is never together with another node => so, it is considered as a member of core part
			#print("a singleton node cluster is added !!!!")
			clu.counter = clu.counter + 1
			partition[i] = clu.counter
		}
	}
	
	return(partition)
}


#################################################################
# This function has the aim of creating a vector in which
#   TRUE values indicate pairs of nodes which are always (or almost always) together with anohter nodes in the same cluster over all partitions.
#   Those nodes are also considered a member of the core part.
#   For instance: c(TRUE, FALSE, FALSE, ..) => this indicates that the first pair of nodes (e.g. nodes 1 and 2) are always together in the same cluster
# Concretely, it does two things for each partition:
#   - it creates an edgelist => a matrix of 2 colums. Then, it constitutes the comembersip vector (i.e. whether two nodes belong to the same cluster or not)
#   An example:
#   [1,]    1    2
#   [2,]    1    3
#   [3,]    1    4
#   [4,]    2    3
#   [5,]    2    4
#   [6,]    3    4
#   In the end, a data frame where each row corresponds to comembership vector
#       part1  part2  part3
#   [1,] FALSE FALSE FALSE
#   [2,] TRUE TRUE TRUE
#   [3,] FALSE TRUE FALSE
#   ..
#   - Second, it creates a final vector, called 'agreement comembership'.
#       As said before, in this vector TRUE values indicate pairs of nodes which are always (or almost always with 'core.part.threshold') 
#       together with anohter nodes in the same cluster over all partitions.
#
#
# partitions: all partitions in the corresponding solution class in list format
# n: the number of nodes in the graph
# core.part.threshold: an integer value representing the threshold value for core part agreement. 1 means that all pairs of nodes is always together over all partitions.
#
#################################################################
determine.final.agreement.comembership = function(partitions, n, core.part.threshold)
{
	m = length(partitions)
	comembership.df = c()
	
	el <- t(combn(x=n, m=2, simplify=TRUE))
	for(partition in partitions){
		comembership <- sapply(1:nrow(el),function(r) partition[el[r,1]]==partition[el[r,2]])
		comembership.df = cbind(comembership.df, comembership)
	}
	
	
	final.comembership = sapply(1:nrow(comembership.df),
			function(r){ 
				true.freq = length(which(comembership.df[r,] == TRUE))
				true.prop = true.freq/m
				if(true.prop >= core.part.threshold)
					return(TRUE)
				else
					return(FALSE)
			})
	
	return(final.comembership)
}



#################################################################
# This function has the aim of creating a vector in which
#   TRUE values indicate pairs of nodes which are never (or almost never) together with anohter nodes in the same cluster over all partitions.
#   Those nodes are also considered a member of the core part.
#   For instance: c(TRUE, FALSE, FALSE, ..) => this indicates that the first pair of nodes (e.g. nodes 1 and 2) are never together in the same cluster
# Concretely, it does two things for each partition:
#   - it creates an edgelist => a matrix of 2 colums. Then, it constitutes the comembersip vector (i.e. whether two nodes belong to the same cluster or not)
#   An example:
#   [1,]    1    2
#   [2,]    1    3
#   [3,]    1    4
#   [4,]    2    3
#   [5,]    2    4
#   [6,]    3    4
#   In the end, a data frame where each row corresponds to comembership vector
#       part1  part2  part3
#   [1,] FALSE FALSE FALSE
#   [2,] TRUE TRUE TRUE
#   [3,] FALSE TRUE FALSE
#   ..
#   - Second, it creates a final vector, called 'disagreement comembership'.
#       As said before, in this vector TRUE values indicate pairs of nodes which are never (or almost never with 'core.part.threshold') 
#       together with anohter nodes in the same cluster over all partitions.
#
#
# partitions: all partitions in the corresponding solution class in list format
# n: the number of nodes n the graph
# core.part.threshold: an integer value representing the threshold value for core part agreement. 1 means that all pairs of nodes is always together over all partitions.
#
#################################################################
determine.final.disagreement.comembership = function(partitions, n, core.part.threshold)
{
	m = length(partitions)
	comembership.df = c()
	
	el <- t(combn(x=n, m=2, simplify=TRUE))
	for(partition in partitions){
		comembership <- sapply(1:nrow(el),function(r) partition[el[r,1]]==partition[el[r,2]])
		comembership.df = cbind(comembership.df, comembership)
	}
	
	
	final.comembership = sapply(1:nrow(comembership.df),
			function(r){
				false.freq = length(which(comembership.df[r,] == FALSE))
				false.prop = false.freq/m
				if(false.prop >= core.part.threshold)
					return(TRUE)
				else
					return(FALSE)
			})
		
	return(final.comembership)
}





#################################################################
# It is a wrapper function to find core part of a set of partitions.
#   It iterates over the set of 'core.part.threshold' values.
#
# partitions: all partitions in the corresponding solution class in list format
# characterization.silh.folder
# g: graph for plotting
# clu.desc: description string for the solution class. if it is equal to "k2", that means the second solution class in kmedoids results. 
#           When it is equal to "all", that means we compute the general core aprt, i.e. without any restriction of solution class => all partitions
# force: whether or not the existing files are overwritten by a fresh call of all corresponding methods (e.g partitioning method)
#
#################################################################
process.core.parts = function(partitions, characterization.silh.folder, g, clu.desc, force){
	
	for(core.part.threshold in CORE.PART.THRESHOLDS){
			characterization.core.folder = file.path(characterization.silh.folder,
				paste0(CORE.PART.FOLDER.NAME,"=",sprintf("%.2f",core.part.threshold)))
			if(!dir.exists(characterization.core.folder))
						dir.create(path=characterization.core.folder, showWarnings=FALSE, recursive=TRUE)
			
			characterization.file.path = file.path(characterization.core.folder, paste0(MBRSHP.FILE.PREFIX,"-",CLUSTER.ID.FOLDER.PREFIX,clu.desc,".txt"))
			if(!file.exists(characterization.file.path) || force){
    			#unlink(clu.characterization.measure.folder, recursive=TRUE)
    	
    			n = length(partitions[[1]]) # all partitions have the same size => n
    			
    			# the edge list information indicating which pairs of nodes are always together in the core part according to 'core.part.threshold'
    			final.agr.comembership = determine.final.agreement.comembership(partitions, n, core.part.threshold)
                # partial core part information: this is not final partition for core part
    			core.part.partition.without.singleton = determine.core.part.agreement.partition(core.part.threshold, final.agr.comembership, n)
    			
    			# the edge list information indicating which pairs of nodes are never together in the core part according to 'core.part.threshold'
    			final.disagr.comembership = determine.final.disagreement.comembership(partitions, n, core.part.threshold)
                # the final core part info
    			core.part.partition = determine.core.part.disagreement.partition(final.disagr.comembership, core.part.partition.without.singleton)
    			
    			
    			write.table(x=core.part.partition, row.names=FALSE, col.names=FALSE, file=characterization.file.path)
    
    			# plot
    			plot.file = file.path(characterization.core.folder, paste0(CLUSTER.ID.FOLDER.PREFIX,clu.desc))
    			plot.title = paste0("cluster ",clu.desc,"\n",
    				"core part threshold=", sprintf("%.2f",core.part.threshold), "\n",
    				"Each white node indicates single node cluster")
    			plot.network(g=g, membership=core.part.partition, plot.file, format="JPEG", method="circle", plot.title=plot.title)
			}
		}
		
}

#################################################################
# It is a wrapper function to find core part of a set of partitions.
#   It iterates over a range of 'k' values indicating the number of clusters, i.e. number of solution classes, in kmedoids analysis.
#   In each 'k', it loads all obtained partitions.
#
#
# g: graph for plotting
# analysis.silh.folder
# characterization.silh.folder
# clu.analysis.mbrshp: the clustering result of kmedoids for a specific number of solution classes, i.e. the variable 'k'
# k: a number of cluster in the kmedoids clustering, i.e. number of solution classes. Usually, the value associated with best silhouette is chosen.
# dist.mtrx: dissimilarity matrix
# part.folder: partition folder in which all partitions are stored
# force: whether or not the existing files are overwritten by a fresh call of all corresponding methods (e.g partitioning method)
#
#################################################################
determine.core.parts = function(g, analysis.silh.folder, characterization.silh.folder, clu.analysis.mbrshp, k, dist.mtrx,
											part.folder, force)
{
	for(clu.no in seq(1,k)){ # for each cluster of similar solutions, i.e. for each solution class
			
		part.ids = which(clu.analysis.mbrshp == clu.no)
		m = length(part.ids) # nb current solution
		
		analysis.clu.folder = file.path(analysis.silh.folder, paste0(CLUSTER.ID.FOLDER.PREFIX,clu.no))

					
		curr.partitions = list()		
		for(i in seq(1,m)){
			part.id = part.ids[i]-1 # minus 1 since ids start from 0
			curr.partitions[[i]] = read.table(file=file.path(analysis.clu.folder,paste0(MBRSHP.FILE.PREFIX,part.id,".txt")))$V1
		}
		#print(curr.partitions)
		process.core.parts(curr.partitions, characterization.silh.folder, g, clu.no, force)
	}
}




#################################################################
# It finds the core part over all solutions, without restriction of solution classes
#
# g: graph for plotting
# analysis.silh.folder
# characterization.silh.folder
# clu.analysis.mbrshp: the clustering result of kmedoids for a specific number of solution classes, i.e. the variable 'k'
# part.folder: partition folder in which all partitions are stored
# force: whether or not the existing files are overwritten by a fresh call of all corresponding methods (e.g partitioning method)
#
#################################################################
determine.general.core.part = function(g, analysis.silh.folder, characterization.silh.folder, clu.analysis.mbrshp, part.folder, force)
{
	partitions = load.membership.files(part.folder)
	m = length(clu.analysis.mbrshp)
	part.ids = seq(1,m)
	
	process.core.parts(partitions, characterization.silh.folder, g, ALL, force)
}



#################################################################
# It is a wrapper function to determine all core parts: 1) core parts of each solution class, 2) core part of all partitions (i.e. general core part)
#
#
# g: graph for plotting
# clu.characterization.measure.folder
# eval.folder
# clu.analysis.measure.folder
# dist.mtrx: dissimilarity matrix
# part.folder: partition folder in which all partitions are stored
# measure: dissimilarity measure which is used in k-mdedoids
# force: whether or not the existing files are overwritten by a fresh call of all corresponding methods (e.g partitioning method)
#
#################################################################
determine.all.core.parts = function(g, clu.characterization.measure.folder, eval.folder, clu.analysis.measure.folder, dist.mtrx, part.folder, measure, force)
{
	table.file = file.path(eval.folder, paste0(EVAL.BEST.K.FOR.KMEDOIDS,"-",measure,".csv"))
	if(file.exists(table.file)){
		df = read.csv(table.file, check.names=FALSE, stringsAsFactors=FALSE)
		k = as.numeric(df[,BEST.K.FOR.SILH.COL.NAME])
		best.silh = df[,BEST.SILH.COL.NAME]
			silh.folder.name = paste0("k=",k,"_silh=", sprintf("%.4f",best.silh))
		clu.analysis.mbrshp = as.numeric(unlist(strsplit(df[,BEST.MEM.FOR.SILH.COL.NAME], ",")))
		analysis.silh.folder = file.path(clu.analysis.measure.folder,silh.folder.name)
		characterization.silh.folder = file.path(clu.characterization.measure.folder,silh.folder.name)
		
		determine.core.parts(g, analysis.silh.folder, characterization.silh.folder, clu.analysis.mbrshp, k, dist.mtrx, part.folder, force)
						
		determine.general.core.part(g, analysis.silh.folder, characterization.silh.folder, clu.analysis.mbrshp, part.folder, force)
	}

}




#################################################################
# We have already computed the dissimilarity score between pairs of partitions.
# In this function, we choose a subset of the whole dissimilarity matrix for each solution class of kmedoids results.
# Then, we compute the mean scores, i.e. a mean score for each partition (over the other partitions).
# Finally, we choose the partition with lowest mean score ==> the partition which is the most similar (i.e. less distance) to the others.
# If there are many candidates fopr the representative partition, i.e. many max scores, then, we choose the first partition.
# In the end, we plot the representative partitions with the corresponding graph.
#
# analysis.silh.folder
# characterization.silh.folder
# clu.analysis.mbrshp
# k: a number of cluster in the kmedoids clustering, i.e. number of solution classes. Usually, the value associated with best silhouette is chosen.
# dist.mtrx: dissimilarity matrix
# net.folder: folder where input graphs are stored
# force: whether or not the existing files are overwritten by a fresh call of all corresponding methods (e.g partitioning method)
#
#################################################################
determine.representative.partitions = function(analysis.silh.folder, characterization.silh.folder, clu.analysis.mbrshp, k, dist.mtrx,
											net.folder, force)
{
   
    characterization.repr.folder = file.path(characterization.silh.folder, REPRESENTATIVE.FOLDER.NAME)
    repr.mbrshp.files = list.files(path=characterization.repr.folder,pattern=paste0("^",MBRSHP.FILE.PREFIX,".*\\.txt$"))

    if(length(repr.mbrshp.files)==0 || force){
        
        # for each cluster (of similar solutions) of k-medoids clustering, i.e. solution classes
    	for(clu.no in seq(1,k)){
    		analysis.clu.folder = file.path(analysis.silh.folder, paste0(CLUSTER.ID.FOLDER.PREFIX,clu.no))
    #		if(!dir.exists(analysis.clu.folder))
    #			dir.create(path=analysis.clu.folder, showWarnings=FALSE, recursive=TRUE)
		
	
		    part.ids = which(clu.analysis.mbrshp == clu.no)
		    m = length(part.ids) # nb current solution
		    curr.dist.mtrx = as.matrix(dist.mtrx[part.ids, part.ids])
		    
    		if(m>1){
    			# ----------------------------------------------------------
    			# mean.sd.result => for each solution, we take an average of distance scores obtained when compared with other solutions
    			# curr.dist.mtrx is a subset of the original dist.mtrx => so need to compute again mean.sd.result
    			mean.sd.result = compute.mean.and.sd.dist.scores(curr.dist.mtrx, is.mtrx.symmetric=TRUE, m, by.row=TRUE, by.col=FALSE)
    			colnames(mean.sd.result) = c(MEAN.COL.NAME, SD.COL.NAME)
    			rownames(mean.sd.result) = paste("sol", seq(0, m-1))
    			table.file = file.path(analysis.clu.folder, paste0(EVAL.MEAN.SD.DIST.SCORES.FILE.PREFIX,".csv"))
    			write.csv(x=mean.sd.result, file=table.file, row.names=TRUE)
    			# -----------------------------------------------------------
    			min.indx = which.min(mean.sd.result[,MEAN.COL.NAME]) # get 1st min occurence
    			repr.part.id = part.ids[min.indx] - 1 # minus 1 since ids start from 0
    	#		print(repr.part.id)
    		} else { # if(m == 1)
    			repr.part.id = part.ids[1]-1 # the unique solution id => # minus 1 since ids start from 0
    		}
    		
    		if(!dir.exists(characterization.repr.folder))
    			dir.create(path=characterization.repr.folder, showWarnings=FALSE, recursive=TRUE)
    		# copy the membership file
    		from.file=file.path(analysis.clu.folder,paste0(MBRSHP.FILE.PREFIX,repr.part.id,".txt"))
    		to.file=file.path(characterization.repr.folder,paste0(MBRSHP.FILE.PREFIX,"-",CLUSTER.ID.FOLDER.PREFIX,clu.no,".txt"))
    		file.copy(from=from.file, to=to.file, overwrite=TRUE)
    		
    		
        	# PLOT
        	repr.mbrshp = read.table(file=to.file)$V1
        	graph.name = paste0(SIGNED.UNWEIGHTED.FILE, ".G") # our reference graph to get optimal and worst imbalance count value
        	network.path = file.path(net.folder, graph.name)
        	g = read.graph.ils(network.path)
        	plot.file = file.path(characterization.repr.folder, paste0(CLUSTER.ID.FOLDER.PREFIX,clu.no))
        	plot.network(g, membership=repr.mbrshp, plot.file, format="JPEG", method="circle")
  
	    }
    }      
	
    # update membership files, i.e. create relative plot mebership files, except the first membership file. And plot them
    nb.mbrshp.file = k
    mbrshp.files = paste0(MBRSHP.FILE.PREFIX,"-",CLUSTER.ID.FOLDER.PREFIX, seq(1,k))
    plot.mbrshp.file = list.files(path=characterization.repr.folder,pattern=paste0("^",PLOT.PREFIX,".*\\.txt$"))
    
    if(length(plot.mbrshp.file)==0 || force){
    	# ------------------------------------------------------------------------------------------------------------------
    	

    	# ----------------------------------------------------------------------
    	# update membershps for plot, write them into different files
    	# ----------------
    	# first, copy files, then update them except the 1st one
    	for(i in seq(1,nb.mbrshp.file)){
    		filename = mbrshp.files[i] # first membership
    		file.copy(from=file.path(characterization.repr.folder,paste0(filename,".txt")),
    			to=file.path(characterization.repr.folder,paste0(PLOT.PREFIX,"-",filename,".txt")), overwrite=TRUE)
    	}
    	# ----------------
    	# note that the first membership remains unchanged for plot
    	if(nb.mbrshp.file > 1){
    		for(i in seq(1,nb.mbrshp.file-1)){ # each time, change only mem2
    			filename1 = mbrshp.files[i]
    			mem1 = read.table(file=file.path(characterization.repr.folder,paste0(PLOT.PREFIX,"-",filename1,".txt")))$V1
    			filename2 = mbrshp.files[i+1]
    			mem2 = read.table(file=file.path(characterization.repr.folder,paste0(PLOT.PREFIX,"-",filename2,".txt")))$V1
    			# create relative membership
    			new.mem2 = create.relative.plot.membership(mem1, mem2)
    			write.table(file=file.path(characterization.repr.folder,paste0(PLOT.PREFIX,"-",filename2,".txt")), x=new.mem2, col.names=F, row.names=F)
    		}
    	}
    	# ----------------------------------------------------------------------
    }
     
     
	# ----------------
	# plot
	graph.name = paste0(SIGNED.UNWEIGHTED.FILE, ".G") # our reference graph to get optimal and worst imbalance count value
	network.path = file.path(net.folder, graph.name)
	g = read.graph.ils(network.path)
	for(i in seq(1,nb.mbrshp.file)){
		filename = mbrshp.files[i]
		plot.file.path = file.path(characterization.repr.folder,filename)
		if(!file.exists(plot.file.path) || force){
    		mem = read.table(file=file.path(characterization.repr.folder,paste0(PLOT.PREFIX,"-",filename,".txt")))$V1
    		plot.network(g, membership=mem, plot.file=plot.file.path, format="JPEG", method="circular", plot.title="")
		}
	}
	# ----------------
	# ------------------------------------------------------------------------------------------------------------------

}




###############################################################
# It is a wrapper function for finding the representative partitions of each solution class in kmedoids clustering result.
# Instead of processing each possible number of cluster values in kmedoids, it chooses the k value associated with the best silhouette score.
#
#
# clu.characterization.measure.folder
# eval.folder
# clu.analysis.measure.folder
# dist.mtrx: dissimilarity matrix
# net.folder: folder where input graphs are stored
# measure: dissimilarity measure which is used in k-mdedoids
# force: whether or not the existing files are overwritten by a fresh call of all corresponding methods (e.g partitioning method)
#
##################################################################
determine.all.representative.partitions = function(clu.characterization.measure.folder, eval.folder, clu.analysis.measure.folder,
				dist.mtrx, net.folder, measure, force)
{
	table.file = file.path(eval.folder, paste0(EVAL.BEST.K.FOR.KMEDOIDS,"-",measure,".csv"))
	if(file.exists(table.file)){
		df = read.csv(table.file, check.names=FALSE,stringsAsFactors=FALSE)
		k = as.numeric(df[,BEST.K.FOR.SILH.COL.NAME])
		best.silh = df[,BEST.SILH.COL.NAME]
		silh.folder.name = paste0("k=",k,"_silh=", sprintf("%.4f",best.silh))
		clu.analysis.mbrshp = as.numeric(unlist(strsplit(df[,BEST.MEM.FOR.SILH.COL.NAME], ",")))
		analysis.silh.folder = file.path(clu.analysis.measure.folder,silh.folder.name)
		characterization.silh.folder = file.path(clu.characterization.measure.folder,silh.folder.name)
		
		determine.representative.partitions(analysis.silh.folder, characterization.silh.folder, clu.analysis.mbrshp, k, dist.mtrx,
											net.folder, force)
	}
}


#################################################################
#
# It characterizes the partitions of a given network based on the considered algorithm name and graph type (weighted or not, etc.).
# It is a wrapper function.
#
# n: graph size
# l0: number of clusters
# d: density
# prop.mispl: proportion of misplaced links
# prop.neg: proportion of negative links
# network.no: network id (the identifiers start from 1)
# cor.clu.exact.algo: the name of correlation clustering algorithm that has been executed for optimal solutions
# comp.measures: dissimilarity measure to compare the optimal solutions
# k.limits: the lower and upper bound limits regarding the the number of cluster values in kmedoids
# force: whether or not the existing files are overwritten by a fresh call of all corresponding methods (e.g partitioning method)
#
##################################################################
perform.cluster.characterization = function(n, l0, d, prop.mispl, prop.neg, network.no,
	 cor.clu.exact.algo, comp.measures, k.limits, force)
{

	#tlog(16, "start to perform cluster characterization with exact algorithms")	
	#for(algo.name in cor.clu.exact.algos){
		tlog(20, "performing cluster characterization => algo.name: ", cor.clu.exact.algo)
		
		for(graph.desc.name in c(SIGNED.UNWEIGHTED.FILE)){ # OPPOSITE.SIGNED.UNWEIGHTED.FILE
			tlog(24, "performing cluster characterization => graph.desc.name: ", graph.desc.name)
			
			net.folder = get.input.network.folder.path(n, l0, d, prop.mispl, prop.neg, network.no)
			if(dir.exists(net.folder)){
				graph.name = paste0(SIGNED.UNWEIGHTED.FILE, ".G") # our reference graph to get optimal and worst imbalance count value
				network.path = file.path(net.folder, graph.name)
				g = read.graph.ils(network.path)
		
				net.folder = get.input.network.folder.path(n, l0, d, prop.mispl, prop.neg, network.no)
				part.folder = get.part.folder.path(n, l0, d, prop.mispl, prop.neg, network.no, cor.clu.exact.algo, graph.desc.name)
				eval.algo.folder = get.eval.folder.path(n, l0, d, prop.mispl, detected.imb=NA, prop.neg, k=ALL, network.no, cor.clu.exact.algo, graph.desc.name)
				clu.analysis.folder = get.clu.analysis.folder.path(n, l0, d, prop.mispl, prop.neg, network.no, cor.clu.exact.algo, graph.desc.name)
				clu.characterization.folder = get.clu.characterization.folder.path(n, l0, d, prop.mispl, prop.neg, network.no, cor.clu.exact.algo, graph.desc.name)
				
				partitions = load.membership.files(part.folder)
				m = length(partitions) # nb partition
				
				for(measure in comp.measures){
					tlog(20, "performing cluster characterization => measure: ", measure)
			
					clu.analysis.measure.folder = file.path(clu.analysis.folder,measure)
					clu.characterization.measure.folder = file.path(clu.characterization.folder,measure)
					if(!dir.exists(clu.characterization.measure.folder))
						dir.create(path=clu.characterization.measure.folder, showWarnings=FALSE, recursive=TRUE)
						
					
					# ---------------------------------------------------
					matrix.file = file.path(eval.algo.folder,paste0(EVAL.DIST.MATRIX.FILE.PREFIX,"-",measure,".csv"))
					dist.mtrx = as.matrix(read.csv(matrix.file, row.names = 1, header= TRUE, check.names=FALSE)) # square matrix
					dist.mtrx = matrix(as.numeric(sprintf("%.4f", dist.mtrx)), m, m)
					# ---------------------------------------------------			
					
					#unlink(clu.characterization.measure.folder, recursive=TRUE)
					determine.all.representative.partitions(clu.characterization.measure.folder, eval.algo.folder, clu.analysis.measure.folder,
						 dist.mtrx, net.folder, force)
									
					determine.all.core.parts(g, clu.characterization.measure.folder, eval.algo.folder, clu.analysis.measure.folder,
						dist.mtrx, part.folder, force)
				}
			}
		
		}
	#}
	
}
	




#################################################################
#
# It is the starting method in the aim of characterizing the partitions of the considered networks thanks to kmedoids analsis (i.e. with solution classes). 
#   It handles all networks by graph.sizes, prop.mispls, my.prop.negs and in.rand.net.folders
#
# graph.sizes: a vector of values regarding graph sizes to be considered
# d: density (it is a single value)
# l0: number of clusters to be considered (it is a single value)
# prop.mispls: a vector of values regarding proportion of misplaced links
# prop.negs: a vector of values regarding proportion of negative links (for now, it is not operational)
# in.rand.net.folders: a vector of values regarding input random graph folders. Sequantial integers (1, .., 10)
# cor.clu.exact.algo: the name of correlation clustering algorithm that has been executed for optimal solutions
# comp.measures: dissimilarity measure to compare the optimal solutions
# k.limits: the lower and upper bound limits regarding the the number of cluster values in kmedoids
# force: whether or not the existing files are overwritten by a fresh call of all corresponding methods (e.g partitioning method)
#
##################################################################
perform.all.cluster.characterization = function(graph.sizes, d, l0, prop.mispls, prop.negs, in.rand.net.folders,
	cor.clu.exact.algo, comp.measures, k.limits, force)
{
	tlog("starts performing cluster characterization")
	for(n in graph.sizes){
		tlog(8, "performing cluster characterization => n: ", n)
		
		for(prop.mispl in prop.mispls){
			tlog(8, "performing cluster characterization => prop.mispl: ", prop.mispl)
			
		    if(is.na(prop.negs) && d == 1){
		        prop.negs = compute.prop.neg(n, d, k, prop.mispl)
		    }
			
			for(prop.neg in prop.negs){
				tlog(12, "performing cluster characterization => prop.neg: ", prop.neg)
				
				for(network.no in in.rand.net.folders){
					tlog(16, "performing cluster characterization => network.no: ", network.no)
						
					perform.cluster.characterization(n, l0, d, prop.mispl, prop.neg, network.no,
					            cor.clu.exact.algo, comp.measures, k.limits, force)
	
				}
				
			}
			
		}
		
	}
	
}
