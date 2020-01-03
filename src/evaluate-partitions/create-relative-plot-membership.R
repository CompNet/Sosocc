

######################################################################
# It basically binds mem1 and mem2 and outputs a data frame. This data frame will show which cluster in 'mem1' is binded to which cluster in 'mem2'.
# The goal is to create a membership vector relatively to mem1. Thus, when we plot 2 membership based on the same graph, we can catch the small differences.
# It is sometimes difficult to catch those little differences when node coloring is not handled well since node coloring  depends only on how membership vector is created
# Even though 2 membership vectors are very similar, they might be colored very differetnly because of the node coloring/cluster numbering.
# To tackle this issue, we are able adjust the mem2 to see color changing for small membership changes.
# An example:
# mem1: 1,1,2,2,2,2
# mem2: 2,2,1,1,1,1
# Based on this example, nodes 1 and 2 will be colored differently in 2 graphs. But if 'mem2' could be in this way:
# mem2: 1,1,2,2,2,2 => hence, they will be colored in the same way.
#
# This method performs the following things:
# 1) use 'purity' as similarity measure. We use 'purity' because it is asymetric and it captures difference and union of 2 set
#    
# 2) apply the purity measure at cluster level
# An example of output:
#       k   nb.item.in.clu purity         associated cluster id
# [1,] "1" "9"            "0.8889"        "1"                  
# [2,] "2" "6"            "1"             "2"                  
# [3,] "3" "2"            "0.5"           "3,5"                
# [4,] "4" "4"            "1"             "4"                  
# [5,] "5" "3"            "1"             "3" 
# We read line 1 in this way: cluster1 in 'mem1', which has 9 elements, is associated to the cluster1 in 'mem2' with 0.889 sim score
# In other words, 8 of 9 elements of 'mem1' is included in 'mem2'
#
#
# 2) Then, we observe that some purity socres are 0.5, and associated cluster id is multiple. To tackle this, we duplicate each row with multpile cluster id:
#   k   nb.item.in.clu purity         associated cluster id
#  "3" "2"            "0.5"           "3,5"   
#  becomes:
#  "3" "2"            "0.5"           "3"  
#  "3" "2"            "0.5"           "5"
#
# Then, remove non-unique associated cluster ids:
# For instance:
#   k   nb.item.in.clu purity         associated cluster id
#  "3" "2"            "0.5"           "3" 
#  "5" "3"            "1"             "3" 
# Here, 3, which is associated cluster id, is non-unique, because there are 2 rows which have 3.
# We need to select one of them: To us, what is most important is nb item in cluster, then purity score.
# in this case, the best row is: 
#  5" "3"            "1"             "3" 
# Even though nb item in cluster is 2 in this case:
#  5" "2"            "1"             "3" 
# this row is still the best. Becase the purity is 1, which higher than the other value (0.5).
# 
# At the end, we obtain:
# index k nb.item.in.clu    purity associated cluster id
# 1     1 1              9 0.8888889                     1
# 2     2 2              6 1.0000000                     2
# 6     4 3              2 0.5000000                     5
# 3     5 4              4 1.0000000                     4
# 4     6 5              3 1.0000000                     3
#
# If needed, for example for node coloring, we need to fill in missing part of the new mem2. If TRUE, do it
#
# Thanks to this data frame, we can easily construct a new 'mem2' vector
#####################################################################
bind.two.membership = function(mem1, mem2, fill.missing.part=TRUE){

	res = compute.purity(partition=mem1, ground.truth=mem2, purity.level="cluster")
  
	# since membership might not be a sequential unique values, purity score might be NaN ==> [1, 2, 2, 1, 5] without 3 and 4
	if(any(is.nan(as.numeric(res[,"purity"])))){
		nan.indx = which(is.nan(as.numeric(res[,"purity"])))
		res = res[-nan.indx,]
#		print("after nan")
#		print(res)
	}
#	print(res)
	
	
	# -------------------------------------------------------
	# remove miltiple cluster ids, make them single cluster id
	mult.corresp.indx = which(grepl(",", res[,"associated cluster id"])) # "associated cluster id"
	upd.res = res
	for(indx in mult.corresp.indx){
		ids = unlist(strsplit(res[indx,"associated cluster id"], split=","))
		duplicated.rows = res[rep(indx, length(ids)),] # duplicate the row by length(ids) times
		duplicated.rows[,"associated cluster id"] = ids
	    upd.res = rbind(upd.res, duplicated.rows)
	}
	# remove the rows containing multiple cluster ids at the end
	mult.corresp.indx = which(grepl(",", upd.res[,"associated cluster id"])) # "associated cluster id"
	if(length(mult.corresp.indx)>0)
		upd.res = upd.res[-mult.corresp.indx,]
	# -------------------------------------------------------
	  
	  
	# -------------------------------------------------------
	# order the data frame
	res2 = data.frame(matrix(as.numeric(upd.res), nrow=nrow(upd.res)))
	colnames(res2) = colnames(res)
	# order the data frame by 'k' in incr order, by cluster size and purity in decr order
	ordered.res = res2[order(-res2[,"k"], res2[,"nb.item.in.clu"], res2[,"purity"], decreasing=TRUE),]
	# add line index as 1st column
	ordered.res = cbind(seq(1,nrow(ordered.res)), ordered.res) 
	colnames(ordered.res) = c("index",colnames(res))
#	print(ordered.res)
#	print("--")
	# -------------------------------------------------------
	
	  
	# -------------------------------------------------------  
	# remove duplicated cluster ids
	freq = table(ordered.res[,"associated cluster id"])
#	print(freq)
	clu.ids = names(freq)[which(freq > 1)]
	for(id in clu.ids){
#		print("for")
#		print(ordered.res)
		indx = which(ordered.res[,"associated cluster id"] == id)
		
		if(length(indx) > 1){ # check if it is still duplicated
#			print("indx")
#			print(indx)
			initial.indx = ordered.res[indx,"index"]
#			print(initial.indx)
			corr.rows = ordered.res[indx,]
			# sort the rows: cluster size is most important criteria, then purity score
			ordered.rows = corr.rows[order(corr.rows[,"nb.item.in.clu"], corr.rows[,"purity"], decreasing=TRUE),]
#			print(ordered.rows)
			# keep only the 1st one in the data frame => considered as 'best'
			best.indx = ordered.rows[1,"index"]
#			print(best.indx)
			excl.indx = indx[-which(initial.indx == best.indx)]
#			print(excl.indx)
			ordered.res = ordered.res[-excl.indx,]
#			print(ordered.res)
			# check if duplicated k' value for best index
#			print("ordered.rows")
#			print(ordered.rows)
			k = ordered.rows[1,"k"] # the best one
			indx = which(ordered.res[,"k"] == k)
#			print("k index")
#			print(indx)
			if(length(indx) > 1){ # if there is other row other than the best one with the same 'k'
		  		excl.indx = indx[-which(indx == best.indx)]
#				print("excl indx")
#				print(excl.indx)
		    	ordered.res = ordered.res[-excl.indx,]
		    }
		}
	}
#	print(ordered.res)
#	print("--")
	# -------------------------------------------------------
	
	  
	# order again for the next operation
	ordered.res = ordered.res[order(-ordered.res[,"k"], ordered.res[,"nb.item.in.clu"], ordered.res[,"purity"], decreasing=TRUE),]
	  
	# -------------------------------------------------------
	# handle duplicated 'k' values if there is not any duplicated cluster ids
	freq = table(ordered.res[,"k"])
	k.vals = names(freq)[which(freq > 1)]
	for(k in k.vals){
		indx = which(ordered.res[,"k"] == k)
		# print(indx)
		if(length(indx) > 1){ # if there is other row other than the best one with the same 'k'
	    	# as it is ordered already, the 1st on is the best
	    	excl.indx = indx[seq(2,length(indx))]
	    	# print(excl.indx)
		    ordered.res = ordered.res[-excl.indx,]
		}
	}
#	print(ordered.res)
#	print("--")
	# -------------------------------------------------------
	 
	
	# for node coloring, this part is needed
	# but, if we want to just bind 2 memberships, we might not need to do this part
	#   => for example: if cluster1 in mem1 is split into 2 unequal-sized clusters in mem2, we bind the bigger part with cluster1
	#   => particularly, we do not need this part when catching different membership files between 2 clusters of 2 different measures
	if(fill.missing.part){
		# -------------------------------------------------------
		# filling the missing part:
		#  => every cluster id must be associated with some value
		# => since we know that values in 'unique(mem2)' should be appear in mem2, calculate missing cluster ids
		missing.clu.ids = setdiff(unique(mem2), ordered.res[,"associated cluster id"])

#		if(length(missing.clu.ids)>0 && length(missing.k.vals)>0){
		if(length(missing.clu.ids)>0){
			# note that ordered.res[,"k"] may contain different clu id than unique(mem2) => due to our create relative membership
			# so it is possible that length(missing.k.vals) != length(missing.k.vals).
			# It is not a bug, it is ok. => juste take as many as you need
			missing.k.vals = setdiff(unique(mem2), ordered.res[,"k"])
			missing.k.vals = missing.k.vals[1:length(missing.clu.ids)]
			
			missing.part = cbind(missing.k.vals,missing.clu.ids)
			# print(missing.part)
			colnames(missing.part) = c("k","associated cluster id")
#			 print(missing.part)
			ordered.res = rbind(ordered.res[,c("k","associated cluster id")], missing.part)
#			 print(ordered.res)
		}
		# -------------------------------------------------------
	}
	
#	print("ordered")
#	print(ordered.res)
#	print("--")
	return(ordered.res)
}








######################################################################
# The problem is the following: it is difficult to distinguish visually if two partitions are similar or when we want to plot them.
#  This is because the membership vectors are created in different times, so they are not supposed to be labelled in the same way.
#  In other words, a cluster labeled as '1' in the first membership (mem1) may be labeled as '2' in the second membership (mem2).
# So, this function creates a new membership vector which will replace 'mem2', because the new vector is built based on 'mem1'.
# So, 'mem2' is supposed to have the most similar form compared to 'mem1'.
# This is performed thanks to the data frame which show the binding inforamation between mem1 and mem2
# At the end, construct the new relative membership based on "k" and "associated cluster id" columns.
# In other words, when mem2 contains 1, put 1. When it contains 6, put 4, and so on.
#####################################################################
create.relative.plot.membership = function(mem1, mem2){
#   print(mem1)
#   print(mem2)
  ordered.res = bind.two.membership(mem1, mem2)

  # print(ordered.res)
  
  # main part
  new.mem2 = rep(NA,length(mem2))
  for(i in seq(1,nrow(ordered.res))){
    k = ordered.res[i,"k"]
    clu.id = ordered.res[i,"associated cluster id"]
    indx = which(mem2 == clu.id)
    new.mem2[indx] = k
  }
  
#  print("mem1")
#  print(mem1)
#  print("old mem2")
#  print(mem2)
#  print("new mem2 in func")
#  print(new.mem2)
  
  return(new.mem2)
}
