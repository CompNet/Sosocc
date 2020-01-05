
library(cluster)


#################################################################
# It finds the centroid partition, i.e. the partition which is most similar to others.
# When finding global centroid, no need to specify 'partition' and 'k'.
#
# when there are multiple candidate for max, it returns the index of the first max value 
#
# dist.mtrx: distance/dissimilarity matrix between all partitions
# partition: the kmedoids clustering result for a specific cluster number value (e.g. k=3)
# k: the specific cluster number value (e.g. k=3) in kmedoids
#
#################################################################
find.centroid = function(dist.mtrx, partition=NA, k=NA){
	centroid.id = NA
	n = nrow(dist.mtrx) # 'n' solutions
	
	if(is.na(k)){ # find the centroid of all points, i.e solutions
		mean.dist.vals = rep(NA, n)
		for(i in seq(1,n)) # note that the distance of a clustering with itself is 0
			mean.dist.vals[i] = mean(dist.mtrx[i,]) 
		centroid.id = which.max(mean.dist.vals)
	}
	else { # find the centroid of a specific cluster
		#print(dist.mtrx)
		indxs = which(partition == k)
		sub.dist.mtrx = dist.mtrx[indxs,indxs]
		if(length(indxs) == 1)
			sub.dist.mtrx = matrix(sub.dist.mtrx, nrow=1, ncol=1)
		
		nsub = length(indxs)
		mean.dist.vals = rep(NA, nsub)
		for(i in seq(1,nsub))
			mean.dist.vals[i] = mean(sub.dist.mtrx[i,])
		max.indx = which.max(mean.dist.vals)
		centroid.id = indxs[max.indx]
	}
	
	return(centroid.id)
}


#################################################################
# evaluating between cluster distances
# see https://infoscience.epfl.ch/record/134302/files/Saitta_et_al_Postprint-2007-IDA-Clust-Index.pdf?version=1
#
# partition: the kmedoids clustering result for a specific cluster number value (e.g. k=3)
# dist.mtrx: distance/dissimilarity matrix between all partitions
# global.centroid.id: the index of partition which is most similar to other partitions (at global level, not at cluster level)
#
#################################################################
compute.bcd = function(partition, dist.mtrx, global.centroid.id){

	n = length(partition)
	k = length(unique(partition))
	
	bcd = 0
#	tot.centroid.id = find.centroid(dist.mtrx)
	for(i in 1:k){
		clu.size = length(which(partition == i))
		clu.centroid.id = find.centroid(dist.mtrx, partition, k=i)
		d = dist.mtrx[global.centroid.id,clu.centroid.id]
		bcd = bcd + (d^2)*clu.size
	}
	bcd = bcd/(n*k)
	return(bcd)
}

#################################################################
# evaluating within cluster distances.
# see https://infoscience.epfl.ch/record/134302/files/Saitta_et_al_Postprint-2007-IDA-Clust-Index.pdf?version=1
#
# partition: the kmedoids clustering result for a specific cluster number value (e.g. k=3)
# dist.mtrx: distance/dissimilarity matrix between all partitions
#
#################################################################
compute.wcd = function(partition, dist.mtrx){
	
	n = length(partition)
	k = length(unique(partition))
	
	wcd = 0
	for(i in 1:k){
		clu.indxs = which(partition == i)
		clu.size = length(clu.indxs)
		clu.centroid.id = find.centroid(dist.mtrx, partition, k=i)
		
		# process for each point of the considered cluster
		sum.d = 0
		for(indx in clu.indxs){
			d = dist.mtrx[indx,clu.centroid.id]
			sum.d = sum.d + d^2
		}		
		wcd = wcd + sqrt(sum.d/clu.size)
	}
	wcd = wcd/k
	return(wcd)
}




#################################################################
# Computes score function (SF), described here: 
    # https://infoscience.epfl.ch/record/134302/files/Saitta_et_al_Postprint-2007-IDA-Clust-Index.pdf?version=1
#
# partition: the kmedoids clustering result for a specific cluster number value (e.g. k=3)
# dist.mtrx: distance/dissimilarity matrix between all partitions
# global.centroid.id: the index of partition which is most similar to other partitions (at global level, not at cluster level)
#
#################################################################
compute.SF = function(partition, dist.mtrx, global.centroid.id){
	bcd = compute.bcd(partition, dist.mtrx, global.centroid.id)
	wcd = compute.wcd(partition, dist.mtrx)
	res = 1 - 1/(exp(exp(bcd-wcd)))
	return(res)
}


#################################################################
# It computes score function (SF) score for each cluster value of kmedoids analysis (e.g. k=2, etc.).
#
# k.values: considered cluster values (e.g. c(2,3,4,5) ..) in k-medoids analysis
# partitions: all the kmedoids clustering results for a range of cluster number values (i.e. k.values)
# dist.mtrx: distance/dissimilarity matrix between all partitions
#
#################################################################
compute.SF.scores = function(k.values, partitions, dist.mtrx){
	result = list()
	result[["k.values"]] = k.values
	result[["partitions"]] = partitions
	
	global.centroid.id = find.centroid(dist.mtrx)
	
	
	i = 0
	scores = c()
	for(partition in partitions){
		i=i+1
		score = compute.SF(partition, dist.mtrx, global.centroid.id)
		scores = c(scores, score)
	}
	
	result[["SF"]] = scores
	return(result)
}

# =======================================================
# =======================================================











#################################################################
# It prepares a list of vectors in which we store the kmedoids results for the considered cluster values 'k.values'
#
# An example of output:
# partitions[["k=1"]] = {1 1 1 1 1 1} 
# partitions[["k=2"]] = {1 2 1 2 1 1}
#
# clu.analysis.measure.folder: folder containing the results of the kmedoids analysis
# folder.names: silhouette folder names (e.g. "silh=0.41-k=2")
# k.values: considered cluster values (e.g. c(2,3,4,5) ..) in k-mdeoids analysis
# 
#################################################################
retreive.clu.result.partitions = function(clu.analysis.measure.folder, folder.names, k.values){
	
	partitions = list()
	for(i in 1:length(k.values)){
		mem.file.path = file.path(clu.analysis.measure.folder, folder.names[i], paste0(MBRSHP.FILE.PREFIX,".txt"))
        # print(mem.file.path)
		k.mem = read.table(mem.file.path)$V1
		k.desc = as.character(k.values[i])
		partitions[[k.desc]] = k.mem
	}
	
	return(partitions)
}














#################################################################
# It shows the difference between partitions via alluvial diagram.
# At first glance, it is a good idea. But, when there are many partitions, it is not easy to read the diagram.
# That is why we propose to create the diagram for some predefined max number of partitions => we call it 'chunk.size'.
# For instance, if there are 20 partitions and if the chunk size is 10, then two sperate diagrams will be created. 
#
# folder: output folder
# mbrshps: partitions
# chunk.size: the max number of partitions that will be considered in a single alluvial diagram
# x.labels: labels associated with the partitions
#
#################################################################
plot.partitions.in.alluvial.diagram = function(folder,mbrshps,chunk.size=10,x.labels){

	nb.mbrshp.file = length(mbrshps)
	if(nb.mbrshp.file == 1) # no need to use alluvial diagram
		return()
	
	counter = 0
	part.counter = 0
	while(counter < nb.mbrshp.file){
		part.counter = part.counter+1
#		print("part")
#		print(part.counter)
		start.indx = counter+1
		end.indx = counter+chunk.size
		if(nb.mbrshp.file < end.indx)
			end.indx = nb.mbrshp.file
		if(end.indx+1 == nb.mbrshp.file){ # do not let have just one file at the end for the next part
			end.indx = nb.mbrshp.file
			counter=counter+1
		}
		
		partitions.df = c()
		for(i in seq(start.indx,end.indx)){
			partitions.df = cbind(partitions.df, mbrshps[[i]])	
		}
		# print(partitions.df)
		
		plot.filename = paste0(ALLUVIAL.DIAGRAM.PLOT.TYPE,"-part",part.counter)
		xlab = x.labels[seq(start.indx,end.indx)]
		plot.alluvial.diagram(folder, plot.filename, partitions.df, xlab, plot.title="", format="JPEG")
			
		counter = counter+chunk.size
	}
	
}
		




#################################################################
# It allows to visualize the partitions locating in the same kmedoids cluster with different graph layouts.
# To ease the perception of the similarity in a sequence of partitions, we update the partition labels in order that the current partition
#   be colored in a similar way compared to the previous partition
#   For instance: c(1,2,2,2,2) may become c(2,1,1,1,1,1).
# Concertely, it first copies the membership files, then updates them (i.e. it changes cluster labels) for visualiztion purposes.
# Finally, it plots the partitions with:
#   1) it plots with different graph layouts
#   2) it plots with alluvial diagram
#
# g: a graph object
# clu.analysis.measure.folder: folder in which the results will be created 
# silh.foldername:
# partitions: all partitions
# clu.analysis.partition: kmedoids clustering result
#
#################################################################
copy.and.plot.corresp.partitions.by.clu.no = function(g, clu.analysis.measure.folder, silh.foldername, partitions,
		clu.analysis.partition)
{
	k = length(unique(clu.analysis.partition)) # number of solution classes
	for(clu in 1:k){ # for each solution class
		indx = which(clu.analysis.partition == clu)
		clu.folder = file.path(clu.analysis.measure.folder, silh.foldername,
				paste0(CLUSTER.ID.FOLDER.PREFIX,clu))
		if(!dir.exists(clu.folder))
			dir.create(clu.folder)

		for(i in indx){
		    # write only memberships into files, do not plot them
			mbrshp = partitions[[i]]
			# since membership files starts from 0, use 'i-1'
			filename = paste0(MBRSHP.FILE.PREFIX,(i-1))
			table.file = file.path(clu.folder, paste0(filename,".txt"))
			write.table(x=mbrshp, file=table.file, row.names=FALSE, col.names=FALSE)

			# #plot
			# plot.file = file.path(clu.folder, filename)
			# plot.network(g, membership=mbrshp, plot.file, format="JPEG", method="circle")
		}

		
		# --------------------------------------------------------------------------------------------------------------
		# update membership files, i.e. create relative plot mebership files, except the first membership file. And plot them
		mbrshp.files = list.files(path=clu.folder,pattern=paste0("^",MBRSHP.FILE.PREFIX,".*\\.txt$"))
		mbrshp.files = gsub(".txt","",mbrshp.files) # remove '.txt'
		nb.mbrshp.file = length(mbrshp.files)
		#print(mbrshp.files)
		
		
		###################################################################
		# 1) update memberships for plot and write them into another files
		###################################################################
		
		# ==================================================================
		# update memberships for plot, write them into different files
		# ----------------
		# first, copy files, then update them except the 1st one
		for(i in seq(1,nb.mbrshp.file)){
			filename = mbrshp.files[i] # first membership
			file.copy(from=file.path(clu.folder,paste0(filename,".txt")),
				 to=file.path(clu.folder,paste0(PLOT.PREFIX,"-",filename,".txt")), overwrite=TRUE)
		}
		# ----------------
		# note that the first membership remains unchanged for plot
		if(nb.mbrshp.file > 1){
			for(i in seq(1,nb.mbrshp.file-1)){ # each time, change only mem2
				filename1 = mbrshp.files[i]
				mem1 = read.table(file=file.path(clu.folder,paste0(PLOT.PREFIX,"-",filename1,".txt")))$V1
				filename2 = mbrshp.files[i+1]
				mem2 = read.table(file=file.path(clu.folder,paste0(PLOT.PREFIX,"-",filename2,".txt")))$V1
			  	# create relative membership
			  	
				new.mem2 = create.relative.plot.membership(mem1, mem2)
				write.table(file=file.path(clu.folder,paste0(PLOT.PREFIX,"-",filename2,".txt")), x=new.mem2, col.names=F, row.names=F)
			}
		}
		# ==================================================================
		

		###################################################################
		# 2) plot with different layouts
		###################################################################
		if(k > 1){
			# ----------------
			# plot graphs
		    # mbrshp.files contains the memberships associated to the cluster 'clu'
			plot.mbrshps = list()
			for(i in seq(1,nb.mbrshp.file)){
				filename = mbrshp.files[i]
				mem = read.table(file=file.path(clu.folder,paste0(PLOT.PREFIX,"-",filename,".txt")))$V1
				plot.mbrshps[[i]] = mem
				# circular
				plot.network(g, membership=mem, plot.file=file.path(clu.folder,filename), format="JPEG", method="circular", plot.title="")


				# unsigned/signed graph based layouts
				g.new = g
				for(layout.method in c(UNSIGNED.GRAPH.LAYOUTS, SIGNED.GRAPH.LAYOUTS)){
					if(!dir.exists(file.path(clu.folder,layout.method)))
						dir.create(path=file.path(clu.folder,layout.method), showWarnings=FALSE, recursive=FALSE)

					x = get.vertex.attribute(graph=g.new,name=paste0("x_",layout.method))
					y = get.vertex.attribute(graph=g.new,name=paste0("y_",layout.method))
					g.new <- set.vertex.attribute(graph=g.new, name="x",value=x)
					g.new <- set.vertex.attribute(graph=g.new, name="y",value=y)
					new.filename = paste0(filename,"-",layout.method)
					# do not specify the layout method name since we alrady integrated this info into the graph object
					plot.network(g.new, membership=mem, plot.file=file.path(clu.folder,layout.method,new.filename), format="JPEG", plot.title="")
				}


			}
			# ----------------
			# ------------------------------------------------------------------------------------------------------------------


			###################################################################
			# 3) plot alluvial diagrams for all partitions (chunk by chunk)
			###################################################################
			# =============================================================================

				plot.mbrshp.files = paste0(PLOT.PREFIX,"-",mbrshp.files)
				x.labels = mbrshp.files
				plot.partitions.in.alluvial.diagram(clu.folder,plot.mbrshps,chunk.size=20,x.labels)

			# =============================================================================

		}

	}
}




#################################################################
# It is a wrapper functon. It iterates over the range of values corresponding to the number of cluster values.
#   One way to do that is to iterate over the silhouette folder names (since each kmedoids clustering result is associated with a silhouette score)
# It first copies the membership files, then updates them (i.e. it changes cluster labels) for visualiztion purposes.
# Finally, it plots the partitions with:
#   1) it plots with different graph layouts
#   2) it plots with alluvial diagram
#
#
# g: a graph object
# clu.analysis.measure.folder: folder in which the results will be created
# silh.folder.names: folder names for each kmedoids clustering result => c("silh=0.42-k=2","silh=0.52-k=3", ..)
# partitions: all partitions
#
#################################################################
copy.and.plot.corresp.partitions = function(g, clu.analysis.measure.folder, silh.folder.names, partitions)
{
    counter = 0
    for(silh.foldername in silh.folder.names){
        counter=counter+1
        clu.analysis.measure.silh.folder =  file.path(clu.analysis.measure.folder, silh.foldername)
        table.file = file.path(clu.analysis.measure.silh.folder, paste0(MBRSHP.FILE.PREFIX,".txt"))
        clu.analysis.partition <- as.numeric(as.matrix(read.table(file=table.file, header=FALSE)))
        
        # apply some threshold, otherwise it takes time
        copy.and.plot.corresp.partitions.by.clu.no(g, clu.analysis.measure.folder, silh.foldername, partitions,
                                          clu.analysis.partition)
    }
}





#################################################################
# It creates kmedoids result folder for each specific number of clusters (i.e. k.values)
#
# An example of output: 
#   c("k=2_silh=0.46", "k=3_silh=0.52", ..)
#
# scores: score values obtained by an internal evaluation measure, such as silhouette
# k.values: a range of values corresponding to the number of clusters (i.e. number of solution classes) in kmedoids clustering analyis (e.g. k=2, k=3, k=4, ..)
# measure.desc: the description string of the dissimilarity measure
#
#################################################################
prepare.cluster.analysis.result.folders = function(scores, k.values, measure.desc){
	k.descs = paste0("k=", k.values)
	descs = paste0(measure.desc,"=", sprintf("%.4f",as.numeric(scores)))
	folder.names = paste(k.descs, descs, sep="_")
	return(folder.names)
}




#################################################################
# In kmeoids clustering, the method is run for  a range of number of cluster values (i.e. k.values).
#   For each kmedoids clustering result (so, for a specific number of clusters), this method allows to record
#    into a folder separately the partitions belonging to the same cluster of a kmedoids clsutering result.
#   For instance: if the is the kmedoids clustering result: c(1,2,1,2,1,1)
#   Then, first, third, fifth and 6.th partitions are recorded into the same folder, and the other partitions are recorded together into another folder.
#
#
# folder: output folder
# partitions: all partitions
# k.values: a range of values corresponding to the number of clusters (i.e. number of solution classes) in kmedoids clustering analyis (e.g. k=2, k=3, k=4, ..)
# avg.sil.list: a list of silhouette values
# measure.desc: the description string of the dissimilarity measure
#
#################################################################
write.cluster.analysis.memberships.into.file = function(folder, partitions, k.values, avg.sil.list, measure.desc){
	folder.names = prepare.cluster.analysis.result.folders(avg.sil.list, k.values, measure.desc)
	
	m = length(folder.names)
	for(k in 1:m){
		k.desc = as.character(k)
		folder.name = folder.names[k]
		res.folder = file.path(folder, folder.name)		
		if(!dir.exists(res.folder))
			dir.create(path=res.folder, showWarnings=FALSE, recursive=FALSE)
		
		table.file = file.path(res.folder, paste0(MBRSHP.FILE.PREFIX,".txt"))
		write.table(x=partitions[[k.desc]], file=table.file, row.names=FALSE, col.names=FALSE)
	}
	
	
}


#################################################################
# It first writes the silhouette scores into file. Then, it plots those values.
#
# silhouette.scores.file: file that will contain silhouette scores (and corresponding k values)
# avg.sil.list: a list of silhouette values
# k.vals: a range of values corresponding to the number of clusters (i.e. number of solution classes) in kmedoids clustering analyis (e.g. k=2, k=3, k=4, ..)
#
#################################################################
record.and.plot.silhouette.scores = function(silhouette.scores.file, avg.sil.list, k.vals){
	
	doc.clu.result.per.k = data.frame(k=k.vals, avg.sil.score=avg.sil.list)
	write.csv(file=paste(silhouette.scores.file,".csv",sep=""), x=doc.clu.result.per.k, row.names=FALSE)
	
	if(all(is.na(avg.sil.list))){
		# This is the case when there is only 1 or 2 optimal solution for a given instance
		# The plot can not be generated for the 1st and last cluster.
		# When there are 2 optimal solutions, this is the case.
		# DO NOTHING
	}
	else {
		pdf(file=paste(silhouette.scores.file,".pdf",sep=""),bg="white",compress=FALSE)
		plot(x=k.vals, y=avg.sil.list, xlab="k (nb cluster)", ylab="silhouette averages scores",
			main="Silhouette averages scores per k", ylim=c(min(min(as.numeric(avg.sil.list), na.rm=TRUE),0),1))
	#	points(x=sub.k.vals, y=sub.avg.sil.list, col="red")
		dev.off()
	}
	return(doc.clu.result.per.k)
}









#################################################################
# It applied the kmedoids clustering method for a specific number of clusters.
#   If user wants to apply another method other than kmdeoids, one needs to change this function.
#   In R, the kmedoids method is run by calling 'pam()' method.
#   For more information about kmedoids, see https://eric.univ-lyon2.fr/~ricco/cours/slides/classif_k_medoides.pdf
#
# method.name: the name of the method to be applied
# diss.mtrx: a dissimilarity matrix
# k: a specific number of clusters (i.e. number of solution classes)
#
#################################################################
apply.clustering.method = function(method.name, diss.mtrx, k){
	
	res = NA
	if(method.name == "kmedoids")
		res = pam(diss.mtrx, k, diss=TRUE) # we specify that we use a dissimilarity matrix
#	else
#		res = ANOTHER.METHOD(diss.mtrx, k)
	
	return(res)
}


#################################################################
# This method is a wrapper function to apply kmedoids clustering method.
#   It iterates over all possible number of clusters in order to apply kmedoids clustering method.
#   In this function, there are two exceptions which are handled:
#       - k=1, i.e. the number of cluster is 1 (i.e. number of solution classes)
#       - k=n, i.e. the number of cluster is n (single node clusters) (i.e. number of solution classes)
#   These two cases are defined for kmedoids clustering. So, we simply put NA for silhouette score.
#   Otherwise, we apply kmedoids clustering method.
#
# dist.mtrx: dissimilarity matrix
# clu.algo.name: a clustering method that can be applied into dissimilarity matrix
# k.limits: the lower and upper bound limits regarding the the number of cluster values (i.e. number of solution classes) in kmedoids
#
# returns a list containing:
#   - avg.sil.list: silhouette scores
#   - partitions: all partitions
#################################################################
perform.clustering = function(dist.mtrx, clu.algo.name="kmedoids", k.limits=c(NA,NA)){
	
	res = list()
	avg.sil.list = c()
	partitions = list()
	
	# -----------------------------------------------------------
	m = nrow(dist.mtrx) # nb solution
	k.lower.bound = 1 # min value by default
	k.upper.bound = m # max value by default

	if(!is.na(k.limits[2]))
		k.upper.bound = min(m, k.limits[2])
    if(k.upper.bound == 1)
        k.lower.bound = 1
	else if(!is.na(k.limits[1]))
		k.lower.bound = k.limits[1]
	# -----------------------------------------------------------
	print(k.lower.bound)
    print(k.upper.bound)
	res[["k.values"]] = seq(k.lower.bound,k.upper.bound, by=1)
	
	# when we start to iterate from 1, we get the following error:
	# Error in sil[, "sil_width"] : nombre de dimensions incorrect
	
	for(k in res[["k.values"]]){ # 'k' must be in {1,2, .., n-1} for pam()
		tlog("k: ", k)
		k.desc = as.character(k)
		
		avg.sil.list[k.desc] = NA
		partitions[[k.desc]] = NA
		
		if(k == 1){
			partitions[[k.desc]] = rep(1, m) # all nodes are in the same cluster
			avg.sil.list[k.desc] = NA
		} else if(k == m){
			partitions[[k.desc]] = seq(1, m, by=1) # all nodes are in different cluster
			avg.sil.list[k.desc] = NA
		} else {
			result = apply.clustering.method(clu.algo.name, dist.mtrx, k)
			partitions[[k.desc]] = result$clustering
			
			sil = silhouette(partitions[[k.desc]], dist.mtrx)
			# 'silhouette' score is computed for each point, take the average of them
			avg.sil.list[k.desc] = sprintf("%.4f", mean(sil[,"sil_width"]))
		}

	}
	# if 'silhouette' score is 1 for point1, that means the point is well placed in its cluster (by means of avg. dist.)
	# if 'silhouette' score is -1 for point1, that is the opposite meaning
	
	res[["avg.sil.list"]] = avg.sil.list
	res[["partitions"]] = partitions
		
	return(res)
}






#################################################################
# It is a wrapper function to perform the kmedoids clustering analysis.
#    It iterates over graph types and dissimilarity measures.
#    It does the following tasks:
#       - call the method kmedoids
#       - record and plot silhouette scores
#       - write cluster analysis.memberships into file
#       - prepare cluster analysis result folders
#       - copy and plot corresp partitions
#
# n: graph size
# l0: number of clusters (i.e. number of solution classes)
# d: density
# prop.mispl: proportion of misplaced links
# prop.neg: proportion of negative links
# network.no: network id (the identifiers start from 1)
# algo.name: the name of correlation clustering algorithm to run
# comp.measures: dissimilarity measure to compare the optimal solutions
# k.limits: the lower and upper bound limits regarding the the number of cluster values (i.e. number of solution classes) in kmedoids
# force: whether or not the existing files are overwritten by a fresh call of all corresponding methods (e.g partitioning method)
#
##################################################################
perform.cluster.analysis = function(n, l0, d, prop.mispl, prop.neg, network.no,
                                    algo.name, comp.measures, k.limits, force)
{

	tlog(16, "start to perform cluster analysis with exact algorithms")	
	tlog(20, "performing cluster analysis => algo.name: ", algo.name)
		
	for(graph.desc.name in c(SIGNED.UNWEIGHTED.FILE)){ # SIGNED.WEIGHTED.FILE
		tlog(24, "performing cluster analysis => graph.desc.name: ", graph.desc.name)
		
		net.folder = get.input.network.folder.path(n, l0, d, prop.mispl, prop.neg, network.no)
#			graph.name = paste0(SIGNED.UNWEIGHTED.FILE, ".G") # our reference graph to get optimal and worst imbalance count value
#			network.path = file.path(net.folder, graph.name)
#			g = read.graph.ils(network.path)

		if(dir.exists(net.folder)){
			graph.name = paste0(GRAPH.FILE.LAYOUT.PREFIX,"-",graph.desc.name, ".graphml")
			network.path = file.path(net.folder, graph.name)
			g = suppressWarnings(read.graph(file=network.path, format="graphml"))
	
			part.folder = get.part.folder.path(n, l0, d, prop.mispl, prop.neg, network.no, algo.name, graph.desc.name)
			eval.algo.folder = get.eval.folder.path(n, l0, d, prop.mispl, detected.imb=NA, prop.neg, k=ALL, network.no, algo.name, graph.desc.name)
			clu.analysis.folder = get.clu.analysis.folder.path(n, l0, d, prop.mispl, prop.neg, network.no, algo.name, graph.desc.name)
			
			partitions = load.membership.files(part.folder)
			m = length(partitions) # nb partition
			
			
			for(measure in comp.measures){
			    matrix.file = file.path(eval.algo.folder,paste0(EVAL.DIST.MATRIX.FILE.PREFIX,"-",measure,".csv"))
			    # print(matrix.file)
			    if(file.exists(matrix.file)){ # if the distance matrix has been already computed
			        tlog(24, "performing cluster analysis => measure: ", measure)
			        
    			    # ---------------------------------------------------
    			    dist.mtrx = as.matrix(read.csv(matrix.file, row.names = 1, header= TRUE, check.names=FALSE)) # square matrix
    			    dist.mtrx = matrix(as.numeric(sprintf("%.4f", dist.mtrx)), m, m)
    			    # ---------------------------------------------------
    			    
        			clu.analysis.measure.folder = file.path(clu.analysis.folder,measure)
        			if(!dir.exists(clu.analysis.measure.folder))
        				dir.create(path=clu.analysis.measure.folder, showWarnings=FALSE, recursive=TRUE)
        			
	
        			# ---------------------------------------------------
        			# ---------------------------------------------------
        			silhouette.scores.file = file.path(clu.analysis.measure.folder,SILH.SCORE.FILENAME)
        			silh.folder.names = NA
        			
        			exists.silh = file.exists(paste0(silhouette.scores.file,".csv"))				
        			clu.result = NA
        			if(!exists.silh)
        			{
        				# ==============================================================
        				# CLUSTERING
        				print("!!!!!! CLUSTERING !!!!!!!!")
        				clu.result = perform.clustering(dist.mtrx, "kmedoids", k.limits)
        				# ==============================================================
        				
        				record.and.plot.silhouette.scores(silhouette.scores.file, clu.result[["avg.sil.list"]],
        				                                  clu.result[["k.values"]])
        				
        				write.cluster.analysis.memberships.into.file(clu.analysis.measure.folder,
        				                    clu.result[["partitions"]], clu.result[["k.values"]], clu.result[["avg.sil.list"]], "silh")	
        				
        				silh.folder.names = prepare.cluster.analysis.result.folders(clu.result[["avg.sil.list"]],
        				                                                            clu.result[["k.values"]], "silh")
        
        			} else if(exists.silh){
        				silh.df = read.csv(file=paste0(silhouette.scores.file,".csv"),check.names=FALSE)
        				silh.folder.names = prepare.cluster.analysis.result.folders(silh.df[["avg.sil.score"]],
        						silh.df[["k"]], "silh")
        				clu.result.partitions = retreive.clu.result.partitions(clu.analysis.measure.folder, silh.folder.names, silh.df[["k"]])
        				
        				clu.result = list()
        				clu.result[["avg.sil.list"]] = silh.df[["avg.sil.score"]]
        				clu.result[["k.values"]] = silh.df[["k"]]
        				clu.result[["partitions"]] = clu.result.partitions
        				
        			}
        		
        							
        			# ---------------------------------------------------
        			max.indx=1 # by default, there is one solution class
        			
        			nb.sol = length(clu.result[["avg.sil.list"]])
        			if(nb.sol == 2){
        				mems = load.membership.files(part.folder)
        				d = compare.partition.pair.as.distance(mems[[1]], mems[[2]], measures=NVI)
        				if(d >= 0.50){
        					max.indx = 2
        				}
        			}
        			#  we do this: -c(1,m), because the first and last silh val are NA, since it is not defined for Silhouette
        			if(nb.sol > 2 && any( (as.numeric(clu.result[["avg.sil.list"]]) >= 0.50)[-c(1,m)] ))
        				max.indx = which.max(as.numeric(clu.result[["avg.sil.list"]]))
        			
        			# if max.index = 1 here, that means that although many sols exists, there is only 1 sol class
        			copy.and.plot.corresp.partitions(g, clu.analysis.measure.folder, silh.folder.names[max.indx], partitions)
        			
        			# if(nb.sol>2 && max.indx == 1){ # although many sols exists, there is only 1 sol class
        			# 	max.indx = which.max(as.numeric(clu.result[["avg.sil.list"]]))
        			# 	copy.and.plot.corresp.partitions(g, clu.analysis.measure.folder, silh.folder.names[max.indx], partitions)
        			# }
        			# ---------------------------------------------------
			
			    } # end if(file.exists(matrix.file))
			
		    } # end comp measures

		} # end if(dir.exists(net.folder))
	
	}

	
}
	





#################################################################
# It is the wrapper function to perform the kmedoids clustering analysis.
#  It iterates over graph sizes, proportion of misplaced links, proportion of negative links and generated different network sizes.
#
# graph.sizes: a vector of values regarding graph sizes to be considered
# d: density (it is a single value)
# l0: number of clusters to be considered (it is a single value)
# prop.mispls: a vector of values regarding proportion of misplaced links
# prop.negs: a vector of values regarding proportion of negative links (for now, it is not operational)
# in.rand.net.folders: a vector of values regarding input random graph folders. Sequantial integers (1, .., 10)
# cor.clu.exact.algo: the name of correlation clustering algorithm to run
# comp.measures: dissimilarity measure to compare the optimal solutions
# k.limits: the lower and upper bound limits regarding the the number of cluster values (i.e. number of solution classes) in kmedoids
# force: whether or not the existing files are overwritten by a fresh call of all corresponding methods (e.g partitioning method)
#
##################################################################
perform.all.cluster.analysis = function(graph.sizes, d, l0, prop.mispls, prop.negs, in.rand.net.folders,
	cor.clu.exact.algo, comp.measures, k.limits, force)
{
	tlog("starts performing cluster analysis")
	for(n in graph.sizes){
		tlog(4, "performing cluster analysis => n: ", n)
		
		for(prop.mispl in prop.mispls){
			tlog(8, "performing cluster analysis => prop.mispl: ", prop.mispl)
			
		    if(is.na(prop.negs) && d == 1){
		        prop.negs = compute.prop.neg(n, d, l0, prop.mispl)
		    }	
			
			for(prop.neg in prop.negs){
				tlog(12, "performing cluster analysis => prop.neg: ", prop.neg)
				
				for(network.no in in.rand.net.folders){
					tlog(16, "performing cluster analysis => network.no: ", network.no)
						
					perform.cluster.analysis(n, l0, d, prop.mispl, prop.neg, network.no,
						 cor.clu.exact.algo, comp.measures, k.limits, force)
	
				}
				
			}
			
		}
		
	}
	
}
