#################################################################
# For a given two dissimilarity measures, it creates alluvial diagram of three layers.
#   Given a cluster of kmedoids, these are:
#   1) the representative partition of the cluster based on first dissimilarity measure
#   2) the partitions detected by only measure1 or only by measure2
#   3) the representative partition of the cluster based on second dissimilarity measure
# If we denote A the number of partitions in the cluster of kmedoids, we generate A-2 alluvial diagrams in the end.
# Note that we do not check if the representative partition of both measures is the same or not.
#
# clu.analysis.difference.folder: folder in which the diagrams are created => output folder
# measure1: first dissimilarity measure
# repr.filename.measure1: representative partition based on first dissimilarity measure
# measure2: second dissimilarity measure
# repr.filename.measure2: representative partition based on second dissimilarity measure
# diff.mbrshp.files: partitions different than those representative ones
# plot.title: plot title
# 
##################################################################
plot.difference.in.alluvial.diagram = function(clu.analysis.difference.folder,measure1,repr.filename.measure1,measure2,
	repr.filename.measure2,diff.mbrshp.files, plot.title=""){
#	mbrshp.files = list.files(path=clu.analysis.difference.folder,pattern=paste0("^",MBRSHP.FILE.PREFIX,".*\\.txt$"))
	diff.mbrshp.files = gsub(".txt","",diff.mbrshp.files) # remove '.txt'
	
	for(mbrshp.file in diff.mbrshp.files){
		# ================================================
		# prepare data
		diff.mem = read.table(file=file.path(clu.analysis.difference.folder,paste0(mbrshp.file,".txt")))$V1
		repr.mem1 = read.table(file=file.path(clu.analysis.difference.folder,repr.filename.measure1))$V1
		repr.mem2 = read.table(file=file.path(clu.analysis.difference.folder,repr.filename.measure2))$V1
		partitions.df = cbind(repr.mem1,diff.mem,repr.mem2)
		x.labels = c(measure1,mbrshp.file,measure2)
		
		# plot - alluvial diagram
		out.folder = clu.analysis.difference.folder
		out.filename = paste0(ALLUVIAL.DIAGRAM.PLOT.TYPE,"-",mbrshp.file)
		plot.alluvial.diagram(out.folder, out.filename, partitions.df, x.labels, plot.title, show.stripes.border=TRUE, format="JPEG")
		# ================================================
	}

}





#################################################################
# For a given cluster value (e.g. 3, 4, etc.): ==> this value can be the best cluster number according to silhouette or any value
#   - It does binding between clusters of measure1 and those of measure2. For instance, clu3 of NMI is associated with clu1 of ARI because the intersection is max
#   - It iterates over the clusters
#   - for each cluster, it identifies the partitions found only by measure2 and those found only by measure1. 
#       This constitutes the partitions differing from one measure to anohter. That is those partitions that we want to analyze and to see the differenc
#   - It shows the difference through alluvial diagram
#
#
# measure1: first dissimilarity measure
# measure2: second dissimilarity measure
# clu.analysis.difference.folder: output folder
# clu.analysis.silh.folder.name.measure1: folder containing the information related to silhouette for measure1
# clu.analysis.silh.folder.name.measure2: folder containing the information related to silhouette for measure1
# clu.characterization.silh.folder.name.measure1: folder in which representative partitions are already computed for measure1
# clu.characterization.silh.folder.name.measure2: folder in which representative partitions are already computed for measure2
# force: whether or not the existing files are overwritten by a fresh call of all corresponding methods (e.g partitioning method)
#
##################################################################
show.difference = function(measure1, measure2, clu.analysis.difference.folder, clu.analysis.silh.folder.name.measure1, clu.analysis.silh.folder.name.measure2,
	 clu.characterization.silh.folder.name.measure1, clu.characterization.silh.folder.name.measure2, force)
{

	cluster.analysis.membership.file1 = file.path(clu.analysis.silh.folder.name.measure1, paste0(MBRSHP.FILE.PREFIX,".txt"))
	cluster.analysis.membership.file2 = file.path(clu.analysis.silh.folder.name.measure2, paste0(MBRSHP.FILE.PREFIX,".txt"))
#	print(cluster.analysis.membership.file1)
#	print(cluster.analysis.membership.file2)
	mem1 = read.table(file=cluster.analysis.membership.file1)$V1
	k.value.measure1 = length(unique(mem1))
	mem2 = read.table(file=cluster.analysis.membership.file2)$V1
	k.value.measure2 = length(unique(mem2))
	binding.result.df = bind.two.membership(mem1, mem2, fill.missing.part=FALSE)
	
	# for each cluster of kmedoids results
	tlog(20, "starting to compare 2 clusters of 2 measures")
	for(i in seq(1,nrow(binding.result.df))){
		k = binding.result.df[i,"k"]
		clu.id = binding.result.df[i,"associated cluster id"]
		tlog(20, "compare 2 clusters: ", paste0(CLUSTER.ID.FOLDER.PREFIX,k), " in ", measure1," vs ",
			 paste0(CLUSTER.ID.FOLDER.PREFIX,clu.id), " in ", measure2)
		
		
		# centroids = representative memberships
		characterization.repr.folder1 = file.path(clu.characterization.silh.folder.name.measure1, REPRESENTATIVE.FOLDER.NAME)
		repr.membership.file1 = file.path(characterization.repr.folder1,paste0(MBRSHP.FILE.PREFIX,"-",CLUSTER.ID.FOLDER.PREFIX,k,".txt"))
		characterization.repr.folder2 = file.path(clu.characterization.silh.folder.name.measure2, REPRESENTATIVE.FOLDER.NAME)
		repr.membership.file2 = file.path(characterization.repr.folder2,paste0(MBRSHP.FILE.PREFIX,"-",CLUSTER.ID.FOLDER.PREFIX,clu.id,".txt"))
		
		clu.folder.measure1 = file.path(clu.analysis.silh.folder.name.measure1, paste0(CLUSTER.ID.FOLDER.PREFIX,k))
		clu.folder.measure2 = file.path(clu.analysis.silh.folder.name.measure2, paste0(CLUSTER.ID.FOLDER.PREFIX,clu.id))
		membership.files.measure1 = list.files(path=clu.folder.measure1,pattern=paste0("^",MBRSHP.FILE.PREFIX,".*\\.txt$"))
#		print("membership files for measure 1")
#		print(membership.files.measure1)
		membership.files.measure2 = list.files(path=clu.folder.measure2,pattern=paste0("^",MBRSHP.FILE.PREFIX,".*\\.txt$"))
#		print("membership files for measure 2")
#		print(membership.files.measure2)
		# difference
#		print(membership.files.measure1)
#		print(membership.files.measure2)

		# setdiff(c(1,5),c(1,3,4)) => {5}
		# setdiff(c(1,3,4),c(1,5)) => {3,4}
		# find the files in measure1 which differs from the files found in measure2
		diff.membership.filenames.measure1 = setdiff(membership.files.measure1, membership.files.measure2)
		# find the files in measure2 which differs from the files found in measure1
		diff.membership.filenames.measure2 = setdiff(membership.files.measure2, membership.files.measure1)
		# then, combine those files
		diff.membership.filenames = c(diff.membership.filenames.measure1,diff.membership.filenames.measure2)
#		print("diff:")
#		print(diff.membership.filenames)
		# identify which membership is located in which folder
		indx.measure1 = which(diff.membership.filenames %in% membership.files.measure1)
		indx.measure2 = which(diff.membership.filenames %in% membership.files.measure2)
		diff.membership.files.measure1 = file.path(clu.folder.measure1, diff.membership.filenames[indx.measure1])
		diff.membership.files.measure2 = file.path(clu.folder.measure2, diff.membership.filenames[indx.measure2])
		diff.membership.files = c(diff.membership.files.measure1, diff.membership.files.measure2)

		# =================================================================
		
		# if there is at least one file which differs from measure1 to measure2 or from measure2 to measure1
		if(length(diff.membership.filenames) > 0){
			# print(" -------------------------- DIFFFFFFFFFFFFFFFFff -------------------------")
			clu.desc.name = paste0(paste0(CLUSTER.ID.FOLDER.PREFIX,k), " vs ", paste0(CLUSTER.ID.FOLDER.PREFIX,clu.id))
			clu.analysis.difference.folder.clu.desc = file.path(clu.analysis.difference.folder, clu.desc.name)
			if(!dir.exists(clu.analysis.difference.folder.clu.desc) || force){ # if the folder exists, then we assume that the diagrams are already created
				dir.create(path=clu.analysis.difference.folder.clu.desc, showWarnings=FALSE, recursive=TRUE)
			
    			repr.filename.measure1 = paste0(REPRESENTATIVE.FILE.PREFIX,"-",MBRSHP.FILE.PREFIX,"-",measure1,".txt")
    #			print(repr.membership.file1)
    #			print(file.path(clu.analysis.difference.folder,repr.filename.measure1))
    			file.copy(from=repr.membership.file1, to=file.path(clu.analysis.difference.folder.clu.desc,repr.filename.measure1), overwrite=TRUE)
    			repr.filename.measure2 = paste0(REPRESENTATIVE.FILE.PREFIX,"-",MBRSHP.FILE.PREFIX,"-",measure2,".txt")
    			file.copy(from=repr.membership.file2, to=file.path(clu.analysis.difference.folder.clu.desc,repr.filename.measure2), overwrite=TRUE)
    			
    			for(diff.membership.file in diff.membership.files){
    				# do not specify target file names: preserve the same names
    				file.copy(from=diff.membership.file, to=clu.analysis.difference.folder.clu.desc, overwrite=TRUE)
    			}
    
    			# plot the diff memberships which is in measure 1, but not in measure 2
    			plot.title = paste0("First and last are representative partitions: (The middle belongs to ",measure1,")")
    			plot.difference.in.alluvial.diagram(clu.analysis.difference.folder.clu.desc,measure1,repr.filename.measure1,
    				measure2,repr.filename.measure2,diff.membership.filenames[indx.measure1],plot.title)
    				
    			# plot the diff memberships which is in measure 2, but not in measure 1
    			plot.title = paste0("First and last are representative partitions: (The middle belongs to ",measure2,")")
    #			plot.title = measure2
    			plot.difference.in.alluvial.diagram(clu.analysis.difference.folder.clu.desc,measure1,repr.filename.measure1,
    				measure2,repr.filename.measure2,diff.membership.filenames[indx.measure2],plot.title)
			
		    }
		}
		
	}

}


	
#################################################################
# It is the wrapper function to show the difference between two dissimilarty measures.
# In this function, we show the difference for only the best cluster number according to silhouette. 
#  
# measure1: first dissimilarity measure
# measure2: second dissimilarity measure
# clu.analysis.difference.folder: output folder
# clu.analysis.measure1.folder: folder containing the information related to kmedoids clustering results for measure1
# clu.analysis.measure2.folder: folder containing the information related to kmedoids clustering results for measure2
# silhouette.scores.file.measure1: the file containing the silhouette scores for measure1 (after k-medoids)
# silhouette.scores.file.measure2: the file containing the silhouette scores for measure2 (after k-medoids)
# clu.characterization.measure1.folder: folder of cluster characterization analysis for measure1
# clu.characterization.measure2.folder: folder of cluster characterization analysis for measure2
# force: whether or not the existing files are overwritten by a fresh call of all corresponding methods (e.g partitioning method)
#
##################################################################
show.difference.for.best.k.in.cluster.analysis = function(measure1, measure2, clu.analysis.difference.folder, clu.analysis.measure1.folder, clu.analysis.measure2.folder, 
		silhouette.scores.file.measure1, silhouette.scores.file.measure2, clu.characterization.measure1.folder, clu.characterization.measure2.folder,
		force)
{
	measures.desc.name = paste0(measure1, " vs ", measure2)
	clu.analysis.difference.folder = file.path(clu.analysis.difference.folder,measures.desc.name,BEST.K.FOLDER.NAME)
	
	
	silh.df1 = read.csv(file=paste0(silhouette.scores.file.measure1,".csv"),check.names=FALSE)
	silh.df2 = read.csv(file=paste0(silhouette.scores.file.measure2,".csv"),check.names=FALSE)
	if(all(is.na(silh.df1[["avg.sil.score"]])) || all(is.na(silh.df2[["avg.sil.score"]]))){
		tlog(20, "warning: there is only NA scores in silhouette score list")
	} else {
	
		silh.folder.names = prepare.cluster.analysis.result.folders(silh.df1[["avg.sil.score"]], silh.df1[["k"]])
		best.indx = which.max(silh.df1[["avg.sil.score"]])
		best.k.measure1 = silh.df1[["k"]][best.indx]
		best.silh.foldername1 = silh.folder.names[best.indx]
		clu.analysis.silh.folder.name.measure1 = file.path(clu.analysis.measure1.folder,best.silh.foldername1)
		clu.characterization.silh.folder.name.measure1 = file.path(clu.characterization.measure1.folder,best.silh.foldername1)
		
		
		silh.folder.names = prepare.cluster.analysis.result.folders(silh.df2[["avg.sil.score"]],silh.df2[["k"]])
		best.indx = which.max(silh.df2[["avg.sil.score"]])
		best.k.measure2 = silh.df2[["k"]][best.indx]
		best.silh.foldername2 = silh.folder.names[best.indx]
		clu.analysis.silh.folder.name.measure2 = file.path(clu.analysis.measure2.folder,best.silh.foldername2)
		clu.characterization.silh.folder.name.measure2 = file.path(clu.characterization.measure2.folder,best.silh.foldername2)
		
		clu.analysis.difference.folder = file.path(clu.analysis.difference.folder, paste0(best.silh.foldername1,"-",best.silh.foldername2))
#		if(!dir.exists(clu.analysis.difference.folder))
#			dir.create(path=clu.analysis.difference.folder, showWarnings=FALSE, recursive=TRUE)
		
		show.difference(measure1, measure2, clu.analysis.difference.folder, clu.analysis.silh.folder.name.measure1, clu.analysis.silh.folder.name.measure2, 
			 clu.characterization.silh.folder.name.measure1, clu.characterization.silh.folder.name.measure2, force)
	}
}




#################################################################
# It is the wrapper function to show the difference between two dissimilarty measures.
# In this function, we show the difference for all cluster numbers (e.g. 1, 2, 3, etc.)
#  
# measure1: first dissimilarity measure
# measure2: second dissimilarity measure
# clu.analysis.difference.folder: output folder
# clu.analysis.measure1.folder: folder containing the information related to kmedoids clustering results for measure1
# clu.analysis.measure2.folder: folder containing the information related to kmedoids clustering results for measure2
# silhouette.scores.file.measure1: the file containing the silhouette scores for measure1 (after k-medoids)
# silhouette.scores.file.measure2: the file containing the silhouette scores for measure2 (after k-medoids)
# clu.characterization.measure1.folder: folder of cluster characterization analysis for measure1
# clu.characterization.measure2.folder: folder of cluster characterization analysis for measure2
# k.limits: the lower and upper bound for the number of cluster values in kmedoids cluster analysis.
#               Normally, we do not really need to know 'k.limits' here, since kmedoids clustering is suppsoed to be performed with k.limits.
# force: whether or not the existing files are overwritten by a fresh call of all corresponding methods (e.g partitioning method)
#
##################################################################
show.difference.for.all.k.in.cluster.analysis = function(measure1, measure2, clu.analysis.difference.folder, clu.analysis.measure1.folder, clu.analysis.measure2.folder, 
		silhouette.scores.file.measure1, silhouette.scores.file.measure2, clu.characterization.measure1.folder, clu.characterization.measure2.folder,
		k.limits, force)
{
	measures.desc.name = paste0(measure1, " vs ", measure2)
	clu.analysis.difference.folder = file.path(clu.analysis.difference.folder,measures.desc.name,ALL.K.FOLDER.NAME)
	
	
	silh.df1 = read.csv(file=paste0(silhouette.scores.file.measure1,".csv"),check.names=FALSE)
	silh.folder.names.measure1 = prepare.cluster.analysis.result.folders(silh.df1[["avg.sil.score"]], silh.df1[["k"]])
	silh.df2 = read.csv(file=paste0(silhouette.scores.file.measure2,".csv"),check.names=FALSE)
	silh.folder.names.measure2 = prepare.cluster.analysis.result.folders(silh.df2[["avg.sil.score"]],silh.df2[["k"]])
	
	k.values = silh.df1[["k"]] # this is the same for silh.df2[["k"]]
	# =================
	# normally, we do not really need to know 'k.limits' here, since kmedoids clustering is suppsoed to be performed with k.limits.
	remove.indx = which((k.values < k.limits[1]) | (k.values > k.limits[2]))
	if(length(remove.indx)>0)
	    k.values = k.values[-remove.indx]
	# =================
	
	for(k in k.values){ # use k values as index
		tlog(20,"k: ",k)
		silh.foldername1 = silh.folder.names.measure1[k]
		silh.foldername2 = silh.folder.names.measure2[k]
		clu.analysis.silh.folder.name.measure1 = file.path(clu.analysis.measure1.folder,silh.foldername1)
		clu.characterization.silh.folder.name.measure1 = file.path(clu.characterization.measure1.folder,silh.foldername1)
		clu.analysis.silh.folder.name.measure2 = file.path(clu.analysis.measure2.folder,silh.foldername2)
		clu.characterization.silh.folder.name.measure2 = file.path(clu.characterization.measure2.folder,silh.foldername2)
		
		diff.folder.name = file.path(clu.analysis.difference.folder, paste0(silh.foldername1,"-",silh.foldername2))
#		if(!dir.exists(diff.folder.name))
#			dir.create(path=diff.folder.name, showWarnings=FALSE, recursive=TRUE)
		
		
		show.difference(measure1, measure2, diff.folder.name, clu.analysis.silh.folder.name.measure1, clu.analysis.silh.folder.name.measure2, 
		 clu.characterization.silh.folder.name.measure1, clu.characterization.silh.folder.name.measure2, force)
	} 
}



#################################################################
# It is the wrapper function to show the difference between two dissimilarty measures.
#   It iterates over exact algorithms used to obtain partitions, graph types and dissimilarity measures.
#
# n: graph size
# l0: number of cluster
# d: density
# prop.mispl: proportion of misplaced links
# prop.neg: proportion of negative links
# network.no: network id (the identifiers start from 1)
# cor.clu.exact.algo: the name of correlation clustering algorithm to run
# comp.measures: dissimilarity measure to compare the optimal solutions
# k.limits: the lower and upper bound for the number of cluster values in kmedoids cluster analysis. 
#            Used only in show.difference.for.all.k.in.cluster.analysis()
# force: whether or not the existing files are overwritten by a fresh call of all corresponding methods (e.g partitioning method)
#
##################################################################
show.cluster.analysis.difference = function(n, l0, d, prop.mispl, prop.neg, network.no,
	 cor.clu.exact.algo, comp.measures, k.limits, force)
{

	#tlog(16, "start to show cluster analysis difference with exact algorithms")	
	#for(algo.name in cor.clu.exact.algos){
		tlog(20, "cluster analysis difference => algo.name: ", cor.clu.exact.algo)
		
		for(graph.desc.name in c(SIGNED.UNWEIGHTED.FILE)){ # OPPOSITE.SIGNED.UNWEIGHTED.FILE
			tlog(24, "cluster analysis difference => graph.desc.name: ", graph.desc.name)
			
			net.folder = get.input.network.folder.path(n, l0, d, prop.mispl, prop.neg, network.no)
			graph.name = paste0(SIGNED.UNWEIGHTED.FILE, ".G") # our reference graph to get optimal and worst imbalance count value
			network.path = file.path(net.folder, graph.name)
			g = read.graph.ils(network.path)
	
			part.folder = get.part.folder.path(n, l0, d, prop.mispl, prop.neg, network.no, cor.clu.exact.algo, graph.desc.name)
			eval.algo.folder = get.eval.folder.path(n, l0, d, prop.mispl, detected.imb=NA, prop.neg, network.no, cor.clu.exact.algo, graph.desc.name)
			clu.analysis.folder = get.clu.analysis.folder.path(n, l0, d, prop.mispl, prop.neg, network.no, cor.clu.exact.algo, graph.desc.name)
			clu.characterization.folder = get.clu.characterization.folder.path(n, l0, d, prop.mispl, prop.neg, network.no, cor.clu.exact.algo, graph.desc.name)
			clu.analysis.difference.folder = get.clu.analysis.difference.folder.path(n, l0, d, prop.mispl, prop.neg, network.no, cor.clu.exact.algo, graph.desc.name)
		
			partitions = load.membership.files(part.folder)
			m = length(partitions) # nb partition
			
#			> combn(x=3, m=2)
#			   		  [,1] [,2] [,3]
#				[1,]    1    1    2
#				[2,]    2    3    3
					
			nb.measure = length(comp.measures)
			comb = combn(x=nb.measure, m=2)
			nb.col = ncol(comb)
			# for each pair of dissimilarity measures
			for(i in seq(1,nb.col)){
				indx1 = comb[1,i]
				indx2 = comb[2,i]
				measure1 = comp.measures[indx1]
				measure2 = comp.measures[indx2]
				tlog(20, "compare measures: ", measure1, " and ", measure2)
				
				
				clu.characterization.measure1.folder = file.path(clu.characterization.folder,measure1)
				clu.characterization.measure2.folder = file.path(clu.characterization.folder,measure2)
				clu.analysis.measure1.folder = file.path(clu.analysis.folder,measure1)
				clu.analysis.measure2.folder = file.path(clu.analysis.folder,measure2)
				silhouette.scores.file.measure1 = file.path(clu.analysis.measure1.folder,SILH.SCORE.FILENAME)
				silhouette.scores.file.measure2 = file.path(clu.analysis.measure2.folder,SILH.SCORE.FILENAME)
				
				show.difference.for.best.k.in.cluster.analysis(measure1, measure2, clu.analysis.difference.folder, clu.analysis.measure1.folder,
					 clu.analysis.measure2.folder, silhouette.scores.file.measure1, silhouette.scores.file.measure2, 
					 clu.characterization.measure1.folder, clu.characterization.measure2.folder, force)
						
				show.difference.for.all.k.in.cluster.analysis(measure1, measure2, clu.analysis.difference.folder, clu.analysis.measure1.folder,
				 	clu.analysis.measure2.folder, silhouette.scores.file.measure1, silhouette.scores.file.measure2,
				 	 clu.characterization.measure1.folder, clu.characterization.measure2.folder, k.limits, force)
				
			}
		
		}
	#}
	
}
	




#################################################################
#
# It is the wrapper function to show the difference between two dissimilarty measures.
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
# k.limits
# force: whether or not the existing files are overwritten by a fresh call of all corresponding methods (e.g partitioning method)
#
#
##################################################################
show.all.cluster.analysis.difference = function(graph.sizes, d, l0, prop.mispls, prop.negs, in.rand.net.folders,
	cor.clu.exact.algo, comp.measures, k.limits, force)
{
	tlog("starts performing cluster analysis")
	for(n in graph.sizes){
		tlog(8, "performing cluster analysis => n: ", n)
		
		for(prop.mispl in prop.mispls){
			tlog(8, "performing cluster analysis => prop.mispl: ", prop.mispl)
			
		    if(is.na(prop.negs) && d == 1){
		        prop.negs = compute.prop.neg(n, d, l0, prop.mispl)
		    }
			
			for(prop.neg in prop.negs){
				tlog(12, "performing cluster analysis => prop.neg: ", prop.neg)
				
				for(network.no in in.rand.net.folders){
					tlog(16, "performing cluster analysis => network.no: ", network.no)
						
					show.cluster.analysis.difference(n, l0, d, prop.mispl, prop.neg, network.no,
						 cor.clu.exact.algo, comp.measures, k.limits, force)
	
				}
				
			}
			
		}
		
	}
	
}