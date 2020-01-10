
###########################################################################
# since we process only distinct heuristic solutions, the size of dist matrix is not 100x100
# Sometimes, we need to obtain that 100x100 matrix, for example for calcualting mean value of dist scores
#
# when process.row=TRUE and  process.col=TRUE, that means a square dist matrix for heursitic sols
# if only process.row=TRUE, that means, the matrix is not sqaure => heuristic sols are only in row side, not column side
###########################################################################
reconstruct.heur.dist.matrix = function(heur.sols.summary.filename, mtrx.file, process.row=TRUE, process.col=TRUE){
	df = read.csv(file=heur.sols.summary.filename, header=TRUE, stringsAsFactors=F)
	heur.partitions.ass = retreive.heuristic.partitions.info(df, HEURISTIC.REPETITIONS)
	nb.heurs = length(heur.partitions.ass)
	
	
	
	# ----------------------
	# turn the recorded matix into bigger matrix in which we duplicate heursitic info
	# the recorded matrix contins only the necessary info
	# but, if a given method finds 2 same sol, we record only one, in this dist matrix, we need to consider all sols
	dist.mtrx = read.csv(mtrx.file, row.names = 1, header= TRUE, check.names=FALSE) # sqaure matrix
	
	
	nb.row = nb.heurs
	if(!process.row)
		nb.row = nrow(dist.mtrx)
	nb.col = nb.heurs
	if(!process.col)
		nb.col = ncol(dist.mtrx)

	
	new.dist.mtrx = matrix(NA,nb.row,nb.col)
	for(i in seq(1,nb.row)){
		sol.id1 = i
		if(process.row)
			sol.id1 = heur.partitions.ass[i]+1 # +1 since sol ids start from 0 => matrix index compability
			
		for(j in seq(1,nb.col)){
			sol.id2 = j
			if(process.col)
				sol.id2 = heur.partitions.ass[j]+1 # +1 since sol ids start from 0 => matrix index compability
			new.dist.mtrx[i,j] = dist.mtrx[sol.id1,sol.id2]
		}
	}
	return(new.dist.mtrx)
}


###########################################################################
#
###########################################################################
is.clusterings.distinct = function(m1, m2){
	# check if they are of the same length
	binding = list()
	for(i in seq(1, length(m1))){
		i1 = as.character(m1[i])
		i2 = as.character(m2[i])
		indx = which(names(binding) %in% i1)
		if(length(indx) == 0)
			binding[[i1]] = i2
		else{
			if(binding[[i1]] != i2)
				return(1)
		}
	}
	return(0)
}




###########################################################################
# this method is used by plots, not here
# warning: for 'df', use 'stringsAsFactors=F'
#
# input:
# id 	freq		association
# 0		3		rep.no=0&id=0;rep.no=2&id=0;rep.no=3&id=0
# 1		1		rep.no=1&id=0
#
# output:
# [0, 1, 0, 0]   ===> new id vs. old id association ==> old ids 0, 2 and 3 are now 0, and old id 1 is now 1
# 
###########################################################################"
retreive.heuristic.partitions.info = function(df, heuristic.repetitions){
	
	# init
	partitions.ass.info = rep(list(NA),length(heuristic.repetitions))
	
	for(line.indx in 1:nrow(df)){
		# freq = df[line.indx, FREQ.COL.NAME]
		partition.new.id = as.numeric(df[line.indx, ID.COL.NAME])
		if(partition.new.id == -1) # TODO: remove it later
			break
			
		ass = df[line.indx, ASSOCIATION.COL.NAME]
		items = unlist(strsplit(ass, split=";")) # rep.no=36&id=0;rep.no=49&id=0
		for(item in items){
			parts = unlist(strsplit(item, split="&"))
			rep.no = as.numeric(unlist(strsplit(parts[1], split="="))[2])
			id = as.numeric(unlist(strsplit(parts[2], split="="))[2])		
			# fill list
			if( !is.na(partitions.ass.info[[rep.no]]) )
				partitions.ass.info[[rep.no]] = c(partitions.ass.info[[rep.no]], id)
			else
				partitions.ass.info[[rep.no]] = c(partition.new.id)
		}
	}
	return(unlist(partitions.ass.info))
}



###########################################################################"
# let us say 'values' for nb cluster is: [4, 6]
# since we have 100 repetitions, there are 100 values which are either 4 or 6
# reconstruct the values based on all heuristic partitions
# output something like: [4, 6, 6, 6, 4, 4 ...]
###########################################################################"
fill.vector.as.heuristic.partitions.info = function(values, heur.partitions.ass){
	#------- construct new values for all heuristic partitions (not for only distinct ones)
	freq = table(heur.partitions.ass)
	new.values = rep(NA, length(heur.partitions.ass))
	new.ids = as.integer(names(freq))
	for(i in new.ids){
		indx = which(heur.partitions.ass == i)
		new.values[indx] = values[i+1] # i+1 since ids start from 0
	}
	return(new.values)
}

#########################################################################################
#########################################################################################


#################################################################
# returns index of the associated heuristic membership which is already added into list
#  if membership candidate already exists. Otherwise, -1
#################################################################
check.if.partition.exists = function(all.mbrshps, mbrshp.candidate){
	indx = -1
	m = length(all.mbrshps)

	if(m>0){
		for(i in seq(1,m)){
			mbrshp = all.mbrshps[[i]]
			
			if(!is.clusterings.distinct(mbrshp, mbrshp.candidate)){
				indx=i
				break
			}
		}
	}
	
	return(indx)
}


#################################################################
#
#################################################################
process.partitions = function(distinct.mbrshps, candidate.mbrshps, heur.freqs, ass.info, rep.no){
	
	
	m = length(distinct.mbrshps) # partition number before updating
	for(i in seq(1, length(candidate.mbrshps))){
		mbrshp.candidate = candidate.mbrshps[[i]]
		indx = check.if.partition.exists(distinct.mbrshps, mbrshp.candidate)
		if(indx != -1){ # already exists
			heur.freqs[indx] = heur.freqs[indx] + 1
			ass.info[indx] = paste(ass.info[indx], paste0("rep.no=",rep.no,"&id=",(i-1)), sep=";")
		} else { # add into list
			m = m +1
			distinct.mbrshps[[m]] = mbrshp.candidate
			heur.freqs[m] = 1
			ass.info[m] = paste0("rep.no=",rep.no,"&id=",(i-1))
		}
	}

	
	result = list()
	result$mbrshps = distinct.mbrshps
	result$freqs = heur.freqs
	result$ass.info = ass.info
	return(result)
}



#################################################################
#
# d: density
# k: nb cluster
#
##################################################################
process.repetitions = function(n, k, d, prop.mispl, prop.neg, network.no,
		cor.clu.heurs.algos, heur.reps, force, plot.formats)
{
	
	for(algo.name in cor.clu.heurs.algos){
		tlog(20, "post processing heuristic partitions => algo.name: ", algo.name)
		
		for(graph.desc.name in c(SIGNED.UNWEIGHTED.FILE)){ # OPPOSITE.SIGNED.UNWEIGHTED.FILE
			tlog(20, "post processing heuristic partitions => graph.desc.name: ", graph.desc.name)
			
			distinct.mbrshps = c()
			heur.freqs = rep(0, length(heur.reps)) # init (in worst case)
			ass.info = rep("", length(heur.reps)) # init (in worst case)
			counter.indx = 1
			part.folder = get.part.folder.path(n, k, d, prop.mispl, prop.neg, network.no, algo.name, graph.desc.name)
			
			if(dir.exists(part.folder)){
				
				# check if already processed
				files = load.membership.files(part.folder)
				if(length(files) == 0 || force){
					# -------------------------------------------------------------
					# maybe, there are existant memberhsip files, but user needs to recreate them: so, we need to delete the old ones
					mbrshp.files = list.files(path=part.folder,pattern=paste0("^",MBRSHP.FILE.PREFIX,"*"))
					unlink(x=file.path(part.folder, mbrshp.files), force=TRUE)
					# -------------------------------------------------------------
					
					# ----------------------------------------------------------------------
					# process exec times
					exec.times = rep(NA, length(HEURISTIC.REPETITIONS))
					for(rep.no in heur.reps){
						part.rep.folder = get.part.folder.path(n, k, d, prop.mispl, prop.neg, network.no, algo.name, graph.desc.name, rep.no)
						exec.times[rep.no] = as.numeric(read.table(file=file.path(part.rep.folder,EXEC.TIME.FILENAME))$V1)
					}
					mtrx = matrix(exec.times, length(HEURISTIC.REPETITIONS), 1) # matrix with 1 column
					colnames(mtrx) = EXEC.TIME.COL.NAME
					sol.ids = HEURISTIC.REPETITIONS-1 # since partitions id start from 0
					rownames(mtrx) = paste("solution", sol.ids)
					table.file = file.path(part.folder,paste0(EVAL.EXEC.TIME.FILENAME,".csv"))
					write.csv(x=mtrx, file=table.file)
					# ----------------------------------------------------------------------
					
					
					# ----------------------------------------------------------------------
					# process memberhsip files
					for(rep.no in heur.reps){
						part.rep.folder = get.part.folder.path(n, k, d, prop.mispl, prop.neg, network.no, algo.name, graph.desc.name, rep.no)
						if(!dir.exists(part.rep.folder))
							dir.create(path=part.rep.folder, showWarnings=FALSE, recursive=TRUE)
						
						# ---------------------------------------------
						# a heuristic method might yield many solutions. But, we get only he first one
						# mbrshps = load.membership.files(part.rep.folder) # DO NOT USE IT
						mbrshps = list()
						id = 0 # the first solution
						table.file = file.path(part.rep.folder, paste0(MBRSHP.FILE.PREFIX,id,".txt"))
						mbrshp <- as.numeric(as.matrix(read.table(file=table.file, header=FALSE)))
						mbrshps[[id+1]] = mbrshp
						# --------------------------------------------
		
						result = process.partitions(distinct.mbrshps, mbrshps, heur.freqs, ass.info, rep.no)
						distinct.mbrshps = result$mbrshps # update
						heur.freqs = result$freqs # update
						ass.info = result$ass.info
					}
					# ----------------------------------------------------------------------
					
					# record distinct membership files into another folder (one level up)
					write.membership.files(part.folder, distinct.mbrshps)
					
					# record association and frequency information
					non.empty.part.indx = which(heur.freqs != 0)
					ids = rep(-1, length(heur.reps))
					ids[non.empty.part.indx] = seq(0, length(non.empty.part.indx)-1)
					df = cbind(ids, heur.freqs, ass.info)	
					colnames(df) = c(ID.COL.NAME, FREQ.COL.NAME, ASSOCIATION.COL.NAME)
					write.csv(x=df, file=file.path(part.folder,paste0(HEURISTIC.SOLS.SUMMARY.FILENAME,".csv")), row.names=FALSE)
				}
			
			}
			
		}
	}
	
}






#################################################################
#
# graph.sizes
# d: density
# k: nb cluster
# prop.mispls
# prop.negs
# cor.clu.heur.algos
# heur.reps: heuristic repetitions. Sequantial integers (1, .., 100)
# cor.clu.exact.algos
# in.rand.g.folders: input random graph folders. Sequantial integers (1, .., 10)
# force
# plot.formats
#
##################################################################
post.processing.heuristic.solutions = function(graph.sizes, d, k, prop.mispls, prop.negs, in.rand.net.folders,
		cor.clu.heur.algos, heur.reps, force, plot.formats)
{
	tlog("starts performing post processing heuristic partitions")
	for(n in graph.sizes){
		tlog(4, "post processing heuristic partitions => n: ", n)
		
		for(prop.mispl in prop.mispls){
			tlog(8, "post processing heuristic partitions => prop.mispl: ", prop.mispl)
			
			my.prop.negs = NA
			if(is.na(prop.negs)){
				prop.neg = compute.prop.neg(n, d, k, prop.mispl)
				my.prop.negs = prop.neg # we show explicitely that we get only 1 value since density=1.0
			}
			
			for(prop.neg in my.prop.negs){
				tlog(12, "post processing heuristic partitions => prop.neg: ", prop.neg)
				
				for(network.no in in.rand.net.folders){
					tlog(16, "post processing heuristic partitions => network.no: ", network.no)
					
					process.repetitions(n, k, d, prop.mispl, prop.neg, network.no,
							cor.clu.heur.algos, heur.reps, force, plot.formats)
					
				}
				
			}
			
		}
		
	}
	
}
