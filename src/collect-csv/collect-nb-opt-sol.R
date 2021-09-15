# TODO: Add comment
# 
# Author: nejat
###############################################################################








#################################################################
#
# It constructs a single row of the csv file in output. It retrieves the statistics 
#   related to the number of optimal solutions for the given values of the input parameters. 
#
# n: graph order
# l0: number of modules
# d: density
# prop.mispl: proportion of misplaced links
# prop.neg: proportion of negative links
# network.no: network id (the identifiers start from 1
# cor.clu.exact.algo: the name of correlation clustering algorithm to run
# force: whether or not the existing files are overwritten by a fresh call of all corresponding methods (e.g partitioning method)
#
##################################################################
retreive.nb.opt.solution = function(n, l0, d, prop.mispl, prop.neg, network.no,
		cor.clu.exact.algo, force)
{
	
	data =c()
	
	#tlog(16, "start to plot transition with exact algorithms")	
	tlog(20, "collecting networks for nb solution => algo.name: ", cor.clu.exact.algo)
	
	for(graph.desc.name in c(SIGNED.UNWEIGHTED.FILE)){ # OPPOSITE.SIGNED.UNWEIGHTED.FILE
		tlog(24, "collecting networks for nb solution => graph.desc.name: ", graph.desc.name)
		#part.folder = part.folder = get.part.folder.path(n, l0, d, prop.mispl, prop.neg, network.no, cor.clu.exact.algo, graph.desc.name)
		eval.folder = get.eval.folder.path(n, l0, d, prop.mispl, NA, prop.neg, k=ALL, network.no, cor.clu.exact.algo, graph.desc.name)
		print(eval.folder)
		curr.line = matrix(NA,1,1)
		if(dir.exists(eval.folder)){
			table.file = file.path(eval.folder, paste0(EVAL.NB.SOL.FILENAME,".csv"))
			#mems = load.membership.files(part.folder)
			curr.line = read.csv(table.file, row.names = 1, header= TRUE, check.names=FALSE, stringsAsFactors=FALSE)
			#colnames(curr.line) = NULL
			
		}
		rownames(curr.line) = paste0("n=", n, ",l0=", l0, ",d=", d, ",prop.mispl=", prop.mispl, ",prop.neg=", prop.neg, ",network.no=", network.no)
		colnames(curr.line) = c(NB.SOL.COL.NAME)
		data = rbind(data, curr.line)
	}
	
	return(data)
}










#################################################################
#
# It is the starting method in the aim of collecting the statistics related to the number of optimal solutions. 
#   It handles all networks by graph.sizes,  prop.mispls, my.prop.negs and in.rand.net.folders
#   For the given value of the parameter l0, it generates a csv file.
#
# graph.sizes: a vector of values regarding graph orders to be considered
# d: density (it is a single value)
# l0: number of modules to be considered (it is a single value)
# prop.mispls: a vector of values regarding proportion of misplaced links
# prop.negs: a vector of values regarding proportion of negative links (for now, it is not operational)
# in.rand.net.folders: a vector of values regarding input random graph folders. Sequantial integers (1, .., 10)
# cor.clu.exact.algo: the name of correlation clustering algorithm to run
# force: whether or not the existing files are overwritten by a fresh call of all corresponding methods (e.g partitioning method)
#
##################################################################
collect.all.nb.opt.solution = function(graph.sizes, d, l0, prop.mispls, prop.negs, in.rand.net.folders,
		cor.clu.exact.algo, comp.measures, force)
{
	
	
	tlog("starts collecting networks")
	for(n in graph.sizes){
		tlog(8, "collecting networks => n: ", n)
		
		all.data = c()
		
		for(prop.mispl in prop.mispls){
			tlog(8, "collecting networks => prop.mispl: ", prop.mispl)
			
		    if(is.na(prop.negs) && d == 1){
		        prop.negs = compute.prop.neg(n, d, l0, prop.mispl)
		    }
			
			for(prop.neg in prop.negs){
				tlog(12, "collecting networks => prop.neg: ", prop.neg)
				
				for(network.no in in.rand.net.folders){
					tlog(16, "partitioning networks => network.no: ", network.no)
					
					data = retreive.nb.opt.solution(n, l0, d, prop.mispl, prop.neg, network.no,
							cor.clu.exact.algo, comp.measures, force)
					if(!is.na(data))
						all.data = rbind(all.data, data)
				}
			}
			
		}
		
		if(length(all.data)>0){
			colnames(all.data) = c(NB.SOL.COL.NAME)
			#rownames(all.conn.comp.data) = NULL
			if(!dir.exists(OUTPUT.CSV.FOLDER))
				dir.create(OUTPUT.CSV.FOLDER, recursive=FALSE, showWarnings=FALSE)
			
			write.csv(file=file.path(OUTPUT.CSV.FOLDER,paste0("nb-sol-l0=",l0,"_",cor.clu.exact.algo,"_n=",n,"_d=",d,".csv")), x=all.data)
		}
		
	}
	
	
}
