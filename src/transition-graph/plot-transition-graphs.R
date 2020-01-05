
#################################################################
#
# d: density
# k: nb cluster
#
##################################################################
get.diam.vals = function(g, edit.dist.matrix){
	nb.sol = nrow(edit.dist.matrix)
	edit.dist.matrix = matrix(as.numeric(edit.dist.matrix),nb.sol,nb.sol)
	comp.mem = components(g)$membership
	k = length(unique(comp.mem))

	diam.vals = c()
	for(i in 1:k){
		indxs = which(comp.mem == i)
		diam.vals = c(diam.vals, max(edit.dist.matrix[indxs,indxs]))
	}
	return(diam.vals)
}
			
			
#################################################################
#
# Note that weight is intentionnaly 0.5, because, in legend we omit all edges whose weight is below 1
# Thus, we distinguish edge related to edit distance and edge between diff connected components
#
##################################################################
indicate.edit.dist.betw.components = function(g, edit.dist.matrix){
				
	nb.sol = nrow(edit.dist.matrix)
	edit.dist.matrix = matrix(as.numeric(edit.dist.matrix),nb.sol,nb.sol)
	comp.mem = components(g)$membership
	k = length(unique(comp.mem))
	
	if(k>1){
		comb = combn(k, 2, simplify = TRUE)
	
		for(i in 1:ncol(comb)){
			c1 = comb[1,i]
			c2 = comb[2,i]
			indxs1 = which(comp.mem == c1)
			nb.row = length(indxs1)
			indxs2 = which(comp.mem == c2)
			nb.col = length(indxs2)
			subm = matrix(edit.dist.matrix[indxs1,indxs2],nb.row,nb.col)
			
			if(nb.row == 1 && nb.col == 1){
				g = g + edge(indxs1,indxs2,label = subm[1,1], color = "black", weight=BETW.COMP.EDGE.WEIGHT)
				
			} else {
				min.indx = which.min(subm)
				d = subm[min.indx]
				col.indx = ceiling(min.indx/nb.row)
				row.indx = min.indx%%nb.row

				if(row.indx == 0)
					row.indx = nb.row
				g = g + edge(indxs1[row.indx],indxs2[col.indx], label = d, color = "black", weight=BETW.COMP.EDGE.WEIGHT)

			}
		}
	}
	
	return(g)
}



#################################################################
#
# d: density
# k: nb cluster
#
##################################################################
build.and.plot.transition.graph = function(n, l0, d, prop.mispl, prop.neg, network.no,
               cor.clu.exact.algo, comp.measures, force, plot.formats)
{
	
	conn.comp.data =c()
	
	net.folder = get.input.network.folder.path(n, l0, d, prop.mispl, prop.neg, network.no)
	
	tlog(16, "start to plot transition with exact algorithms")	
	#for(algo.name in cor.clu.exact.algos){
		tlog(20, "plot transition networks => algo.name: ", cor.clu.exact.algo)
		
		for(graph.desc.name in c(SIGNED.UNWEIGHTED.FILE)){ # OPPOSITE.SIGNED.UNWEIGHTED.FILE
			tlog(24, "plot transition networks => graph.desc.name: ", graph.desc.name)
			part.folder = part.folder = get.part.folder.path(n, l0, d, prop.mispl, prop.neg, network.no, cor.clu.exact.algo, graph.desc.name)
			plot.folder = get.transition.graph.plot.folder.path(n, l0, d, prop.mispl, prop.neg, network.no, cor.clu.exact.algo, graph.desc.name)
			clu.analysis.folder = get.clu.analysis.folder.path(n, l0, d, prop.mispl, prop.neg, network.no, cor.clu.exact.algo, graph.desc.name)
			eval.folder = get.eval.folder.path(n, l0, d, prop.mispl, NA, prop.neg, k=ALL, network.no, cor.clu.exact.algo, graph.desc.name)
			
			mems = load.membership.files(part.folder)
			nb.sol = length(mems)
			print(nb.sol)		
			if(nb.sol>1){
						
				for(measure in comp.measures){
					tlog(20, "performing cluster characterization => measure: ", measure)

				    table.file = file.path(eval.folder, paste0(EVAL.BEST.K.FOR.KMEDOIDS,"-",measure,".csv"))
				    #print(table.file)		
				    
				    if(file.exists(table.file)){
				        #print("exist")
				        df = read.csv(table.file, check.names=FALSE, stringsAsFactors=FALSE)
				        best.k = as.numeric(df[,BEST.K.FOR.SILH.COL.NAME])
				        best.silh.score = df[,BEST.SILH.COL.NAME]
				        #silh.folder.name = paste0("k=",k,"_silh=", sprintf("%.4f",best.silh))
				        kmedoids.mem = NA
				        if(best.k<=10)
				            kmedoids.mem = as.numeric(unlist(strsplit(df[,BEST.MEM.FOR.SILH.COL.NAME], ",")))
				    
					
    			        if(dir.exists(net.folder)){
    						if(!dir.exists(plot.folder))
    							dir.create(path=plot.folder, showWarnings=FALSE, recursive=TRUE)
    										
    						tlog(24, "plotting networks in ", plot.folder)		
    						graph.name = paste0(graph.desc.name, ".G")
    						plot.filename = file.path(plot.folder,paste0(graph.desc.name, ".pdf"))
    						
    						
    						#g.name = paste0(SIGNED.UNWEIGHTED.FILE,".G")
    						network.path = file.path(net.folder,graph.name)
    						g = read.graph.ils(network.path)
    						g.title = ""
    						if(!is.na(best.silh.score))
    							g.title = paste0("Best silh score = ",best.silh.score)
    						
    						
    						# if one of the graphs already exists, do not generate all
    						if(!file.exists(file.path(plot.folder,paste0(MBRSHP.FILE.PREFIX,"0.png"))) || FORCE){
    							tlog("generate graph images")
    							save.graph.images.with.bg.color(g, part.folder, nb.sol, plot.folder, bg.color.membership=kmedoids.mem)
    						} else {
    							tlog("no need to generate graph images")
    						}
    	
    						# -----------------------------------------------
    						# -----------------------------------------------
    			
    						edit.dist.matrix = NA
    						measure = EDIT
    						mtrx.file = file.path(eval.folder,paste0(EVAL.DIST.MATRIX.FILE.PREFIX,"-",measure,".csv"))
    						if(!file.exists(mtrx.file)){
    							tlog("matrix does not exist")
    							edit.dist.matrix = build.edit.dist.matrix(part.folder, nb.sol)		
    							rownames(edit.dist.matrix) = paste("sol",0:(nb.sol-1))
    							colnames(edit.dist.matrix) = paste("sol",0:(nb.sol-1))
    							write.csv(file=mtrx.file,x=edit.dist.matrix)
    						}
    						else{
    							tlog("matrix exists")
    							#print(read.csv(file=mtrx.file))
    							mtrx = as.matrix(read.csv(file=mtrx.file))
    							edit.dist.matrix = mtrx[,2:ncol(mtrx)]
    							
    							rownames(edit.dist.matrix) = mtrx[,1]
    						}
    						
    						
    						# ---------------------------------------------------------------------
    						#  contains only edit distance 1
    						adj.matrix1 = build.adj.matrix.from.edit.dist.matrix(edit.dist.matrix, nb.edit=1)
    						g1 = graph_from_adjacency_matrix(adj.matrix1, mode="undirected", weighted=NULL, diag=FALSE)
    #						write.graph(graph=g1,file=file.path(plot.folder,"transition-graph1.graphml"),"graphml")
    						el1 = get.edgelist(g1)
    						g.title1 = g.title
    						plot.filename1 = "edit1-graph.pdf"
    						# ---------------------------------------------------------------------
    						
    						
    						#  contains only edit distance 2
    						adj.matrix2 = build.adj.matrix.from.edit.dist.matrix(edit.dist.matrix, nb.edit=2)
    						g2 = graph_from_adjacency_matrix(adj.matrix2, mode="undirected", weighted=NULL, diag=FALSE)
    #						write.graph(graph=g2,file=file.path(plot.folder,"transition-graph2.graphml"),"graphml")
    						el2 = get.edgelist(g2)
    						
    						#  contains only edit distance 3
    						adj.matrix3 = build.adj.matrix.from.edit.dist.matrix(edit.dist.matrix, nb.edit=3)
    						g3 = graph_from_adjacency_matrix(adj.matrix3, mode="undirected", weighted=NULL, diag=FALSE)
    #						write.graph(graph=g3,file=file.path(plot.folder,"transition-graph3.graphml"),"graphml")
    						el3 = get.edgelist(g3)						
    							
    						#  contains only edit distance 4
    						adj.matrix4 = build.adj.matrix.from.edit.dist.matrix(edit.dist.matrix, nb.edit=4)
    						g4 = graph_from_adjacency_matrix(adj.matrix4, mode="undirected", weighted=NULL, diag=FALSE)
    #						write.graph(graph=g4,file=file.path(plot.folder,"transition-graph4.graphml"),"graphml")
    						el4 = get.edgelist(g4)	
    						
    						
    						# ---------------------------------------------------------------------
    						adj.matrix12 = adj.matrix1 + adj.matrix2
    						g12 = graph_from_adjacency_matrix(adj.matrix12, mode="undirected", weighted=NULL, diag=FALSE)
    						el12 = get.edgelist(g12)
    						#write.graph(graph=g12,file=file.path(plot.folder,"transition-graph12.graphml"),"graphml")
    						el12 = get.edgelist(g12)
    						g.title12 = g.title
    						plot.filename12 = "edit1_2-graph.pdf"
    						e1.indx = c()
    						e2.indx = c()
    						if(length(el1)>0)
    							e1.indx = sapply(1:nrow(el1), function(i){ get.edge.ids(g12,el1[i,])})
    						if(length(el2)>0)
    							e2.indx = sapply(1:nrow(el2), function(i){ get.edge.ids(g12,el2[i,])})
    						cols = rep(NA,nrow(el12))
    
    						cols[e1.indx] = EDIT1.EDGE.COLOR
    						cols[e2.indx] = EDIT2.EDGE.COLOR
    						E(g12)$color = cols
    						# ---------------------------------------------------------------------
    
    						# ---------------------------------------------------------------------
    						adj.matrix123 = adj.matrix1 + adj.matrix2 + adj.matrix3
    						g123 = graph_from_adjacency_matrix(adj.matrix123, mode="undirected", weighted=NULL, diag=FALSE)
    						el123 = get.edgelist(g123)
    						#write.graph(graph=g123,file=file.path(plot.folder,"transition-graph123.graphml"),"graphml")
    						el123 = get.edgelist(g123)
    						g.title123 = g.title
    						plot.filename123 = "edit1_2_3-graph.pdf"
    						e1.indx = c()
    						e2.indx = c()
    						e3.indx = c()
    						if(length(el1)>0)
    							e1.indx = sapply(1:nrow(el1), function(i){ get.edge.ids(g123,el1[i,])})
    						if(length(el2)>0)
    							e2.indx = sapply(1:nrow(el2), function(i){ get.edge.ids(g123,el2[i,])})
    						if(length(el3)>0)
    							e3.indx = sapply(1:nrow(el3), function(i){ get.edge.ids(g123,el3[i,])})
    
    						cols[e1.indx] = EDIT1.EDGE.COLOR
    						cols[e2.indx] = EDIT2.EDGE.COLOR
    						cols[e3.indx] = EDIT3.EDGE.COLOR
    						E(g123)$color = cols
    						# ---------------------------------------------------------------------
    
    
    						# ---------------------------------------------------------------------
    						adj.matrix1234 = adj.matrix1 + adj.matrix2 + adj.matrix3 + adj.matrix4
    						g1234 = graph_from_adjacency_matrix(adj.matrix1234, mode="undirected", weighted=NULL, diag=FALSE)
    						el1234 = get.edgelist(g1234)
    						#write.graph(graph=g1234,file=file.path(plot.folder,"transition-graph1234.graphml"),"graphml")
    						#print(components(g12)$membership)
    						el1234 = get.edgelist(g1234)
    						g.title1234 = g.title
    						plot.filename1234 = "edit1_2_3_4-graph.pdf"
    						e1.indx = c()
    						e2.indx = c()
    						e3.indx = c()
    						e4.indx = c()
    						if(length(el1)>0)
    							e1.indx = sapply(1:nrow(el1), function(i){ get.edge.ids(g1234,el1[i,])})
    						if(length(el2)>0)
    							e2.indx = sapply(1:nrow(el2), function(i){ get.edge.ids(g1234,el2[i,])})
    						if(length(el3)>0)
    							e3.indx = sapply(1:nrow(el3), function(i){ get.edge.ids(g1234,el3[i,])})
    						if(length(el4)>0)
    							e4.indx = sapply(1:nrow(el4), function(i){ get.edge.ids(g1234,el4[i,])})
    						cols = rep(NA,nrow(el1234))
    
    						cols[e1.indx] = EDIT1.EDGE.COLOR
    						cols[e2.indx] = EDIT2.EDGE.COLOR
    						cols[e3.indx] = EDIT3.EDGE.COLOR
    						cols[e4.indx] = EDIT4.EDGE.COLOR
    						E(g1234)$color = cols
    						# ---------------------------------------------------------------------
    
    
    						# ------------------------------------
    						# prepare layout: higher edit distance is, lower edge weight is.
    						weighted.adj.matrix1234 = 4*adj.matrix1 + 3*adj.matrix2 + EDIT3.EDGE.WEIGHT*adj.matrix3 + EDIT4.EDGE.WEIGHT*adj.matrix4
    						#print(weighted.adj.matrix1234)
    						weighted.g1234 = graph_from_adjacency_matrix(weighted.adj.matrix1234, mode="undirected", weighted=TRUE, diag=FALSE)
    						lyt <- layout.fruchterman.reingold(graph=weighted.g1234)
    						# ------------------------------------
    
    
    
    #						if(nb.sol<100)
    						V(g1)$x <- lyt[,1]
    						V(g1)$y <- lyt[,2]
    						E(g1)$color = EDIT1.EDGE.COLOR
    						E(g1)$weight = EDIT1.EDGE.WEIGHT
    						comp.mem1 = components(g1)$membership
    						diam.vals1 = get.diam.vals(g1, edit.dist.matrix)
    						g1 = indicate.edit.dist.betw.components(g1, edit.dist.matrix)
    						lgnd = list(Edit1=list(weight=EDIT1.EDGE.WEIGHT,color=EDIT1.EDGE.COLOR))
    						plot.transition.graph(g1, g.title1, plot.filename1, part.folder, plot.folder, nb.sol,
    							diam.vals1, comp.mem1, lgnd)


    						weighted.adj.matrix12 = EDIT1.EDGE.WEIGHT*adj.matrix1 + EDIT2.EDGE.WEIGHT*adj.matrix2
    						weighted.g12 = graph_from_adjacency_matrix(weighted.adj.matrix12, mode="undirected", weighted=TRUE, diag=FALSE)
    						V(weighted.g12)$x <- lyt[,1]
    						V(weighted.g12)$y <- lyt[,2]
    						if(!is.null(E(g12)$color))
    							E(weighted.g12)$color = E(g12)$color
    						comp.mem12 = components(weighted.g12)$membership
    						diam.vals12 = get.diam.vals(weighted.g12, edit.dist.matrix)
    						weighted.g12 = indicate.edit.dist.betw.components(weighted.g12, edit.dist.matrix)
    						lgnd = list(Edit1=list(weight=EDIT1.EDGE.WEIGHT,color=EDIT1.EDGE.COLOR),
    							Edit2=list(weight=EDIT2.EDGE.WEIGHT,color=EDIT2.EDGE.COLOR))
    						plot.transition.graph(weighted.g12, g.title12, plot.filename12, part.folder, plot.folder, nb.sol,
    							diam.vals12, comp.mem12, lgnd)


    						weighted.adj.matrix123 = EDIT1.EDGE.WEIGHT*adj.matrix1 + EDIT2.EDGE.WEIGHT*adj.matrix2 + EDIT3.EDGE.WEIGHT*adj.matrix3
    						weighted.g123 = graph_from_adjacency_matrix(weighted.adj.matrix123, mode="undirected", weighted=TRUE, diag=FALSE)
    						V(weighted.g123)$x <- lyt[,1]
    						V(weighted.g123)$y <- lyt[,2]
    						if(!is.null(E(g123)$color))
    							E(weighted.g123)$color = E(g123)$color
    						comp.mem123 = components(weighted.g123)$membership
    						diam.vals123 = get.diam.vals(weighted.g123, edit.dist.matrix)
    						weighted.g123 = indicate.edit.dist.betw.components(weighted.g123, edit.dist.matrix)
    						lgnd = list(Edit1=list(weight=EDIT1.EDGE.WEIGHT,color=EDIT1.EDGE.COLOR),
    							Edit2=list(weight=EDIT2.EDGE.WEIGHT,color=EDIT2.EDGE.COLOR),
    							Edit3=list(weight=EDIT3.EDGE.WEIGHT,color=EDIT3.EDGE.COLOR))
    						plot.transition.graph(weighted.g123, g.title123, plot.filename123, part.folder, plot.folder, nb.sol,
    							diam.vals123, comp.mem123, lgnd)



    						V(weighted.g1234)$x <- lyt[,1]
    						V(weighted.g1234)$y <- lyt[,2]
    						if(!is.null(E(g1234)$color))
    							E(weighted.g1234)$color = E(g1234)$color
    						comp.mem1234 = components(weighted.g1234)$membership
    						diam.vals1234 = get.diam.vals(weighted.g1234, edit.dist.matrix)
    						weighted.g1234 = indicate.edit.dist.betw.components(weighted.g1234, edit.dist.matrix)
    #						if(nb.sol<100)
    						weighted.g1234 = indicate.edit.dist.betw.components(weighted.g1234, edit.dist.matrix)
    						lgnd = list(Edit1=list(weight=EDIT1.EDGE.WEIGHT,color=EDIT1.EDGE.COLOR),
    							Edit2=list(weight=EDIT2.EDGE.WEIGHT,color=EDIT2.EDGE.COLOR),
    							Edit3=list(weight=EDIT3.EDGE.WEIGHT,color=EDIT3.EDGE.COLOR),
    							Edit4=list(weight=EDIT4.EDGE.WEIGHT,color=EDIT4.EDGE.COLOR))
    						plot.transition.graph(weighted.g1234, g.title1234, plot.filename1234, part.folder, plot.folder, nb.sol,
    							diam.vals1234, comp.mem1234, lgnd)

    							
    							
    							
    #						curr.line = c(n, k, d, prop.mispl, prop.neg, network.no, 
    #										paste(components(g1)$membership,collapse=","),
    #										paste(components(g2)$membership,collapse=","),
    #										paste(components(g12)$membership,collapse=",")
    #										)
    #						conn.comp.data = rbind(conn.comp.data, curr.line)
    	
    						# -----------------------------------------------
    						# -----------------------------------------------
    			        }
				        
				    }
					
				}
			}
		}
	#}
	
#	return(conn.comp.data)
}




#################################################################
#
# graph.sizes
# d: density
# l0: nb cluster
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
plot.all.transition.graphs = function(graph.sizes, d, l0, prop.mispls, prop.negs, in.rand.net.folders,
		cor.clu.exact.algo, comp.measures, force, plot.formats)
{
#	all.conn.comp.data = c()
	
	tlog("starts partitioning networks")
	for(n in graph.sizes){
		tlog(8, "partitioning networks => n: ", n)
		
		for(prop.mispl in prop.mispls){
			tlog(8, "partitioning networks => prop.mispl: ", prop.mispl)
			
		    if(is.na(prop.negs) && d == 1){
		        prop.negs = compute.prop.neg(n, d, l0, prop.mispl)
		    }
			
			for(prop.neg in prop.negs){
				tlog(12, "partitioning networks => prop.neg: ", prop.neg)
				
				for(network.no in in.rand.net.folders){
					tlog(16, "partitioning networks => network.no: ", network.no)
					
					build.and.plot.transition.graph(n, l0, d, prop.mispl, prop.neg, network.no,
					                                cor.clu.exact.algo, comp.measures, force, plot.formats)
#					conn.comp.data = build.and.plot.transition.graph(n, l0, d, prop.mispl, prop.neg, network.no,
#							cor.clu.exact.algo, comp.measures, force, plot.formats)
#						all.conn.comp.data = rbind(all.conn.comp.data, conn.comp.data)
				}
				
			}
			
		}
		
	}
	
#	if(length(all.conn.comp.data)>0){
#		colnames(all.conn.comp.data) = c("n", "k", "d", "prop.mispl", "prop.neg", "network.no", "compsEdit1", "compsEdit2", "compsEdit1&2")
#		rownames(all.conn.comp.data) = NULL
#		write.csv(file="transition-graph-conn-compenents-stats-k4.csv", x=all.conn.comp.data)
#	}
}
