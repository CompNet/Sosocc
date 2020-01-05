#source("utils.R")
#source("plot-network.R")
#source("create-relative-plot-membership.R")
#source("define-purity.R")
#library(igraph) 
library(png)

# be aware that there is a warning about hamming.distance() method
# link: https://stackoverflow.com/a/16571702
library(e1071)


BG.COLORS = c(
		adjustcolor("seagreen2", alpha.f = 0.3),
		adjustcolor("darkgoldenrod1", alpha.f = 0.3),
		adjustcolor("darkmagenta", alpha.f = 0.3),
		adjustcolor("red", alpha.f = 0.3),
		adjustcolor("blue", alpha.f = 0.2),
		adjustcolor("green", alpha.f = 0.3),
		adjustcolor("yellow", alpha.f = 0.3),
		adjustcolor("gray47", alpha.f = 0.3),
		adjustcolor("black", alpha.f = 0.4),
		adjustcolor("olivedrab1", alpha.f = 0.5)
)



###################################################################
#
###################################################################
save.graph.images.with.bg.color = function(g, in.dir, nb.sol, out.dir, bg.color.membership=NA){
	for(i in 0:(nb.sol-1)){
		filename = paste0(MBRSHP.FILE.PREFIX,i)
		mem = read.table(file=file.path(in.dir,paste0(filename,".txt")))$V1
		plot.file = file.path(out.dir,paste0(filename,".png"))
		
		bg.col = "transparent"
		if(!any(is.na(bg.color.membership)))
			bg.col=BG.COLORS[bg.color.membership[i+1]]
		plot.network.for.transition.graph(g, mem, plot.file, format="PNG", method="circular", plot.title="", bg.color=bg.col)
	}
}

###################################################################
#
# source: https://stackoverflow.com/questions/16932442/r-display-image-at-exact-node-in-igraph
###################################################################
plot.transition.graph = function(g, g.title, plot.filename, in.dir, out.dir, nb.sol, diam.vals=NA, comp.mem = NA, lgnd = NA){
	
	vertex.labels = rep(NA,nb.sol)
	for(i in 0:(nb.sol-1)){
		filename = paste0(MBRSHP.FILE.PREFIX,i)
		mem = read.table(file=file.path(in.dir,paste0(filename,".txt")))$V1
		vertex.labels[i+1] = paste0("sol",i,"\nk",length(unique(mem)))
	}
	
	rasters = list()
	for(i in 0:(nb.sol-1)){
		filename = paste0(MBRSHP.FILE.PREFIX,i)
		plot.file = file.path(out.dir,paste0(filename,".png"))
		rasters[[as.character(i)]] <- readPNG(plot.file, native=TRUE)
	}

	V(g)$raster <- rasters
	
	lyt = NA
	v.att <- list.vertex.attributes(graph=g)
	if(!("x" %in% v.att) | !("y" %in% v.att)){
		lyt <- layout.fruchterman.reingold(graph=g)
	} else
		lyt = cbind(V(g)$x,V(g)$y)
	
	l.cex = 0.6
	v.size = 15
	if(nb.sol<50){
		v.size = 25
		l.cex = 0.8
	}
	
	V(g)$label.cex = l.cex
	pdf(file.path(out.dir,plot.filename)) # width=15,height=10
	#par(mar=c(0,0,0,0)+.1)
	plot(g, main=g.title, layout=lyt, vertex.shape="raster", vertex.label=vertex.labels, edge.width=E(g)$weight, margin=.1,
			vertex.size=v.size, vertex.size2=v.size, vertex.label.dist=0, vertex.label.degree=0)
		
	if(!is.na(diam.vals)){
		norm.lyt = layout.norm(lyt)
		k = length(unique(comp.mem))
		for(i in 1:k){
			indxs = which(comp.mem == i)
			text(x = min(norm.lyt[indxs,1])-0.1, y = min(norm.lyt[indxs,2])-0.1, labels=paste0("Diam=",diam.vals[i]), cex=0.8, col="red")
#			text(x = mean(norm.lyt[indxs,1]), y = mean(norm.lyt[indxs,2]), labels=paste0("Diam=",diam.vals[i]), cex=0.8, col="red")
		}
	}
	
	# -----------------------------------------------------------
	# legend
	if(!is.na(lgnd)){ # if there is edit distance edges (other than edge between connected components)
		weights = sapply(lgnd, function(item) item$weight)
		colors = sapply(lgnd, function(item) item$col)
		legend("topright", legend=names(lgnd), lwd=weights,col=colors, cex=0.7)
	}
	# -----------------------------------------------------------
	dev.off()
}

