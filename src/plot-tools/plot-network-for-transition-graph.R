#############################################################################################
# Functions used to plot signed networks, including their node partition if any.
# 
# 11/2015 Vincent Labatut
#############################################################################################
library("igraph")
library("expm")




#############################################################################################
# Colors used in the plots.
# Taken from http://colorbrewer2.org/#type=qualitative&scheme=Set1&n=9
#############################################################################################
COLORS <- c(
		rgb(228,26,28,maxColorValue=255),
		rgb(55,126,184,maxColorValue=255),
		rgb(77,175,74,maxColorValue=255),
		rgb(152,78,163,maxColorValue=255),
		rgb(255,127,0,maxColorValue=255),
		rgb(255,255,51,maxColorValue=255),
		rgb(166,86,40,maxColorValue=255),
		rgb(247,129,191,maxColorValue=255),
		rgb(153,153,153,maxColorValue=255),
		rgb(0,0,0,maxColorValue=255),
		rgb(50,50,50,maxColorValue=255)
)




#############################################################################################
# Plot the specified signed graph, generating files of the specified formats.
#
# g: signed graph to plot.
# plot.file: base name (and path) of the generated files.
# format: format(s) to handle (among "PDF", "PNG", and NA for screen).
# method: layout used. One supported by layout.signed.laplacian, or one for 
#		  unsigned graphs, such as kamada.kawai, fruchterman.reingold.
#
# returns: the same graph, with the spatial positions stored as nodal attributes x and y.
#############################################################################################
plot.network.for.transition.graph <- function(g, membership=NA, plot.file, format=c("PDF","PNG",NA), method="bn_zheng", plot.title="", bg.color=NA)
{	
	
	# ----------------------------------------
	# preprocessing: replace negative values by 0 in membership vector
	indx = which(membership < 0)
	membership[indx] = NA
	# ----------------------------------------
	
	# setup node parameters
	vertex.sizes <- 30
	if(vcount(g)>20 && vcount(g)<=28)
		vertex.sizes <- 25
	if(vcount(g)>28 && vcount(g)<=36)
		vertex.sizes <- 20
	
#	vertex.label <- NA
	vertex.label <- V(g)$id
	
	V(g)$label.cex = 8
	if(vcount(g)>24 && vcount(g)<=28)
		V(g)$label.cex = 6
	if(vcount(g)>28 && vcount(g)<=36)
		V(g)$label.cex = 4
	
	# setup link parameters
	if(ecount(g)>0)
	{	edge.colors <- rep(NA,ecount(g))
		edge.colors[E(g)$weight<0] <- adjustcolor("RED", alpha.f=0.5)
		edge.colors[E(g)$weight>0] <- adjustcolor("GREEN", alpha.f=0.5)
		edge.widths <- abs(E(g)$weight)*2
	}
	
	# set up node colors
	if(all(is.na(membership)))
#		vertex.colors <- "SkyBlue2" # default igraph color
		vertex.colors <- "white" # default igraph color
	else
	{	# palette <- rainbow(length(unique(membership)))
		vertex.colors <- COLORS[membership]
	}
	
	# setup layout (only if the graph doesn't already have one)
	v.att <- list.vertex.attributes(graph=g)
	if(!("x" %in% v.att) | !("y" %in% v.att)) 
	{	# layouts for unsigned graphs
		if(method=="kamada.kawai")
		{	gpos <- delete.edges(graph=g,edges=which(E(g)$weight<0))
			lyt <- layout.kamada.kawai(graph=gpos)
		}
		else if(method=="fruchterman.reingold")
		{	gpos <- delete.edges(graph=g,edges=which(E(g)$weight<0))
			lyt <- layout.fruchterman.reingold(graph=gpos)
		}
		# layouts for signed graphs
		else if(method=="bn_zheng")
		{	lyt <- layout.signed.laplacian(g, method)
		}
		else
			lyt <- layout.circle(graph=g)
		
		# store spatial positions as nodal attributes
		V(g)$x <- lyt[,1]
		V(g)$y <- lyt[,2]
	}
	else{
		lyt <- cbind(V(g)$x,V(g)$y)
	}
	
	
	# to save memory, remove links
	g2 <- delete.edges(graph=g,edges=E(g))
	
	# process each specified format
	for(frmt in format)
	{	# create the file
		if(!is.na(frmt))
		{	# set plot file name
			plot.filename <- plot.file
			if(toupper(substr(plot.filename, nchar(plot.filename)-2, nchar(plot.filename)))!=toupper(frmt))
				plot.filename <- paste0(plot.filename ,".",frmt)
			# handle format
			if(frmt=="PNG") # USE ONLY THIS, it is better in terms of image quality
			{	bg = "transparent"
				if(!is.na(bg.color))
					bg = bg.color
				png(filename=plot.filename,width=300,height=300,units="px",pointsize=5,bg=bg)
				par(mar=c(0,0,0,0),xpd = NA)
			}
			if(frmt=="JPEG")
			{	print("jpg")
#				png(filename=plot.filename,width=800,height=800,units="px",pointsize=20,bg="white")
				bg = bg.color
				jpeg(filename=plot.filename,width=500,height=500,units="px",pointsize=20,bg=bg,quality="70")
				par(mar=c(0,0,0,0),xpd = NA)
			}
			else if(frmt=="PDF")
			{	pdf(file=plot.filename,bg="white")
			}
		}
		
		# create the plot
		plot(g2, layout=lyt,
				vertex.size=vertex.sizes, vertex.label=vertex.label, vertex.color=vertex.colors,
		)
		
		title(plot.title, cex.main=0.5)
		
		# finalize plot file
		if(!is.na(frmt))
			dev.off()
	}
	
	return(g)
}




#############################################################################################
### Test using the article toy graph
#g2 <-graph.empty(6,directed=FALSE)
#g2 <- add_edges(g2,edges=c(1,2,2,4,2,6,3,4,4,6,5,6))
#E(g2)$weight <- c(1,-1,-1,1,-1,1)
#plot.network(g2, membership=membership, plot.file="sqdq", format=NA, method="bn_zheng")