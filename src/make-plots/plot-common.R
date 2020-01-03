
# we have 2 color lists:
# 1) distinct color list (that we can distinguish clearly). There are 14 colors. This is main color list
# 2) yellow-orange-red colors in continuity. There are 10 colors. This is mainly used for frequency indication in scatter plot




####################
# 1st color list
####################
 
palette <- as.list(rainbow(14))
palette[1:4] = rainbow(4, v=0.7)
# ==============
palette[[5]] = "#E3CF57" # banana yellow
palette[[6]] = "maroon1" #
palette[[7]] = "darkviolet" # 
palette[[8]] = "plum1" #
palette[[9]] = "seagreen4" #
palette[[10]] = "black" # 
palette[[11]] = "burlywood4" # 
palette[[12]] = "yellow" #
palette[[13]] = "lightseagreen"
palette[[14]] = "gray87" # light gray
# ==============




####################
# 2nd color list
####################
library("RColorBrewer") # we will use this color palette: "OrRd" (there are 9 colors)
nb.color = 10
# color order: blue, green, yellow, orange, red
heat.palette = c("#2E64FE", "#2ECCFA", "#2EFEC8", "#2EFE64", "#3ADF00", "#F7FE2E", "#FACC2E", "#FE9A2E", "#FF4000", "#B40404")
#heat.palette = c("#F4FA58", "#F7D358", brewer.pal(n = nb.color, name = "OrRd")[2:9]) # add




# ------------------------------------------



################################################################
# It shortens the name of a given result type if it is too long
#       (so, depending on the type, not depending on the length, since we know which one is lengthy in advance)
#
# An example: from 'Optimality−analysis−Nmi−ExCC−all_vs_GRASP−CC_l1_a0.8_g0_t3600_i400_n1' to 'Optimality−analysis'
#
# result.type: result type
#
################################################################
shorthen.result.type.name.if.possible = function(result.type){
	if(result.type == EXEC.TIME)
		return(result.type)
	else if(result.type == NB.CLUSTER)
		return(result.type)
	else if(result.type == NB.SOLUTION)
		return(result.type)
	else if(result.type == IMB.OPT)
		return(result.type)
	else if( startsWith(result.type, DIST.SCORE) )
		return(DIST.SCORE)
}


################################################################
# It allows to preapre axis label based on result type.
#   The possible values of the variable 'result.type' are defined in 'src/define-constants.R'
#
# result.type
#
################################################################
prepare.axis.label = function(result.type){
	if(result.type == EXEC.TIME)
		return(EXEC.TIME.LABEL.NAME)
	else if(result.type == NB.CLUSTER)
		return(NB.CLUSTER.LABEL.NAME)
	else if(result.type == NB.SOLUTION)
		return(NB.SOLUTION.LABEL.NAME)
	else if(result.type == IMB.OPT)
		return(IMB.OPT.LABEL.NAME)
	else if( startsWith(result.type, DIST.SCORE) ){
		measure = gsub(paste0(DIST.SCORE,"-"),"", result.type) # ex: from 'Distance-score-Nmi' to 'Nmi'
		return(paste0(DIST.SCORE.LABEL.NAME, " with ", measure))
	}
	else if(startsWith(result.type, MEAN.DIST.SCORE)){
		measure = gsub(paste0(MEAN.DIST.SCORE,"-"),"", result.type) # ex: from 'Distance-score-Nmi' to 'Nmi'
		return(paste0(MEAN.DIST.SCORE.LABEL.NAME, " with ", measure))
	} else
		return(result.type)
}


# 
# ################################################################
# # It allows to prepare a plot title.
# #
# # n
# # prop.mispl
# # detected.imb.interval
# # prop.neg
# # network.no
# # algo.name
# # graph.desc.name
# #
# ################################################################
# prepare.plot.title = function(n=NA, prop.mispl=NA, detected.imb.interval=NA, prop.neg=NA, network.no=NA, algo.name=NA, graph.desc.name=NA){
# 	
# 	title = ""
# 	if(!is.na(graph.desc.name))
# 		title = paste0(title, "type=", graph.desc.name)
# 	if(!is.na(n))
# 		title = paste0(title, " n=", n)
# 	
# 	if(!is.na(prop.mispl) && is.character(prop.mispl))
# 		title = paste0(title, " prop.mispl=", prop.mispl)
# 	else if(!is.na(prop.mispl) && is.numeric(prop.mispl))
# 		title = paste0(title, " prop.mispl=", sprintf("%.3f",prop.mispl))
# 	
# 	if(!is.na(detected.imb.interval) && is.character(detected.imb.interval))
# 	    title = paste0(title, " detected.imb.interval=", detected.imb.interval)
# 	
# 	if(!is.na(prop.neg) && is.character(prop.neg))
# 		title = paste0(title, " prop.neg=", prop.neg)
# 	else if(!is.na(prop.neg) && is.numeric(prop.neg))
# 		title = paste0(title, " prop.neg=", sprintf("%.4f",prop.neg))
# 	title = paste0(title, "\n")
# 	
# 	if(!is.na(network.no))
# 		title = paste0(title, " network.no=", network.no)
# 	if(!is.na(algo.name))
# 		title = paste0(title, " method=", algo.name)
# 
# 
# 	return(title)
# }






################################################################
# It allows to draw a scatter plot with two different use cases.
#
# Use case 1: Basic scatterplot. What is important is 'y.values'. So, 'x.values' is only sequence values. 
#           Here, we do not really need an automatic legend
#   make.scatter.plot(y.values=c(1,3,1,4,1,6), x.values=c(1,2,3,4,5,6), title="title", xlabel="xlab", ylabel="ylab",
#           NA, NA, JUST.PLOT, ylimit=NA, xlimit=NA, max.nb.elt=NA, enable.legend=FALSE) ==> the point (3,1) appears as green dot
#
# Use case 2: x and y values are separate vectors. This is related to: https://stackoverflow.com/questions/15828058/add-legend-with-color-and-range-in-r
#               So, if there are repeated points, they will appear in different color.
#   make.scatter.plot(y.values=c(1,1,1,4,1,6), x.values=c(3,3,3,9,3,2), title="title", xlabel="xlab", ylabel="ylab",
#           NA, NA, JUST.PLOT, ylimit=NA, xlimit=NA, max.nb.elt=NA, enable.legend=FALSE) ==> the point (3,1) appears as green dot
#
#
# y.values
# x.values
# title
# xlabel
# ylabel
# plot.folder
# plot.filename
# plot.format
# ylimit
# xlimit
# max.nb.elt: I am not sure why I introduced this parameter. User can give NA.
# enable.legend
#
################################################################
make.scatter.plot = function(y.values, x.values, title, xlabel, ylabel, plot.folder, plot.filename, plot.format, ylimit=NA, xlimit=NA,
                             max.nb.elt, enable.legend){
    
    # color order for 'heat.palette': blue, green, yellow, orange, red
    
    if(is.na(ylimit))
        ylimit = NULL
    if(is.na(xlimit))
        xlimit = NULL
    
    #print("make.scatter.plot")
    
    #	# =================================================
    #	# hex binning ==> intreresting but it does not work with par() for layout plots, because the plot type is glob/gtable etc.
    #	library(gridExtra)
    #	library(hexbin)
    #	library("RColorBrewer")
    #	plotList <- lapply(1:4, function(i) {
    #				VAR1 = rnorm(i*1000)  
    #				VAR2 = rnorm(i*1000)  
    #				hexbinplot(VAR1 ~ VAR2, colorkey=TRUE, mincnt=1, maxcnt=50, colorcut = seq(0, 1, length = 8), colramp = function(n) { brewer.pal(n = n, name = "OrRd") })
    #			})
    #	
    #	do.call(grid.arrange, c(plotList, ncol=2))
    #	# =================================================
    
    if(plot.format == PLOT.AS.PDF){
        if(!dir.exists(plot.folder))
            dir.create(path=plot.folder, showWarnings=FALSE, recursive=TRUE)
        pdf(file=file.path(plot.folder, plot.filename), compress=FALSE)
    }
    
    nb.elt = length(y.values)
    
    # interesting: https://stackoverflow.com/questions/15828058/add-legend-with-color-and-range-in-r	
    if(!is.na(max.nb.elt))
        nb.elt = max.nb.elt
    
    # -------------------------------------------------
    # compute frequency intervals based on nb elt
    nb.color = length(heat.palette)
    
    freq.intervals = NA
    if(nb.elt == 1)
        freq.intervals = factor(nb.elt)
    else if(nb.elt < nb.color)
        freq.intervals = cut(seq(1, nb.elt), breaks = nb.elt) # example: (0.995,1.83] (1.83,2.67] (2.67,3.5] (3.5,4.33] ... (5.17,6]
    else # if(length(nb.elt) > 1)
        freq.intervals = cut(seq(1, nb.elt), breaks = nb.color, dig.lab=2)
    # -------------------------------------------------
    
    freqs = table(paste0(x.values, "-", y.values)) # if the point (3,7) frequently appears, it will be colored differently than the others.
    # example:
    # 2-2     2-6      3-1         9-4 
    # 1        1        3           1
    unique.pairs = names(freqs)
    
    interval.vals = freq.intervals[freqs] # existing intervals. Example: (0.995,1.83] (0.995,1.83] (2.67,3.5]   (0.995,1.83]
    
    
    # sort curr intervals in increasing order, thus strong
    sort.indx = order(interval.vals) 
    color.indx = as.numeric(interval.vals) # levels
    b <- heat.palette[color.indx]
    
    legend.intervals = sort(unique(interval.vals))
    legend.col = heat.palette[as.numeric(legend.intervals)]
    
    
    
    # reconstruct x and y values based on only distinct (x,y) pairs ==> "12-23", "4-12" then new.x.values=[12,4] and new.y.values=[23,12]
    new.x.values = as.numeric(sapply(unique.pairs, function(str.pair) unlist(strsplit(str.pair, split="-"))[1]))
    new.y.values = as.numeric(sapply(unique.pairs, function(str.pair) unlist(strsplit(str.pair, split="-"))[2]))	
    
    
    print(new.x.values)
    print(new.y.values)
    plot(x=new.x.values, y=new.y.values, col=b, pch=16, xlim=xlimit, ylim=ylimit, las=2, main=title, xlab=xlabel, ylab=ylabel)
    
    if(enable.legend)
        legend("bottomright",legend=legend.intervals,col=legend.col,pch=16,cex=0.4)
    
    
    if(plot.format == PLOT.AS.PDF)
        dev.off()
}




################################################################
# It allows to draw a single plot.
#
# Example:
# make.line.plot(y.values=c(1,2,4), x.values=c(2,6,9), title="title", xlabel="x lab", ylabel="y lab", NA, NA, plot.format=JUST.PLOT, increasing.order=T)
#
# y.values
# x.values
# title
# xlabel
# ylabel
# plot.folder
# plot.filename
# plot.format
# increasing.order
# y.axis.log.scale
# ylimit
#
################################################################
make.line.plot = function(y.values, x.values, title, xlabel, ylabel, plot.folder, plot.filename, plot.format, increasing.order, 
                          y.axis.log.scale=FALSE, ylimit=NA){
    
    # print("make.line.plot")
    log.scale = ""
    if(y.axis.log.scale)
        log.scale = "y"
    
    if(is.null(x.values))
        x.values = seq(1, length(y.values))
    
    if(increasing.order){
        ordr = order(y.values)
        y.values = y.values[ordr]
        xlabel = paste0(xlabel, " in increasing order")
    }
    
    if(plot.format == PLOT.AS.PDF){
        if(!dir.exists(plot.folder))
            dir.create(path=plot.folder, showWarnings=FALSE, recursive=TRUE)
        pdf(file=file.path(plot.folder, plot.filename), compress=FALSE)
    }
    
    if(is.na(ylimit))
        plot(x=x.values, y=y.values, log=log.scale, main=title, col=palette[[1]], type="o", lwd=2, xlab=xlabel, ylab=ylabel)
    else
        plot(x=x.values, y=y.values, log=log.scale, main=title, ylim=ylimit, col=palette[[1]], type="o", lwd=2, xlab=xlabel, ylab=ylabel)
    
    if(plot.format == PLOT.AS.PDF)
        dev.off()
}


################################################################
# It allows to draw multi lines in the same plot.
#
# Example:
# vector.list = list(c(1,2,3), c(2,3,4))
# make.multi.line.plot(vector.list, xvalues=c(2,6,9), xlabel="x lab", ylabel="y lab",
#           ylimit=NA, title="sdf", NA, NA, JUST.PLOT, TRUE, legend.lines=c("A","B"), legend.title="legend title")
#
# vector.list: a list of vectors. Each vector, i.e. y, corresponds to a line plot
# xvalues: x values in the plot. It should be of the same length as the length of y (i.e. length of each vector)
# xlabel: x label
# ylabel: y label
# ylimit: y limit. A vector of two values consisting of lower and upper value of y-axis
# title: title of the plot
# plot.folder: folder where the plot wil be recorded. If plot.format=JUST.PLOT, then user can give NA
# plot.filename: plot filename. If plot.format=JUST.PLOT, then user can give NA
# plot.format: one of the plot formats. PLOT.AS.PDF, PLOT.AS.JPEG, PLOT.AS.PNG or JUST.PLOT (only displaying)
# enable.legend
# legend.lines: a string vector of A values where A is the number of lines used in the plot => desciption of each line
# legend.title
#
################################################################
make.multi.line.plot = function(vector.list, xvalues, xlabel, ylabel, ylimit=NA, title, plot.folder, plot.filename,
                                plot.format, enable.legend, legend.lines=NA, legend.title=NA)
{
    
    # find y limits
    if(is.na(ylimit)){
        all.y.values = unlist(vector.list)
        ylimit = c(min(all.y.values), max(all.y.values))
    }
    
    if(is.null(xvalues))
        xvalues = seq(1,  length(vector.list[[1]])) # seq(1, nb.x.elt)
    
    nb.y.elt = length(vector.list)
    
    if(plot.format == PLOT.AS.PDF){
        if(!dir.exists(plot.folder))
            dir.create(path=plot.folder, showWarnings=FALSE, recursive=TRUE)
        pdf(file=file.path(plot.folder, plot.filename), compress=FALSE)
    }
    
    if(enable.legend)
        # Add extra space to right of plot area; change clipping to figure
        par(mar=c(5.1, 4.1, 4.1, 9.1), xpd=TRUE)
    
    plot(x=xvalues, y=vector.list[[1]], main=title, ylim=ylimit, col=palette[[1]], type="o", lwd=2, xlab=xlabel, ylab=ylabel)
    if(nb.y.elt > 1){ # we need to put this statement just in case
        for(i in 2:nb.y.elt){
            points(x=xvalues, y=vector.list[[i]], col=palette[[i]], type="o", lwd=2)
        }
    }
    
    if(enable.legend)
        legend("topright", inset=c(-0.4,0), title=legend.title, legend=legend.lines, cex=0.8, col=unlist(palette)[1:nb.y.elt], pch=21, lty=1)
    
    if(plot.format == PLOT.AS.PDF)
        dev.off()
}



#################################################################
# It allows to draw an histogram.
#
# values
# title
# xlabel
# plot.folder
# plot.filename
# plot.format
#
#################################################################
make.histogram.plot = function(values, title, xlabel, plot.folder, plot.filename, plot.format)
{ 
    if(plot.format == PLOT.AS.PDF){
        if(!dir.exists(plot.folder))
            dir.create(path=plot.folder, showWarnings=FALSE, recursive=TRUE)
        pdf(file=file.path(plot.folder, plot.filename), compress=FALSE)
    }
    
    hist(x=values, main=title, xlab=xlabel)
    
    if(plot.format == PLOT.AS.PDF)
        dev.off()
}



################################################################################
# ==============================================================================
# ==============================================================================
################################################################################







################################################################################
################################################################################
# BEGIN BAR PLOT ----> we first used all plots through  bar plots. But, then we switched to violin plots. 
#                     Nevertheless, we plot some plots (e.g. single class proportion) with bar plot, though



# source of bar plots: https://www.statmethods.net/graphs/bar.html

## an example of stacked bar plot
## a (and b) is a vector of size 100, contains imb optimality values
## convert opt values into 3 categories: "x<90", "90<=x<100", "x=100"
## each x axis value corresponds to a network
## y axis corresponds to optimality frequency
#a = c(2, 2, 1, 1, 2, 0, 2, 1, 1, 2) # network 1
#b = c(2, 1, 1, 1, 2, 2, 1, 2, 0, 1) # network 2
#counts = cbind(table(a), table(b))
## 				network1 network2
## x<90             1        1
## 90<=x<100        4        5
## x=100            5        4
#rownames(counts) = c("x<90", "90<=x<100", "x=100")
#colnames(counts) = c("network1", "network2")
## colors for row names
#barplot(counts, col=c("red", "orange", "darkblue"), legend=rownames(counts))

# an example of grouped bar plot
# use the previous example, just set the following parameter to TRUE: 'beside'
# barplot(counts, col=c("red", "orange", "darkblue"), legend=rownames(counts), beside=TRUE)





################################################################
# 
#            0.3       0.2
# [1,] 0.2222222 0.1100000
# [2,] 0.2222222 0.6666667
# [3,] 0.6666667 0.6666667
# [4,] 0.6666667 1.0000000
#
# combine.strategy: "cbind"
# ################################################################
make.bar.plot.with.proportion = function(vector.list, ctgry.names, param.values, ylimit=c(0,1),
		plot.folder, plot.filename, plot.format, title, enable.legend=TRUE, legend.title="")
{
	#mtrx = do.call(cbind, vector.list)
	mtrx = t(do.call(rbind, vector.list)) # transpose() is workaround

	rownames(mtrx) = ctgry.names
	colnames(mtrx) = param.values # this is used as x tick labels

	colors = unlist(palette[1:nrow(mtrx)])
	colors <- adjustcolor(colors, alpha.f = 0.7)

	if(plot.format == PLOT.AS.PDF){
		if(!dir.exists(plot.folder))
			dir.create(path=plot.folder, showWarnings=FALSE, recursive=TRUE)
		pdf(file=file.path(plot.folder, plot.filename), compress=FALSE)
	}

	# las=2: makes x labels vertical
	if(enable.legend){
		# Add extra space to right of plot area; change clipping to figure
		par(mar=c(5.1, 4.1, 4.1, 9.1), xpd=TRUE)

		barplot(mtrx, ylim=ylimit, col=colors, legend=ctgry.names, beside=TRUE, las=2,
				main=title, args.legend=list(title=legend.title, x="topright", inset=c(-0.4,0)), border=NA) # omit bar border
	}
	else
		barplot(mtrx, ylim=ylimit, col=colors, beside=TRUE, las=2, main=title, border=NA) # omit bar border

	if(plot.format == PLOT.AS.PDF)
		dev.off()
}



################################################################
#
# count.levels can be provided to use always the same levels in layout plots
# combine.strategy: "rbind" or "cbind"
################################################################
make.bar.plot = function(vector.list, count.levels=NA, row.names, col.names, combine.strategy, is.beside, ylimit=NA,
		title, plot.folder, plot.filename, plot.format, enable.legend=TRUE, legend.title="")
{
	print("make.bar.plot")

	# find count levels
	if(is.na(count.levels)){
		level.values = lapply(vector.list, function(v) names(table(v)))
		count.levels = sort(unique(as.numeric(do.call("c", level.values ))))
	}

	if(is.na(ylimit)){
		max.freqs = sapply(vector.list, function(vals) max(table(vals)))

		# I think it is not meaningful to compute lower bound of ylimit, because if for min value it will show a line (i.e. without height), not a bar
#		min.freqs = sapply(vector.list, function(vals) min(table(vals)))
		ylimit = c(0, max(max.freqs))
	}

	#  compute counts matrix based on count levels. we use factors() for factor levels (i.e. to get the same level set)
	# TODO maybe we can handle more properly legend arg of legend()  ==== i.e. not from rownames(count.matrix)
	count.mtrx = do.call(combine.strategy, lapply(vector.list, function(v) table( factor(v, levels=count.levels))))


#	ind = which(count.mtrx[,2] == 0)
#	count.mtrx[ind,2] = NA

	if(combine.strategy == "rbind"){
		rownames(count.mtrx) = row.names
		colnames(count.mtrx) = count.levels
	} else { # combine.strategy == "cbind"
		rownames(count.mtrx) = count.levels
		colnames(count.mtrx) = col.names
	}

	colors = unlist(palette[1:nrow(count.mtrx)])
	colors <- adjustcolor(colors, alpha.f = 0.7)

	if(plot.format == PLOT.AS.PDF){
		if(!dir.exists(plot.folder))
			dir.create(path=plot.folder, showWarnings=FALSE, recursive=TRUE)
		pdf(file=file.path(plot.folder, plot.filename), compress=FALSE)
	}


	# las=2: makes x labels vertical
	if(enable.legend){
		# Add extra space to right of plot area; change clipping to figure
		par(mar=c(5.1, 4.1, 4.1, 9.1), xpd=TRUE)

		barplot(count.mtrx, ylim=ylimit, col=colors, legend=rownames(count.mtrx), beside=is.beside, las=2,
				main=title, args.legend=list(title=legend.title, x="topright", inset=c(-0.4,0)), border=NA) # omit bar border
	}
	else
		barplot(count.mtrx, ylim=ylimit, col=colors, beside=is.beside, las=2, main=title, border=NA) # omit bar border

	if(plot.format == PLOT.AS.PDF)
		dev.off()
}


################################################################
#
################################################################
make.grouped.bar.plot = function(vector.list, count.levels=NA, ctgry.names, is.ctgry.names.on.x.axis=TRUE, ylimit=NA,
		plot.folder, plot.filename, plot.format, title, enable.legend=TRUE, legend.title)
{
#	test
#	a = c(2, 3, 2, 4, 4, 5, 2, 3, 4, 4)
#	b = c(3, 10, 3, 3, 3, 4, 2, 3, 3, 6)
#	my.l = list(a, b)
#	make.grouped.bar.plot(my.l,c("net1","net2"))

	if(is.ctgry.names.on.x.axis)
		make.bar.plot(vector.list, count.levels, row.names=NA, col.names=ctgry.names, combine.strategy="cbind", is.beside=TRUE,
				ylimit, title, plot.folder, plot.filename, plot.format, enable.legend, legend.title)
	else # if(is.ctgry.names.on.legend)
		make.bar.plot(vector.list, count.levels, row.names=ctgry.names, col.names=NA, combine.strategy="rbind", is.beside=TRUE,
				ylimit, title, plot.folder, plot.filename, plot.format, enable.legend, legend.title)
}



################################################################
#
################################################################
make.stacked.bar.plot = function(vector.list, count.levels=NA, ctgry.names, is.ctgry.names.on.x.axis=TRUE, ylimit=NA,
		plot.folder, plot.filename, plot.format, title, enable.legend=TRUE, legend.title){
#	test
#	a = c(0, 1, 2, 2, 2, 1, 0, 0, 2, 2)
#	b = c(0, 2, 1, 1, 2, 1, 2, 1, 1, 1)
#	my.l = list(a, b)
#	make.stacked.bar.plot(my.l,c("net1","net2"))

	make.bar.plot(vector.list, count.levels, row.names=NA, col.names=ctgry.names, combine.strategy="cbind", is.beside=FALSE,
			ylimit, title, plot.folder, plot.filename, plot.format, enable.legend, legend.title)
}

# END BAR PLOT 
################################################################################
################################################################################











################################################################################
################################################################################

#################################################################
# It is the function 'perform.layout.vio.plot.by()' which uses this function.
#
# Example 1: 
# matrix.list = list( matrix(c(6)), matrix(c(1)) )
# make.vioplot(matrix.list, plot.folder=NA, plot.filename=NA, plot.format=JUST.PLOT, title="title", x.axis.names=c("a","b"), ylabel='ylab', ylimit=c(1,10))
#
#
# Example 2:
# matrix.list = list( matrix(c(3,4,6,1)), matrix(c(6,3,4,1)) )
# make.vioplot(matrix.list, plot.folder=NA, plot.filename=NA, plot.format=JUST.PLOT, title="title", x.axis.names=c("a","b"), ylabel='ylab', ylimit=c(1,10))
#
#
# matrix.list
# plot.folder
# plot.filename
# plot.format
# title
# x.axis.names
# ylabel
# ylimit
#
#################################################################
make.vioplot = function(matrix.list, plot.folder, plot.filename, plot.format, title, x.axis.names, ylabel, ylimit){
    #### print("HERE1")
    
    library(vioplot)
	
	# find y limits
	if(all(is.na(ylimit))){
		all.y.values = unlist(vector.list)
		ylimit = c(min(all.y.values), max(all.y.values))
	} 

	if(plot.format == PLOT.AS.PDF){
		if(!dir.exists(plot.folder))
			dir.create(path=plot.folder, showWarnings=FALSE, recursive=TRUE)
		pdf(file=file.path(plot.folder, plot.filename), compress=FALSE)
	}

    matrix.list2 = lapply(matrix.list, function(vec) if(all(is.na(vec))) return(-50) else vec)
    matrix.list3 = lapply(matrix.list2, function(vec) if(length(unique(vec))==1) return(vec[1]) else vec)

    mean.values = sapply(matrix.list3, function(v) mean(v, na.rm=TRUE))
    
    
	# workaround: names() does not work if there is only 1 boxplot
	if(length(x.axis.names) == 1){ # single boxplot
#		vioplot(vector.list, outline=TRUE, xlab=x.axis.names, main=title, ylab=ylabel, ylim=ylimit)
#		points(mean.values,col="red",pch=18)
	    
		vioplot(matrix.list3, xlab=x.axis.names, las=2, colMed="blue", pchMed=16, cex.lab=1.5, cex.axis=1.5, cex.main=1.5, # las=2 makes labels vertical on x axis
				main=title, ylab=ylabel, ylim=ylimit, rectCol="palevioletred", lineCol="violetred", areaEqual = FALSE)
        points(mean.values,col="green",pch=24)
        

	} else { # multiple boxplots
#		boxplot(vector.list, outline=TRUE, names=x.axis.names, las=2, # las=2 makes labels vertical on x axis
#				main=title, ylab=ylabel, ylim=ylimit)
#		points(mean.values,col="red",pch=18)

		vioplot(matrix.list3, names=x.axis.names, las=2, colMed="blue", pchMed=16, cex.lab=1.5, cex.axis=1.5, cex.main=1.5, # las=2 makes labels vertical on x axis
				main=title, ylab=ylabel, ylim=ylimit, rectCol="palevioletred", lineCol="violetred", areaEqual = FALSE)
        points(mean.values,col="green",pch=24)
	}
		
	if(plot.format == PLOT.AS.PDF)
		dev.off()
	
}



#################################################################
# it is the function 'perform.layout.vio.plot.by.with.sol.classes()' which uses this function.
# For a quick visualization, run the small example below:
## Small example:
# matrix.list.of.list2 = list(list(), list("1"=matrix(c(2,5)),"2"=matrix(c(2,7))), list("1"=matrix(c(3,4)),"2"=matrix(c(1,2))))
# make.grouped.vioplot2(matrix.list.of.list2, plot.folder=NA, plot.filename=NA,
#   plot.format=JUST.PLOT, param.values=c("n=20","n=24"), ctgry.names=c("1","2"), title="title", ylabel="ylab", ylimit=c(1,10))
# 
#
# matrix.list.of.list2
# plot.folder
# plot.filename
# plot.format
# param.values
# ctgry.names
# title
# x.axis.names
# ylabel
# ylimit
#
# 
#
# # Another example of 'matrix.list.of.list2':
# # [[1]]
# #     list()
# # 
# # [[2]]
# #     [[2]]$`1` ====> results with 1 solution class
# #     network6 k=1  network8 k=1 network10 k=1 
# #     0.8750        0.9167        1.0000 
# #     
# #     [[2]]$`2`  ====> results with 2 solution classes
# #     network1 k=1 network1 k=2 network2 k=1 network2 k=2 network4 k=1 network4 k=2 
# #     1.0000       1.0000       1.0000       1.0000       1.0000       1.0000 
# #     network5 k=1 network5 k=2 network7 k=1 network7 k=2 network9 k=1 network9 k=2 
# #     1.0000       1.0000       0.9583       0.9167       1.0000       1.0000 
# # 
# # [[3]]
# #     [[3]]$`1` ====> results with 1 solution class
# #     network3 k=1 network4 k=1 network7 k=1 network8 k=1 
# #     0.8750       0.8333       0.7500       0.8750 
# #     
# #     [[3]]$`2` ====> results with 2 solution classes
# #     network1 k=1  network1 k=2  network5 k=1  network5 k=2  network9 k=1 
# #     0.8333        0.7917        1.0000        1.0000        0.9583 
# #     network9 k=2 network10 k=1 network10 k=2 
# #     0.9583        1.0000        0.7500 
# # 
# # [[4]]
# #     [[4]]$`1` ====> results with 1 solution class
# #     network2 k=1 network3 k=1 network4 k=1 network5 k=1 network7 k=1 network9 k=1 
# #     0.7917       0.8750       0.7500       0.3750       0.9167       0.8750 
# #     
# #     [[4]]$`2` ====> results with 2 solution classes
# #     network8 k=1  network8 k=2 network10 k=1 network10 k=2 
# #     1             1             1             1 
# # 
# # [[5]]
# #     [[5]]$`1` ====> results with 1 solution class
# #     network7 k=1 
# #     0.5833 
# #     
# #     [[5]]$`2` ====> results with 2 solution classes
# #     network1 k=1  network1 k=2  network3 k=1  network3 k=2  network4 k=1 
# #     1.0000        1.0000        0.9583        1.0000        0.5000 
# #     network4 k=2  network5 k=1  network5 k=2  network6 k=1  network6 k=2 
# #     0.4583        0.9583        0.8750        0.4583        0.7083 
# #     network8 k=1  network8 k=2 network10 k=1 network10 k=2 
# #     0.5833        0.6250        0.5833        0.6250 
# #     
# #     [[5]]$`4` ====> results with 4 solution classes
# #     network9 k=1 network9 k=2 network9 k=3 network9 k=4 
# #     0.8750       0.9167       0.8750       0.9167 
# # 
# #     [[5]]$`6` ====> results with 6 solution classes
# #     network2 k=1 network2 k=2 network2 k=3 network2 k=4 network2 k=5 network2 k=6 
# #     0.5417       0.7083       0.5833       0.7083       0.5833       0.8333 
#
#################################################################
make.grouped.vioplot2 = function(matrix.list.of.list2, plot.folder, plot.filename, plot.format,
		param.values, ctgry.names, title, ylabel, ylimit){
    #### print("HERE2")

    
	library(vioplot)
	
	count = 1
	v.list = list()
	last.x.indx = 1
	x.at = c()
	x2.at = c()
	x.axis.names = c()
	
	
	for(i in 1:length(param.values)){
		l = matrix.list.of.list2[[i]] # list
		if(length(l) == 0){
			l = list("NA"=c(-50)) # workaround: if there is no info, put negative value, it will be invisible in the plot
		}
		#print(l)
		
		first = last.x.indx
		last = last.x.indx + length(l)-1
		x.at = c(x.at, seq(first, last))
		x2.at = c(x2.at, ceiling((first+last)/2) )
		last.x.indx = last.x.indx + length(l)+ 2 # 2 is for spaces between diff x axis param values
		x.axis.names = c(x.axis.names, names(l))
		
		for(j in 1:length(l)){
			v.list[[count]] = l[[j]]
			count = count + 1
		}
		
	}
	
	mean.values = sapply(1:length(v.list), function(i) mean(v.list[[i]], na.rm=TRUE))


	# source: https://tomizonor.wordpress.com/2013/04/18/color-boxplot/
	#c1 = unlist(palette[1:nb.categry])
	c1 = "red"
	c2 <- adjustcolor(c1, alpha.f = 0.2)
	# palette
#	boxplot(v.list, names = x.axis.names, xaxs = FALSE, at = x.at, horizontal=FALSE,  las=2, ylim=ylimit,
#		col=c2, medcol=c1, whiskcol=c1, staplecol=c1, boxcol=c1, outcol=c1, cex.axis=0.7)
#	points(x=x.at, y=mean.values,col="red",pch=18)

	vioplot(v.list, names = x.axis.names, xaxs = FALSE, at = x.at, horizontal=FALSE,  las=2, ylim=ylimit, colMed="blue", pchMed=16, rectCol="palevioletred", lineCol="violetred",
	col=c2, medcol=c1, whiskcol=c1, staplecol=c1, boxcol=c1, outcol=c1, cex.lab=1.5, cex.axis=0.6, cex.main=1.5) # cex.axis=0.7,
    points(x.at, mean.values,col="green",pch=24)

	axis(side=1,at=x2.at,labels=param.values, mgp=c(3,2.5,1), lwd.ticks = FALSE, line=FALSE, cex.axis=1.25) # cex.axis=0.7
}



# #################################################################
# # DELETE IT !!
# #################################################################
# make.grouped.vioplot = function(vector.list, plot.folder, plot.filename, plot.format,
# 		param.values, ctgry.names, title, x.axis.names, ylabel, ylimit){
#     # print("in make.grouped.boxplot")
#     
#     library(vioplot)
# 	
# 	vector.list2 = lapply(vector.list, function(vec) if(all(is.na(vec))) return(0) else vec)
# 	
# 	v.list = list()
# 	if(length(param.values) == length(vector.list2)){
# 		count = 1
# 		for(i in 1:length(param.values)){
# 			if(is.matrix(vector.list2[[i]])){
# 				m = vector.list2[[i]]
# 				for(j in 1:ncol(m)){
# 					v.list[[count]] = m[,j]
# 					count = count + 1
# 				}
# 			} else { # vector
# 				v.list[[count]] = vector.list2[[i]]
# 				count = count + 1
# 			}
# 		}
# 	} else {
# 		v.list = vector.list2
# 	}
# 	
# 	nb.categry = length(ctgry.names)
# 	
# 	mean.values = sapply(1:length(v.list), function(i) mean(v.list[[i]], na.rm=TRUE))
# 	
# 	x.axis.names = rep("",length(v.list))
# 	indx = 1
# 	indx2 = ceiling(nb.categry/2)
# 	x.at = c()
# 	for(i in 1:length(param.values)){
# 		x.at = c(x.at, seq(indx,indx+nb.categry-1))
# 		x.axis.names[indx2] = param.values[i]
# 		indx = indx + nb.categry + 2 # +2 for space
# 		indx2 = indx2 + nb.categry
# 	}
# 
# 	
# 	# source: https://tomizonor.wordpress.com/2013/04/18/color-boxplot/
# 	c1 = unlist(palette[1:nb.categry])
# 	c2 <- adjustcolor(c1, alpha.f = 0.2)
# 	# palette
# #	boxplot(v.list, names = x.axis.names, xaxs = FALSE, at = x.at, horizontal=FALSE,  las=2, ylim=ylimit,
# #		col=c2, medcol=c1, whiskcol=c1, staplecol=c1, boxcol=c1, outcol=c1)
# #	points(x=x.at, y=mean.values,col="red",pch=18)
# 	
# 	vioplot(v.list, names = x.axis.names, xaxs = FALSE, at = x.at, horizontal=FALSE,  las=2, ylim=ylimit, colMed="blue", pchMed=16, rectCol="palevioletred", lineCol="violetred",
# 	col=c2, medcol=c1, whiskcol=c1, staplecol=c1, boxcol=c1, outcol=c1, cex.lab=1.5, cex.axis=1.5, cex.main=1.5)
#     points(x.at, mean.values,col="green",pch=24)
# }




# ================================================================================================
# ================================================================================================

#################################################################
# It controls the choice of plot function depending on plot.args$plot.type.
#
# vector.list: data to be plotted
# plot.format: PLOT.AS.PDF or PLOT.AS.PLOT.OBJECT
# plot.types: list of plot type (DENS, HIST, MULTI.BAR, LINE, MULTI.LINE, BOX)
# ctgry.names: mainly x axis tick names (for barplot and boxplot)
#
#################################################################
make.plot = function(vector.list, plot.args, legend.args)
{
	

	
	# ----------------------------------------------------------------------------------------
	# density
	if(plot.args$plot.type == DENS){
		# TODO
	}
	# ----------------------------------------------------------------------------------------
	
	# ----------------------------------------------------------------------------------------
	# hist
	if(plot.args$plot.type == HIST){
		values = vector.list[[1]]
		make.histogram.plot(values, plot.args$plot.title, plot.args$xlabel, plot.args$plot.folder,
				plot.args$plot.filename, plot.args$plot.format)
	}
	# ----------------------------------------------------------------------------------------
	
	# ----------------------------------------------------------------------------------------
	# multi bar
	if(plot.args$plot.type == MULTI.BAR){
		# x axis will be shown "nb cluster", and each x value is grouped by networks, so category.names will be on x axis
		if(plot.args$use.stacked.bars)
			make.stacked.bar.plot(vector.list, plot.args$count.levels, plot.args$ctgry.names, plot.args$is.ctgry.names.on.x.axis, 
					plot.args$ylimit, plot.args$plot.folder, plot.args$plot.filename, plot.args$plot.format, plot.args$plot.title,
					legend.args$enabled, legend.args$title)
		else
			make.grouped.bar.plot(vector.list, plot.args$count.levels, plot.args$ctgry.names, plot.args$is.ctgry.names.on.x.axis, 
					plot.args$ylimit, plot.args$plot.folder, plot.args$plot.filename, plot.args$plot.format, plot.args$plot.title,
					legend.args$enabled, legend.args$title)
	}
	
	
	if(plot.args$plot.type == MULTI.BAR.WITH.PROPORTION){
        print("MULTI.BAR.WITH.PROPORTION")
		make.bar.plot.with.proportion(vector.list, plot.args$ctgry.names, plot.args$param.values, ylimit=plot.args$ylimit,
			plot.args$plot.folder, plot.args$plot.filename, plot.args$plot.format, plot.args$plot.title, legend.args$enabled, legend.args$title)
	}
	# ----------------------------------------------------------------------------------------
	
	# ----------------------------------------------------------------------------------------
	# box
	if(plot.args$plot.type == VIO){
		make.vioplot(vector.list, plot.args$plot.folder, plot.args$plot.filename, plot.args$plot.format, 
				plot.args$plot.title, plot.args$x.axis.names, plot.args$ylabel, plot.args$ylimit)
	}
	
	# box2
	if(plot.args$plot.type == VIO2){
		make.grouped.vioplot2(vector.list, plot.args$plot.folder, plot.args$plot.filename, plot.args$plot.format, 
				plot.args$param.values, plot.args$ctgry.names, plot.args$plot.title, plot.args$ylabel, plot.args$ylimit)

	}
	
	# if(plot.args$plot.type == GROUPED.VIO){
	# 	make.grouped.vioplot(vector.list, plot.args$plot.folder, plot.args$plot.filename, plot.args$plot.format, 
	# 			plot.args$param.values, plot.args$ctgry.names, plot.args$plot.title, plot.args$x.axis.names, plot.args$ylabel, plot.args$ylimit)
	# }
	# ----------------------------------------------------------------------------------------

	# ----------------------------------------------------------------------------------------
	# line
	if(plot.args$plot.type == LINE){
		y.values = vector.list[[1]]	
		make.line.plot(y.values, plot.args$xvalues, plot.args$plot.title, plot.args$xlabel, plot.args$ylabel, plot.args$plot.folder,
				plot.args$plot.filename, plot.args$plot.format, plot.args$increasing.order, plot.args$y.axis.log.scale, plot.args$ylimit)
	}
		

	
	# ----------------------------------------------------------------------------------------

	# ----------------------------------------------------------------------------------------
	# multi line
	if(plot.args$plot.type == MULTI.LINE){
		make.multi.line.plot(vector.list, plot.args$xvalues, plot.args$xlabel, plot.args$ylabel, plot.args$ylimit, 
				plot.args$plot.title, plot.args$plot.folder, plot.args$plot.filename, plot.args$plot.format,
				legend.args$enabled, legend.args$legend, legend.args$title)
	}
	# ----------------------------------------------------------------------------------------


	if(plot.args$plot.type == SCATTER){
		y.values = vector.list
		make.scatter.plot(y.values, plot.args$xvalues, plot.args$plot.title, plot.args$xlabel, plot.args$ylabel, plot.args$plot.folder,
				plot.args$plot.filename, plot.args$plot.format, plot.args$ylimit, plot.args$xlimit,
				 plot.args$max.nb.elt, legend.args$enabled)
	}

}
