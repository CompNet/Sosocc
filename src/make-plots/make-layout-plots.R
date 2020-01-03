


##################################################################
# It performs two things: 
#   1) to adjust max value for two parameters
#   2) to determine the used values of prop.negs. We can not know them in advance ..
#
##################################################################
check.and.adjust.upper.bounds = function(g.params, by, param.values, subplot.by, subplot.param.values){

    max.index.for.param.values = 1
    max.indx.for.subplot.param.values = 1

    initial.g.params = g.params
    for(i in 1:length(subplot.param.values)){
        subplot.value = subplot.param.values[i]
        g.params = update.g.params.by(g.params=initial.g.params, by=subplot.by, value=subplot.value)
        
        for(j in 1:length(param.values)){
            param.value = param.values[j]
            g.params = update.g.params.by(g.params=g.params, by=by, value=param.value)
            
            # ===========================
            # up to now, all params are provided, except prop.negs, now complete this
            eval.upper.folder = get.eval.folder.path(g.params$n, g.params$l0, g.params$d, g.params$prop.mispl, g.params$detected.imb.interval,
                                                     NA, NA, NA, NA, NA)
            existing.folders = list.dirs(path = eval.upper.folder, full.names = FALSE, recursive = FALSE)
            my.prop.negs = as.numeric(gsub("propNeg=","",existing.folders))  # it might give 'numeric(0)' if there is not any folder (which is equivalent to do "dir.exists(eval.folder)")
            # ===========================
            

            for(prop.neg in my.prop.negs){
                eval.folder = get.eval.folder.path(g.params$n, g.params$l0, g.params$d, g.params$prop.mispl, g.params$detected.imb.interval,
                                                   prop.neg, g.params$k, g.params$network.no, g.params$algo.name, g.params$graph.desc.name)

                if(dir.exists(eval.folder)){
                    # 'subplot.param.values' is in the outer loop
                    # as soon as 'i' increases, 'max.index.for.param.values' should be, as well
                    max.indx.for.subplot.param.values = i
                    
                    # 'param.values' is in the inner loop
                    # we will retreive the max value used over all 'subplot.param.values' values
                    if(j>max.index.for.param.values)
                        max.index.for.param.values = j
                }
                
            }
            
            
        }
    }
    
    param.values = param.values[1:max.index.for.param.values]
    subplot.param.values = subplot.param.values[1:max.indx.for.subplot.param.values]

    result = list(param.values=param.values, subplot.param.values=subplot.param.values)
    return(result)
}



#################################################################
#
# g.params
# by
# values
# subplot.by
# subplot.param.values
#

# Example of output: a list of a list of matrix values
#                       Inner list: graph sizes (here: n \in {20,24})
#                       Outer list: prop.mispl values (here: props.mispl \in {0.05, 0.10, 0.15})
# RESULT TYPE = NB.SOLUTION
# $`0.05` ========> prop.mispl
#   $`0.05`[[1]] ========> n=20
#                    [,1]
#   network1 sol0     1
#   network2 sol0     2
#   network3 sol0     2
#   network4 sol0     1
#   network5 sol0     2
#   network6 sol0     5
#   network7 sol0     2
#   network8 sol0     1
#   network9 sol0     4
#   network10 sol0    1
# 
#   $`0.05`[[2]] ========> n=24
#                    [,1]
#   network1 sol0     1
#   network2 sol0     1
#   network3 sol0     1
#   network4 sol0     1
#   network5 sol0     1
#   network6 sol0     1
#   network7 sol0     1
#   network8 sol0     1
#   network9 sol0     1
#   network10 sol0    1
# 
# $`0.1` ========> prop.mispl
#   $`0.1`[[1]] ========> n=20
#                    [,1]
#   network1 sol0    10
#   network2 sol0    40
#   network3 sol0     6
#   network4 sol0    12
#   network5 sol0     6
#   network6 sol0     2
#   network7 sol0     2
#   network8 sol0     4
#   network9 sol0     2
#   network10 sol0    4
# 
#   $`0.1`[[2]] ========> n=24
#                   [,1]
#   network1 sol0     2
#   network2 sol0     2
#   network3 sol0     1
#   network4 sol0     2
#   network5 sol0     2
#   network6 sol0     6
#   network7 sol0     5
#   network8 sol0     3
#   network9 sol0     2
#   network10 sol0    4
# 
# $`0.15` ========> prop.mispl
#   $`0.15`[[1]] ========> n=20
#                  [,1]
#   network1 sol0    35
#   network2 sol0    53
#   network3 sol0    31
#   network4 sol0    14
#   network5 sol0     1
#   network6 sol0    15
#   network7 sol0     6
#   network8 sol0     2
#   network9 sol0     1
#   network10 sol0    6
# 
#   $`0.15`[[2]] ========> n=24
#                   [,1]
#   network1 sol0    10
#   network2 sol0     1
#   network3 sol0     3
#   network4 sol0     4
#   network5 sol0    12
#   network6 sol0     1
#   network7 sol0     4
#   network8 sol0     8
#   network9 sol0     8
#   network10 sol0    7



##################################################################
#
#
##################################################################
retreive.layout.data.by = function(result.filename, colname, g.params, by, param.values, subplot.by, subplot.param.values)
{
    #	g.params = update.g.params.by(g.params=g.params, by=by, value=ALL.PARAM.NAME)
    #	g.params = update.g.params.by(g.params=g.params, by=subplot.by, value=ALL.PARAM.NAME)

    data = list()
    for(subplot.value in subplot.param.values){
        g.params = update.g.params.by(g.params=g.params, by=subplot.by, value=subplot.value)
        vector.list = retreive.csv.results(g.params, result.filename, colname, by, param.values)
        data[[as.character(subplot.value)]] = vector.list
    }
    return(data)
}



##################################################################
#
#
##################################################################
retreive.layout.data.by.with.sol.classes = function(result.filename, colname, g.params, measure, by, param.values, subplot.by, subplot.param.values)
{
    #	g.params = update.g.params.by(g.params=g.params, by=by, value=ALL.PARAM.NAME)
    #	g.params = update.g.params.by(g.params=g.params, by=subplot.by, value=ALL.PARAM.NAME)
    
    data = list()
    for(subplot.value in subplot.param.values){
        g.params = update.g.params.by(g.params=g.params, by=subplot.by, value=subplot.value)
        vector.list = retreive.csv.results.by.sol.classes(g.params, result.filename, colname, measure, by, param.values)
        data[[as.character(subplot.value)]] = vector.list
    }
    
    return(data)
}



#################################################################
#
# g.params
# by
# values
# subplot.by
# subplot.param.values
#
##################################################################
determine.ylimit.as.value.for.layout.by = function(layout.data, subplot.param.values){
	max.vals = c()
	min.vals = c()
	for(subplot.value in subplot.param.values){
		vector.list = layout.data[[as.character(subplot.value)]]
		if(all(is.na(unlist(vector.list)))){
			max.vals = c(max.vals, NA)
			min.vals = c(min.vals, NA)
		} else {
			max.vals = c(max.vals, max(unlist(vector.list),na.rm=TRUE))
			min.vals = c(min.vals, min(unlist(vector.list),na.rm=TRUE))
		}
	}

	
	if(all(is.na(min.vals)) || all(is.na(max.vals))){
		return(NULL)
	}
	
	ylimit = c(min(min.vals,na.rm=TRUE), max(max.vals,na.rm=TRUE))
	return(ylimit)
}


#################################################################
# This method is used for barplot (also histogram if we use it)
# g.params
# by
# values
# subplot.by
# subplot.param.values
#
##################################################################
determine.ylimit.as.freq.for.layout.by = function(layout.data, subplot.param.values){
	max.vals = c()
	min.vals = c()
	for(subplot.value in subplot.param.values){
		vector.list = layout.data[[as.character(subplot.value)]]
		max.freqs = sapply(vector.list, function(vals) max(table(vals)))
		min.freqs = sapply(vector.list, function(vals) min(table(vals)))
		max.vals = c(max.vals, max(max.freqs))
		min.vals = c(min.vals, min(min.freqs))
	}
#	ylimit = c(min(min.vals), max(max.vals))
	# I think it is not meaningful to compute lower bound of ylimit, because if for min value it will show a line (i.e. without height), not a bar
	ylimit = c(0, max(max.vals))
	return(ylimit)
}


#############################################
#############################################



#################################################################
#
# g.params
# by
# values
# subplot.by
# subplot.param.values
#
##################################################################
perform.layout.multi.bar.with.prop.plot.by = function(result.type, layout.data, g.params, by, param.values, subplot.by, subplot.param.values,
		is.ctgry.names.on.x.axis, ctgry.names=NA, xlabel=NA, ylabel=NA, ylimit=NA)
{

	if(is.na(xlabel))
		xlabel = by
	
	# layout.data = retreive.layout.data.by(result.type, g.params, by, param.values, subplot.by, subplot.param.values)	
	#print(layout.data)

# 	# ----------------------------------------------------
# 	# workaround
# 		if(by != DETECTED.IMB.PROP.PARAM.NAME){
# 		param.values = c()
# 		prop.seq.vals = seq(MIN.PROP.VAL, MAX.PROP.VAL, 0.05)
# 		for(i in 1:(length(prop.seq.vals)-1)){
# 			param.values  = c(param.values, paste0(prop.seq.vals[i],"-",prop.seq.vals[i+1]) )
# 		}
#         # ---- ATTENTION
# 		#param.values = c(param.values, paste0(prop.seq.vals[length(prop.seq.vals)],">"))
#         # ---- ATTENTION
# 	}
# 	# ----------------------------------------------------


	#shorthen.result.type.name.if.possible(result.type)
	layout.plot.args = list(xlabel=xlabel, ylabel="Proportion")
	layout.legend.args = list(title=result.type, legend=ctgry.names, col=unlist(palette))
	plot.args = list(ylimit=ylimit, plot.type=MULTI.BAR.WITH.PROPORTION,
			ctgry.names=ctgry.names, param.values=param.values, xlabel=NA, ylabel=NA, enable.plot.title=FALSE)
	legend.args = list(enabled=FALSE, title=NA)
	

	perform.generic.layout.plot.by(layout.data, result.type, g.params, by, param.values,
			subplot.by, subplot.param.values, layout.plot.args, layout.legend.args, plot.args, legend.args)
}




# #################################################################
# # OLD
#
# # g.params
# # by
# # values
# # subplot.by
# # subplot.param.values
# #
# ##################################################################
# perform.layout.multi.bar.plot.by = function(result.type, layout.data, g.params, by, param.values, subplot.by, subplot.param.values,
# 		is.ctgry.names.on.x.axis, ctgry.names=NA, xlabel=NA, ylabel=NA, ylimit=NA, use.stacked.bars=FALSE)
# {
# 	print("perform.layout.multi.bar.plot.by")
# 	# --------------------------------------------
# 	count.levels = c()
# 	tmp.g.params = g.params
# 	for(subplot.value in subplot.param.values){
# 		tmp.g.params = update.g.params.by(tmp.g.params, by=subplot.by, value=subplot.value)
# 		vector.list = retreive.results.by.result.type.and.param.name(result.type, tmp.g.params, by, param.values)
# 		level.values = lapply(vector.list, function(v) names(table(v)))
# 		curr.count.levels = sort(unique(as.numeric(do.call("c", level.values ))))
# 		count.levels = c(count.levels, curr.count.levels)
# 	}
# 	count.levels = sort(unique(count.levels))
# 	# --------------------------------------------
# 	
# 	if(is.na(xlabel))
# 		xlabel = by
# 	
# 	layout.data = retreive.layout.data.by(result.type, g.params, by, param.values, subplot.by, subplot.param.values)
# 	
# 	if(is.na(ylimit))
# 		ylimit = determine.ylimit.as.freq.for.layout.by(layout.data, subplot.param.values)
# 	
# 	# TODO gets error bc of 'retreive.nb.partition()'
# 	# ----------------------------------
# 	nb.sol.list = list()
# 	if(is.na(ctgry.names)){
# 		tmp.g.params = g.params
# 		for(subplot.value in subplot.param.values){
# 			tmp.g.params = update.g.params.by(tmp.g.params, by=subplot.by, value=subplot.value)
# 			nb.sol.vector = unlist(retreive.nb.partition(tmp.g.params, by, folders=param.values))
# 			nb.sol.list[[as.character(subplot.value)]] = nb.sol.vector
# 		}
# 		
# #		ctgry.names = paste0(param.values, "\n(", nb.tot.solution.vector, " sols)")
# #		if(by == ALGO.NAME && all(param.values == COR.CLU.HEURISTIC.ALGOS)) # if we process plots by heuristic algos
# #			ctgry.names = shorthen.cor.clu.heur.algo.names(param.values)
# 	}
# 	# ----------------------------------
# 	
# 	ctgry.names = param.values
# 	if(by == ALGO.NAME && all(param.values == COR.CLU.HEURISTIC.ALGOS)) # if we process plots by heuristic algos
# 		ctgry.names = shorthen.cor.clu.heur.algo.names(param.values)
# 
# 
# 
# 	
# 	layout.plot.args = list(xlabel=xlabel, ylabel="Frequency")
# 	layout.legend.args = list(title=shorthen.result.type.name.if.possible(result.type), legend=count.levels, col=unlist(palette))
# 	plot.args = list(ylimit=ylimit, plot.type=MULTI.BAR, is.ctgry.names.on.x.axis=is.ctgry.names.on.x.axis,
# 			count.levels=count.levels, ctgry.names=ctgry.names, xlabel=NA, ylabel=NA, enable.plot.title=FALSE, 
# 			use.stacked.bars=use.stacked.bars)
# 	legend.args = list(enabled=FALSE, title=NA)
# 	
# 	perform.generic.layout.plot.by(layout.data, result.type, g.params, by, param.values,
# 			subplot.by, subplot.param.values, layout.plot.args, layout.legend.args, plot.args, legend.args)
# }





#################################################################
#
# g.params
# by
# param.values
#
##################################################################
perform.layout.multi.stacked.bar.plot.by = function(result.type, g.params, by, param.values, subplot.by, subplot.param.values, 
		is.ctgry.names.on.x.axis, ctgry.names, xlabel=NA, ylabel=NA, ylimit=NA)
{	
	print("perform.layout.multi.stacked.bar.plot.by")
	
	tmp.g.params = g.params
	heights = c()
	for(subplot.value in subplot.param.values){
		tmp.g.params = update.g.params.by(tmp.g.params, by=subplot.by, value=subplot.value)
		vector.list = retreive.results.by.result.type.and.param.name(result.type, tmp.g.params, by, param.values)
		curr.heights = sapply(vector.list, function(vec) length(unique(vec[which(!is.na(vec))])) )
		heights = c(heights, max(curr.heights))
	}
	
	if(is.na(ylimit)){
		ylimit = c(0, max(heights))
	}
	
	
	perform.layout.multi.bar.plot.by(result.type, g.params, by, param.values, subplot.by, subplot.param.values, 
			is.ctgry.names.on.x.axis=TRUE, ctgry.names=ctgry.names, use.stacked.bars=TRUE, ylimit=ylimit)
}



#################################################################
#
##################################################################
perform.layout.multi.line.plot.by = function(result.type, g.params, by, param.values, subplot.by, subplot.param.values,
		xvalues=NA, xlabel, ylabel, ylimit=NA)
{	
	print("perform.layout.multi.line.plot.by")
	
	layout.data = retreive.layout.data.by(result.type, g.params, by, param.values, subplot.by, subplot.param.values)	
	
	if(is.na(ylimit))
		ylimit = determine.ylimit.as.value.for.layout.by(layout.data, subplot.param.values)
	
	lines.names = param.values
	if(by == ALGO.NAME && all(param.values == COR.CLU.HEURISTIC.ALGOS)) # if we process plots by heuristic algos
		lines.names = shorthen.cor.clu.heur.algo.names(param.values)
	
	layout.plot.args = list(xlabel=xlabel, ylabel=ylabel)
	layout.legend.args = list(title=by, legend=lines.names, col=unlist(palette))
	plot.args = list(ylimit=ylimit, plot.type=MULTI.LINE, xvalues=xvalues, ylabel=NA, xlabel=NA, enable.plot.title=FALSE)
	legend.args = list(enabled=FALSE, title=NA, legend=NA)
	
	
	perform.generic.layout.plot.by(layout.data, result.type, g.params, by, param.values,
			subplot.by, subplot.param.values, layout.plot.args, layout.legend.args, plot.args, legend.args)
}







#################################################################
#
# g.params
# by
# values
# subplot.by
# subplot.param.values
#
##################################################################
perform.layout.vio.plot.by = function(result.type, layout.data, g.params, by, param.values, subplot.by, subplot.param.values,
		xlabel=NA, ylabel=NA, ylimit=NA)
{	
	print("perform.layout.vio.plot.by")
	
	if(is.na(xlabel))
		xlabel = by
	if(is.na(ylabel))
		ylabel=prepare.axis.label(result.type)

	if(is.na(ylimit))
		ylimit = determine.ylimit.as.value.for.layout.by(layout.data, subplot.param.values)


	x.axis.names = param.values
	if(by == ALGO.NAME) # if we process plots by heuristic algos
		x.axis.names = shorthen.cor.clu.heur.algo.names(param.values)

	layout.plot.args = list(xlabel=xlabel, ylabel=ylabel)
	layout.legend.args = NA
	plot.args = list(ylimit=ylimit, plot.type=VIO, x.axis.names=x.axis.names, xlabel=NA, ylabel=NA, enable.plot.title=FALSE)
	legend.args = NA

	perform.generic.layout.plot.by(layout.data, result.type, g.params, by, param.values,
			subplot.by, subplot.param.values, layout.plot.args, layout.legend.args, plot.args, legend.args)
}



#################################################################
#
# g.params
# by
# values
# subplot.by
# subplot.param.values
#
##################################################################
perform.layout.vio.plot.by.with.sol.classes = function(result.type, layout.data, g.params, by, param.values, subplot.by, subplot.param.values,
		xlabel=NA, ylabel=NA, ylimit=NA)
{	
	# print("perform.layout.box.plot.by.with.sol.classes")

	if(is.na(xlabel))
		xlabel = by
	if(is.na(ylabel))
		ylabel=prepare.axis.label(result.type)

	if(is.na(ylimit))
		ylimit = determine.ylimit.as.value.for.layout.by(layout.data, subplot.param.values)


	layout.plot.args = list(xlabel=xlabel, ylabel=ylabel)
	layout.legend.args = NA
	plot.args = list(ylimit=ylimit, plot.type=VIO2, param.values=param.values, xlabel=NA, ylabel=NA, enable.plot.title=FALSE)
	legend.args = NA
	
	perform.generic.layout.plot.by(layout.data, result.type, g.params, by, param.values,
			subplot.by, subplot.param.values, layout.plot.args, layout.legend.args, plot.args, legend.args)
}





# #################################################################
# #
# # g.params
# # by
# # values
# # subplot.by
# # subplot.param.values
# #
# ##################################################################
# perform.layout.grouped.vio.plot.by = function(result.type, g.params, by, param.values, subplot.by, subplot.param.values,
# 		ctgry.names=NA, xlabel=NA, ylabel=NA, ylimit=NA)
# {	
# 	print("perform.layout.grouped.box.plot.by")
# 	layout.data = retreive.layout.data.by(result.type, g.params, by, param.values, subplot.by, subplot.param.values)
# 	#print(layout.data)
# 	
# 	
# 	if(is.na(ylimit))
# 		ylimit = NULL
# 	# it may not be a good solution because of outliers
# 	#ylimit = determine.ylimit.as.value.for.layout.by(layout.data, subplot.param.values)
# 	
# 	#shorthen.result.type.name.if.possible(result.type)
# 	layout.plot.args = list(xlabel=xlabel, ylabel=ylabel)
# 	layout.legend.args = list(title=result.type, legend=ctgry.names, col=unlist(palette))
# 	plot.args = list(ylimit=ylimit, plot.type=GROUPED.VIO,
# 			ctgry.names=ctgry.names, param.values=param.values, xlabel=NA, ylabel=NA, enable.plot.title=FALSE)
# 	legend.args = list(enabled=FALSE, title=NA)
# 	
# 	perform.generic.layout.plot.by(layout.data, result.type, g.params, by, param.values,
# 			subplot.by, subplot.param.values, layout.plot.args, layout.legend.args, plot.args, legend.args)
# }





#################################################################
#
# g.params
# by
# values
# subplot.by
# subplot.param.values
#
##################################################################
perform.layout.line.plot.by = function(result.type, layout.data, g.params, by, param.values, subplot.by, subplot.param.values,
		xvalues=NA, xlabel=NA, ylabel=NA, ylimit=NA, increasing.order=FALSE, y.axis.log.scale=FALSE)
{	

	# --------------------------------------
	if(result.type == NB.SOLUTION || result.type == EXEC.TIME){
		# TODO WARNING: we need to transform 'vector.list' for line plot beforehand
		
		for(subplot.value in subplot.param.values){
			vector.list = layout.data[[as.character(subplot.value)]]
			vector.values = unlist(vector.list) # if there are x network instances, there will be x values	
			names(vector.values) = param.values
			vector.list = list(vector.values) # list with 1 element
			layout.data[[as.character(subplot.value)]] = vector.list
		}
	}
	# --------------------------------------
	print(layout.data)
	
	
	if(is.na(xlabel))
		xlabel = by
	if(is.na(ylabel))
		ylabel=prepare.axis.label(result.type)
	
	if(y.axis.log.scale)
		ylabel = paste0(ylabel, " (log scale)")
	
	if(is.na(ylimit))
		ylimit = determine.ylimit.as.value.for.layout.by(layout.data, subplot.param.values)
	
	layout.plot.args = list(xlabel=xlabel, ylabel=ylabel)
	layout.legend.args = NA
	plot.args = list(ylimit=ylimit, plot.type=LINE, xvalues=xvalues, xlabel=NA, ylabel=NA, x.values=param.values,
			increasing.order=increasing.order, enable.plot.title=FALSE, y.axis.log.scale=y.axis.log.scale)
	legend.args = NA
	
	print("!!")
	print(plot.args)
	perform.generic.layout.plot.by(layout.data, result.type, g.params, by, param.values,
			subplot.by, subplot.param.values, layout.plot.args, layout.legend.args, plot.args, legend.args)
}



# #################################################################
# IT IS NOT USED ANYMORE. AND, IT NEEDS TO BE UPDATED.
#
# # g.params
# # subplot.by
# # subplot.param.values
# #
# ##################################################################
# perform.layout.scatter.plot.by = function(main.result.type, y.result.type, x.result.type, g.params, subplot.by, subplot.param.values,
# 		ylabel=NA, xlabel=NA, ylimit=NA, xlimit=NA)
# {
# 	print("perform.layout.scatter.plot.by")
# 	x.values.list = retreive.results.by.result.type.and.param.name(x.result.type, g.params, subplot.by, subplot.param.values)
# 	names(x.values.list) = as.character(subplot.param.values)
# 	y.values.list = retreive.results.by.result.type.and.param.name(y.result.type, g.params, subplot.by, subplot.param.values)
# 	names(y.values.list) = as.character(subplot.param.values)
# 	
# 		
# 	if(is.na(xlabel))
# 		xlabel = prepare.axis.label(x.result.type)
# 	if(is.na(ylabel))
# 		ylabel=prepare.axis.label(y.result.type)
# 	
# 	
# 	if(is.na(xlimit)){
# 		min.vals = sapply(x.values.list, function(x.values) min(x.values,na.rm=TRUE))
# 		max.vals = sapply(x.values.list, function(x.values) max(x.values,na.rm=TRUE))
# 		xlimit = c(min(min.vals), max(max.vals))
# 	}
# 	if(is.na(ylimit)){
# 		min.vals = sapply(y.values.list, function(y.values) min(y.values,na.rm=TRUE))
# 		max.vals = sapply(y.values.list, function(y.values) max(y.values,na.rm=TRUE))
# 		# when there is only 1 opt sol, this causes NA values since we cannot obtain any dist score
# 		ylimit = c(min(min.vals,na.rm=TRUE), max(max.vals,na.rm=TRUE))
# 	}
# 
# 
# 	
# #	print(y.values.list)
# #	max.nb.elt = max(sapply(y.values.list, function(y.values){ max(table(y.values)) }))
# 	max.nb.elt = max(
# 		sapply(seq(1,length(y.values.list)),
# 		 function(i){
# 		 	y.values = y.values.list[[i]]
# 		 	x.values = x.values.list[[i]]
# 		 	if(all(is.na(x.values)) || all(is.na(y.values)))
# 		 		return(NA)
# 		 		
# 		 	# find the most frequent value
# 		 	freq = table(y.values)
# 		 	j = which.max(freq)
# 		 	freq.val = names(freq)[j]
# 		 	# look at the x criteria for the solutions that contain the most frequent value
# 		 	# for example: when dist score 0.1 is most frequent, which nb cluster is most frequent among them?
# 		 	indx = which(y.values == freq.val)
# 		 	return(max(table(x.values[indx]), na.rm=TRUE))
# 		 	
# 	}), na.rm=TRUE)
# 	
# 	# -------------------------------------------------
# 	# compute frequency intervals based on nb elt
# 	nb.color = length(heat.palette)
# 	
# 
# 	if(is.integer(max.nb.elt)){ # if NA, it becomes -Inf
# 		if(max.nb.elt == 1)
# 			freq.intervals = factor(max.nb.elt)
# 		else if(max.nb.elt < nb.color)
# 			freq.intervals = cut(seq(1, max.nb.elt), breaks = max.nb.elt)
# 		else # if(length(nb.elt) > 1)
# 			freq.intervals = cut(seq(1, max.nb.elt), breaks = nb.color, dig.lab=2)
# 		
# 		legend.intervals = sort(unique(freq.intervals))
# 	
# 		# -------------------------------------------------
# 	
# 
# 		layout.plot.args = list(xlabel=xlabel, ylabel=ylabel)
# 		layout.legend.args = list(title="Frequency intervals", legend=legend.intervals, col=heat.palette)
# 		plot.args = list(ylimit=ylimit, xlimit=xlimit, plot.type=SCATTER, x.values.list=x.values.list, xvalues=NA, ylabel=NA, xlabel=NA,
# 				 enable.plot.title=FALSE, max.nb.elt=max.nb.elt)
# 		legend.args = list(enabled=FALSE, title=NA, legend=NA)
# 		
# 		perform.generic.layout.plot.by(y.values.list, main.result.type, g.params, by=NA, param.values=NA,
# 				subplot.by, subplot.param.values, layout.plot.args, layout.legend.args, plot.args, legend.args)
# 	}
# }




#################################################################
#
# g.params
# by
# values
# subplot.by
# subplot.param.values
#
##################################################################
perform.generic.layout.plot.by = function(layout.data, result.type, g.params, by, param.values,
		subplot.by, subplot.param.values, layout.plot.args, layout.legend.args=NA, plot.args, legend.args=NA)
{	
	# -----------------------------------------------------
	# check if all values are NA, if so, do not plot
	na.res = c()
	#for(subplot.value in subplot.param.values){
	for(subplot.value in names(layout.data)){
		vector.list = layout.data[[as.character(subplot.value)]]
		na.res = c(na.res, all(is.na(vector.list)))
	}
	if(all(na.res))
		return()
	# -----------------------------------------------------

	# =============================
	g.params$k = NA # the folder related to solution class is not needed to show in the hierarchy of plot folders  
	# =============================
	
	
	plot.folder = get.plot.folder.path.by(g.params, by, subplot.by)
	print(plot.folder)
	if(!dir.exists(plot.folder))
		dir.create(path=plot.folder, showWarnings=FALSE, recursive=TRUE)
		
	plot.args$plot.title = NA
	plot.args$plot.folder = NA
	plot.args$plot.format = NA
	plot.args$plot.format = JUST.PLOT
	
	plot.filename = paste0(result.type, "-", by, "_vs_", subplot.by,"-", plot.args$plot.type, ".PDF")
	
	
	# a nice source about margins in plot: https://www.r-graph-gallery.com/74-margin-and-oma-cheatsheet/
	
	# ----------------------------------------------------
	# we need to adjust plot horizontal length when there is a lot of x values for multi bar and box plots
	# Let's say, there there are 4 plots in a line, and each plot has 10 x values. It is too much. Need to reduce nb plot per line
	# TODO we should prioritiz plot.type=MULTI.BAR
	nrow = 4
	ncol = 4
	if(length(subplot.param.values) <= 6){
		#nrow = 3
		#ncol = 2
		nrow = 2
		ncol = 3
	}
	else if(length(subplot.param.values) <= 9){
		nrow = 3
		ncol = 3
	}
	else if(length(subplot.param.values) <= 12 && length(param.values) > 9){
		nrow = 3
		ncol = 4
	}
	else if(length(subplot.param.values) <= 12 && length(param.values) < 9){
		nrow = 4
		ncol = 3
	}
	# ----------------------------------------------------
	
	#print(file.path(plot.folder, plot.filename))
	#pdf(file=file.path(plot.folder, plot.filename), compress=FALSE, width=5*ncol, height=2*nrow) # for k_init=2
	pdf(file=file.path(plot.folder, plot.filename), compress=FALSE, width=4*ncol, height=2*nrow) # for k_init=3
	#pdf(file=file.path(plot.folder, plot.filename), compress=FALSE, width=2.5*ncol, height=2*nrow) # for k_init=4
	#pdf(file=file.path(plot.folder, plot.filename), compress=FALSE, width=3*ncol, height=2*nrow) # for k_init=2
	#pdf(file=file.path(plot.folder, plot.filename), compress=FALSE, width=2*ncol, height=2*nrow) # for k_init=3

	par(mfrow = c(nrow, ncol)) # 2-by-2 grid of plots
#	par(oma = c(4, 4, 5, 0)) # make room (i.e. the 4's) for the overall x and y axis titles
#	par(mar = c(4, 3, 2.5, 1)) # make the plots be closer together

#	par(oma = c(2, 1, 0, 0)) # make room (i.e. the 4's) for the overall x and y axis titles
#	par(mar = c(1.5,3, 3, 1)) # make the plots be closer together
	par(oma = c(0, 0, 0, 0)) # make room (i.e. the 4's) for the overall x and y axis titles
	par(mar = c(2.5,3, 3, 1)) # make the plots be closer together
	
	#for(subplot.value in subplot.param.values){
	for(subplot.value in names(layout.data)){
		vector.list = layout.data[[as.character(subplot.value)]]
		#g.params = update.g.params.by(g.params=g.params, by=subplot.by, value=subplot.value)
		g.params = NA

		curr.plot.args = plot.args # reset at every iteration
				
		if(plot.args$plot.type == VIO){
			# ---------------------------------------------------------------
			# This fact is actually is not important, but when reading the plot, it is easy to understand
			# It specifies the values used in the boxplot is 'nb optimal solution', or something generic like 'value'
			vio.size.desc = ""
			if(startsWith(result.type, DIST.SCORE) || startsWith(result.type, ALL.DIST.SCORE) || startsWith(result.type, MEAN.DIST.SCORE))
				vio.size.desc = "val" # value
			else
				vio.size.desc = "sol" # nb solution
			# ---------------------------------------------------------------
			
#			nb.sol.vector = unlist(retreive.nb.partition(g.params, by, folders=param.values))
			nb.sol.vector = sapply(vector.list, function(vec) length(vec[which(!is.na(vec))]))
#			curr.plot.args$x.axis.names = paste0(curr.plot.args$x.axis.names, " (", nb.sol.vector, box.size.desc, ")") 
			curr.plot.args$x.axis.names = paste0(curr.plot.args$x.axis.names)
		} 

#		if(startsWith(result.type, RESULT.CORRESPONDANCE)){
#			curr.plot.args$xvalues = plot.args$x.values.list[[as.character(subplot.value)]]
#		}
		
#		print(vector.list)
#		print(file.path(plot.folder, plot.filename))
#		if(all(is.na(vector.list)))
#			plot.new()
#		else
		make.plot(vector.list, curr.plot.args, legend.args)
		#print(vector.list)
		
		short.plot.title = paste0(subplot.by,"=",subplot.value)
#		title(short.plot.title, cex.main=0.75, line = 0)
		title(short.plot.title, cex.main=1.1, line = 0.5)	

	}
	
#	# add an overall title
#	g.params = update.g.params.by(g.params=g.params, by=by, value=NA)
#	g.params = update.g.params.by(g.params=g.params, by=subplot.by, value=NA)
#	plot.title = create.plot.title(g.params)
#	title(plot.title, line = 0, outer = TRUE)
#	
#	# print the overall labels
#	mtext(paste0("x axis: ", layout.plot.args$xlabel), side = 1, outer = TRUE, line = 2) # x axis label
#	mtext(paste0("y axis: ", layout.plot.args$ylabel), side = 2, outer = TRUE, line = 2) # y axis label
#
#	# add an overall legend
#	if(!is.na(layout.legend.args)){
##		plot(1, type = "n", axes=FALSE, xlab="", ylab="")
#	
#
#		# source: https://stackoverflow.com/questions/24082485/r-plot-a-centered-legend-at-outer-margins-of-multiple-plots
#		reset <- function() {
#			par(mfrow=c(1, 1), oma=rep(0, 4), mar=rep(0, 4), new=TRUE)
#			plot(0:1, 0:1, type="n", xlab="", ylab="", axes=FALSE)
#		}
#		
#		reset()
#	
#		if(length(layout.legend.args$legend) < 20){
#			legend(x = "top",inset = 0, title = layout.legend.args$title,
#					legend = layout.legend.args$legend, 
#					col=layout.legend.args$col, lwd=2, cex=.35, horiz = TRUE) # do not change this setttings
#		} else {
#			legend(x = "top",inset = 0, title = layout.legend.args$title,
#					legend = layout.legend.args$legend, 
#					col=layout.legend.args$col, lwd=2, cex=.35, horiz = FALSE, ncol=length(layout.legend.args$legend)/2) # do not change this setttings
#		}
#	}

	dev.off()
	
	
}













#################################################################
#
# d: density
# k: nb cluster
#
##################################################################
perform.single.network.layout.plots.for.exact.algos = function(g.params, by, param.values, subplot.by, subplot.param.values,
		comp.measures, force, plot.format)
{
	tlog(20, "perform.layout.plots.for.exact.algos")
	
	#  ---- UPDATE g.params -----------------------------------------------------
	# TODO handle it in a cleaner way
	# param
	folder.value = ALL.PARAM.NAME
	if(by == ALGO.NAME)
		folder.value = ALL.COR.CLU.EXACT.ALGOS.PARAM.NAME
	g.params = update.g.params.by(g.params=g.params, by=by, value=folder.value)
	
	# subplot param
	folder.value = ALL.PARAM.NAME
	if(subplot.by == ALGO.NAME)
		folder.value = ALL.COR.CLU.EXACT.ALGOS.PARAM.NAME
	g.params = update.g.params.by(g.params=g.params, by=subplot.by, value=folder.value)
	# --------------------------------------------------------------------------
	
	result.type=NB.SOLUTION
	result.filename = paste0(csv.result.files[[result.type]],".csv")
	colname = csv.related.colnames[[result.type]]
	layout.data = retreive.layout.data.by(result.filename, colname, g.params, by, param.values, subplot.by, subplot.param.values)
	#print(layout.data)
	# perform.layout.line.plot.by(result.type, layout.data, g.params, by, param.values, subplot.by, subplot.param.values,
	# 		xvalues=param.values, xlabel=by, ylabel=prepare.axis.label(result.type), ylimit=NA)
	perform.layout.vio.plot.by(result.type, layout.data, g.params, by, param.values, subplot.by, subplot.param.values, ylimit=NA)


	result.type=EXEC.TIME
	result.filename = paste0(csv.result.files[[result.type]],".csv")
	colname = csv.related.colnames[[result.type]]
	layout.data = retreive.layout.data.by(result.filename, colname, g.params, by, param.values, subplot.by, subplot.param.values)
	#print(layout.data)
	# perform.layout.line.plot.by(result.type, layout.data, g.params, by, param.values, subplot.by, subplot.param.values,
	#                             xvalues=param.values, xlabel=by, ylabel=prepare.axis.label(result.type), ylimit=NA, y.axis.log.scale=TRUE)
	perform.layout.vio.plot.by(result.type, layout.data, g.params, by, param.values, subplot.by, subplot.param.values, ylimit=NA)


	result.type=IMB.PERC
	result.filename = paste0(csv.result.files[[result.type]],".csv")
	colname = csv.related.colnames[[result.type]]
	layout.data = retreive.layout.data.by(result.filename, colname, g.params, by, param.values, subplot.by, subplot.param.values)
	#print(layout.data)
	#perform.layout.line.plot.by(result.type, layout.data, g.params, by, param.values, subplot.by, subplot.param.values,
	#                            xvalues=param.values, ylabel=NA, ylimit=c(0,38)) ==> TODO bug
	perform.layout.vio.plot.by(result.type, layout.data, g.params, by, param.values, subplot.by, subplot.param.values, ylimit=c(0,38))


	result.type=IMB.COUNT
	result.filename = paste0(csv.result.files[[result.type]],".csv")
	colname = csv.related.colnames[[result.type]]
	layout.data = retreive.layout.data.by(result.filename, colname, g.params, by, param.values, subplot.by, subplot.param.values)
	#print(layout.data)
	#perform.layout.line.plot.by(result.type, layout.data, g.params, by, param.values, subplot.by, subplot.param.values,
	#                            xvalues=param.values, ylabel=NA, ylimit=c(0,38)) ==> TODO bug: x values are not integer
	perform.layout.vio.plot.by(result.type, layout.data, g.params, by, param.values, subplot.by, subplot.param.values, ylimit=NA)


	result.type=NB.CLUSTER
	result.filename = paste0(csv.result.files[[result.type]],".csv")
	colname = csv.related.colnames[[result.type]]
	layout.data = retreive.layout.data.by(result.filename, colname, g.params, by, param.values, subplot.by, subplot.param.values)
	#print(layout.data)
	#perform.layout.line.plot.by(result.type, layout.data, g.params, by, param.values, subplot.by, subplot.param.values,
	#                            xvalues=param.values, ylabel=NA, ylimit=c(0,38)) ==> TODO bug
	perform.layout.vio.plot.by(result.type, layout.data, g.params, by, param.values, subplot.by, subplot.param.values, ylimit=c(2,14))


    # TODO
	# perform.layout.multi.bar.plot.by(NB.CLUSTER, g.params,by, param.values, subplot.by, subplot.param.values,
	#                                  is.ctgry.names.on.x.axis=TRUE, ctgry.names=NA, xlabel=NA, ylabel=NA, use.stacked.bars=FALSE)
	#


	# new:  class core part size in core part analysis
	for(core.part.threshold in CORE.PART.THRESHOLDS){
		result.type=paste0(CLASS.CORE.PART.SIZE,"-tresh=",sprintf("%.2f",core.part.threshold))

		result.filename = paste0(csv.result.files[[CLASS.CORE.PART.SIZE]],"-tresh=",sprintf("%.2f",core.part.threshold),".csv")
		colname = csv.related.colnames[[CLASS.CORE.PART.SIZE]]
		layout.data = retreive.layout.data.by(result.filename, colname, g.params, by, param.values, subplot.by, subplot.param.values)
		#print(layout.data)
		perform.layout.vio.plot.by(result.type, layout.data, g.params, by, param.values, subplot.by, subplot.param.values,
				xlabel=by, ylabel=prepare.axis.label(result.type), ylimit=c(0,1))
	}


	for(measure in comp.measures){
	# do not use the for statement below, because it changes g.params values
		tlog("measure: ",measure)

		result.type=paste0(ALL.DIST.SCORE,"-",measure)
		result.filename = paste0(csv.result.files[[ALL.DIST.SCORE]],"-",measure,".csv")
		colname = paste0(csv.related.colnames[[ALL.DIST.SCORE]],"-",measure)
		layout.data = retreive.layout.data.by(result.filename, colname, g.params, by, param.values, subplot.by, subplot.param.values)
		#print(layout.data)
		perform.layout.vio.plot.by(result.type, layout.data, g.params, by, param.values, subplot.by, subplot.param.values,
			 #ylimit=NA) # ylimit=DIST.SCORE.INTERVAL
		 	ylimit=c(0,1)) # TODO this is a y limit for a normalized dissimilarity measure. For VI [0,1] does not fit !!!
 	}

	 
# 	
# 	# TODO: I think it is not meaningful
# ##	for(measure in comp.measures){
# ##		tlog("measure: ",measure)	
# ##		# x and y axis will be result type (i.e. nb cluster vs dist scores)
# ##		for(param.value in param.values){
# ##			tlog("by: ", by, " value: ", param.value)
# ##			g.params = update.g.params.by(g.params=g.params, by=by, value=param.value)
# ##			
# ##			# in order to perform 'result correspondance', x and y values should be of the same length
# ##			main.result.type = paste0(RESULT.CORRESPONDANCE,"-",NB.CLUSTER,"_vs_",paste0(MEAN.DIST.SCORE,"-",measure))	
# ##			perform.layout.scatter.plot.by(main.result.type, y.result.type=paste0(MEAN.DIST.SCORE,"-",measure),
# ##					 x.result.type=NB.CLUSTER, g.params, subplot.by, subplot.param.values,
# ##					ylabel=NA, xlabel=NA, ylimit=NA, xlimit=NA) # ylimit=DIST.SCORE.INTERVAL
# ##		}
# ##	}

}




#################################################################
#
# d: density
# k: nb cluster
#
##################################################################
perform.summary.network.layout.plots.for.exact.algos = function(g.params, by, param.values, subplot.by, subplot.param.values,
		comp.measures, force, plot.format)
{
	tlog(20, "perform.layout.plots.for.exact.algos")
	
	#  ---- UPDATE g.params -----------------------------------------------------
	# TODO handle it in a cleaner way
	# param
	folder.value = ALL.PARAM.NAME
	if(by == ALGO.NAME)
		folder.value = ALL.COR.CLU.EXACT.ALGOS.PARAM.NAME
	g.params = update.g.params.by(g.params=g.params, by=by, value=folder.value)
	
	# subplot param
	folder.value = ALL.PARAM.NAME
	if(subplot.by == ALGO.NAME)
		folder.value = ALL.COR.CLU.EXACT.ALGOS.PARAM.NAME
	g.params = update.g.params.by(g.params=g.params, by=subplot.by, value=folder.value)
	# --------------------------------------------------------------------------




	result.type=NB.SOLUTION # usually used for plot file name or shorthen description name for axis label
	result.filename = paste0(csv.result.files[[result.type]],".csv")
	colname = csv.related.colnames[[result.type]]

	layout.data = retreive.layout.data.by(result.filename, colname, g.params, by, param.values, subplot.by, subplot.param.values)
	print(layout.data)

	perform.layout.vio.plot.by(result.type, layout.data, g.params, by, param.values, subplot.by, subplot.param.values,
	                           xlabel=by, ylabel=prepare.axis.label(result.type), ylimit=c(0,500)) # ylimit=c(0,175)


	for(measure in comp.measures){
	    tlog(NB.SOLUTION.BY.SOL.CLASSES," => measure: ",measure)
	    
    	result.type=paste0(NB.SOLUTION.BY.SOL.CLASSES,"-",measure) # usually used for plot file name or shorthen description name for axis label
    	result.filename = paste0(csv.result.files[[NB.SOLUTION.BY.SOL.CLASSES]],".csv")
    	colname = csv.related.colnames[[NB.SOLUTION.BY.SOL.CLASSES]]
    	layout.data = retreive.layout.data.by.with.sol.classes(result.filename, colname, g.params, measure, by, param.values, subplot.by, subplot.param.values)
    	#print(layout.data)
    	# perform.layout.vio.plot.by(result.type, layout.data, g.params, by, param.values, subplot.by, subplot.param.values,
    	#                            xlabel=by, ylabel=prepare.axis.label(result.type), ylimit=c(0,500)) # ylimit=c(0,175)
    	perform.layout.vio.plot.by.with.sol.classes(result.type, layout.data, g.params, by, param.values, subplot.by, subplot.param.values,
    	                                             			xlabel=by, ylabel=prepare.axis.label(result.type), ylimit=NA)
    }

	result.type=EXEC.TIME # usually used for plot file name or shorthen description name for axis label
	result.filename = paste0(csv.result.files[[result.type]],".csv")
	colname = csv.related.colnames[[result.type]]
	layout.data = retreive.layout.data.by(result.filename, colname, g.params, by, param.values, subplot.by, subplot.param.values)
	#print(layout.data)
	perform.layout.vio.plot.by(result.type, layout.data, g.params, by, param.values, subplot.by, subplot.param.values,
	                           xlabel=by, ylabel=prepare.axis.label(result.type), ylimit=NA)




#     # TODO
# 	#result.type=NB.COMP.TRANSITION.GRAPH
# 	#perform.layout.multi.bar.with.prop.plot.by(result.type, g.params,by, param.values, subplot.by, subplot.param.values,
# 	#		is.ctgry.names.on.x.axis=FALSE,
# 	#		#ctgry.names=c(EDIT1.COMP.MEM.COL.NAME, EDIT12.COMP.MEM.COL.NAME, EDIT123.COMP.MEM.COL.NAME, EDIT1234.COMP.MEM.COL.NAME),
# 	#		ctgry.names=c(EDIT1.COMP.MEM.COL.NAME, EDIT1234.COMP.MEM.COL.NAME),
# 	#		xlabel=NA, ylabel=NA, ylimit=c(0,1))



	result.type=SINGLE.CLU.PROP.FOR.KMEDOIDS # usually used for plot file name or shorthen description name for axis label
	result.filename = paste0(csv.result.files[[result.type]],".csv")
	colname = csv.related.colnames[[result.type]]
	layout.data = retreive.layout.data.by(result.filename, colname, g.params, by, param.values, subplot.by, subplot.param.values)
	#print(layout.data)
	perform.layout.multi.bar.with.prop.plot.by(result.type, layout.data, g.params,by, param.values, subplot.by, subplot.param.values,
			is.ctgry.names.on.x.axis=FALSE,
			#ctgry.names=c(PROP.SINGLE.CLU.FOR.SILH.COL.NAME,PROP.SINGLE.CLU.FOR.SF.COL.NAME),
			ctgry.names=c(PROP.SINGLE.CLU.FOR.SILH.COL.NAME),
			xlabel=NA, ylabel=NA, ylimit=c(0,1))



	if(by != DETECTED.IMB.PROP.PARAM.NAME){
		result.type=IMB.PERC # usually used for plot file name or shorthen description name for axis label
		result.filename = paste0(csv.result.files[[result.type]],".csv")
		colname = csv.related.colnames[[result.type]]
		layout.data = retreive.layout.data.by(result.filename, colname, g.params, by, param.values, subplot.by, subplot.param.values)
		#print(layout.data)
		perform.layout.vio.plot.by(result.type, layout.data, g.params, by, param.values, subplot.by, subplot.param.values, ylimit=c(0,40))

		result.type=IMB.COUNT # usually used for plot file name or shorthen description name for axis label
		result.filename = paste0(csv.result.files[[result.type]],".csv")
		colname = csv.related.colnames[[result.type]]
		layout.data = retreive.layout.data.by(result.filename, colname, g.params, by, param.values, subplot.by, subplot.param.values)
		#print(layout.data)
		perform.layout.vio.plot.by(result.type, layout.data, g.params, by, param.values, subplot.by, subplot.param.values, ylimit=NA)
	}




	result.type=NB.CLUSTER # usually used for plot file name or shorthen description name for axis label
	result.filename = paste0(csv.result.files[[result.type]],".csv")
	colname = csv.related.colnames[[result.type]]
	layout.data = retreive.layout.data.by(result.filename, colname, g.params, by, param.values, subplot.by, subplot.param.values)
	#print(layout.data)
	perform.layout.vio.plot.by(result.type, layout.data, g.params, by, param.values, subplot.by, subplot.param.values, ylimit=c(2,14))




# 	# TODO: bug
# 	# 		if(by != DETECTED.IMB.PROP.PARAM.NAME){
# 	# 		perform.layout.multi.bar.plot.by(NB.CLUSTER, g.params,by, param.values, subplot.by, subplot.param.values,
# 	# 				is.ctgry.names.on.x.axis=TRUE, ctgry.names=NA, xlabel=NA, ylabel=NA, use.stacked.bars=FALSE)
# 	# 	}



	# new:  class core part size in core part analysis
	for(core.part.threshold in CORE.PART.THRESHOLDS){
		result.type=paste0(CLASS.CORE.PART.SIZE,"-tresh=",sprintf("%.2f",core.part.threshold)) # usually used for plot file name or shorthen description name for axis label
		result.filename = paste0(csv.result.files[[CLASS.CORE.PART.SIZE]],"-tresh=",sprintf("%.2f",core.part.threshold),".csv")
		colname = csv.related.colnames[[CLASS.CORE.PART.SIZE]]
		layout.data = retreive.layout.data.by(result.filename, colname, g.params, by, param.values, subplot.by, subplot.param.values)
		#print(layout.data)
		perform.layout.vio.plot.by(result.type, layout.data, g.params, by, param.values, subplot.by, subplot.param.values,
				xlabel=by, ylabel=prepare.axis.label(result.type), ylimit=c(0,1))
	}



	for(measure in comp.measures){
	    tlog("measure: ",measure)

    	# new: class core part size by sol classes in core part analysis
    	for(core.part.threshold in CORE.PART.THRESHOLDS){
    		result.type=paste0(CLASS.CORE.PART.SIZE.BY.SOL.CLASSES,"-tresh=",sprintf("%.2f",core.part.threshold)) # usually used for plot file name or shorthen description name for axis label
    		result.filename = paste0(csv.result.files[[CLASS.CORE.PART.SIZE.BY.SOL.CLASSES]],"-tresh=",sprintf("%.2f",core.part.threshold),".csv")
    		colname = csv.related.colnames[[CLASS.CORE.PART.SIZE.BY.SOL.CLASSES]]
    		layout.data = retreive.layout.data.by.with.sol.classes(result.filename, colname, g.params, measure, by, param.values, subplot.by, subplot.param.values)
    		#print(layout.data)
    		perform.layout.vio.plot.by.with.sol.classes(result.type, layout.data, g.params, by, param.values, subplot.by, subplot.param.values,
    		                                            xlabel=by, ylabel=prepare.axis.label(result.type), ylimit=c(0,1))
    	}
	}



	for(measure in comp.measures){
		tlog("measure: ",measure)

		result.type=paste0(ALL.DIST.SCORE,"-",measure) # usually used for plot file name or shorthen description name for axis label
		result.filename = paste0(csv.result.files[[ALL.DIST.SCORE]],"-",measure,".csv")
		print(result.filename)
		colname = paste0(csv.related.colnames[[ALL.DIST.SCORE]],"-",measure)
		layout.data = retreive.layout.data.by(result.filename, colname, g.params, by, param.values, subplot.by, subplot.param.values)
		#print(layout.data)
		perform.layout.vio.plot.by(result.type, layout.data, g.params, by, param.values, subplot.by, subplot.param.values,
			 #ylimit=NA) # ylimit=DIST.SCORE.INTERVAL
		 	ylimit=c(0,1)) # TODO this is a y limit for a normalized dissimilarity measure. For VI [0,1] does not fit !!!
	}

}




#################################################################
#
#################################################################
inner.make.layout.plots = function(initial.g.params, by, param.values, subplot.by, subplot.param.values, comp.measures,	force, plot.format){
    
    result = check.and.adjust.upper.bounds(initial.g.params, by, param.values, subplot.by, subplot.param.values)
    print(result)


    if( !(by == NETWORK.NO || subplot.by == NETWORK.NO) && initial.g.params$network.no == SUMMARY.FOLDER.NAME){

        # by vs. subplot.by ==> ex: GRAPH.SIZE vs PROP.MISPL
        perform.summary.network.layout.plots.for.exact.algos(initial.g.params, by=by, param.values=result$param.values,
                                subplot.by=subplot.by, subplot.param.values=result$subplot.param.values, comp.measures,	force, plot.format)

        # # subplot.by vs. by ==> ex: PROP.MISPL vs GRAPH.SIZE
        perform.summary.network.layout.plots.for.exact.algos(initial.g.params, by=subplot.by, param.values=result$subplot.param.values,
                                subplot.by=by, subplot.param.values=result$param.values, comp.measures,	force, plot.format)

    } else {

        # by vs. subplot.by ==> ex: GRAPH.SIZE vs PROP.MISPL
        perform.single.network.layout.plots.for.exact.algos(initial.g.params, by=by, param.values=result$param.values,
                                subplot.by=subplot.by, subplot.param.values=result$subplot.param.values, comp.measures,	force, plot.format)

        # # subplot.by vs. by ==> ex: PROP.MISPL vs GRAPH.SIZE
        perform.single.network.layout.plots.for.exact.algos(initial.g.params, by=subplot.by, param.values=result$subplot.param.values,
                                subplot.by=by, subplot.param.values=result$param.values, comp.measures,	force, plot.format)

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
# cor.clu.exact.algos
# in.rand.g.folders: input random graph folders. Sequantial integers (1, .., 10)
# force
# plot.formats
#
##################################################################
make.layout.plots.by.graph.size.and.imbalance = function(graph.sizes, d, l0, prop.mispls, detected.imb.intervals, prop.negs, in.rand.net.folders,
		cor.clu.exact.algo, comp.measures, force, plot.format)
{		
	tlog("starts making plots at network level")	
    
    # note that if density=1, prop.negs is supposed to be NA. When it is NA, we assign an appropriate value to prop.negs in 'retreive.csv.results()' in 'retreive-results.R'
    for(prop.neg in prop.negs){
        tlog(12, "making plots at network level => prop.neg: ", prop.neg)
        
    	for(network.no in in.rand.net.folders){
    		tlog(16, "making plots at instance level => network.no: ", network.no)
    		
    		for(graph.desc.name in c(SIGNED.UNWEIGHTED.FILE)){ # OPPOSITE.SIGNED.UNWEIGHTED.FILE
    			tlog(16, "evaluating partitions => graph.desc.name: ", graph.desc.name)
    			tlog(20, "making layout plots => algo.name: ", cor.clu.exact.algo)
    			
    			initial.g.params = build.g.params(NA, l0, d, NA, NA, prop.neg, k=ALL, network.no, cor.clu.exact.algo, graph.desc.name)
    			
    			inner.make.layout.plots(initial.g.params, by=PROP.MISPL, param.values=prop.mispls,
    			                        subplot.by=GRAPH.SIZE, subplot.param.values=graph.sizes, comp.measures, force, plot.format)
    			
    			inner.make.layout.plots(initial.g.params, by=DETECTED.IMB.PROP.PARAM.NAME, param.values=detected.imb.intervals,
    			                        subplot.by=GRAPH.SIZE, subplot.param.values=graph.sizes, comp.measures, force, plot.format)
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
# cor.clu.exact.algos
# in.rand.g.folders: input random graph folders. Sequantial integers (1, .., 10)
# force
# plot.formats
#
##################################################################
make.layout.plots.by.graph.size.and.network.no = function(graph.sizes, d, l0, prop.mispls, detected.imb.intervals, prop.negs, in.rand.net.folders,
		cor.clu.exact.algo, comp.measures, force, plot.format)
{		
	tlog("starts making plots at network level")	
	
	for(prop.neg in prop.negs){
		tlog(12, "making plots at network level => prop.neg: ", prop.neg)
		
		for(prop.mispl in prop.mispls){
			tlog(16, "making plots at instance level => prop.mispl: ", prop.mispl)

			for(graph.desc.name in c(SIGNED.UNWEIGHTED.FILE)){ # OPPOSITE.SIGNED.UNWEIGHTED.FILE
				tlog(16, "evaluating partitions => graph.desc.name: ", graph.desc.name)

				# only exact solutions
				tlog(16, "making layout plots for exact solutions")
				#for(algo.name in cor.clu.exact.algos){
					tlog(20, "making layout plots => algo.name: ", cor.clu.exact.algo)

					initial.g.params = build.g.params(NA, l0, d, prop.mispl, NA, prop.neg, k=ALL, NA, cor.clu.exact.algo, graph.desc.name)

					inner.make.layout.plots(initial.g.params, by=NETWORK.NO, param.values=in.rand.net.folders,
					                        subplot.by=GRAPH.SIZE, subplot.param.values=graph.sizes, comp.measures, force, plot.format)
			}
		}
	    
	    # =======
	    
	    for(detected.imb.interval in detected.imb.intervals){
	        tlog(16, "making plots at instance level => detected.imb.interval: ", detected.imb.interval)
	        
	        for(graph.desc.name in c(SIGNED.UNWEIGHTED.FILE)){ # OPPOSITE.SIGNED.UNWEIGHTED.FILE
	            tlog(16, "evaluating partitions => graph.desc.name: ", graph.desc.name)
	            
	            # only exact solutions
	            tlog(16, "making layout plots for exact solutions")
	            #for(algo.name in cor.clu.exact.algos){
	            tlog(20, "making layout plots => algo.name: ", cor.clu.exact.algo)
	            
	            initial.g.params = build.g.params(NA, l0, d, NA, detected.imb.interval, prop.neg, k=ALL, NA, cor.clu.exact.algo, graph.desc.name)
	            
	            # initial 'in.rand.net.folders' makes no sense regarding 'detected.imb.intervals': there might be more than 10 even though the limit was 10.
	            my.in.rand.net.folders = seq(1,1000) # put something large, it will be adjusted later
	            inner.make.layout.plots(initial.g.params, by=NETWORK.NO, param.values=my.in.rand.net.folders,
	                                    subplot.by=GRAPH.SIZE, subplot.param.values=graph.sizes, comp.measures, force, plot.format)
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
# cor.clu.exact.algos
# in.rand.g.folders: input random graph folders. Sequantial integers (1, .., 10)
# force
# plot.formats
#
##################################################################
make.layout.plots.by.imbalance.and.network.no = function(graph.sizes, d, l0, prop.mispls, detected.imb.intervals, prop.negs, in.rand.net.folders,
		cor.clu.exact.algos, comp.measures, force, plot.format)
{		
	tlog("starts making layout plots at network level")	
	for(n in graph.sizes){
		tlog(4, "making plots at instance level => n: ", n)
		
		for(prop.neg in prop.negs){
			tlog(12, "making layout plots at network level => prop.neg: ", prop.neg)
			
			for(graph.desc.name in c(SIGNED.UNWEIGHTED.FILE)){ # OPPOSITE.SIGNED.UNWEIGHTED.FILE
				tlog(16, "evaluating partitions => graph.desc.name: ", graph.desc.name)
				
				# only exact solutions
				tlog(16, "making layout plots for exact solutions")
				for(algo.name in cor.clu.exact.algos){
					tlog(20, "making layout plots => algo.name: ", algo.name)
					
					initial.g.params = build.g.params(n, l0, d, NA, NA, prop.neg, k=ALL, NA, algo.name, graph.desc.name)
					
					# initial 'in.rand.net.folders' makes no sense regarding 'detected.imb.intervals': there might be more than 10 even though the limit was 10.
					my.in.rand.net.folders = seq(1,1000) # put something large, it will be adjusted later
					
					inner.make.layout.plots(initial.g.params, by=NETWORK.NO, param.values=my.in.rand.net.folders,
					                        subplot.by=PROP.MISPL, subplot.param.values=prop.mispls, comp.measures, force, plot.format)
					
					inner.make.layout.plots(initial.g.params, by=NETWORK.NO, param.values=my.in.rand.net.folders,
					                        subplot.by=DETECTED.IMB.PROP.PARAM.NAME, subplot.param.values=detected.imb.intervals,
					                        comp.measures, force, plot.format)
				}
			
				
				
			}
		}
	}
	
}










