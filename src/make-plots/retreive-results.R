





#################################################################
# It reads the values regarding solution classes of a set of networks from the file EVAL.BEST.K.FOR.KMEDOIDS.
#
# An example of file: EVAL.BEST.K.FOR.KMEDOIDS
#             Best k for Silhouette
# network1 sol0	    3
# network2 sol0	    NA ====> NA means there is only 1 partition. So, ignored.
# network3 sol0	    1 ====> this means there are at least 2 partitions, and the found nb solution class is 1
# network4 sol0	    1
# network5 sol0	    NA
# network6 sol0	    1
# network7 sol0	    2
# network8 sol0	    2
# network9 sol0	    2
# network10 sol0	3
#
#################################################################
retreive.sol.class.values = function(g.params, measure){
    
    # upper.eval.folder = get.eval.folder.path(g.params$n, g.params$k, g.params$d, g.params$prop.mispl, g.params$prop.neg, NA, NA, NA, NA)
    # existing.dirs = list.dirs(path = upper.eval.folder, full.names = FALSE, recursive = FALSE)
    # indx = which(startsWith(existing.dirs,"k="))
    # sol.class.values = as.integer(gsub("k=","",existing.dirs[indx]))
    
    result.filename = paste0(EVAL.BEST.K.FOR.KMEDOIDS,"-",measure,".csv")
    colname = BEST.K.FOR.SILH.COL.NAME
    raw.sol.class.values = retreive.csv.result(g.params, result.filename, colname)
    sol.class.values = unique(raw.sol.class.values) # it is a matrix with 1 column and several rows
    sol.class.values = as.vector(sol.class.values)
    return(sol.class.values)
}

#################################################################
# It reads a column vector from csv result file. There is normally at least one column. Each column is normally named.
#   If file does not exist, then it returns matrix(NA). Otherwise, it converts a vector to matrix of 1 column, then returns it.
#   Sinec the number of solution class is not specifie, k=ALL here.
#
# g.params 
# result.filename
# colname
#
##################################################################
retreive.csv.result = function(g.params, result.filename, colname)
{
    eval.folder = get.eval.folder.path(g.params$n, g.params$l0, g.params$d, g.params$prop.mispl, g.params$detected.imb.interval, g.params$prop.neg, g.params$k,
                                       g.params$network.no, g.params$algo.name, g.params$graph.desc.name)
    #print(eval.folder)
    #print(result.filename)
    table.file = file.path(eval.folder, result.filename)
    #print(table.file)
    
    # -------------------------------
    values = matrix(NA) # before, this was used
    # -------------------------------
    
    if(file.exists(table.file)){
        df = as.matrix(read.csv(table.file, row.names = 1, header= TRUE, check.names=FALSE, stringsAsFactors=FALSE))
        values = df[,colname]
        values = matrix(as.numeric(values),ncol=length(colname))
        rownames(values) = rownames(df)
    }
    
    return(values)
}



#################################################################
# It is similar to the function "retreive.csv.result()". 
# The only difference is that it reads the file from the folder specific to number of solution classes.
# Example: g.params$k = 2
#
# g.params
# result.filename
# colname
#
##################################################################
retreive.csv.result.by.sol.classes = function(g.params, result.filename, colname, measure)
{
    vector.list = list()
    sol.class.values = retreive.sol.class.values(g.params, measure)
    
    # ===========
    # remove NA value if it exists
    na.indx = which(is.na(sol.class.values))
    if(length(na.indx)>0)
        sol.class.values = sol.class.values[-na.indx]
    if(length(sol.class.values)>0)
        sol.class.values = sol.class.values[order(as.integer(sol.class.values))] # sort the values related to the number of solution classes
    # ===========

    if(length(sol.class.values)>0){
        for(k in sol.class.values){
            g.params$k = paste0("k=",k,"-",measure)
            result = retreive.csv.result(g.params, result.filename, colname)
            vector.list[[g.params$k]] = result
        } 
    }else {
        return(list())
    }
    
    print(vector.list)
    
    return(vector.list)
}


#################################################################
# It iterates over the variable 'folders' (e.g. graph sizes) and then retreives the results and writes into a list.
# Since we dont use any number of solution classes, g.params$k = ALL here.
#
# Examples of output: a list of matrix values (a matrix of 1 column)
# Note that we use a workaround ehre for 'prop.neg' is NA (i.e. density=1).
#   When we plot a layout plot with graph sizes and/or prop.mispl, etc.
#   This blocks us to compute 'prop.neg' through 'compute.prop.neg()'.
#   Since we have access to all graph params (except prop.neg), we can compue here.

#
#
# g.params
# result.filename, 
# colname
# by=NA
# folders=NA
#
#
#
# result.type = NB.SOLUTION
# Example 1
# by:  Prop mispl 
# param.values 0.05 0.1 0.15 0.20

# [[1]] ===> Prop mispl = 0.05
#                     [,1]
#     network1 sol0     1
#     network2 sol0     2
#     network3 sol0     2
#     network4 sol0     1
#     network5 sol0     2
#     network6 sol0     5
#     network7 sol0     2
#     network8 sol0     1
#     network9 sol0     4
#     network10 sol0    1
# 
# [[2]] ===> Prop mispl = 0.10
#                     [,1]
#     network1 sol0    10
#     network2 sol0    40
#     network3 sol0     6
#     network4 sol0    12
#     network5 sol0     6
#     network6 sol0     2
#     network7 sol0     2
#     network8 sol0     4
#     network9 sol0     2
#     network10 sol0    4
# 
# [[3]] ===> Prop mispl = 0.15
#                     [,1]
#     network1 sol0    35
#     network2 sol0    53
#     network3 sol0    31
#     network4 sol0    14
#     network5 sol0     1
#     network6 sol0    15
#     network7 sol0     6
#     network8 sol0     2
#     network9 sol0     1
#     network10 sol0    6
# 
# [[4]] ===> Prop mispl = 0.20
#                     [,1]
#     network1 sol0     1
#     network2 sol0    16
#     network3 sol0     8
#     network4 sol0    32
#     network5 sol0     1
#     network6 sol0    48
#     network7 sol0    13
#     network8 sol0     1
#     network9 sol0     1
#     network10 sol0  319
#
##################################################################
retreive.csv.results = function(g.params, result.filename, colname, by=NA, folders=NA)
{	
    # nb.result=1
    # if(!is.na(folders))
    #     nb.result = length(folders)
    
    vector.list = list()
    for(i in seq(1,length(folders))){ # if folders=NA, then by is also NA at the same time
        eval.folder = NA
        if(!is.na(folders)){
            folder.name = folders[i] # in in case of class type 'vector' and 'list'
            if(is.list(folder.name)) # if 'folders' is a list object, then 'folder' is a list with 1 element 
                folder.name = unlist(folder.name)[1]
            g.params = update.g.params.by(g.params, by, folder.name)
        }
        
        # -------------------------------------------------------
        # workaround: since when density=1, prop.neg=NA.
        # So, we need to handle this exceptional case here since prop.neg is changed in terms of 'prop.mispl' (and 'n')
        if(g.params$d == 1)
            g.params$prop.neg = compute.prop.neg(g.params$n, g.params$d, g.params$l0, g.params$prop.mispl)
        # -------------------------------------------------------
        
        
        values = retreive.csv.result(g.params, result.filename, colname)
        vector.list[[i]] = values
    }
    return(vector.list)
}



#################################################################
# It is similar to 'retreive.csv.results()'. The only exception is that we retreive all the values of solution classes.
#   So, we have in the end a list of a list of matrix values.
#
#
# g.params 
# result.filename 
# colname 
# by 
# folders
#
#
#
# Example of output:
#
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
#
##################################################################
retreive.csv.results.by.sol.classes = function(g.params, result.filename, colname, measure, by=NA, folders=NA)
{	
    # nb.result=1
    # if(!is.na(folders))
    #     nb.result = length(folders)
    
    vector.list = list()
    for(i in seq(1,length(folders))){ # if folders=NA, then by is also NA at the same time
        eval.folder = NA
        if(!is.na(folders)){
            folder.name = folders[i] # in in case of class type 'vector' and 'list'
            if(is.list(folder.name)) # if 'folders' is a list object, then 'folder' is a list with 1 element 
                folder.name = unlist(folder.name)[1]
            g.params = update.g.params.by(g.params, by, folder.name)
        }
        
        # -------------------------------------------------------
		# workaround: since when density=1, prop.neg=NA.
		# So, we need to handle this exceptional case here since prop.neg is changed in terms of 'prop.mispl' (and 'n')
        if(g.params$d == 1)
            g.params$prop.neg = compute.prop.neg(g.params$n, g.params$d, g.params$l0, g.params$prop.mispl)
		# -------------------------------------------------------
        
        values = retreive.csv.result.by.sol.classes(g.params, result.filename, colname, measure)
        vector.list[[i]] = values
    }
    return(vector.list)
}


