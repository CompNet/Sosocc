#################################################################
# Correlation Clustering (CC) problem
#################################################################

# ===============================================================
# Heuristic Approach: ILS and GRASP
# ===============================================================

CODE.COR.CLU.ILS = "ILS"
COR.CLU.ILS <- "ILS"
COR.CLU.ILS.CC <- paste(COR.CLU.ILS,"CC",sep="-")
CODE.COR.CLU.GRASP = "GRA"
COR.CLU.GRASP <- "GRASP"
COR.CLU.GRASP.CC <- paste(COR.CLU.GRASP,"CC",sep="-")

CODE.COR.CLU.VNS = "VNS"
COR.CLU.VNS = "VNS"
COR.CLU.VNS.CC = paste(COR.CLU.VNS,"CC",sep="-")
COR.CLU.VNS.RCC = paste(COR.CLU.VNS,"RCC",sep="-")

CODE.COR.CLU.VOTE.BOEM = "V-B"
COR.CLU.VOTE.BOEM = "VOTE-BOEM"
COR.CLU.VNS.VOTE.BOEM.CC = paste(COR.CLU.VOTE.BOEM,"CC",sep="-")
COR.CLU.VNS.VOTE.BOEM.RCC = paste(COR.CLU.VOTE.BOEM,"RCC",sep="-")

ILS.GRASP.LIB.FOLDER = file.path(LIB.FOLDER,paste(COR.CLU.ILS,COR.CLU.GRASP,sep="-"))

# ===============================================================
# Exact Approach: ExCC
# ===============================================================

COR.CLU.ExCC <- "ExCC"
COR.CLU.ExCC.ENUM.ALL <- "ExCC-all" # TODO
ExCC.LIB.FOLDER = file.path(LIB.FOLDER,COR.CLU.ExCC)
ExCC.JAR.PATH = paste(ExCC.LIB.FOLDER,"cplex-partition.jar",sep="/") # gaia cluster - CERI
#CPLEX.BIN.PATH = "/users/narinik/Cplex/ILOG/CPLEX_Studio128/cplex/bin/x86-64_linux/"
CPLEX.BIN.PATH = "/opt/ibm/ILOG/CPLEX_Studio128/cplex/bin/x86-64_linux/" # in my personel computer



# ---------------------------------------------------------------
# ---------------------------------------------------------------


#############################################################################################
#
#############################################################################################
get.k.from.algo.name = function(token, part.folder){
    # example1: kFrom=(ALGO.NAME):k+1
    # example2: k=?
    k = NA
    
    tmp2 <- strsplit(x=token, split="=", fixed=TRUE)[[1]] # c("kFrom=ILS-RCC:k+1", "..")
    if(tmp2[1] == "k"){ # the 'k' is already provided by user
        k = tmp2[2]
    } else if(tmp2[1] == "kFrom"){ # the 'k' is not provided by user, but we get this info through the result of an algo
        tmp3 <- strsplit(x=tmp2[2], split=":", fixed=TRUE)[[1]] # c("ILS-RCC","k+1")
        list.chars = strsplit(tmp3[1],'')[[1]]
        # remove the paranthesis located in the begining and the end
        k.from.algo.name <- paste(list.chars[2:(length(list.chars)-1)], collapse="") # except the 1st and the last character which are paranthesis
        
        partition.file = file.path(part.folder,paste0(k.from.algo.name,"-membership.txt"))
        
        
        if(!file.exists(partition.file)){
            tlog("........Partition file ",partition.file," not found")
            return(-1)
        }
        else{
            partition <- as.matrix(read.table(partition.file))
            k = length(unique(partition))
        }
        
        if(grepl('\\+', tmp3[2])){ # if it contains the '+' sign ==> 2 possiblities: "k" or "k+<NUMBER>" like "k+2"
            tmp4 <- strsplit(x=tmp3[2], split="+", fixed=TRUE)[[1]][2]
            k.incr.nbr = as.integer(tmp4)
            k = k + k.incr.nbr
        }
    }
    
    return(k)
}


#############################################################################################
# This method is especially designed when the "k" parameter is needed for execution of the algo.
# But, it handles with the other case where there is not any "k" value needed
# There are 2 input options when user needs to provide "k":
# 1): kFrom=(ALGO.NAME):k+1			====> take the k value of an algo
# 2) k=<NUMBER>                     ====> take the provided k value
#############################################################################################
break.down.algo.code = function(algo.name){
    # break down the specified code (short name)
    # ================================================================================================================
    tmp=NA
    if(grepl('\\(', algo.name) && grepl('\\)', algo.name)){ # if "(" and ")" paranthesis symbols are used
        list.chars = strsplit(algo.name, "")[[1]] # get a list of characters of the string 'alog.name'
        start.paranthesis = which(list.chars=="(")
        end.paranthesis = which(list.chars==")")
        # substiture the content of the paranthesis part with the word 'REPLACE'.
        # Thus, it will be easy to split the 'algo.name' by the symbol "_"
        algo.name2 = paste(c(list.chars[1:(start.paranthesis-1)], "REPLACE", list.chars[(end.paranthesis+1):length(list.chars)]), collapse="")
        tmp <- strsplit(x=algo.name2, split="_", fixed=TRUE)[[1]]
        
        # find the index whose item contains "REPLACE"
        indx = which(grepl("REPLACE", tmp))
        tmp[indx] = gsub("REPLACE", paste(list.chars[start.paranthesis:end.paranthesis], collapse=""), tmp[indx]) 
        
    } else {
        tmp <- strsplit(x=algo.name, split="_", fixed=TRUE)[[1]]
    }
    
    # At the end, the result of the example "ILS-RCC_kFrom=(ILS-CC_l1_a1_g0_p3_t3_i10):k+1_l1_a1_g0_p3_t3_i10"
    #	[1] "ILS-RCC"                              
    #	[2] "kFrom=(ILS-CC_l1_a1_g0_p3_t3_i10):k+1"
    #	[3] "l1"                                   
    #	[4] "a1"                                   
    #	[5] "g0"                                   
    #	[6] "p3"                                   
    #	[7] "t3"                                   
    #	[8] "i10"
    # ================================================================================================================
    
    return(tmp) # list of parameters
}






############################################################################
# It reads a .G graph file, and returns the contents as a data frame object.
#
# network.path: the file path which stores the .G graph file
#
############################################################################
read.graph.ils.file.as.df = function(network.path){
	# skip the first line bc it does not contain graph info
	df = read.table(
			file=network.path, 
			header=FALSE, 
			sep="\t", 
			skip=1, 
			check.names=FALSE
	)
	# df$V1: vertex1
	# df$V2: vertex2
	# df$V3: weight
	return(df)
}



############################################################################
#  It reads a .G graph file, and returns the contents as a igraph graph object.
#  To handle isolated nodes, first we had to find the max vertex id.
#  Then, indicate explicitely vertices ids in graph.data.frame()
#
# network.path: the file path which stores the .G graph file
#
############################################################################
read.graph.ils = function(network.path){
	df = read.graph.ils.file.as.df(network.path)
	
	edg.list = df[,c(1, 2)]
	max.v.id = max(unique(c(edg.list[,1], edg.list[,2])))
	
	g <- graph.data.frame(edg.list, vertices=seq(0,max.v.id), directed=FALSE)
	cat("max id: ",max.v.id, "\n")
	E(g)$weight = df[, 3]
	# V(g)$id = seq(0,max.v.id)
	# V(g)$id = seq(1,max.v.id+1)
	
	return(g)
}


############################################################################
# It writes the graph object into a file
#
# graph: the igraph graph object
# file path: the file path which will store the graph content in the format of .G graph file
############################################################################
write.graph.ils = function(graph, file.path){
# export using a format compatible with pILS
	t <-get.edgelist(graph=graph)
	t =  matrix(as.integer(t), nrow(t), ncol(t))
	if(t[1,1] == 1)
		t <- t - 1	# start numbering nodes at zero
	
	t <- cbind(t,E(graph)$weight)		# add weights as the third column
	
	write.table(data.frame(vcount(graph),ecount(graph)), file=file.path, append=FALSE, sep="\t", row.names=FALSE, col.names=FALSE) # write header
	write.table(t, file=file.path, append=TRUE, sep="\t", row.names=FALSE, col.names=FALSE) # write proper graph
}








#############################################################################################
# Returns the code (short name) for the Iterated Local Search (ILS) partioning method. See 
# the algorithm documentation for more details.
#
# rcc: whether to solve the correlation clustering (FALSE) or relaxed CC problem (TRUE).
# l: neighborhood size (during the local search).
# k: number of clusters (max number for RCC)
# k.from: the algo name from which we use the 'k' value. for RCC problem
# rel.k.val: the relative 'k' value. It is used with "k.from" parameter. Ex: "k", "k+2", etc.
# alpha: randomness factor.
# gain: 0 for min imbalance, 
#       1 for max modularity gain function, 
#       2 for max negative modularity gain function, 
#       3 for max positive-negative modularity gain function, 
#       4 for max pos-neg mod gain function II, 
#       5 for max pos-neg mod gain function III
# perturbation: maximal level of perturbation.
# vns: enbale VNS metaheuristic
#
# returns: the short name corresponding to the ILS method with the specified parameters.
#############################################################################################
get.ils.code <- function(l, alpha, rcc, gain, perturbation, time.limit=3600, iter.nbr=10,
                         vns=FALSE, k=NA, k.from=NA, rel.k.val=NA)
{	result <- NA

if(rcc){
    # example1: ILS-RCC_kFrom=ILS-RCC:k+1_....
    # example2: ILS-RCC_k=2_....
    
    if(vns){
        result <- COR.CLU.VNS.RCC
    } else {
        
        result <- COR.CLU.ILS.RCC
        if(is.na(k))
            result <- paste0(result,"_kFrom=(",k.from,"):",rel.k.val)
        else
            result <- paste0(result,"_k=",k)
    }
}
else{
    if(vns){
        result <- COR.CLU.VNS.CC
    } else {
        result <- COR.CLU.ILS.CC
    }
    
}


result <- paste0(result,"_l",l)
result <- paste0(result,"_a",alpha)
result <- paste0(result,"_g",gain)
result <- paste0(result,"_p",perturbation)
result <- paste0(result,"_t",time.limit)
result <- paste0(result,"_i",iter.nbr)

return(result)
}






#############################################################################################
# Returns the inline command for the Iterated Local Search (ILS) partioning method. See 
# the algorithm documentation for more details.
#
# algo.name: short code associated to the algorithm.
# input.folder: relative path to the folder containing targeted graph file.
# out.folder: relative path to the folder in which the produced files will be placed.
# k: nb cluster to be detected for RCC (when RCC is enbaled)
#
# returns: the command allowing to invoke the program externally.
#############################################################################################
get.ils.command <- function(algo.name, input.folder, out.folder, graph.name)
{	
    # break down the specified code (short name)
    tmp = break.down.algo.code(algo.name)
    
    
    params <- c()
    k=NA
    next.param.indx=NA
    # =======================================================================================
    # rcc flag
    # example1: ILS-RCC_kFrom=(ILS-CC_l1_a1_g0_p3_t3_i10):k+1_l1_a1_g0_p3_t3_i10
    # example2: ILS-RCC_k=2_l1_a1_g0_p3_t3_i10
    algo.name <- strsplit(x=tmp[1], split="-", fixed=TRUE)[[1]][1]
    rcc.flag <- strsplit(x=tmp[1], split="-", fixed=TRUE)[[1]][2]
    if(rcc.flag=="RCC"){
        params["rcc"] <- 1
        k = get.k.from.algo.name(tmp[2], out.folder)
        next.param.indx = 3
    }
    else{
        params["rcc"] <- 0
        next.param.indx = 2
    }
    
    # VNS meta heuristic
    if(algo.name == COR.CLU.VNS)
        params["vns"]=1
    else if(algo.name == COR.CLU.ILS)
        params["vns"]=0
    
    # =======================================================================================
    
    
    
    for(s in tmp[next.param.indx:length(tmp)])
    {	params <- c(params,substr(s,2,nchar(s)))
    names(params)[length(params)] <- substr(s,1,1)
    }
    
    
    # init
    input.file <- file.path(input.folder, graph.name)
    output.folder <- file.path(out.folder)
    command.folder <- file.path(ILS.GRASP.LIB.FOLDER)
    result <- file.path(command.folder, "graspcc")
    #result <- paste0("mpirun -n 1 ",result)
    
    
    
    # build command
    result <- paste0(result, " --vns ",params["vns"])
    result <- paste0(result, " --rcc ",params["rcc"])
    if(rcc.flag=="RCC")
        result <- paste0(result, " --k ",k)
    result <- paste0(result, " --alpha ",params["a"])
    result <- paste0(result, " --iterations ",params["i"])
    result <- paste0(result, " --neighborhood_size ",params["l"])
    result <- paste0(result, " --time-limit ",params["t"])
    result <- paste0(result, " --input-file \"",input.file,"\"")
    result <- paste0(result, " --output-folder \"",output.folder,"\"")
    result <- paste0(result, " --gain-function-type ",params["g"])
    result <- paste0(result, " --strategy ","ILS") 
    result <- paste0(result, " --perturbationLevelMax ",params["p"])
    
    return(result)
}


#############################################################################################
# Returns the code (short name) for the Grasp partioning method. See the algorithm documentation
# for more details.
#
# rcc: whether to solve the correlation clustering (FALSE) or relaxed CC problem (TRUE).
# l: neighborhood size (during the local search).
# k: number of clusters (max number for RCC)
# alpha: randomness factor.
# gain: 0 for min imbalance, 
#       1 for max modularity gain function, 
#       2 for max negative modularity gain function, 
#       3 for max positive-negative modularity gain function, 
#       4 for max pos-neg mod gain function II, 
#       5 for max pos-neg mod gain function III
# perturbation: maximal level of perturbation.
#
# returns: the short name corresponding to the Grasp method with the specified parameters.
#############################################################################################
get.grasp.code <- function(rcc, l, k, alpha, gain, time.limit, iter.nbr, oneOptNeig)
{	
    # TODO k can be specified for GRASP as well?
    
    result <- COR.CLU.GRASP
    if(!oneOptNeig && iter.nbr == -1)
        result <- COR.CLU.VOTE.BOEM
    
    if(rcc)
        result <- paste0(result,"-RCC")
    else
        result <- paste0(result,"-CC")
    
    
    result <- paste0(result,"_l",l)
    result <- paste0(result,"_a",alpha)
    result <- paste0(result,"_g",gain)
    result <- paste0(result,"_t",time.limit)
    result <- paste0(result,"_i",iter.nbr)
    result <- paste0(result,"_n", oneOptNeig)
    
    return(result)
}


#############################################################################################
# Returns the inline command for the Greedy Randomized Adaptive Search Procedure (Grasp) partioning 
# method. See the algorithm documentation for more details.
#
# algo.name: short code associated to the algorithm.
# input.folder: relative path to the folder containing targeted graph file.
# out.folder: relative path to the folder in which the produced files will be placed.
# time.limit: maximum duration of the processing.
# iter.nbr: maximum number of iterations of the processing.
#
# returns: the command allowing to invoke the program externally.
#############################################################################################
get.grasp.command <- function(algo.name, input.folder, out.folder, graph.name)
{	
    # break down the specified code (short name)
    tmp = break.down.algo.code(algo.name)
    
    params <- c()
    k=NA
    next.param.indx=NA
    # =======================================================================================
    # rcc flag
    # example1: GRASP-RCC_kFrom=(ILS-CC_l1_a1_g0_p3_t3_i10):k+1_l1_a1_g0_p3_t3_i10
    # example2: GRASP-RCC_k=2_l1_a1_g0_p3_t3_i10
    algo.name <- strsplit(x=tmp[1], split="-", fixed=TRUE)[[1]][1]
    rcc.flag <- strsplit(x=tmp[1], split="-", fixed=TRUE)[[1]][2]
    if(rcc.flag=="RCC"){
        params["rcc"] <- 1
        k = get.k.from.algo.name(tmp[2], out.folder)
        next.param.indx = 3
    }
    else{
        params["rcc"] <- 0
        next.param.indx = 2
    }
    # =======================================================================================
    
    for(s in tmp[next.param.indx:length(tmp)])
    {	params <- c(params,substr(s,2,nchar(s)))
    names(params)[length(params)] <- substr(s,1,1)
    }
    
    
    # init
    input.file <- file.path(input.folder, graph.name)
    output.folder <- file.path(out.folder)
    command.folder <- file.path(ILS.GRASP.LIB.FOLDER)
    result <- file.path(command.folder, "graspcc")
    # result <- paste0("mpirun -n 1 ",result)
    
    
    # build command
    result <- paste0(result, " --alpha ",params["a"])
    result <- paste0(result, " --iterations ",params["i"])
    result <- paste0(result, " --neighborhood_size ",params["l"])
    result <- paste0(result, " --rcc ",params["rcc"])
    if(rcc.flag=="RCC")
        result <- paste0(result, " --k ",k)
    
    result <- paste0(result, " --time-limit ",params["t"])
    result <- paste0(result, " --input-file \"",input.file,"\"")
    result <- paste0(result, " --output-folder \"",output.folder,"\"")
    result <- paste0(result, " --gain-function-type ",params["g"])
    result <- paste0(result, " --strategy ","GRASP")
    result <- paste0(result, " --firstImprovementOnOneNeig ",params["n"])
    
    return(result)
}




#############################################################################################
# 
#############################################################################################
get.ExCC.code <- function(enum.all)
{
	result <- COR.CLU.ExCC
	if(enum.all)
		result <- COR.CLU.ExCC.ENUM.ALL
	return(result)
}


#############################################################################################
# 
#############################################################################################
get.ExCC.command <- function(algo.name, input.folder, out.folder, graph.name)
{
	is.cp = "true" # in any case, use cutting plane approach
	is.enum.all = "false"
	# tilim = 3600 # 1 hour
	tilim = -1 #  no time limit
	if(algo.name == COR.CLU.ExCC.ENUM.ALL){
		is.enum.all = "true"
		tilim = -1
	}
		
	print(graph.name)
	input.file = paste("'", input.folder, "/", graph.name, "'", sep="")
	
	# An example:
	# java -Djava.library.path=/users/narinik/Cplex/ILOG/CPLEX_Studio128/cplex/bin/x86-64_linux/
	# -DinFile=data/test.G -DoutDir=. -Dcp=false -DenumAll=false -Dtilim=-1 -DlazyInBB=false
	# -DusercutInBB=false -DMaxTimeForRelaxationImprovement=-1 -jar exe/cplex-partition.jar
	
	cmd = 
		paste(
				"java",		
				paste("-Djava.library.path=", CPLEX.BIN.PATH, sep=""),
				paste0("-DinFile=", input.file),
				paste0("-DoutDir=", out.folder),
				paste0("-Dcp=",is.cp),
				paste0("-DenumAll=",is.enum.all),
				paste0("-Dtilim=",tilim),
				paste0("-DMaxTimeForRelaxationImprovement=","-1"), # no specific time limit, use the default one
				"-DlazyInBB=false",
				"-DuserCutInBB=false",
				"-jar",
				ExCC.JAR.PATH,
				sep=" "
		)

	print(cmd)

	return(cmd)
}







#############################################################################################
#
#############################################################################################
prepare.algo.output.filename = function(part.folder, algo.name, g.name){
    
    if(startsWith(algo.name,COR.CLU.ILS) || startsWith(algo.name,COR.CLU.GRASP)
       || startsWith(algo.name,COR.CLU.VNS) || startsWith(algo.name,COR.CLU.VOTE.BOEM))
    {
        # identify the latest folder
        temp.folder <- file.path(part.folder, g.name)
        details <- file.info(list.dirs(temp.folder))
        if(nrow(details)==0)
            stop("load.mestrado.partition: No subfolder was found in folder \"", temp.folder ,"\", cannot load the partition.")
        details <- details[with(details, order(as.POSIXct(mtime))), ]
        temp.folders <- rownames(details)
        folder.name <- file.path(temp.folders[length(temp.folders)])
        
        # set up the file name
        tmp <- strsplit(x=algo.name, split="_", fixed=TRUE)[[1]]
        rcc.flag <- strsplit(x=tmp[1], split="-", fixed=TRUE)[[1]][2]
        f.name=NA
        if(rcc.flag=="RCC")
            f.name <- "rcc-result.txt"
        else
            f.name <- "cc-result.txt"
        file.name <- file.path(folder.name, f.name)
        
        # retain only the partiton info file by changing its name
        id=0
        file.copy(from=file.name, to=file.path(part.folder, paste0(ALGO.RESULT.FILE.PREFIX,id,".txt")))
        
        # possibly remove the original algorithm files
        unlink(x=temp.folder, recursive=TRUE)
    }
    
    else if(algo.name == COR.CLU.ExCC)
    {
        ExCC.output.file <- file.path(part.folder, "ExCC-result.txt")
        id=0
        file.rename(from=ExCC.output.file, to=file.path(part.folder, paste0(ALGO.RESULT.FILE.PREFIX,id,".txt")))
    }
    else if(algo.name == COR.CLU.ExCC.ENUM.ALL){
        # do nothing
    }

    
}




#############################################################################################
# Returns the full name based on the normalized (short) name. Note that for parameterized 
# algorithms, this will just return a clean version of the short name, since it contains 
# the parameter values.
#
# algo.names: short names of the considered algorithms.
#
# returns: the corresponding full names, to be used in plots for instance.
#############################################################################################
get.algo.names <- function(algo.names)
{	result <- c()
	
	for(algo.name in algo.names)
	{
		# parameters
		result <- c(result, gsub(pattern="_", replacement=" ", x=algo.name, fixed=TRUE))
	}
	
	return(result)
}



#############################################################################################
# Returns the inline command for the specified algorithm. The "..." parameters are fetched
# to the algorithm-specific function.
#
# algo.name: short code associated to the algorithm, including its parameter values.
#
# returns: the command allowing to invoke the program externally.
#############################################################################################
get.algo.commands <- function(algo.names, ...)
{	result <- c()

    # substring(x, 1, nchar(prefix)) == prefix
    
    for(algo.name in algo.names)
    {	
        if(startsWith(algo.name,COR.CLU.ILS) || startsWith(algo.name,COR.CLU.VNS))
            result <- c(result, get.ils.command(algo.name, ...))
        else if(startsWith(algo.name,COR.CLU.GRASP) || startsWith(algo.name,COR.CLU.VOTE.BOEM))
            result <- c(result, get.grasp.command(algo.name, ...))
        else if(startsWith(algo.name,COR.CLU.ExCC))
            result <- c(result, get.ExCC.command(algo.name, ...))
    }
    
    return(result)
}
