
# ===============================================================
# Exact Approach: ExCC
# ===============================================================

NB.THREAD = 4

COR.CLU.ExCC <- "ExCC"
COR.CLU.ExCC.ENUM.ALL <- "OneTreeCC"
ExCC.LIB.FOLDER = file.path(LIB.FOLDER,COR.CLU.ExCC)
ExCC.JAR.PATH = paste(ExCC.LIB.FOLDER,paste0(COR.CLU.ExCC,".jar"),sep="/")
#CPLEX.BIN.PATH = "/users/narinik/Cplex/ILOG/CPLEX_Studio128/cplex/bin/x86-64_linux/"
CPLEX.BIN.PATH = "/opt/ibm/ILOG/CPLEX_Studio201/cplex/bin/x86-64_linux/"

ExCCAll.LIB.FOLDER = file.path(LIB.FOLDER,COR.CLU.ExCC)
ExCCAll.JAR.PATH = paste(ExCC.LIB.FOLDER,paste0(COR.CLU.ExCC,".jar"),sep="/")
ExCCAll.MAX.NB.SOLS = 1000
ExCCAll.MAX.TIME.LIMIT = 3600*12 # 12 hours

# ---------------------------------------------------------------
# ---------------------------------------------------------------

ENUMCC = "EnumCC"
ENUMCC.LIB.FOLDER = file.path(LIB.FOLDER,ENUMCC)
ENUMCC.FOLDER = file.path(LIB.FOLDER,ENUMCC)
ENUMCC.JAR.PATH = paste(ENUMCC.LIB.FOLDER,paste0(ENUMCC,".jar"),sep="/")

ENUMCC.MAX.NB.SOLS = 1000 # we know that our method is not efficient for very large number of solutions
ENUMCC.MAX.TIME.LIMIT = 3600*12 # 12 hours





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
# It returns the description code of the method ExCC or OneTreeCC.
#############################################################################################
get.ExCC.code <- function(enum.all)
{
	result <- COR.CLU.ExCC
	if(enum.all)
		result <- COR.CLU.ExCC.ENUM.ALL
	return(result)
}


#############################################################################################
# It returns the command of the method ExCC. Note that OneTreeCC is also called through this function.
#############################################################################################
get.ExCC.command <- function(algo.name, input.folder, out.folder, graph.name)
{
	is.cp = "true" # in any case, use cutting plane approach
	is.enum.all = "false"
	# tilim = 3600 # 1 hour
	if(algo.name == COR.CLU.ExCC.ENUM.ALL){
		is.enum.all = "true"
	}
	
	print(graph.name)
	input.file = paste("'", input.folder, "/", graph.name, "'", sep="")
	input.file.for.g = file.path(input.folder, graph.name)
	g = read.graph.ils(input.file.for.g)
	
	initSolutionFilePath = file.path(out.folder,"..","..",COR.CLU.ExCC,SIGNED.UNWEIGHTED.FILE,"membership0.txt")
	if(!file.exists(initSolutionFilePath))
		initSolutionFilePath="''"
	
	cmd = "NONE"
	if(is.enum.all == "false"){
		# An example:
		# java -Djava.library.path=/users/narinik/Cplex/ILOG/CPLEX_Studio128/cplex/bin/x86-64_linux/
		# -DinFile="in/""$name" -DoutDir="out/""$modifiedName" -DenumAll=false -Dcp=true -DMaxTimeForRelaxationImprovement=20
		# -DuserCutInBB=false -DinitSolutionFilePath="$initSolutionFilePath" -DLPFilePath="$LPFilePath"
		# -DonlyFractionalSolution=false -DfractionalSolutionGapPropValue=-1.0 -DnbThread=2 -Dverbose=true -Dtilim=200 -jar exe/ExCC.jar
		
		# -------------------------------------------------------------------------
		# TODO handle this in a better way, for instance, use small values for sparse networks
		# TODO write a function called 'estimate.max.time.for.relaxation.improvement(..)'
		maxTimeForRelaxationImprovement = "600"
		if(vcount(g)>39 && vcount(g)<=50)
			maxTimeForRelaxationImprovement = "4500" # 1h15m
		else if(vcount(g)>50)
			maxTimeForRelaxationImprovement = "10000" # 166 mins 
		# -------------------------------------------------------------------------
		
		cmd = 
			paste(
				"java",		
				paste("-Djava.library.path=", CPLEX.BIN.PATH, sep=""),
				paste0("-DinFile=", input.file),
				paste0("-DoutDir=", out.folder),
				paste0("-Dcp=",is.cp),
				paste0("-DenumAll=",is.enum.all),
				paste0("-Dtilim=",-1),
				paste0("-DtilimForEnumAll=",-1),
				paste0("-DsolLim=",1),
				paste0("-DMaxTimeForRelaxationImprovement=",maxTimeForRelaxationImprovement),
				"-DlazyInBB=false",
				"-DuserCutInBB=false",
				paste0("-DinitSolutionFilePath=",initSolutionFilePath),
				"-Dverbose=true",
				paste0("-DnbThread=", NB.THREAD),
				"-DLPFilePath=''",
				"-DonlyFractionalSolution=false",
				"-DfractionalSolutionGapPropValue=-1",
				"-jar",
				ExCC.JAR.PATH,
				sep=" "
			)
	} else { # if(is.enum.all == "true") >> OneTreeCC
		# An example:
		# java -DinFile="in/""$name" -DoutDir="out/""$modifiedName" -DenumAll=true
		# -Dcp=false -DinitSolutionFilePath="$initSolutionFilePath" -DLPFilePath="$LPFilePath"
		# -DnbThread=2 -Dverbose=true -Dtilim=-1 -DtilimForEnumAll=60 -DsolLim=100 -jar exe/ExCC.jar
		
		LP.filepath = file.path(out.folder,"..","..",COR.CLU.ExCC,SIGNED.UNWEIGHTED.FILE,"strengthedModelAfterRootRelaxation.lp")
		if(file.exists(LP.filepath)) {
			cmd = 
				paste(
					"java",		
					paste("-Djava.library.path=", CPLEX.BIN.PATH, sep=""),
					paste0("-DinFile=", input.file),
					paste0("-DoutDir=", out.folder),
					paste0("-Dcp=","false"),
					paste0("-DenumAll=","true"),
					paste0("-Dtilim=",-1),
					paste0("-DtilimForEnumAll=",ExCCAll.MAX.TIME.LIMIT),
					paste0("-DsolLim=",ExCCAll.MAX.NB.SOLS),
					paste0("-DMaxTimeForRelaxationImprovement=","-1"), # no specific time limit, use the default one
					"-DlazyInBB=false",
					"-DuserCutInBB=false",
					paste0("-DinitSolutionFilePath=",initSolutionFilePath),
					"-Dverbose=true",
					paste0("-DnbThread=", NB.THREAD),
					paste0("-DLPFilePath='",LP.filepath,"'"),
					"-DonlyFractionalSolution=false",
					"-DfractionalSolutionGapPropValue=-1",
					"-jar",
					ExCC.JAR.PATH,
					sep=" "
				)
		}
	}

	print(cmd)
	return(cmd)
}







#############################################################################################
# It returns the description code of the method EnumCC.
#############################################################################################
get.EnumCC.code <- function(maxNbEdit)
{
	result <- paste0(ENUMCC,"-maxNbEdit",maxNbEdit)
	return(result)
}


#############################################################################################
# It returns the command of the method EnumCC.
#############################################################################################
get.EnumCC.command <- function(algo.name, input.folder, out.folder, graph.name)
{
	print(algo.name)
	base.algo.name <- strsplit(x=algo.name, split="-", fixed=TRUE)[[1]][1]
	params.str <- gsub(paste0(base.algo.name,"-"),"",algo.name)
	print(params.str)
	
	print(graph.name)
	input.file = paste("'", input.folder, "/", graph.name, "'", sep="")
	
	maxNbEdit = as.integer(gsub("maxNbEdit","",params.str))
	
	cmd = "NONE"
	LP.filepath = file.path(out.folder,"..","..",COR.CLU.ExCC,SIGNED.UNWEIGHTED.FILE,"strengthedModelAfterRootRelaxation.lp")
	
	if(file.exists(LP.filepath)) {
		cmd = 
			paste(
				"java",		
				paste("-Djava.library.path=", CPLEX.BIN.PATH, sep=""),
				paste0("-DinFile=", input.file),
				paste0("-DoutDir=", out.folder),
				paste0("-DLPFilePath=", LP.filepath),
				paste0("-DinitMembershipFilePath=", file.path(out.folder,"..","..",COR.CLU.ExCC,SIGNED.UNWEIGHTED.FILE,"membership0.txt")),
				paste0("-DnbThread=",NB.THREAD),
				paste0("-DmaxNbEdit=",maxNbEdit),
				paste0("-DsolLim=",ENUMCC.MAX.NB.SOLS),
				paste0("-Dtilim=",ENUMCC.MAX.TIME.LIMIT),
				paste0("-DJAR_filepath_RNSCC=",paste(ENUMCC.LIB.FOLDER,paste0("RNSCC.jar"),sep="/")),
				"-jar",
				ENUMCC.JAR.PATH,
				sep=" "
			)
	}
	
	print(cmd)
	return(cmd)
}





#############################################################################################
# It prepares the output of a partitioning algorithm. It can be:
# - renaming the output result files
# - creating an output file called 'allResults.txt', which indicates the absolute paths of the solutions obtained.
#############################################################################################
prepare.algo.output = function(part.folder, algo.name, g.name){

    if(algo.name == COR.CLU.ExCC)
    {
        ExCC.output.file <- file.path(part.folder, "ExCC-result.txt")
        id=0
        file.rename(from=ExCC.output.file, to=file.path(part.folder, paste0(ALGO.RESULT.FILE.PREFIX,id,".txt")))
		########################
		sol.paths = list.files(path = part.folder,
				pattern = paste0("^", ALGO.RESULT.FILE.PREFIX, ".*\\.txt$"), full.names = TRUE, recursive = FALSE)
		#write.table(x=sol.paths, file=file.path(part.folder,"allResults.txt"), row.names=F, col.names=F, quote=F)
		sol.paths.ordered = paste0(part.folder,"/",MBRSHP.FILE.PREFIX,seq(0,length(sol.paths)-1),".txt")
		write.table(x=sol.paths.ordered, file=file.path(part.folder,"allResults.txt"), row.names=F, col.names=F, quote=F)
		
	}
    else if(algo.name == COR.CLU.ExCC.ENUM.ALL){
		sol.paths = list.files(path = part.folder,
				pattern = paste0("^", ALGO.RESULT.FILE.PREFIX, ".*\\.txt$"), full.names = TRUE, recursive = FALSE)
		#write.table(x=sol.paths, file=file.path(part.folder,"allResults.txt"), row.names=F, col.names=F, quote=F)
		sol.paths.ordered = paste0(part.folder,"/",MBRSHP.FILE.PREFIX,seq(0,length(sol.paths)-1),".txt")
		write.table(x=sol.paths.ordered, file=file.path(part.folder,"allResults.txt"), row.names=F, col.names=F, quote=F)
	}
	else if(startsWith(algo.name,ENUMCC)){
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
        if(startsWith(algo.name,COR.CLU.ExCC))
            result <- c(result, get.ExCC.command(algo.name, ...))
		else if(startsWith(algo.name,COR.CLU.ExCC.ENUM.ALL))
			result <- c(result, get.ExCC.command(algo.name, ...))
        else if(startsWith(algo.name,ENUMCC))
            result <- c(result, get.EnumCC.command(algo.name, ...))
	}
	
	return(result)
}
