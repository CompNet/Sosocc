
#################################################################
# Exact Approach of Correlation Clustering (CC) problem
#################################################################

COR.CLU.ExCC <- "ExCC" # only for 1 optimal solution
COR.CLU.ExCC.ENUM.ALL <- "ExCC-all" # for all optimal solutions
ExCC.LIB.FOLDER = file.path(LIB.FOLDER,COR.CLU.ExCC)
ExCC.JAR.PATH = paste(ExCC.LIB.FOLDER,"cplex-partition.jar",sep="/") 
#CPLEX.BIN.PATH = "/users/narinik/Cplex/ILOG/CPLEX_Studio128/cplex/bin/x86-64_linux/" # in gaia cluster - CERI
CPLEX.BIN.PATH = "/opt/ibm/ILOG/CPLEX_Studio128/cplex/bin/x86-64_linux/" # in my personel computer





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
	tilim = 3600 # 1 hour
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

	if(algo.name == COR.CLU.ExCC)
	{
	    # since the ExCC program outputs the result by writing into a file whose the name has a suffix 'O'. 
	    # This is because the program can also handle enumerate all optimal solutions. 
	    # And in this case, it outputs like that: "sol0.txt", "sol1.txt", and so on
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
		if(startsWith(algo.name,COR.CLU.ExCC))
			result <- c(result, get.ExCC.command(algo.name, ...))
	}
	
	return(result)
}
