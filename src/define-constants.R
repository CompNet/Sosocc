#############################################################################################
# This script allows defining all the file constants used by the other scripts and functions.
# This include both file names and folder names.

#############################################################################################


ALL = "All"



#############################################################################################
# Folders
#############################################################################################
# main folder
MAIN.FOLDER <- "."
#MAIN.FOLDER <- "/home/vlabatut/eclipse/remworkspaces/Networks/NetVotes"
# external libraries folder
LIB.FOLDER <- file.path(MAIN.FOLDER,"lib")
# general input folder
IN.FOLDER <- file.path(MAIN.FOLDER,"in")
#IN.FOLDER <- file.path(MAIN.FOLDER,"in-k=5-n20_50")
#IN.FOLDER <- file.path(MAIN.FOLDER,"in-k=2-n20_50")
	# folder containing signed random networks
 	RANDOM.NETWORKS.FOLDER <- file.path(IN.FOLDER,"random-networks")
	REAL.NETWORKS.FOLDER <- file.path(IN.FOLDER,"real-networks")
# general ouput folder
OUT.FOLDER <- file.path(MAIN.FOLDER,"out")   # ====> TODO "out-random-networks", "out-real-networks"
#OUT.FOLDER <- file.path(MAIN.FOLDER,"out-k=5-n20_50")
#OUT.FOLDER <- file.path(MAIN.FOLDER,"out-k=2-n20_50")

	# folder containing the partitions corresponding to the document-wise networks
#	PARTITIONS.FOLDER <- file.path(OUT.FOLDER,"partitions")
	PARTITIONS.FOLDER <- file.path(OUT.FOLDER,"partitions")
		OPPOSITE.NETWORKS.FOLDER = file.path(PARTITIONS.FOLDER,"partitions")
		
	# folder containing the evluation of partitions
	EVALUATION.FOLDER <- file.path(OUT.FOLDER,"evaluate-partitions")
	PLOT.FOLDER <- file.path(OUT.FOLDER,"layout-plots")
	CLU.ANALYSIS.FOLDER <- file.path(OUT.FOLDER,"cluster-analysis")
	CLUSTER.CHARACTERIZATION.FOLDER <- file.path(OUT.FOLDER,"cluster-characterization")
	CLU.ANALYSIS.DIFFERENCE.FOLDER <- file.path(OUT.FOLDER,"cluster-analysis-difference")
	
	# NETWORKS.STAT.FOLDER = file.path(OUT.FOLDER,"networks-stats")
	OUTPUT.CSV.FOLDER = file.path(OUT.FOLDER,"output-csv-data")
	
EXACT.SUMMARY.FOLDER.NAME = "exact-summary"
SUMMARY.FOLDER.NAME = "summary"
PLOTS.FOLDER.NAME = "plots"

TRANSITION.GRAPH.PLOT.FOLDER = file.path(OUT.FOLDER,"transition-graph")

ALL.PARAM.NAME = "All"

NOT.AVAILABLE.FOLDER.NAME = "NotAvailable"

# for getting plot folder name
ALL.COR.CLU.EXACT.ALGOS.PARAM.NAME = "All-Exact"

# cluster analysis
REPRESENTATIVE.FILE.PREFIX = "representative"
REPRESENTATIVE.FOLDER.NAME = "representative-partitions"
CORE.PART.FOLDER.NAME = "core-part"

# cluster analysis difference
BEST.K.FOLDER.NAME = "best-cluster"
ALL.K.FOLDER.NAME = "all-clusters"



#############################################################################################
# Some parameters used for evaluating the partition results
#############################################################################################
MAX.NB.SOLS.FOR.PROCESSING = 5000
EVAL.NB.ENUM.FILENAME = "nb-enum"
NB.ENUM.COL.NAME = "nb enumeration"
EVAL.EXEC.TIME.FILENAME = "exec-time"
EVAL.DELAY.EXEC.TIME.FILENAME = "delay-exec-time"
EVAL.ENUM.EXEC.TIME.FILENAME = "enum-exec-time"
DELAY.EXEC.TIME.COL.NAME = "delay exec time (s)"
ENUM.STEP.NO.COL.NAME="Enum Step No"


#############################################################################################
# Clustering Distance Measure Names
#############################################################################################
#NMI = "Nmi"
#RAND = "Rand"
#ADJUSTED.RAND = "AdjustedRand"
#ADJUSTED.RAND2 = "AdjustedRand2"
#HARMONIC.MEAN.PURITY = "HarmonicMeanPurity"
NMI = "Nmi"
RI = "RI" # Rand index
HA.ARI = "HA's ARI" # Hubert and Ara-bie's adjusted Rand index
MA.ARI = "MA's ARI" # Morey and Agresti’s RI
JACCARD = "Jaccard"
FM = "FM" # Fowlkes and Mallows’s index
HARMONIC.MEAN.PURITY = "HarmonicMeanPurity" # Harmonic mean purity (F-measure)
MI = "MI"
VI = "VI" # Variation of Information
NVI = "NVI" # Normalized Variation of Information => see Vinh et al. 2010 (not the seminal paper of VI)
NVI.by.k = "NVI by k" # Normalized version of Variation of Information (nomalized by 2*log(k))
NVI.by.n = "NVI by n" # Normalized version of Variation of Information (nomalized by log(n))
EDIT = "Edit"
EDIT.NORM = "EditNorm"

##################################################################################
# network/graph parameter names
##################################################################################

GRAPH.SIZE = "Graph size"
PROP.MISPL = "Prop mispl"
PROP.NEG = "Prop neg"
NETWORK.NO = "Network no"
ALGO.NAME = "Algo name"
GRAPH.DESC.NAME = "Graph desc"

# for I(P) interval ay y-axis
DETECTED.IMB.LOWER.BOUND = 0.0 # ----> for "Detected Imbalance percentage"
DETECTED.IMB.UPPER.BOUND = 1.0 # ----> for "Detected Imbalance percentage"
DETECTED.IMB.PROP.LABEL.NAME = "Detected Imbalance proportion"
DETECTED.IMB.PROP.PARAM.NAME = "detectedImbalance"



##################################################################################
# PLOT TYPES
##################################################################################

DENS = "density"
HIST = "hist"
MULTI.BAR = "multi-bar"
LINE = "line"
MULTI.LINE = "multi-line"
VIO = "vio"
VIO2 = "vio2"
GROUPED.VIO = "grouped-vio"
SCATTER = "scatter"
MULTI.BAR.WITH.PROPORTION = "multi-bar-prop"

##################################################################################
# PLOT TYPES
##################################################################################

ALLUVIAL.DIAGRAM.PLOT.TYPE = "alluvial-diagram"
PLOT.PREFIX = "plot"
PLOT.AS.PDF = "PDF"
PLOT.AS.JPEG = "JPEG"
PLOT.AS.PNG = "PNG"
JUST.PLOT = "NA"


##################################################################################
# LAYOUT TYPES
##################################################################################


FRUCHTERMAN.REINGOLD = "fruchterman.reingold"
KAMADA.KAWAI = "kamada.kawai"
KUNEGIS = "kunegis"
RW.KUNEGIS = "rw_kunegis"
SN.KUNEGIS = "sn_kunegis"
SN.ZHENG = "sn_zheng"
BN.ZHENG = "bn_zheng"
#SIGNED.GRAPH.LAYOUTS = c(KUNEGIS, RW.KUNEGIS, SN.KUNEGIS, SN.ZHENG, BN.ZHENG)
SIGNED.GRAPH.LAYOUTS = c(KUNEGIS, BN.ZHENG)
#UNSIGNED.GRAPH.LAYOUTS = c(FRUCHTERMAN.REINGOLD, KAMADA.KAWAI)
UNSIGNED.GRAPH.LAYOUTS = c()

##################################################################################
# RESULT TYPES - key words
##################################################################################

EXEC.TIME = "Exec-time"
NB.CLUSTER = "Nb-cluster"
NB.SOLUTION = "Nb-solution"
NB.SOLUTION.BY.SOL.CLASSES = "Nb-solution-by-sol-classes"
NB.OPT.SOLUTION = "Nb-optimal-solution"
IMB.COUNT = "Imbalance-count"
IMB.PERC = "Imbalance-percentage"
IMB.OPT = "Imbalance-optimality"
DIST.SCORE = "Distance-score"
ALL.DIST.SCORE = "All-distance-scores"

MEAN.DIST.SCORE = "Mean-of-distance-scores"
NB.COMP.TRANSITION.GRAPH = "Nb-component-transition-graph"
BEST.K.FOR.KMEDOIDS = "Best-k-for-kmedoids"
SINGLE.CLU.PROP.FOR.KMEDOIDS = "single-cluster-proportion-for-kmedoids"

CLASS.CORE.PART.SIZE = "class-core-part-size"
CLASS.CORE.PART.SIZE.BY.SOL.CLASSES = "class-core-part-size-by-sol-classes"

##################################################################################
# LABEL NAMES
##################################################################################

EXEC.TIME.LABEL.NAME = "Execution time (sec)"
NB.CLUSTER.LABEL.NAME = "Number of cluster"
NB.SOLUTION.LABEL.NAME = "Number of solution"
IMB.COUNT.LABEL.NAME = "Imbalance count"
IMB.PERC.LABEL.NAME = "Imbalance percentage (%)"
IMB.OPT.LABEL.NAME = "Imbalance optimality (%)"
DIST.SCORE.LABEL.NAME = "Distance score"
#RESULT.CORRESPONDANCE.LABEL.NAME = "Result correspondance" # with scatterplot on x and y axes
MEAN.DIST.SCORE.LABEL.NAME = "Mean of distance scores"




##################################################################################
# INTERVALS
##################################################################################

IMB.OPT.INTERVAL = c(0,100)


#############################################################################################
# Files
#############################################################################################
# script filename
RECORD.MEM.INFO.SCRIPTNAME = "record-mem.sh"

USED.MEM.VAL.FILENAME = "used-ram-memory.txt" #warning: be consistent with 'record-mem.sh' for the filename

# graph files
GRAPH.FILE.LAYOUT.PREFIX = "layout"
SIGNED.UNWEIGHTED.FILE = "signed-unweighted" # undirected as well
SIGNED.WEIGHTED.FILE = "signed-weighted" # undirected as well

# partitioning
ALGO.RESULT.FILE.PREFIX = "sol"
MBRSHP.FILE.PREFIX = "membership"
#PLOT.MBRSHP.FILE.PREFIX = paste0("plot-",MBRSHP.FILE.PREFIX)
EXEC.TIME.FILENAME = "exec-time.txt"

# k-medoids cluster analysis 
CLU.ANALYSIS.FILENAME = "cluster-analysis"


# Partition evaluation result files
EVAL.EXEC.TIME.FILENAME = "exec-time"
#EVAL.MEAN.EXEC.TIME.FILENAME = "mean-exec-time"
EVAL.NB.SOL.FILENAME = "nb-solution"
EVAL.NB.OPT.SOL.FILENAME = "nb-opt-solution"
EVAL.DIST.MATRIX.FILE.PREFIX = "dist-matrix"
EVAL.CLU.INFO.TABLE.FILENAME = "clusters"
EVAL.IMB.INFO.TABLE.FILENAME = "imbalance" # general name => something whihc covers imb count/perc etc.
EVAL.IMB.COUNT.TABLE.FILENAME = "imbalance-count"
EVAL.IMB.PERC.TABLE.FILENAME = "imbalance-perc"
EVAL.IMB.OPT.TABLE.FILENAME = "imbalance-optimality"
EVAL.MEAN.SD.DIST.SCORES.FILE.PREFIX = "mean-sd-dist-scores"
EVAL.ALL.DIST.SCORES.FILE.PREFIX = "dist-scores"
BIG.MEMORY.MATRIX.DESCRIPTOR.FILENAME = "dist-mtrx.desc"
BIG.MEMORY.MATRIX.BACKING.FILENAME = "dist-mtrx.bin"

EVAL.EDIT.DIST.COMPS = "Edit-Dist-Components"
EVAL.EDIT.DIST.NB.COMP = "Edit-Dist-NbComponents"
EVAL.EDIT1.DIST.NB.COMP.PROP = "Edit1-Dist-NbComponents-proportion"

EVAL.BEST.K.FOR.KMEDOIDS = "Best-k-for-kmedoids"
#EVAL.BEST.K.FOR.KMEDOIDS.PROP = "Best-k-for-kmedoids-proportion"
EVAL.SINGLE.CLU.PROP.FOR.KMEDOIDS = "Single-cluster-prop-for-kmedoids"



# cluster analysis
SILH.SCORE.FILENAME = "silhouette-scores"
SF.SCORE.FILENAME = "SF-scores"
CLUSTER.ID.FOLDER.PREFIX = "clu"

EDIT.DIST.MATRIX.FILENAME = "edit-dist-mtrx"



#############################################################################################
# Transition graph
#############################################################################################
EDIT1.EDGE.WEIGHT = 4
EDIT1.EDGE.COLOR = "green"
EDIT2.EDGE.WEIGHT = 3
EDIT2.EDGE.COLOR = "blue"
EDIT3.EDGE.WEIGHT = 1
EDIT3.EDGE.COLOR = "orange"
EDIT4.EDGE.WEIGHT = 0.75
EDIT4.EDGE.COLOR = "red"
BETW.COMP.EDGE.WEIGHT = 0.5




#############################################################################################
# CSV Column names
#############################################################################################

# colnames used for each network
EXEC.TIME.COL.NAME = "exec time (s)"
NB.SOL.COL.NAME = "nb solution"
NB.OPT.SOL.COL.NAME = "nb optimal solution"
NB.CLU.COL.NAME = "nb cluster"
CLUSTER.COL.NAME.PREFIX = "clu"
IMB.COUNT.COL.NAME = "imbalance count"
IMB.PERC.COL.NAME = "imbalance percentage"
IMB.OPTIMALITY.PERC.COL.NAME = "imbalance optimality (%)"
MEAN.COL.NAME = "mean"
SD.COL.NAME = "sd"
ID.COL.NAME = "id"
FREQ.COL.NAME = "frequency"
ASSOCIATION.COL.NAME = "association"


ALL.COL.NAME = "All"

# colnames used for each prop.neg (over networks)
MEAN.PREFIX = "Mean of"
SD.PREFIX = "Std of"
MEAN.NB.SOL.COL.NAME = paste(MEAN.PREFIX, NB.SOL.COL.NAME)
SD.NB.SOL.COL.NAME = paste(SD.PREFIX, NB.SOL.COL.NAME)
MEAN.NB.CLU.COL.NAME = paste(MEAN.PREFIX, NB.CLU.COL.NAME)
SD.NB.CLU.COL.NAME = paste(SD.PREFIX, NB.CLU.COL.NAME)
MEAN.IMB.COUNT.COL.NAME = paste(MEAN.PREFIX, IMB.COUNT.COL.NAME)
SD.IMB.COUNT.COL.NAME = paste(SD.PREFIX, IMB.COUNT.COL.NAME)
MEAN.IMB.PERC.COL.NAME = paste(MEAN.PREFIX, IMB.PERC.COL.NAME)
SD.IMB.PERC.COL.NAME = paste(SD.PREFIX, IMB.PERC.COL.NAME)
MEAN.IMB.OPT.COL.NAME = paste(MEAN.PREFIX, IMB.OPTIMALITY.PERC.COL.NAME)
SD.IMB.OPT.COL.NAME = paste(SD.PREFIX, IMB.OPTIMALITY.PERC.COL.NAME)
# for dist scores  for each prop.neg (over networks)
MEAN.MEAN.COL.NAME = paste(MEAN.PREFIX, MEAN.COL.NAME)
SD.MEAN.COL.NAME = paste(SD.PREFIX, MEAN.COL.NAME)
MEAN.SD.COL.NAME = paste(MEAN.PREFIX, SD.COL.NAME)
SD.SD.COL.NAME = paste(SD.PREFIX, SD.COL.NAME)


EDIT1.COMP.MEM.COL.NAME = "Edit1-Components"
EDIT12.COMP.MEM.COL.NAME = "Edit12-Components"
EDIT123.COMP.MEM.COL.NAME = "Edit123-Components"
EDIT1234.COMP.MEM.COL.NAME = "Edit1234-Components"

BEST.K.FOR.SILH.COL.NAME = "Best k for Silhouette"
BEST.SILH.COL.NAME = "Best Silhouette"
BEST.MEM.FOR.SILH.COL.NAME = "Best clustering for Silhouette"

PROP.SINGLE.CLU.FOR.SILH.COL.NAME = "Proportion of k=1 for Silhouette"


#############################################################################################
# kmedoids-specific config
#############################################################################################
KMEDOIDS.PLOT.THRESHOLD = 5 # 5 means that we plot the graph with partition info only for 'k=1', 'k=2', 'k=3', 'k=4' and 'k=5'
SILH.THRESH = 0.50 # a threshold used to get a reasonable cluster structure

#############################################################################################
# Core part analysis
#############################################################################################
CORE.PART.THRESHOLDS = c(1.00)
# CORE.PART.THRESHOLDS = c(0.80, 0.90, 1.00)

PROP.CLASS.CORE.PART.SIZE.COL.NAME = "Proportion of class core part size in core part analysis"
EVAL.CLASS.CORE.PART.SIZE.FILE.NAME = "class-core-part-size"
#PROP.GENERAL.CORE.PART.SIZE.COL.NAME = "Proportion of general core part size in core part analysis"
#EVAL.GENERAL.CORE.PART.SIZE.FILE.NAME = "general-core-part-size"
#PROP.CORE.PART.ALL.UNSTABLE.SIZE.COL.NAME = "Proportion of unstable node size in core part analysis"
#EVAL.CORE.PART.ALL.UNSTABLE.SIZE.FILE.NAME = "core-part-all-unstable-size"
#PROP.CORE.PART.SIM.BETW.SOL.CLASSES.COL.NAME = "Proportion of similarity between solution classes in core part analysis"
#EVAL.CORE.PART.SIM.BETW.SOL.CLASSES.FILE.NAME = "core-part-sim-betw-sol-classes"

#############################################################################################
# Some Inner parameters in rder to create distance matrix in parallel mode
#############################################################################################
# parameters used for creating distance matrix in parallel mode
PAR.MODE.THRESH.NB.MEM.FILE = 600
PAR.MODE.NB.CORE.DEFAULT = 4
PAR.MODE.CHUNK.SIZE.DEFAULT = 200






#############################################################################################
# MAPPING
#############################################################################################
csv.result.files = list()

csv.result.files[[NB.SOLUTION]] = EVAL.NB.SOL.FILENAME
csv.result.files[[NB.SOLUTION.BY.SOL.CLASSES]] = EVAL.NB.SOL.FILENAME
csv.result.files[[NB.OPT.SOLUTION]] = EVAL.NB.SOL.FILENAME
csv.result.files[[EXEC.TIME]] = EVAL.EXEC.TIME.FILENAME
csv.result.files[[NB.CLUSTER]] = EVAL.CLU.INFO.TABLE.FILENAME
csv.result.files[[NB.COMP.TRANSITION.GRAPH]] = EVAL.EDIT1.DIST.NB.COMP.PROP
#csv.result.files[[BEST.K.FOR.KMEDOIDS]] = EVAL.BEST.K.FOR.KMEDOIDS
csv.result.files[[SINGLE.CLU.PROP.FOR.KMEDOIDS]] = EVAL.SINGLE.CLU.PROP.FOR.KMEDOIDS
csv.result.files[[IMB.COUNT]] = EVAL.IMB.INFO.TABLE.FILENAME
csv.result.files[[IMB.PERC]] = EVAL.IMB.INFO.TABLE.FILENAME
csv.result.files[[ALL.DIST.SCORE]] = EVAL.ALL.DIST.SCORES.FILE.PREFIX # prefix
csv.result.files[[CLASS.CORE.PART.SIZE.BY.SOL.CLASSES]] = EVAL.CLASS.CORE.PART.SIZE.FILE.NAME # prefix
csv.result.files[[CLASS.CORE.PART.SIZE]] = EVAL.CLASS.CORE.PART.SIZE.FILE.NAME # prefix


csv.related.colnames = list()
csv.related.colnames[[NB.SOLUTION]] = NB.SOL.COL.NAME
csv.related.colnames[[NB.SOLUTION.BY.SOL.CLASSES]] = NB.SOL.COL.NAME
csv.related.colnames[[NB.OPT.SOLUTION]] = NB.OPT.SOL.COL.NAME
csv.related.colnames[[EXEC.TIME]] = EXEC.TIME.COL.NAME
csv.related.colnames[[NB.CLUSTER]] = NB.CLU.COL.NAME
csv.related.colnames[[NB.COMP.TRANSITION.GRAPH]] =  c(EDIT1.COMP.MEM.COL.NAME, EDIT1234.COMP.MEM.COL.NAME)
#csv.related.colnames[[BEST.K.FOR.KMEDOIDS]] = BEST.K.FOR.SILH.COL.NAME,
csv.related.colnames[[SINGLE.CLU.PROP.FOR.KMEDOIDS]] = PROP.SINGLE.CLU.FOR.SILH.COL.NAME
csv.related.colnames[[IMB.COUNT]] = IMB.COUNT.COL.NAME
csv.related.colnames[[IMB.PERC]] = IMB.PERC.COL.NAME
csv.related.colnames[[ALL.DIST.SCORE]] = DIST.SCORE # prefix
csv.related.colnames[[CLASS.CORE.PART.SIZE.BY.SOL.CLASSES]] = PROP.CLASS.CORE.PART.SIZE.COL.NAME
csv.related.colnames[[CLASS.CORE.PART.SIZE]] = PROP.CLASS.CORE.PART.SIZE.COL.NAME



