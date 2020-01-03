

# ===============================================================
#
# set CPLEX.BIN.PATH in src/define-algos.R
# ===============================================================


## libraries for parallel processing
#library(foreach)
#library(doParallel)

source("src/define-imports.R")


plot.layout <- c(
    # "kamada.kawai"
    # "fruchterman.reingold"
    # "bn_zheng"
    "circle"
)

#####################
plot.format <- c(
    #PLOT.AS.PDF
    #PLOT.AS.JPEG
    #PLOT.AS.PNG
    JUST.PLOT
)
#COMPRESS=FALSE
FORCE = FALSE

keep.algo.log.files = FALSE

# COR.CLU.EXACT.ALGO = get.ExCC.code(enum.all=FALSE)  # only for 1 optimal solution ==> it should be used only for debug or to know the optimal cost
COR.CLU.EXACT.ALGO = get.ExCC.code(enum.all=TRUE)   # for all optimal solutions


# not all the measures are distance measure, 
# but we will convert all similarity measure into distance in the function 'compare.partition.pair'
# IMPORTANT: if you want to calculate VI and EDIT distances, which are not normalized measures,
#           and if you want to compute both the unnormalized and normalized versions:
#           So, write first the unnormalized version, then add the normalized version. For instance, "VI", "NVI"
#           if you want to use other distance measures, ordering is not important.
#           We do this, because calculating EDIT distance might be time consuming, so its normalization should not be
COMP.MEASURES <- c(
		#VI # variation of information
		NVI # according to the paper: Xuan Vinh et al. (2010)
#		EDIT # edit distance
#		EDIT.NORM # normalized edit distance
)
# TODO: in compare-partitions.R, handle better the global variables in clusters

# Note that we use VI to get the closest solution association between ExCC optimal solutions and heuristic solutions. 
# Actually, what we want is just to know if a heuristic solution is optimal or not. so, distance is 0 or not..

#
#K.LIMITS = c(3,3) # LIMIT K-MEDOIDS RESULTS
K.LIMITS = c(NA,NA) # LIMIT K-MEDOIDS RESULTS


# #GRAPH.SIZES = seq(16, 36, by=4)
# GRAPH.SIZES = c(20, 24)
# DENSITY = 1
# #K = 2
# #K = 3
# K = 4

GRAPH.SIZES = c(20,24)

#K = 2
#K = 3
K = 4

PROP.MISPLS = c( seq(0.05, 0.20, by=0.05))

#PROP.MISPLS = c( seq(0.05, 1.00, by=0.05))

#PROP.MISPLS = c( seq(0.05, 0.40, by=0.05)) # k_init=4
#PROP.MISPLS = c( seq(0.05, 0.60, by=0.05)) # k_init=3
#PROP.MISPLS = c( seq(0.05, 0.85, by=0.05) ) # k_init=2

DETECTED.IMB.INTERVALS = c()
DETECTED.IMB.INTERVAL.SEQ.VALS = seq(0.00, 1.00, 0.05)
for(i in 1:(length(DETECTED.IMB.INTERVAL.SEQ.VALS)-1)){
    lower.bound = DETECTED.IMB.INTERVAL.SEQ.VALS[i]
    upper.bound = DETECTED.IMB.INTERVAL.SEQ.VALS[i+1]
    desc = paste0("[",lower.bound,",",upper.bound,")")
    DETECTED.IMB.INTERVALS = c(DETECTED.IMB.INTERVALS, desc)
}


# graph density
DENSITY = 0.125 # 0.250 ====> we can only give a single value (not a vector !!)


# ----------------------------------------------
PROP.NEGS = c(0.1, 0.2, 0.3)

##PROP.NEGS = NA # DO NOT DEFINE FOR NOW
# Note that PROP.NEGS should be 'NA' when DENSITY=1.0
# In the point of view of a developer, if PROP.NEGS = NA, we should know that prop.neg can be only 1 value and it can be computed from (graph.size, prop.mispl)
# Otherwise, user should be able to give a common range of prop.negs for each pair of (graph.size, prop.mispl)
# ----------------------------------------------

INPUT.RANDOM.NETWORKS = seq(1, 2, by=1)
# INPUT.RANDOM.NETWORKS = seq(1, 10, by=1)
# --------------------------------------------------------------




# ====================================
# TODO
# Likewise, the code works for the signed graphs with integer weights ?
# TODO: for(graph.desc.name in c(SIGNED.UNWEIGHTED.FILE)){ # SIGNED.WEIGHTED.FILE
#           ..... ==> boyle olmasin, main.R'de bunlari variable'a al
#
# - it the responsability of the user who will ensure if the RAM requirement of his/her system is ok for large graphs,
#     since Cplex may require a lot of RAM for graphs whose size is larger than 28
# - scale of y-axis should not be normalized for VI (or for another unnormalized measure)
# - transition graph ? plots?
# - opposite networks ? I think, it is another subject/work
# - the choice of internal evaluation measure is not possible for now. We use onyl Silhouette measure.
#       Later, we can use other measures, if needed. But, the code needs to be changed a bit.
# - tracking the used memory during Cplex is probably system dependent, i.e. only Linux ..
# - manage the output of ExCC:
#   - TODO: for Cplex output, we have the following Java warning: In readLines(con) : incomplete final line found on './out/partitions/../sol49.txt'
# - the whole scripts are run for all types of input graphs ? like SIGNED.WEIGHTED.FILE ? TODO
#       But, it is sure that directed signed graphs are not supported
# TODO: in compare-partitions.R, handle better the global variables in clusters
# TODO: in 'evaluate-partitions.R', you can also put the used RAM statistic as a output csv file
# TODO: Actually, prop.neg can be computed without knowing prop.mispl when d=1. You may simplify some process in layout plot prcoess?
#
# Remark: CORE.PART.THRESHOLD should be set to 1, not less than 1. Because, when it is less than 1, there might be multiple way of showing core part.
# ====================================



#################################################################
# POST-PROCESSING INPUT NETWORKS ==> since layout algos are stochastic, the obtained layout is used whenever plot is used
#################################################################
add.layouts.for.all.networks(GRAPH.SIZES, DENSITY, K, PROP.MISPLS, PROP.NEGS, INPUT.RANDOM.NETWORKS,
                    UNSIGNED.GRAPH.LAYOUTS, SIGNED.GRAPH.LAYOUTS)



##################################################################
# PARTITION NETWORKS
##################################################################
partition.networks(GRAPH.SIZES, DENSITY, K, PROP.MISPLS, PROP.NEGS, INPUT.RANDOM.NETWORKS,
	                COR.CLU.EXACT.ALGO, keep.algo.log.files, plot.format, plot.layout, FORCE)


#################################################################
# CREATE GEPHI NETWORKS with partition info for all solutions ==> it is a bit slower (depending on the number of obtained partitions)
#################################################################
create.gephi.networks(GRAPH.SIZES, DENSITY, K, PROP.MISPLS, PROP.NEGS, INPUT.RANDOM.NETWORKS, COR.CLU.EXACT.ALGO, FORCE)






#################################################################
# EVALUATE PARTITIONS
#################################################################
evaluate.partitions(GRAPH.SIZES, DENSITY, K, PROP.MISPLS, PROP.NEGS, INPUT.RANDOM.NETWORKS,
                    COR.CLU.EXACT.ALGO, COMP.MEASURES, FORCE)




#################################################################
# CLUSTER ANALYSIS
#################################################################
K.LIMITS = c(1,20)
perform.all.cluster.analysis(GRAPH.SIZES, DENSITY, K, PROP.MISPLS, PROP.NEGS, INPUT.RANDOM.NETWORKS,
        COR.CLU.EXACT.ALGO, COMP.MEASURES, K.LIMITS, FORCE)


# # ################################################################
# # # EVALUATE PARTITIONS WITH KMEDOID
# # ################################################################
evaluate.partitions.with.kmedoid.results(GRAPH.SIZES, DENSITY, K, PROP.MISPLS, PROP.NEGS, INPUT.RANDOM.NETWORKS,
        COR.CLU.EXACT.ALGO,
        #VI,
        COMP.MEASURES,
        FORCE)


# ##############################################################
# # CLUSTER CHARACTERIZATION
# ##############################################################
perform.all.cluster.characterization(GRAPH.SIZES, DENSITY, K, PROP.MISPLS, PROP.NEGS, INPUT.RANDOM.NETWORKS,
            COR.CLU.EXACT.ALGO, COMP.MEASURES, K.LIMITS, FORCE)


#################################################################
# OLD ? - CLUSTER ANALYSIS DIFFERENCES - MEASURES ===> needs multiple comp measures ...
#################################################################
## COMP.MEASURES <- c(
##     NMI,
##     VI
## )
# show.all.cluster.analysis.difference(GRAPH.SIZES, DENSITY, K, PROP.MISPLS, PROP.NEGS, INPUT.RANDOM.NETWORKS,
#         COR.CLU.EXACT.ALGO, COMP.MEASURES, K.LIMITS, FORCE)


# ##################################################################
# ## EVALUATE PARTITIONS WITH CORE PART ANALYSIS AND REPRESENTATIVES
# ##################################################################
evaluate.partitions.with.cluster.characterization(GRAPH.SIZES, DENSITY, K, PROP.MISPLS, PROP.NEGS, INPUT.RANDOM.NETWORKS,
                COR.CLU.EXACT.ALGO,
                #VI,
                COMP.MEASURES,
                FORCE)







##################################################################
## REORGANIZE
##################################################################
reorganize.all.csv.results.by.detected.imbalance(GRAPH.SIZES, DENSITY, K, PROP.MISPLS, PROP.NEGS, INPUT.RANDOM.NETWORKS,
                                              COR.CLU.EXACT.ALGO, FORCE)

# # for proportion of misplaced AND detected imbalance intervals
reorganize.all.csv.results.by.sol.class(GRAPH.SIZES, DENSITY, K, PROP.MISPLS, DETECTED.IMB.INTERVALS, PROP.NEGS, INPUT.RANDOM.NETWORKS,
                                             COR.CLU.EXACT.ALGO, COMP.MEASURES, FORCE)




###################################################################
### COMBINE EVALUATED PARTITIONS
###################################################################
## for proportion of misplaced  AND  detected imbalance intervals
take.average.over.networks(GRAPH.SIZES, DENSITY, K, PROP.MISPLS, DETECTED.IMB.INTERVALS, PROP.NEGS, INPUT.RANDOM.NETWORKS,
        COR.CLU.EXACT.ALGO, COMP.MEASURES, FORCE)






#-----------------------------------------------------------------------------------------------------------
##################################################################
## LAYOUT PLOTS
##################################################################
make.layout.plots.by.graph.size.and.imbalance(GRAPH.SIZES, DENSITY, K, PROP.MISPLS, DETECTED.IMB.INTERVALS, PROP.NEGS, INPUT.RANDOM.NETWORKS,
	COR.CLU.EXACT.ALGO, COMP.MEASURES, FORCE, plot.format)

make.layout.plots.by.graph.size.and.network.no(GRAPH.SIZES, DENSITY, K, PROP.MISPLS, DETECTED.IMB.INTERVALS, PROP.NEGS, INPUT.RANDOM.NETWORKS,
	COR.CLU.EXACT.ALGO, COMP.MEASURES, FORCE, plot.format)

make.layout.plots.by.imbalance.and.network.no(GRAPH.SIZES, DENSITY, K, PROP.MISPLS, DETECTED.IMB.INTERVALS, PROP.NEGS, INPUT.RANDOM.NETWORKS,
	COR.CLU.EXACT.ALGO, COMP.MEASURES, FORCE, plot.format)


###################################################################
### SUMMARY FOLDER - LAYOUT PLOTS
###################################################################
make.layout.plots.by.graph.size.and.imbalance(GRAPH.SIZES, DENSITY, K, PROP.MISPLS, DETECTED.IMB.INTERVALS, PROP.NEGS, c(SUMMARY.FOLDER.NAME),
	 COR.CLU.EXACT.ALGO, COMP.MEASURES, FORCE, plot.format)





#################################################################
# WRITE SOME RESULTS INTO CSV FILES
#################################################################

collect.all.best.k.for.kmedoids(GRAPH.SIZES, DENSITY, K, PROP.MISPLS, PROP.NEGS, INPUT.RANDOM.NETWORKS,
		COR.CLU.EXACT.ALGO, COMP.MEASURES, FORCE, plot.format)

collect.all.nb.opt.solution(GRAPH.SIZES, DENSITY, K, PROP.MISPLS, PROP.NEGS, INPUT.RANDOM.NETWORKS,
		COR.CLU.EXACT.ALGO, COMP.MEASURES, FORCE, plot.format)

collect.core.parts(GRAPH.SIZES, DENSITY, K, PROP.MISPLS, PROP.NEGS, INPUT.RANDOM.NETWORKS,
	COR.CLU.EXACT.ALGO, COMP.MEASURES, FORCE, plot.format)

