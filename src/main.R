
# =============================================================================
# VARIABLES
#   ==> REmark 1: Do not forget to set CPLEX.BIN.PATH correclty in src/define-algos.R
#   ==> Remark 2: PROP.NEGS should be set to 'NA' when DENSITY=1.0
#   ==> Remark 3: CORE.PART.THRESHOLD should be set to 1. Because, when it is less than 1,
#                  there might be multiple way of building core part.
#   ==> Remark 4: It is the responsability of the user who will ensure if the RAM requirement
#                  of his/her system is ok for large graphs, because Cplex may require 
#                  lots of RAM for graphs whose size is larger than 28.
# =============================================================================


## libraries for parallel processing
#library(foreach)
#library(doParallel)

source("src/define-imports.R")


plot.layout <- c( # ==========> it is not taken into account everywhere !! TODO
    # "kamada.kawai"
    # "fruchterman.reingold"
    # "bn_zheng"
    "circle"
)

#####################
plot.format <- c( # ==========> it is not taken into account everywhere !! TODO
    #PLOT.AS.PDF
    #PLOT.AS.JPEG
    #PLOT.AS.PNG
    JUST.PLOT
)

FORCE = FALSE

keep.algo.log.files = FALSE

# COR.CLU.EXACT.ALGO = get.ExCC.code(enum.all=FALSE)  # only for 1 optimal solution ==> it should be used only for debug or to know the optimal cost
COR.CLU.EXACT.ALGO = get.ExCC.code(enum.all=TRUE)   # for all optimal solutions


# not all the measures are distance measure, 
# but we will convert all similarity measure into distance in the function 'compare.partition.pair'
# IMPORTANT: if you want to calculate EDIT distance, which is not normalized measures,
#           and if you want to compute both the unnormalized and normalized versions:
#           So, write first the unnormalized version, then add the normalized version.
#           if you want to use other distance measures, ordering is not important.
#           We do this, because calculating EDIT distance might be time consuming, so its normalization should not be
COMP.MEASURES <- c(
		#VI # variation of information
		NVI # according to the paper: Xuan Vinh et al. (2010)
#		EDIT # edit distance
#		EDIT.NORM # normalized edit distance
)


#K.LIMITS = c(3,3) # LIMIT K-MEDOIDS RESULTS
K.LIMITS = c(NA,NA) # LIMIT K-MEDOIDS RESULTS


GRAPH.SIZES = c(20,24)
L0 = 4
PROP.MISPLS = c( seq(0.05, 0.20, by=0.05))
DENSITY = 0.125 # graph density: we can only give a single value (not a vector !!)
INPUT.RANDOM.NETWORKS = seq(1, 2, by=1)
# ----------------------------------------------
PROP.NEGS = c(0.2, 0.3)
# Note that PROP.NEGS should be 'NA' when DENSITY=1.0
# In the point of view of a developer, if PROP.NEGS = NA, we should know that prop.neg can be only 1 value and it can be computed from (graph.size, prop.mispl)
# Otherwise, user should be able to give a common range of prop.negs for each pair of (graph.size, prop.mispl)
# ----------------------------------------------

DETECTED.IMB.INTERVALS = c()
DETECTED.IMB.INTERVAL.SEQ.VALS = seq(0.00, 1.00, 0.05) # tis vector does not matter. In any case, we check if the folder exists or not
for(i in 1:(length(DETECTED.IMB.INTERVAL.SEQ.VALS)-1)){
    lower.bound = DETECTED.IMB.INTERVAL.SEQ.VALS[i]
    upper.bound = DETECTED.IMB.INTERVAL.SEQ.VALS[i+1]
    desc = paste0("[",lower.bound,",",upper.bound,")")
    DETECTED.IMB.INTERVALS = c(DETECTED.IMB.INTERVALS, desc)
}



##########################################################################
# Another scenario: This is used in the article
##########################################################################
#GRAPH.SIZES = seq(16, 36, by=4)
#PROP.MISPLS = c( seq(0.05, 0.40, by=0.05)) # when l0=4 is selected
#PROP.MISPLS = c( seq(0.05, 0.60, by=0.05)) # when l0=3 is selected
#PROP.MISPLS = c( seq(0.05, 0.85, by=0.05) ) # when l0=2 is selected
#DENSITY = 1 # we can only give a single value (not a vector !!)
#L0 = 2 # OR L0 = 3 OR L0 = 4
#PROP.NEGS = NA
#INPUT.RANDOM.NETWORKS = seq(1, 100, by=1)
# =======================================================================







# =============================================================================
# FUNCTIONS
# =============================================================================

#################################################################
# POST-PROCESSING INPUT NETWORKS ==> since layout algos are stochastic,
#       the obtained layout is used whenever plot is used
#################################################################
add.layouts.for.all.networks(GRAPH.SIZES, DENSITY, L0, PROP.MISPLS, PROP.NEGS, INPUT.RANDOM.NETWORKS,
                    UNSIGNED.GRAPH.LAYOUTS, SIGNED.GRAPH.LAYOUTS)


##################################################################
# PARTITION NETWORKS
##################################################################
partition.networks(GRAPH.SIZES, DENSITY, L0, PROP.MISPLS, PROP.NEGS, INPUT.RANDOM.NETWORKS,
	                COR.CLU.EXACT.ALGO, keep.algo.log.files, plot.format, plot.layout, FORCE)


#################################################################
# CREATE GEPHI NETWORKS with partition info for all solutions
#   ==> it is a bit slower (depending on the number of obtained partitions)
#################################################################
create.gephi.networks(GRAPH.SIZES, DENSITY, L0, PROP.MISPLS, PROP.NEGS, INPUT.RANDOM.NETWORKS, COR.CLU.EXACT.ALGO, FORCE)



#################################################################
# EVALUATE PARTITIONS
#################################################################
evaluate.partitions(GRAPH.SIZES, DENSITY, L0, PROP.MISPLS, PROP.NEGS, INPUT.RANDOM.NETWORKS,
                    COR.CLU.EXACT.ALGO, COMP.MEASURES, FORCE)


#################################################################
# CLUSTER ANALYSIS
#################################################################
K.LIMITS = c(1,20)
perform.all.cluster.analysis(GRAPH.SIZES, DENSITY, L0, PROP.MISPLS, PROP.NEGS, INPUT.RANDOM.NETWORKS,
        COR.CLU.EXACT.ALGO, COMP.MEASURES, K.LIMITS, FORCE)


# # ################################################################
# # # EVALUATE PARTITIONS WITH KMEDOID
# # ################################################################
evaluate.partitions.with.kmedoid.results(GRAPH.SIZES, DENSITY, L0, PROP.MISPLS, PROP.NEGS, INPUT.RANDOM.NETWORKS,
        COR.CLU.EXACT.ALGO,
        #VI,
        COMP.MEASURES,
        FORCE)


# ##############################################################
# # CLUSTER CHARACTERIZATION
# ##############################################################
perform.all.cluster.characterization(GRAPH.SIZES, DENSITY, L0, PROP.MISPLS, PROP.NEGS, INPUT.RANDOM.NETWORKS,
            COR.CLU.EXACT.ALGO, COMP.MEASURES, K.LIMITS, FORCE)


#################################################################
# OLD ? - CLUSTER ANALYSIS DIFFERENCES - MEASURES ===> needs multiple comp measures ...
#################################################################
## COMP.MEASURES <- c(
##     NMI,
##     VI
## )
# show.all.cluster.analysis.difference(GRAPH.SIZES, DENSITY, L0, PROP.MISPLS, PROP.NEGS, INPUT.RANDOM.NETWORKS,
#         COR.CLU.EXACT.ALGO, COMP.MEASURES, K.LIMITS, FORCE)


# ##################################################################
# ## EVALUATE PARTITIONS WITH CORE PART ANALYSIS AND REPRESENTATIVES
# ##################################################################
evaluate.partitions.with.cluster.characterization(GRAPH.SIZES, DENSITY, L0, PROP.MISPLS, PROP.NEGS, INPUT.RANDOM.NETWORKS,
                COR.CLU.EXACT.ALGO,
                #VI,
                COMP.MEASURES,
                FORCE)



##################################################################
## REORGANIZE
##################################################################
reorganize.all.csv.results.by.detected.imbalance(GRAPH.SIZES, DENSITY, L0, PROP.MISPLS, PROP.NEGS, INPUT.RANDOM.NETWORKS,
                                              COR.CLU.EXACT.ALGO, FORCE)

# # for proportion of misplaced AND detected imbalance intervals
reorganize.all.csv.results.by.sol.class(GRAPH.SIZES, DENSITY, L0, PROP.MISPLS, DETECTED.IMB.INTERVALS, PROP.NEGS, INPUT.RANDOM.NETWORKS,
                                             COR.CLU.EXACT.ALGO, COMP.MEASURES, FORCE)



###################################################################
### COMBINE EVALUATED PARTITIONS
###################################################################
## for proportion of misplaced  AND  detected imbalance intervals
take.average.over.networks(GRAPH.SIZES, DENSITY, L0, PROP.MISPLS, DETECTED.IMB.INTERVALS, PROP.NEGS, INPUT.RANDOM.NETWORKS,
        COR.CLU.EXACT.ALGO, COMP.MEASURES, FORCE)




#-----------------------------------------------------------------------------------------------------------
##################################################################
## LAYOUT PLOTS
##################################################################
make.layout.plots.by.graph.size.and.imbalance(GRAPH.SIZES, DENSITY, L0, PROP.MISPLS, DETECTED.IMB.INTERVALS, PROP.NEGS, INPUT.RANDOM.NETWORKS,
	COR.CLU.EXACT.ALGO, COMP.MEASURES, FORCE, plot.format)

make.layout.plots.by.graph.size.and.network.no(GRAPH.SIZES, DENSITY, L0, PROP.MISPLS, DETECTED.IMB.INTERVALS, PROP.NEGS, INPUT.RANDOM.NETWORKS,
	COR.CLU.EXACT.ALGO, COMP.MEASURES, FORCE, plot.format)

make.layout.plots.by.imbalance.and.network.no(GRAPH.SIZES, DENSITY, L0, PROP.MISPLS, DETECTED.IMB.INTERVALS, PROP.NEGS, INPUT.RANDOM.NETWORKS,
	COR.CLU.EXACT.ALGO, COMP.MEASURES, FORCE, plot.format)


###################################################################
### SUMMARY FOLDER - LAYOUT PLOTS
###################################################################
make.layout.plots.by.graph.size.and.imbalance(GRAPH.SIZES, DENSITY, L0, PROP.MISPLS, DETECTED.IMB.INTERVALS, PROP.NEGS, c(SUMMARY.FOLDER.NAME),
	 COR.CLU.EXACT.ALGO, COMP.MEASURES, FORCE, plot.format)





#################################################################
# WRITE SOME RESULTS INTO CSV FILES
#################################################################

collect.all.best.k.for.kmedoids(GRAPH.SIZES, DENSITY, L0, PROP.MISPLS, PROP.NEGS, INPUT.RANDOM.NETWORKS,
		COR.CLU.EXACT.ALGO, COMP.MEASURES, FORCE, plot.format)

collect.all.nb.opt.solution(GRAPH.SIZES, DENSITY, L0, PROP.MISPLS, PROP.NEGS, INPUT.RANDOM.NETWORKS,
		COR.CLU.EXACT.ALGO, COMP.MEASURES, FORCE, plot.format)

collect.core.parts(GRAPH.SIZES, DENSITY, L0, PROP.MISPLS, PROP.NEGS, INPUT.RANDOM.NETWORKS,
	COR.CLU.EXACT.ALGO, COMP.MEASURES, FORCE, plot.format)







# ====================================
# TODOs
# ====================================
#
# - the code works for the signed graphs with integer weights ?
# - for(graph.desc.name in c(SIGNED.UNWEIGHTED.FILE)){ # SIGNED.WEIGHTED.FILE
#       * do not do in that way, put the potetial values into global list in main.R
#       * In folder hierarchy, we need to add another child folder called "signed-unweighted" in in/random-networks/...
#
# - scale of y-axis should not be normalized for VI (or for another unnormalized measure)
# - the choice of internal evaluation measure is not possible for now. We use only Silhouette measure.
#       * Later, we can use other measures, if needed. But, the code needs to be changed a bit.
# - tracking the used memory during Cplex is probably system dependent, i.e. only Linux ..
# - manage the output of ExCC:
# - for Cplex output, we have the following Java warning: In readLines(con) : incomplete final line found on './out/partitions/../sol49.txt'
# - the whole scripts are run for all types of input graphs ? like SIGNED.WEIGHTED.FILE ?
#       * But, it is sure that directed signed graphs are not supported
# - in compare-partitions.R, handle better the global variables in clusters
# - in 'evaluate-partitions.R', you can also put the used RAM statistic as a output csv file
# - Actually, prop.neg can be computed without knowing prop.mispl when d=1. You may simplify some process in layout plot prcoess?
#
# # ====================================

