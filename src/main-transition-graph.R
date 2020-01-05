
# =============================================================================
#   ==> Remark 1: It is supposed that partitions are already obtained.
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
		#NVI # according to the paper: Xuan Vinh et al. (2010)
	EDIT, # edit distance
	EDIT.NORM # normalized edit distance
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
# EVALUATE PARTITIONS - it requires to exist eval files related to edit distance
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


##################################################################
## EDIT DISTANCE BASED CONNECTED COMP
##################################################################
create.all.edit.dist.connected.comps(GRAPH.SIZES, DENSITY, L0, PROP.MISPLS, PROP.NEGS, INPUT.RANDOM.NETWORKS,
                    COR.CLU.EXACT.ALGO, COMP.MEASURES, FORCE)


###################################################################
### COMBINE EVALUATED PARTITIONS
###################################################################
take.average.over.networks(GRAPH.SIZES, DENSITY, L0, PROP.MISPLS, PROP.NEGS, INPUT.RANDOM.NETWORKS,
        COR.CLU.EXACT.ALGO, COMP.MEASURES, FORCE)



################################################################
# PLOT TRANSITION GRAPHS
################################################################
plot.all.transition.graphs(GRAPH.SIZES, DENSITY, L0, PROP.MISPLS, PROP.NEGS, INPUT.RANDOM.NETWORKS,
		COR.CLU.EXACT.ALGO, COMP.MEASURES, FORCE, plot.format)



#################################################################
# COLLECT INTO CSV
#################################################################
collect.nb.comp.in.all.transition.graphs(GRAPH.SIZES, DENSITY, L0, PROP.MISPLS, PROP.NEGS, INPUT.RANDOM.NETWORKS,
                            COR.CLU.EXACT.ALGO, COMP.MEASURES, FORCE)






