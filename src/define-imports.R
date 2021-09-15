#############################################################################################
# Just loads the scripts necessary to process the dataset. This script is needed when
# using foreach (parallel processing), since each worker must load all these dependencies.
# 
#############################################################################################

source("src/define-constants.R")
source("src/define-algos.R")
source("src/common.R")
source("src/graph-common.R")
source("src/define-paths.R")
source("src/define-random-generator.R")

source("src/plot-tools/plot-network.R")
source("src/plot-tools/define-alluvial-diagram.R")
source("src/plot-tools/plot-delay-exec-times.R")

source("src/post-processing-input-networks/add-layouts-into-input-networks.R")

source("src/partition-networks/partition-networks.R")
source("src/partition-networks/load-membership.R")
source("src/partition-networks/create-gephi-networks.R")

source("src/evaluate-partitions/common.R")
source("src/evaluate-partitions/evaluate-partitions.R")
source("src/evaluate-partitions/evaluate-EnumCC.R")
source("src/evaluate-partitions/evaluate-imbalance.R")
source("src/evaluate-partitions/compare-partitions.R")
source("src/evaluate-partitions/take-average-over-networks.R")

source("src/reorganize/reorganize-csv-results-by-solution-class.R")
source("src/reorganize/reorganize-csv-results-by-detected-imbalance.R")


source("src/evaluate-partitions/create-relative-plot-membership.R")
source("src/evaluate-partitions/evaluate-partitions-with-kmedoid-results.R")
source("src/evaluate-partitions/evaluate-partitions-with-cluster-characterization.R")

source("src/make-plots/retreive-results.R")
source("src/make-plots/plot-common.R")
source("src/make-plots/make-layout-plots.R")
source("src/make-plots/define-custom-plot-functions.R")

source("src/cluster-analysis/define-purity.R")
source("src/cluster-analysis/perform-cluster-analysis.R")
source("src/cluster-analysis/show-cluster-analysis-difference.R")

source("src/cluster-characterization/perform-cluster-characterization.R")

source("src/collect-csv/collect-best-k-for-kmedoids.R")
source("src/collect-csv/collect-nb-opt-sol.R")
source("src/collect-csv/collect-core-parts.R")
source("src/collect-csv/collect-exec-time.R")
source("src/collect-csv/collect-delay-exec-time.R")
source("src/collect-csv/collect-Edit-dist-connected-comps.R")

