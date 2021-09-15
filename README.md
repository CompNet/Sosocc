# Sosocc
Space of Optimal Solutions of the Correlation Clustering Problem

* Copyright 2020-21 Nejat Arınık

Sosocc is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation. For source availability and license information see the file `LICENCE`

* Lab site: http://lia.univ-avignon.fr/
* GitHub repo: https://github.com/CompNet/Sosocc
* Contact: Nejat Arınık <arinik9@gmail.com>, Vincent Labatut <vincent.labatut@univ-avignon.fr>


-----------------------------------------------------------------------

# Description
This set of `R` scripts was designed to analyze the space of Optimal Solutions of the *Correlation Clustering Problem*. When solving an instance of such problem, several or even many optimal solutions (i.e. partitions) may coexist. If multiple optimal partitions coexist, one can then wonder how different/diverse they are. Put differently, we want to know what we loose when considering only one solution, while there might be multiple ones. In order to answer these questions, one should ideally enumerate completely the space of optimal solutions, and perform its analysis. To this end, we propose a new efficient solution space enumeration method and a cluster analysis-based framework in order to first enumerate the space of optimal partitions and then empirically study such space.

In this repository, we are able to run three exact partitioning methods for the CC problem:
* *ExCC*: It is a method aiming to obtain a single optimal solution. Its source code is found [here](https://github.com/CompNet/ExCC). It first employs a root relaxation phase, where we add violated valid inequalities through Cutting Plane method, then it proceeds to the branching phase. Use `get.ExCC.code(enum.all=FALSE)` in the code.
* *OneTreeCC* (also called *ExCC-all*): It is a solution space enumeration method incorprated in the commercial solver CPLEX.  Its source code is found [here](https://github.com/CompNet/ExCC). This two-step method first build and explore the search tree, i.e. B&B tree, in order to find efficiently the first optimal solution, then enumerate all the other optimal solutions based on the same tree. Use `get.ExCC.code(enum.all=TRUE)` in the code.
* *EnumCC*: It is a solution space enumeration method that we propose in *[Arınık'21]*. Its source code is found [here](https://github.com/CompNet/EnumCC). It takes in input a distance parameter `r_{max}` to explore the neighborhoof of an optimal solution. In our experiments, we have shown that it is more convenient to set `r_{max}=3` in general. Use `get.EnumCC.code(maxNbEdit=3)` in the code. Note that you need to run `ExCC` before calling `EnumCC`, since the latter needs an already-found optimal solution.


# Data
Our tool is applied to a set of signed networks generated thanks to our [signed graph generator](https://github.com/CompNet/SignedBenchmark). The details about the generator are explained [here](https://www.overleaf.com/read/pdggvqbsmrch). You may consider downloading our generated signed graphs associated with our article *[Arınık'20]*, which are in the folder `Input Signed Networks.tar.gz` on [FigShare](https://doi.org/10.6084/m9.figshare.8233340). The advantage of dowloading our data is that we also provide you with the optimal solutions for each signed network (in the folder `All Partitions Results.tar.gz` on [FigShare](https://doi.org/10.6084/m9.figshare.8233340). 

Moreover, You may also consider downloading our second set of generated signed graphs (both complete and incomplete signed graphs) associated with our submitted article *[Arınık'21]* on [FigShare](https://doi.org/10.6084/m9.figshare.15043911).

To show explicitly the folder structure used in the signed graph generation and for a quick test, we have already put some generated networks in `in/random-networks` and *some* corresponding optimal partitions in `out/partitions`. 


# Organization
Here are the folders composing the project:
* Folder `src`: contains the source code (R scripts).
* Folder `in`: contains the generated signed networks. 
* Folder `lib`: contains executable files related to the used external partitioning methods.
  * Folder `ExCC`: Executable file of the method `ExCC` whose the name will be `cplex-partition.jar`. See the *Installation* section for more details.
* Folder `out`: contains the folders and files produced by our scripts. See the *Use* section for more details.


# Installation
1. Install the [`R` language](https://www.r-project.org/)
2. Install the following R packages:
   * [`igraph`](http://igraph.org/r/) Tested with the version 1.2.6.
   * [`XML`](https://cran.r-project.org/web/packages/XML/index.html)
   * [`alluvial`](https://cran.r-project.org/web/packages/alluvial/)
   * [`cluster`](https://cran.r-project.org/web/packages/cluster/)
   * [`stringr`](https://cran.r-project.org/web/packages/stringr/)
   * [`plotrix`](https://cran.r-project.org/web/packages/plotrix/)
   * [`RColorBrewer`](https://cran.r-project.org/web/packages/RColorBrewer/)
   * [`vioplot`](https://cran.r-project.org/web/packages/vioplot/)
   * [`clues`](https://cran.r-project.org/web/packages/clues/)
   * [`NMF`](https://cran.r-project.org/web/packages/NMF/)
   * [`entropy`](https://cran.r-project.org/web/packages/entropy/)
   * [`e1071`](https://cran.r-project.org/web/packages/e1071/)
   * [`parallel`](https://cran.r-project.org/web/packages/parallel/)
   * [`doParallel`](https://cran.r-project.org/web/packages/doParallel/)
   * [`iterators`](https://cran.r-project.org/web/packages/iterators/)
   * [`bigmemory`](https://cran.r-project.org/web/packages/bigmemory/)
   * [`expm`](https://cran.r-project.org/web/packages/expm/)
3. Install [`IBM CPlex`](https://www.ibm.com/developerworks/community/blogs/jfp/entry/CPLEX_Is_Free_For_Students?lang=en). Tested with the versions 12.8 and 20.1. Set correctly the variable `CPLEX.BIN.PATH` in `define-algos.R` (e.g. `/opt/ibm/ILOG/CPLEX_Studio128/cplex/bin/x86-64_linux/`).
   * For ubuntu, type the following command:
     * `sudo ./cplex_studio<YOUR_VERSION>.linux-x86-64.bin` 
       * The default installation location for education version is: `/opt/ibm/ILOG/CPLEX_Studio<YOUR_VERSION`.
       * The default installation location for trial version is: `/opt/ibm/ILOG/CPLEX_Studio_Community<YOUR_VERSION/cplex/bin/x86-64_linux/`.
4. Download the project of `ExCC` on [github](https://github.com/arinik9/ExCC). First, configure and then compile it. To test it, you can run the file `run.sh`.If everything works (i.e. if a file `sol0.txt` created in the output folder), move the executable file `ExCC.jar`, which is in `exe`, into the `lib/ExCC` folder in this project.
5. Download the project of `ClusteringEditDistance` on [github](https://github.com/arinik9/ClusteringEditDistance). Move the executable file `ClusteringEditDist.jar` into the lib `folder ` in this project. This jar file allows to compute the edit distance between membership vectors.
6. Download the signed networks on [`figshare`](https://doi.org/10.6084/m9.figshare.8233340) or generate your own signed networks based on our [signed graph generator](https://github.com/CompNet/SignedBenchmark).


# Use
1. Set correctly the variables `CPLEX.BIN.PATH`.
2. Open the `R` console.
3. Set the current directory as the working directory, using `setwd("<my directory>")`.
4. Run the main script `src/main.R`.


The script will produce the following subfolders in the folder `out`:
* `partitions`: Folder containing all obtained partitions. Note that we prefer to use the term *module* in the same sense of *cluster* in this level in order to distinguish it from the term *cluster* used in k-medoids clustering analysis.
  * `n=<GRAPH_SIZE>_l0=<NB_INIT_MODULES>_dens=<DENSITY>`: A subfolder for specific values of graph size, number of initial modules in the network and density.
    * `propMispl=<PROP_MISPL>`: A subfolder for a specific value of proportion of misplaced links, which is used for signed network generation.
      * `propNeg=<PROP_NEG>`: A subfolder for a specific value of proportion of negative links
        * `network=<NETWORK_ID>`: A subfolder for a specific value of network id (for replication pruposes).
          * `<CC_ALGO_NAME>`: A subfolder for a specific partitioning method solving the CC problem.
        	* `<GRAPH_TYPE>`: A subfolder for a specific network type, e.g. `signed_unweighted`. The distinction of this level is actually useless for now.
* `evaluate-partitions`: Folder containing all results in csv format.
  * `n=<GRAPH_SIZE>_l0=<NB_INIT_MODULES>_dens=<DENSITY>`: A subfolder for specific values of graph size, number of initial modules in the network and density.
    * `propMispl=<PROP_MISPL>`: A subfolder for a specific value of proportion of misplaced links, which is used for signed network generation.
      * `propNeg=<PROP_NEG>`: A subfolder for a specific value of proportion of negative links
		* `<K_DESC>`: This folder is used to distinguish two cases: 1) *All*, 2) `k=<K_VALUE>-<COMP_MEASURE>`. The former exists for the case where the notion of solution classs is not considered. In the latter, the networks are reorganized based on the detected number of solution classes. This is convenient to see and plot the results.
          * `network=<NETWORK_ID>`: A subfolder for a specific value of network id (for replication pruposes).
            * `<CC_ALGO_NAME>`: A subfolder for a specific partitioning method solving the CC problem.
        	  * `<GRAPH_TYPE>`: A subfolder for a specific network type, e.g. `signed_unweighted`. The distinction of this level is actually useless for now.
    * `detectedImbalance=<IMB_PROP_INTERVAL>`: This folder allows to regroup the networks for a specific interval of detected imbalance (e.g. *0.05-0.1*). This is different from `propMispl`. In other words, this folder does not contain new information compared to `propMispl=<PROP_MISPL>`, it is just a reorganization of results by a different aspect. This is convenient to see and plot the results.
      * `propNeg=<PROP_NEG>`: A subfolder for a specific value of proportion of negative links.
		* `<K_DESC>`: This folder is used to distinguish two cases: 1) *All*, 2) `k=<K_VALUE>-<COMP_MEASURE>`. The former exists for the case where the notion of solution classs is not considered. In the latter, the networks are reorganized based on the detected number of solution classes. This is convenient to see and plot the results.
          * `network=<NETWORK_ID>`: A subfolder for a specific value of network id (for replication pruposes). Note that, network ids might be differant than those in `propMispl=<PROP_MISPL>`. So, they are renumbered.
            * `<CC_ALGO_NAME>`: A subfolder for a specific partitioning method solving the CC problem.
        	  * `<GRAPH_TYPE>`: A subfolder for a specific network type, e.g. `signed_unweighted`. The distinction of this level is actually useless for now.
* `cluster-analysis`: Folder containing the results related to k-medoids clustering analysis.
  * `n=<GRAPH_SIZE>_l0=<NB_INIT_MODULES>_dens=<DENSITY>`: A subfolder for specific values of graph size, number of initial modules in the network and density.
    * `propMispl=<PROP_MISPL>`: A subfolder for a specific value of proportion of misplaced links, which is used for signed network generation.
      * `propNeg=<PROP_NEG>`: A subfolder for a specific value of proportion of negative links
        * `network=<NETWORK_ID>`: A subfolder for a specific value of network id (for replication pruposes).
          * `<CC_ALGO_NAME>`: A subfolder for a specific partitioning method solving the CC problem.
        	* `<GRAPH_TYPE>`: A subfolder for a specific network type, e.g. `signed_unweighted`. The distinction of this level is actually useless for now.
        	  * `<COMP_MEASURE>`: A subfolder for a specific value of external comparison measure, e.g. `NMI`
        		* `k=<K_VALUE>-sil=<SILHOUETTE_SCORE>`: A subfolder containing the results of the k-medoids clustering result for a specific value of `k`, where `k` is the desired number of solution classes. 
				  * `clu<SOLUTION_CLASS_ID>`: In this subfolder(s) the partition files related to a given signed network are divided into `k` folders. We create those subfolder(s) only for the value of `k` which is associated with the best silhouette score.
* `cluster-characterization`: Folder containing the results related tocluster characterization process. 1) core part, 2) representative
  * `n=<GRAPH_SIZE>_l0=<NB_INIT_MODULES>_dens=<DENSITY>`: A subfolder for specific values of graph size, number of initial modules in the network and density.
    * `propMispl=<PROP_MISPL>`: A subfolder for a specific value of proportion of misplaced links, which is used for signed network generation.
      * `propNeg=<PROP_NEG>`: A subfolder for a specific value of proportion of negative links
        * `network=<NETWORK_ID>`: A subfolder for a specific value of network id (for replication pruposes).
          * `<CC_ALGO_NAME>`: A subfolder for a specific partitioning method solving the CC problem.
        	* `<GRAPH_TYPE>`: A subfolder for a specific network type, e.g. `signed_unweighted`. The distinctin of this level is actually useless for now.
        	  * `<COMP_MEASURE>`: A subfolder for a specific value of external comparison measure, e.g. `NMI`
        		* `k=<K_VALUE>-sil=<SILHOUETTE_SCORE>`: A subfolder containing the results of the k-medoids clustering result for a specific value of `k`, where `k` is the desired number of solution classes. Note that Silhouette is not defined for `k=1` and `k=n`.
				  * `core-part=<CORE_PART_THRESHOLD>`: In this subfolder, the core part of each solution class, as well as all the partitions, are shown. Our definition of *core part* depends on `CORE_PART_THRESHOLD`. If this value equals 1, this means that the core part of a solution class consists of being groups of nodes that are always together over all considered partitions. Actually, this value should be always 1 if one wants to have only a single core part structure.
				  * `representative-partitions`: In this subfolder, the optimal partitions that are considered *representative* for each solution class are found. By *representative*, we mean the partition which is most similar to others in the same solution class. Additionnaly, we place in this folder the presentative partition of all the partitions. Note that, even though there are multiple candidates to satisfy our definition of *representative*, we keep only a single one (but this is rare)
* `layout-plots`: The subfolder hierarchy is similar to the previous folders. This folder will contain the plots of type *layout*, i.e. each plot contains several subplots corresponding to an input parameter, and each subplot is divided into multiple parts corresponding to different values of a second input parameter. Therefore, in each plot file, two input parameters are investigated, the other parameters are fixed. Note that we also make distinction of two cnocepts: 1) single network layout plot, and 2) summary network layout plot. In the former, each network is treated separately. Although this is good to plot, what is most interesting is the latter, which corresponds to a combination of all networks under the same input parameter set.
* `output-csv-data`: Some results are recorded in csv format (e.g. number of optimal solutions per network).


# References
* **[Arınık'20]** N. Arınık & R. Figueiredo & V. Labatut. *Multiplicity and Diversity: Analyzing the Optimal Solution Space of the Correlation Clustering Problem*, Journal of Complex Networks, 8(6) 2020. [doi: https://doi.org/10.1093/comnet/cnaa025](https://github.com/CompNet/Sosocc)
* **[Arınık'21]** N. Arınık & R. Figueiredo & V. Labatut. *Efficient Enumeration of Correlation Clustering Optimal Solution Space (submitted)*, Journal of Global Optmization, 2021.
