
source("core-part-analysis-utils.R")
#library(igraph)

POSNEG.DESC.NAME = "pos/neg"
NEG.TRI.DESC.NAME = "neg triangle"
QUASI.NEG.TRI.DESC.NAME = "quasi neg triangle"
POT.NEG.TRI.DESC.NAME = "potential neg triangle"
POS.TRI.DESC.NAME = "pos triangle"
QUASI.POS.TRI.DESC.NAME = "quasi pos triangle"
POT.POS.TRI.DESC.NAME = "potential pos triangle"
EGO.DESC.NAME = "ego2"
DIST.DES.NAME = "mean dist"
EIGEN.DESC.NAME = "eigen"
AUTH.DESC.NAME = "authority"
ALPHA.DESC.NAME = "alpha"
POWER.DESC.NAME = "power"
PRANK.DESC.NAME = "page rank"
SUBG.DESC.NAME = "subgraph"
BETW.DESC.NAME = "betweenness"
CLOSE.DESC.NAME = "closeness"
ECCEN.DESC.NAME = "eccentricity"
ARTPOI.DESC.NAME = "articulation point"
CORE.DESC.NAME = "coreness"
DIV.DESC.NAME = "diversity"
DEGR.DESC.NAME = "degree"
#STR.DESC.NAME = "strength"
LOC.TRANS.DESC.NAME = "local transitivity"
BURT.DESC.NAME = "Burst's constraint"



make.stats.for.core.part.analysis = function(graph) {
    result = c()

    A = as_adj(g, attr=NULL)

    result = rbind(result, sprintf("%.4f",ego_size(graph, order=2, nodes = V(graph))) )
    dists.mtrx = distances(graph, weights = NULL, algorithm = "unweighted")
    result = rbind(result, round(apply(dists.mtrx,2,mean),digits=4) )
    result = rbind(result, sprintf("%.4f",eigen_centrality(graph, directed = FALSE, scale = TRUE, weights = NULL)$vector) )

    result = rbind(result, sprintf("%.4f",authority_score(graph, scale = TRUE, weights = NULL)$vector) )
    # result = rbind(result, sprintf("%.4f",hub_score(graph, scale = TRUE, weights = NULL)$vector) )

    result = rbind(result, sprintf("%.4f",alpha_centrality(graph, nodes = V(graph), alpha = 1, loops = FALSE, exo = 1, sparse = FALSE)) )
    result = rbind(result, sprintf("%.4f",power_centrality(graph, nodes = V(graph), loops = TRUE, exponent = 1, rescale = FALSE, sparse = FALSE)) )
    result = rbind(result, sprintf("%.4f",page_rank(graph, algo = "prpack", vids = V(graph), directed = FALSE, damping = 0.85, weights = NULL)$vector) )

    result = rbind(result, sprintf("%.4f",subgraph_centrality(graph, diag = FALSE)) )
    result = rbind(result, sprintf("%.4f",betweenness(graph, v = V(graph), directed = FALSE, weights = NULL, nobigint = TRUE, normalized = FALSE)) )
    result = rbind(result, sprintf("%.4f",closeness(graph, vids = V(graph), weights = NULL, normalized = FALSE)) )
    result = rbind(result, sprintf("%.4f",eccentricity(graph, vids = V(graph))) )

    art.poi.ids =  articulation_points(graph)
    art.poi.res = rep(0, vcount(graph))
    art.poi.res[art.poi.ids] = 1
    result = rbind(result, art.poi.res )
    # --------------
    result = rbind(result, sprintf("%.4f",coreness(graph)) )
    #?? Participation
    #?? Intern Intensity
    # ?? Extern Intensity
    result = rbind(result, sprintf("%.4f",diversity(graph, weights = NULL, vids = V(graph))) )
    # ?? heterogeneity
    # --------------
    result = rbind(result, sprintf("%.4f",degree(graph, v = V(graph), loops = FALSE, normalized = FALSE)) )
    #result = rbind(result, sprintf("%.4f",strength(graph, vids = V(graph), loops = FALSE, weights = NULL)) )
    result = rbind(result, sprintf("%.4f",transitivity(graph, type = "local", vids = V(graph), weights = NULL)) )
    result = rbind(result, sprintf("%.4f",constraint(graph, nodes = V(graph), weights = NULL)) )

    rownames(result) = c(
                        EGO.DESC.NAME ,DIST.DES.NAME, EIGEN.DESC.NAME, AUTH.DESC.NAME, ALPHA.DESC.NAME,
                        POWER.DESC.NAME, PRANK.DESC.NAME, SUBG.DESC.NAME, BETW.DESC.NAME, CLOSE.DESC.NAME, ECCEN.DESC.NAME, ARTPOI.DESC.NAME,
                        CORE.DESC.NAME, DIV.DESC.NAME, DEGR.DESC.NAME, LOC.TRANS.DESC.NAME, BURT.DESC.NAME 
                    )
    colnames(result) = paste0("V",seq(0,(vcount(graph)-1)))

    return(result)
}

######################################################################
######################################################################
######################################################################


g = read.graph("signed-unweighted.graphml","graphml")
A = as_adj(g, attr="weight")
gpos = delete.edges(graph=g, edges=which(E(g)$weight<0))
Apos = as_adj(gpos)
gneg = delete.edges(graph=g, edges=which(E(g)$weight>0))
E(gneg)$weight = E(gneg)$weight * (-1) # make them positive
Aneg = as_adj(gneg)


signed.result = c()
pos.deg = degree(gpos, v = V(gpos), loops = FALSE, normalized = FALSE)
neg.deg = degree(gneg, v = V(gneg), loops = FALSE, normalized = FALSE)
ratio.pos.neg = sprintf("%.4f",pos.deg/neg.deg)
pos.tri.res = sapply(1:vcount(g), function(v.id) count.pos.triangle.at.node(A, v.id))
quasi.pos.tri.res = sapply(1:vcount(g), function(v.id) count.quasi.pos.triangle.at.node(A, v.id))
potential.pos.tri.res = sapply(1:vcount(g), function(v.id) count.potential.pos.triangle.at.node(A, v.id))
neg.tri.res = sapply(1:vcount(g), function(v.id) count.neg.triangle.at.node(A, v.id))
quasi.neg.tri.res = sapply(1:vcount(g), function(v.id) count.quasi.neg.triangle.at.node(A, v.id))
potential.neg.tri.res = sapply(1:vcount(g), function(v.id) count.potential.neg.triangle.at.node(A, v.id))
signed.result = rbind(signed.result, ratio.pos.neg, pos.tri.res, quasi.pos.tri.res, potential.pos.tri.res, neg.tri.res, quasi.neg.tri.res, potential.neg.tri.res)
rownames(signed.result) = c(POSNEG.DESC.NAME, POS.TRI.DESC.NAME ,QUASI.POS.TRI.DESC.NAME, POT.POS.TRI.DESC.NAME, NEG.TRI.DESC.NAME ,QUASI.NEG.TRI.DESC.NAME ,POT.NEG.TRI.DESC.NAME)
colnames(signed.result) = paste0("V",seq(0,(vcount(g)-1)))

print("obtaining results regarding positive graph")
pos.result = make.stats.for.core.part.analysis(gpos)
print("obtaining results regarding negative graph")
neg.result = make.stats.for.core.part.analysis(gneg)
print("writing into files")
write.csv(file="pos-result.csv", x=pos.result)
write.csv(file="neg-result.csv", x=neg.result)
write.csv(file="signed-result.csv", x=signed.result)
