
source("core-part-analysis-utils2.R")
#library(igraph)




QUASI.POS.A1.TRI.DESC.NAME = "quasi pos triangle A1"
QUASI.POS.A2.TRI.DESC.NAME = "quasi pos triangle A2"
QUASI.NEG.TRI.A1.DESC.NAME = "quasi neg triangle A1"
QUASI.NEG.TRI.A2.DESC.NAME = "quasi neg triangle A2"
POS.TRI.A.DESC.NAME = "pos triangle A"
NEG.TRI.A.DESC.NAME = "neg triangle A"
QUASI.POS.B1.TRI.DESC.NAME = "quasi pos triangle B1"
QUASI.POS.B2.TRI.DESC.NAME = "quasi pos triangle B2"
QUASI.POS.B3.TRI.DESC.NAME = "quasi pos triangle B3"
QUASI.NEG.TRI.B1.DESC.NAME = "quasi neg triangle B1"
QUASI.NEG.TRI.B2.DESC.NAME = "quasi neg triangle B2"
QUASI.NEG.TRI.B3.DESC.NAME = "quasi neg triangle B3"
POS.TRI.B.DESC.NAME = "pos triangle B"
NEG.TRI.B.DESC.NAME = "neg triangle B"
QUASI.POS.C1.TRI.DESC.NAME = "quasi pos triangle C1"
QUASI.POS.C2.TRI.DESC.NAME = "quasi pos triangle C2"
QUASI.NEG.TRI.C1.DESC.NAME = "quasi neg triangle C1"
QUASI.NEG.TRI.C2.DESC.NAME = "quasi neg triangle C2"
POS.TRI.C.DESC.NAME = "pos triangle C"
NEG.TRI.C.DESC.NAME = "neg triangle C"
QUASI.POS.D1.TRI.DESC.NAME = "quasi pos triangle D1"
QUASI.POS.D2.TRI.DESC.NAME = "quasi pos triangle D2"
QUASI.NEG.TRI.D1.DESC.NAME = "quasi neg triangle D1"
QUASI.NEG.TRI.D2.DESC.NAME = "quasi neg triangle D2"
POS.TRI.D.DESC.NAME = "pos triangle D"
NEG.TRI.D.DESC.NAME = "neg triangle D"


######################################################################
######################################################################
######################################################################

aggrega.mod = "sum"
g = read.graph("signed-unweighted.graphml","graphml")
A = as_adj(g, attr="weight")
#membership = read.table(file="membership-repr.txt")$V1
membership = read.table(file="membership0.txt")$V1
result = c()


tri.A.res = c()
for(v.id in 1:vcount(g)){
    tri.A.res = cbind(tri.A.res, count.triangle.A.at.node(A, membership, v.id, aggrega.mod))
}
rownames(tri.A.res) = c(QUASI.POS.A1.TRI.DESC.NAME, QUASI.POS.A2.TRI.DESC.NAME, QUASI.NEG.TRI.A1.DESC.NAME, QUASI.NEG.TRI.A2.DESC.NAME, POS.TRI.A.DESC.NAME, NEG.TRI.A.DESC.NAME)


tri.B.res = c()
for(v.id in 1:vcount(g)){
    tri.B.res = cbind(tri.B.res, count.triangle.B.at.node(A, membership, v.id, aggrega.mod))
}
rownames(tri.B.res) = c(QUASI.POS.B1.TRI.DESC.NAME, QUASI.POS.B2.TRI.DESC.NAME, QUASI.POS.B3.TRI.DESC.NAME, QUASI.NEG.TRI.B1.DESC.NAME, 
                             QUASI.NEG.TRI.B2.DESC.NAME, QUASI.NEG.TRI.B3.DESC.NAME, POS.TRI.B.DESC.NAME, NEG.TRI.B.DESC.NAME)

  
tri.C.res = c()
for(v.id in 1:vcount(g)){
    tri.C.res = cbind(tri.C.res, count.triangle.C.at.node(A, membership, v.id))
}
rownames(tri.C.res) = c(QUASI.POS.C1.TRI.DESC.NAME, QUASI.POS.C2.TRI.DESC.NAME, QUASI.NEG.TRI.C1.DESC.NAME, QUASI.NEG.TRI.C2.DESC.NAME, POS.TRI.C.DESC.NAME, NEG.TRI.C.DESC.NAME)
  
tri.D.res = c()
for(v.id in 1:vcount(g)){
    tri.D.res = cbind(tri.D.res, count.triangle.D.at.node(A, membership, v.id, aggrega.mod))
}
rownames(tri.D.res) = c(QUASI.POS.D1.TRI.DESC.NAME, QUASI.POS.D2.TRI.DESC.NAME, QUASI.NEG.TRI.D1.DESC.NAME, QUASI.NEG.TRI.D2.DESC.NAME, POS.TRI.D.DESC.NAME, NEG.TRI.D.DESC.NAME)


result = rbind(result, tri.A.res, tri.B.res, tri.C.res, tri.D.res)
colnames(result) = paste0("V",seq(0,(vcount(g)-1)))

print("writing into files")
write.csv(file=paste0("partition-based-result-mem0-",aggrega.mod,".csv"), x=result)
