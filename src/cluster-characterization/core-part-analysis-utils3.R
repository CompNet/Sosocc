


##############################################
#
#
##############################################
compute.internal.degree <- function(A, membership, v.id) {
    same.clu.ids = which(membership == membership[v.id])
    # suppose that there is no self loop
    return(sum(A[v.id,same.clu.ids])) 
}



##############################################
#
#
##############################################
compute.participation.coeff = function(A, membership, v.id){
    nb.clu = length(unique(membership))
    deg = sum(A[v.id,])

    res = 0
    for(k in 1:nb.clu){
        target.clu.neig.ids = which(membership==k & A[v.id,]!=0)
        ext.deg = sum(A[v.id,target.clu.neig.ids])
        res = res + (ext.deg/deg)^2
    }
    return(1-res)
}




##############################################
#
#
##############################################
compute.z.score = function(A, membership, v.id, FUN){
    #print("-------------")
    #cat("node: ", v.id, "\n")

    same.clu.ids = which(membership == membership[v.id])
    v.score  = FUN(A, membership, v.id)
    #print(v.score)
    clu.scores = sapply(same.clu.ids, function(v) FUN(A, membership, v) )
    #print(clu.scores)
    #print(mean(clu.scores))
    #print(sd(clu.scores))

    if(length(same.clu.ids) == 1)
        return(0)
    if(sd(clu.scores) == 0)
        return(0)
    else{
        z.score = (v.score - mean(clu.scores))/sd(clu.scores) 
        return(z.score) 
    }
}



##############################################
#
#
##############################################

g = read.graph("signed-unweighted.graphml","graphml")
membership = read.table(file="membership0.txt")$V1
A = as_adj(g, attr="weight")
gpos = delete.edges(graph=g, edges=which(E(g)$weight<0))
Apos = as_adj(gpos)
v.id = 1
int.z.scores = sapply(V(gpos), function(v.id) compute.z.score(Apos, membership, v.id, compute.internal.degree) )
part.scores = sapply(V(gpos), function(v.id) compute.participation.coeff(Apos, membership, v.id) )
print(int.z.scores)
print(part.scores)
plot(x=part.scores, y=int.z.scores)
text(x=part.scores, y=int.z.scores, labels=(V(g)-1), cex= 0.7,  pos=3)
