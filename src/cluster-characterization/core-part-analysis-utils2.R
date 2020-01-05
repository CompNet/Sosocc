


##############################################
# node ids should start from 1, as opposed to .G format. 
# Then, you can convert it (minus 1)
#
##############################################
count.triangle.A.at.node = function(A, membership, v.id, aggrega.mod="sum"){
    #print("-------------")
    #cat("node: ", v.id, "\n")
    nb.clu = length(unique(membership))

    quasi.pos.A = c()
    quasi.neg.A = c()
    pos.A = c()
    neg.A = c()
    for(k in 1:nb.clu){
        #cat("cluster: ",k,"\n")
        if(k != membership[v.id]){ # diff than its own cluster
            target.clu.neig.ids = which(membership==k & A[v.id,]!=0)

            if(length(target.clu.neig.ids)>1){
                #print("quasi.pos.A")
                curr.quasi.pos.A = count.quasi.pos.A.and.C(A, membership, v.id, target.clu.neig.ids)
                quasi.pos.A = rbind(quasi.pos.A, curr.quasi.pos.A)
                #print("quasi.neg.A")
                curr.quasi.neg.A = count.quasi.neg.A.and.C(A, membership, v.id, target.clu.neig.ids)
                quasi.neg.A = rbind(quasi.neg.A, curr.quasi.neg.A)
                #print("pos.A")
                curr.pos.A = count.pos.A.and.C(A, membership, v.id, target.clu.neig.ids)
                pos.A = rbind(pos.A, curr.pos.A)
                #print("neg.A")
                curr.neg.A = count.neg.A.and.C(A, membership, v.id, target.clu.neig.ids)
                neg.A = rbind(neg.A, curr.neg.A)
            }
        }
    }

    # apply(X,2,sum) # sum by column
    if(aggrega.mod == "sum")
        return( c(apply(quasi.pos.A,2,sum),apply(quasi.neg.A,2,sum),apply(pos.A,2,sum),apply(neg.A,2,sum)) )
    else if(aggrega.mod == "mean")
        return( c(round(apply(quasi.pos.A,2,mean),digits=2),round(apply(quasi.neg.A,2,mean),digits=2),round(apply(pos.A,2,mean),digits=2),round(apply(neg.A,2,mean),digits=2)) )
    else if(aggrega.mod == "sd")
        return( c(round(apply(quasi.pos.A,2,sd),digits=2),round(apply(quasi.neg.A,2,sd),digits=2),round(apply(pos.A,2,sd),digits=2),round(apply(neg.A,2,sd),digits=2)) )
}



##############################################
#
#
##############################################
count.quasi.pos.A.and.C = function(A, membership, v.id, neig.ids){
    nb.neig = length(neig.ids)
    count.A1 = 0
    count.A2 = 0
    for(i in 1:(nb.neig-1)){
        neig.id1 = neig.ids[i]

            for(j in (i+1):nb.neig){
                 neig.id2 = neig.ids[j]
                 if(A[v.id,neig.id1]<0 && A[v.id,neig.id2]>0 && A[neig.id1,neig.id2]>0)
                    count.A1 = count.A1 + 1
                 else if(A[v.id,neig.id1]>0 && A[v.id,neig.id2]<0 && A[neig.id1,neig.id2]>0)
                    count.A1 = count.A1 + 1
                 else if(A[v.id,neig.id1]>0 && A[v.id,neig.id2]>0 && A[neig.id1,neig.id2]<0)
                    count.A2 = count.A2 + 1
            }

    }
    return(c(count.A1,count.A2))
}


##############################################
#
#
##############################################
count.quasi.neg.A.and.C = function(A, membership, v.id, neig.ids){
    nb.neig = length(neig.ids)
    count.A1 = 0
    count.A2 = 0
    for(i in 1:(nb.neig-1)){
        neig.id1 = neig.ids[i]

        for(j in (i+1):nb.neig){
             neig.id2 = neig.ids[j]
             if(A[v.id,neig.id1]>0 && A[v.id,neig.id2]<0 && A[neig.id1,neig.id2]<0)
                count.A1 = count.A1 + 1
             else if(A[v.id,neig.id1]<0 && A[v.id,neig.id2]>0 && A[neig.id1,neig.id2]<0)
                count.A1 = count.A1 + 1
             else if(A[v.id,neig.id1]<0 && A[v.id,neig.id2]<0 && A[neig.id1,neig.id2]>0)
                count.A2 = count.A2 + 1
        }

    } 

    return(c(count.A1,count.A2))
}


##############################################
#
#
##############################################
count.neg.A.and.C = function(A, membership, v.id, neig.ids){
    nb.neig = length(neig.ids)
    count.A1 = 0
    for(i in 1:(nb.neig-1)){
        neig.id1 = neig.ids[i]

        if(A[v.id,neig.id1]<0){
            for(j in (i+1):nb.neig){
                 neig.id2 = neig.ids[j]
                 if(A[v.id,neig.id2]<0 && A[neig.id1,neig.id2]<0)
                    count.A1 = count.A1 + 1
            }

        }

    }
    return(c(count.A1))
}


##############################################
#
#
##############################################
count.pos.A.and.C = function(A, membership, v.id, neig.ids){
    nb.neig = length(neig.ids)
    count.A1 = 0
    for(i in 1:(nb.neig-1)){
        neig.id1 = neig.ids[i]

        if(A[v.id,neig.id1]>0){
            for(j in (i+1):nb.neig){
                 neig.id2 = neig.ids[j]
                 if(A[v.id,neig.id2]>0 && A[neig.id1,neig.id2]>0)
                    count.A1 = count.A1 + 1
            }

        }

    }
    return(c(count.A1))
}



# ========================================================================



##############################################
# node ids should start from 1, as opposed to .G format. 
# Then, you can convert it (minus 1)
#
##############################################
count.triangle.B.at.node = function(A, membership, v.id, aggrega.mod="sum"){
    # cat("node: ", v.id, "\n")
    nb.clu = length(unique(membership))
    source.clu.neig.ids = which(membership[v.id] == membership)
    if(length(source.clu.neig.ids)==1)
        return(c(0,0,0, 0,0,0, 0, 0))


    quasi.pos.B = c()
    quasi.neg.B = c()
    pos.B = c()
    neg.B = c()
    for(k in 1:nb.clu){
        #cat("cluster: ",k,"\n")
        
        if(k!=membership[v.id]){ # k is diff than its own cluster
            target.clu.neig.ids = which(membership==k & A[v.id,]!=0)

            if(length(target.clu.neig.ids)>1){
                #print("quasi.pos.B")
                curr.quasi.pos.B = count.quasi.pos.B(A, membership, v.id, source.clu.neig.ids, target.clu.neig.ids)
                quasi.pos.B = rbind(quasi.pos.B, curr.quasi.pos.B)
                #print("quasi.neg.B")
                curr.quasi.neg.B = count.quasi.neg.B(A, membership, v.id, source.clu.neig.ids, target.clu.neig.ids)
                quasi.neg.B = rbind(quasi.neg.B, curr.quasi.neg.B)
                #print("pos.B")
                curr.pos.B = count.pos.B(A, membership, v.id, source.clu.neig.ids, target.clu.neig.ids)
                pos.B = rbind(pos.B, curr.pos.B)
                #print("neg.B")
                curr.neg.B = count.neg.B(A, membership, v.id, source.clu.neig.ids, target.clu.neig.ids)
                neg.B = rbind(neg.B, curr.neg.B)
            }
        }
    }

    # apply(X,2,sum) # sum by column
    if(aggrega.mod == "sum")
        return( c(apply(quasi.pos.B,2,sum),apply(quasi.neg.B,2,sum),apply(pos.B,2,sum),apply(neg.B,2,sum)) )
    else if(aggrega.mod == "mean")
        return( c(round(apply(quasi.pos.B,2,mean),digits=2),round(apply(quasi.neg.B,2,mean),digits=2),round(apply(pos.B,2,mean),digits=2),round(apply(neg.B,2,mean),digits=2)) )
    else if(aggrega.mod == "sd")
        return( c(round(apply(quasi.pos.B,2,sd),digits=2),round(apply(quasi.neg.B,2,sd),digits=2),round(apply(pos.B,2,sd),digits=2),round(apply(neg.B,2,sd),digits=2)) )
}


##############################################
#
#
##############################################
count.quasi.pos.B = function(A, membership, v.id, source.neig.ids, target.neig.ids){
    count.B1 = 0
    count.B2 = 0
    count.B3 = 0
    for(neig.id1 in source.neig.ids){
        for(neig.id2 in target.neig.ids){
             if(A[v.id,neig.id1]>0 && A[v.id,neig.id2]<0 && A[neig.id1,neig.id2]>0)
                count.B1 = count.B1 + 1
            else if(A[v.id,neig.id1]>0 && A[v.id,neig.id2]>0 && A[neig.id1,neig.id2]<0)
                count.B2 = count.B2 + 1
            else if(A[v.id,neig.id1]<0 && A[v.id,neig.id2]>0 && A[neig.id1,neig.id2]>0)
                count.B3 = count.B3 + 1
        }
    }
    return(c(count.B1,count.B2,count.B3))
}


##############################################
#
#
##############################################
count.quasi.neg.B = function(A, membership, v.id, source.neig.ids, target.neig.ids){
    count.B1 = 0
    count.B2 = 0
    count.B3 = 0
    for(neig.id1 in source.neig.ids){
        for(neig.id2 in target.neig.ids){
             if(A[v.id,neig.id1]<0 && A[v.id,neig.id2]>0 && A[neig.id1,neig.id2]<0)
                count.B1 = count.B1 + 1
            else if(A[v.id,neig.id1]<0 && A[v.id,neig.id2]<0 && A[neig.id1,neig.id2]>0)
                count.B2 = count.B2 + 1
            else if(A[v.id,neig.id1]>0 && A[v.id,neig.id2]<0 && A[neig.id1,neig.id2]<0)
                count.B3 = count.B3 + 1
        }
    }
    return(c(count.B1,count.B2,count.B3))
}


##############################################
#
#
##############################################
count.pos.B = function(A, membership, v.id, source.neig.ids, target.neig.ids){
    count.B1 = 0
    for(neig.id1 in source.neig.ids){
        for(neig.id2 in target.neig.ids){
             if(A[v.id,neig.id1]>0 && A[v.id,neig.id2]>0 && A[neig.id1,neig.id2]>0)
                count.B1 = count.B1 + 1
        }
    }
    return(c(count.B1))
}


##############################################
#
#
##############################################
count.neg.B = function(A, membership, v.id, source.neig.ids, target.neig.ids){
    count.B1 = 0
    for(neig.id1 in source.neig.ids){
        for(neig.id2 in target.neig.ids){
             if(A[v.id,neig.id1]<0 && A[v.id,neig.id2]<0 && A[neig.id1,neig.id2]<0)
                count.B1 = count.B1 + 1
        }
    }
    return(c(count.B1))
}


# ========================================================================



##############################################
# node ids should start from 1, as opposed to .G format. 
# Then, you can convert it (minus 1)
#
##############################################
count.triangle.C.at.node = function(A, membership, v.id){
    # cat("node: ", v.id, "\n")
    nb.clu = length(unique(membership))
    source.clu.neig.ids = which(membership[v.id] == membership)
    v.indx = which(source.clu.neig.ids == v.id)
    source.clu.neig.ids = source.clu.neig.ids[-v.indx]
    if(length(source.clu.neig.ids)<2)
        return(c(0,0, 0,0, 0, 0))

    #print("quasi.pos.C")
    quasi.pos.C = count.quasi.pos.A.and.C(A, membership, v.id, source.clu.neig.ids)
    #print("quasi.neg.C")
    quasi.neg.C = count.quasi.neg.A.and.C(A, membership, v.id, source.clu.neig.ids)
    #print("pos.C")
    pos.C = count.pos.A.and.C(A, membership, v.id, source.clu.neig.ids)
    #print("neg.C")
    neg.C = count.neg.A.and.C(A, membership, v.id, source.clu.neig.ids)
       
    return( c(quasi.pos.C,quasi.neg.C,pos.C,neg.C) )
}




# ========================================================================



##############################################
# node ids should start from 1, as opposed to .G format. 
# Then, you can convert it (minus 1)
#
##############################################
count.triangle.D.at.node = function(A, membership, v.id, aggrega.mod="sum"){
    # cat("node: ", v.id, "\n")
    nb.clu = length(unique(membership))
    if(nb.clu<3)
        return(c(0,0, 0,0, 0, 0))


    quasi.pos.D = c()
    quasi.neg.D = c()
    pos.D = c()
    neg.D = c()
    for(k1 in 1:(nb.clu-1)){
        #cat("cluster: ",k1,"\n")
        target1.clu.neig.ids = which(membership==k1 & A[v.id,]!=0)

        for(k2 in (k1+1):nb.clu){
            #cat("cluster: ",k2,"\n")
            target2.clu.neig.ids = which(membership==k2 & A[v.id,]!=0)
           
             if(k1!=membership[v.id] && k2!=membership[v.id]){ # k is diff than its own cluster                
                #print("quasi.pos.D")
                curr.quasi.pos.D = count.quasi.pos.D(A, membership, v.id, target1.clu.neig.ids, target2.clu.neig.ids)
                quasi.pos.D = rbind(quasi.pos.D, curr.quasi.pos.D)
                #print("quasi.neg.D")
                curr.quasi.neg.D = count.quasi.neg.D(A, membership, v.id, target1.clu.neig.ids, target2.clu.neig.ids)
                quasi.neg.D = rbind(quasi.neg.D, curr.quasi.neg.D)
                #print("pos.D")
                curr.pos.D = count.pos.D(A, membership, v.id, target1.clu.neig.ids, target2.clu.neig.ids)
                pos.D = rbind(pos.D, curr.pos.D)
                #print("neg.D")
                curr.neg.D = count.neg.D(A, membership, v.id, target1.clu.neig.ids, target2.clu.neig.ids)
                neg.D = rbind(neg.D, curr.neg.D)
            }
        }
    }

    # apply(X,2,sum) # sum by column
    if(aggrega.mod == "sum")
        return( c(apply(quasi.pos.D,2,sum),apply(quasi.neg.D,2,sum),apply(pos.D,2,sum),apply(neg.D,2,sum)) )
    else if(aggrega.mod == "mean")
        return( c(round(apply(quasi.pos.D,2,mean),digits=2),round(apply(quasi.neg.D,2,mean),digits=2),round(apply(pos.D,2,mean),digits=2),round(apply(neg.D,2,mean),digits=2)) )
    else if(aggrega.mod == "sd")
        return( c(round(apply(quasi.pos.D,2,sd),digits=2),round(apply(quasi.neg.D,2,sd),digits=2),round(apply(pos.D,2,sd),digits=2),round(apply(neg.D,2,sd),digits=2)) )
}


##############################################
#
#
##############################################
count.quasi.pos.D = function(A, membership, v.id, target1.neig.ids, target2.neig.ids){
    count.D1 = 0
    count.D2 = 0
    for(neig.id1 in target1.neig.ids){
        for(neig.id2 in target2.neig.ids){
            if(A[v.id,neig.id1]<0 && A[v.id,neig.id2]>0 && A[neig.id1,neig.id2]>0)
                count.D1 = count.D1 + 1
            else if(A[v.id,neig.id1]>0 && A[v.id,neig.id2]<0 && A[neig.id1,neig.id2]>0)
                count.D1 = count.D1 + 1
            else if(A[v.id,neig.id1]>0 && A[v.id,neig.id2]>0 && A[neig.id1,neig.id2]<0)
                count.D2 = count.D2 + 1
        }
    }
    return(c(count.D1,count.D2))
}


##############################################
#
#
##############################################
count.quasi.neg.D = function(A, membership, v.id, target1.neig.ids, target2.neig.ids){
    count.D1 = 0
    count.D2 = 0
    for(neig.id1 in target1.neig.ids){
        for(neig.id2 in target2.neig.ids){
            if(A[v.id,neig.id1]>0 && A[v.id,neig.id2]<0 && A[neig.id1,neig.id2]<0)
                count.D1 = count.D1 + 1
            else if(A[v.id,neig.id1]<0 && A[v.id,neig.id2]>0 && A[neig.id1,neig.id2]<0)
                count.D1 = count.D1 + 1
            else if(A[v.id,neig.id1]<0 && A[v.id,neig.id2]<0 && A[neig.id1,neig.id2]>0)
                count.D2 = count.D2 + 1
        }
    }
    return(c(count.D1,count.D2))
}


##############################################
#
#
##############################################
count.pos.D = function(A, membership, v.id, target1.neig.ids, target2.neig.ids){
    count.D1 = 0
    for(neig.id1 in target1.neig.ids){
        for(neig.id2 in target2.neig.ids){
             if(A[v.id,neig.id1]>0 && A[v.id,neig.id2]>0 && A[neig.id1,neig.id2]>0)
                count.D1 = count.D1 + 1
        }
    }
    return(c(count.D1))
}


##############################################
#
#
##############################################
count.neg.D = function(A, membership, v.id, target1.neig.ids, target2.neig.ids){
    count.D1 = 0
    for(neig.id1 in target1.neig.ids){
        for(neig.id2 in target2.neig.ids){
             if(A[v.id,neig.id1]<0 && A[v.id,neig.id2]<0 && A[neig.id1,neig.id2]<0)
                count.D1 = count.D1 + 1
        }
    }
    return(c(count.D1))
}
