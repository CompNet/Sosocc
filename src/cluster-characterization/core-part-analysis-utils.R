


##############################################
# node ids should start from 1, as opposed to .G format. 
# Then, you can convert it (minus 1)
#
##############################################
count.neg.triangle.at.node = function(A, v.id){
    # cat("node: ", v.id, "\n")
    neig.ids = which(A[v.id,] != 0)
    nb.neig = length(neig.ids)
    count = 0
    for(i in 1:(nb.neig-1)){
        neig.id1 = neig.ids[i]

        if(A[v.id,neig.id1]<0){
            for(j in (i+1):nb.neig){
                 neig.id2 = neig.ids[j]
                 if(A[v.id,neig.id2]<0 && A[neig.id1,neig.id2]<0)
                    count = count + 1
            }
        }
    }

    return(count)
}



##############################################
#
#
##############################################
count.quasi.neg.triangle.at.node = function(A, v.id){
    # cat("node: ", v.id, "\n")
    neig.ids = which(A[v.id,] != 0)
    nb.neig = length(neig.ids)
    count = 0
    for(i in 1:(nb.neig-1)){
        neig.id1 = neig.ids[i]

        if(A[v.id,neig.id1]<0){
            for(j in (i+1):nb.neig){
                 neig.id2 = neig.ids[j]
                 if(A[v.id,neig.id2]<0 && A[neig.id1,neig.id2]>0)
                    count = count + 1
            }
        }
    }

    return(count)
}



##############################################
#
#
##############################################
count.potential.neg.triangle.at.node = function(A, v.id){
    # cat("node: ", v.id, "\n")
    neig.ids = which(A[v.id,] != 0)
    nb.neig = length(neig.ids)
    count = 0
    for(i in 1:(nb.neig-1)){
        neig.id1 = neig.ids[i]

        for(j in (i+1):nb.neig){
             neig.id2 = neig.ids[j]
             if(A[v.id,neig.id1]<0 && A[v.id,neig.id2]>0 && A[neig.id1,neig.id2]<0)
                count = count + 1
            else if(A[v.id,neig.id1]>0 && A[v.id,neig.id2]<0 && A[neig.id1,neig.id2]<0)
                count = count + 1
        }
    }

    return(count)
}



##############################################
#
#
##############################################
count.pos.triangle.at.node = function(A, v.id){
    # cat("node: ", v.id, "\n")
    neig.ids = which(A[v.id,] != 0)
    nb.neig = length(neig.ids)
    count = 0
    for(i in 1:(nb.neig-1)){
        neig.id1 = neig.ids[i]

        if(A[v.id,neig.id1]>0){
            for(j in (i+1):nb.neig){
                 neig.id2 = neig.ids[j]
                 if(A[v.id,neig.id2]>0 && A[neig.id1,neig.id2]>0)
                    count = count + 1
            }
        }
    }

    return(count)
}


##############################################
#
#
##############################################
count.quasi.pos.triangle.at.node = function(A, v.id){
    # cat("node: ", v.id, "\n")
    neig.ids = which(A[v.id,] != 0)
    nb.neig = length(neig.ids)
    count = 0
    for(i in 1:(nb.neig-1)){
        neig.id1 = neig.ids[i]

        if(A[v.id,neig.id1]>0){
            for(j in (i+1):nb.neig){
                 neig.id2 = neig.ids[j]
                 if(A[v.id,neig.id2]>0 && A[neig.id1,neig.id2]<0)
                    count = count + 1
            }
        }
    }

    return(count)
}


##############################################
#
#
##############################################
count.potential.pos.triangle.at.node = function(A, v.id){
    # cat("node: ", v.id, "\n")
    neig.ids = which(A[v.id,] != 0)
    nb.neig = length(neig.ids)
    count = 0
    for(i in 1:(nb.neig-1)){
        neig.id1 = neig.ids[i]

        for(j in (i+1):nb.neig){
             neig.id2 = neig.ids[j]
             if(A[v.id,neig.id1]>0 && A[v.id,neig.id2]<0 && A[neig.id1,neig.id2]>0)
                count = count + 1
            else if(A[v.id,neig.id1]<0 && A[v.id,neig.id2]>0 && A[neig.id1,neig.id2]>0)
                count = count + 1
        }
    }

    return(count)
}

