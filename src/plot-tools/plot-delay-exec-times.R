

# source: https://stackoverflow.com/questions/20924705/plot-negative-values-in-logarithmic-scale-with-ggplot-2

library(ggplot2)
library(scales)
library(ggallin)


palette <- as.list(rainbow(14))
palette[1:4] = rainbow(4, v=0.7)
# ==============
palette[[5]] = "#E3CF57" # banana yellow
palette[[6]] = "maroon1" #
palette[[7]] = "darkviolet" # 
palette[[8]] = "plum1" #
palette[[9]] = "seagreen4" #
palette[[10]] = "black" # 
palette[[11]] = "burlywood4" # 
palette[[12]] = "yellow" #
palette[[13]] = "lightseagreen"
palette[[14]] = "gray87" # light gray
# ==============




###############################################################################
# For each value of the parameters n, l0 and qneg, it plots the difference of execution times between
#   EnumCC(3) and OneTreeCC, EnumCC(3) minus OneTreeCC, represented on the log-scaled y-axis of the plots.
# When such difference takes a negative value, this means EnumCC(3) runs faster than OneTreeCC(). 
#   The set of 20 graphs generated for each parameter set are indexed, in the x-axis. Finally, to guide our discussion, 
#   we show for each graph the maximal number of solutions found by the method(s) within a time limit and the number of jumps 
#   related to EnumCC(3), where the latter is shown in parentheses.
#
###############################################################################
plot.delay.exec.times = function(graph.sizes, d, l0.values, prop.negs, in.rand.net.folders, algo.name, force){
	
	for(n in graph.sizes){
		for(l0 in l0.values){
			for(qneg in prop.negs){
				cat("plotting for n=",n," and qneg=",qneg,"\n")
				
				df.edit = read.csv(file.path(OUTPUT.CSV.FOLDER,paste0("delay-exec-time-",algo.name,"-l0=",l0,"_n=",n,"_d=",d,".csv")), header=1, check.names=F, row.names=1)
				nb.sols = read.csv(file.path(OUTPUT.CSV.FOLDER,paste0("nb-sol-l0=",l0,"_",algo.name,"_n=",n,"_d=",d,".csv")), header=1, check.names=F, row.names=1)
				df.edit[,"nbSols"] = nb.sols
				nb.comp = read.csv(file.path(OUTPUT.CSV.FOLDER,paste0("Edit-Dist-NbComponents-l0=",l0,"_n=",n,"_d=",d,".csv")), header=1, check.names=F, row.names=1)
				df.edit[,"nbComps"] = nb.comp[,algo.name]
				#print(df.edit)
				
				csv.n.values.str = sapply(row.names(df.edit), function(x) unlist(strsplit(x,","))[1])
				csv.n.values = as.integer(sapply(csv.n.values.str, function(x) unlist(strsplit(x,"="))[2]))
				csv.prop.mispl.values.str = sapply(row.names(df.edit), function(x) unlist(strsplit(x,","))[4])
				csv.prop.mispl.values = as.numeric(sapply(csv.prop.mispl.values.str, function(x) unlist(strsplit(x,"="))[2]))
				
				csv.prop.neg.values.str = sapply(row.names(df.edit), function(x) unlist(strsplit(x,","))[5])
				csv.prop.neg.values = as.numeric(sapply(csv.prop.neg.values.str, function(x) unlist(strsplit(x,"="))[2]))
				
				nb.instances.per.param = length(in.rand.net.folders)
				nb.columns = length(which((csv.prop.mispl.values*100)%%2 == 0)) # nb columns
				
				indxs = which(csv.n.values == n & ((csv.prop.mispl.values*100)%%2 == 0)) # even values
				if(!is.na(qneg))
				    indxs = which(csv.n.values == n & ((csv.prop.mispl.values*100)%%2 == 0) & csv.prop.neg.values==qneg)
				names(indxs) = NULL
				data.mtrx = rbind(cbind(df.edit[indxs,1], algo.name,df.edit[indxs,"nbSols"], df.edit[indxs,"nbComps"]))
				#print(data.mtrx)
				df <- data.frame(x=rep(seq(1,nb.instances.per.param),nb.columns), y=as.numeric(data.mtrx[,1]), method=data.mtrx[,2], nbSols=data.mtrx[,3], nbComps=data.mtrx[,4], propMispl=csv.prop.mispl.values[indxs])
				
				ggplot(df, aes(x = x,y = y))+
				     geom_line(aes(color = method, linetype = method), size=1.5)+
				     geom_point(aes(color = method))+ 
				     geom_text(aes(label=ifelse(method==algo.name,paste0(nbSols,"\n(",nbComps,")"),"")),hjust=0, vjust=0, size=3.5)+
				     scale_color_manual(values = c("darkred", "steelblue", "#E3CF57", "palevioletred"))+
				     scale_y_continuous(trans = pseudolog10_trans,breaks=c(-5000,-1000,-100,-10,-1,0,1,10,100,1000,5000))+
				     #scale_x_continuous(sec.axis=sec_axis( name="nbSols"))+
				     geom_hline(yintercept=0, linetype="dashed", color="black") +
				     facet_wrap(~ propMispl, ncol=5)+
				     theme_light(base_size = 21)
				
				 ggsave(filename=file.path(OUTPUT.CSV.FOLDER,paste0("delay-exact-time_",algo.name,"_d=",d,"_n",n,"_l0=",l0,"_qneg=",qneg,".pdf")),
						 width=16)
			 
				#sp = sp + geom_hline(yintercept=0, linetype="dashed", color="black")
				#sp = sp + geom_vline(xintercept=20, linetype="dashed", color="black")
				#sp = sp + geom_vline(xintercept=40, linetype="dashed", color="black")
				#sp = sp + geom_vline(xintercept=60, linetype="dashed", color="black")
			}
		}
	}

}

