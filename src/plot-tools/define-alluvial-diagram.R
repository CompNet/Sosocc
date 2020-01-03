
library(alluvial)



###############################################################
#
# n: The number of vertical layer (i.e partition) desired. 
# 	  Not that it migh be lower than 'n' at the end if the best k value is low (i.e. 2 or 3, etc.)
# themes: Themes covered in a given instance (i.e. France, Agri, 12-13).
# 		   Note that themes will differ from other instances (i.e. for another period)
# partitions.df: A data frame object in which each column corresponds to a partition
# x.labels: Labels in x axis, for each vertical layer
#
################################################################
plot.alluvial.diagram = function(out.folder, out.filename, partitions.df, x.labels, plot.title, show.stripes.border=FALSE, format=c("PDF","PNG",NA))
{
	border=NA
	
	if(show.stripes.border){
#		border="black"
		border="white"	
#		out.filename=paste0("bordered-",out.filename)
	}
	filename=file.path(out.folder,out.filename)
	print(filename)
	
	# process each specified format
	for(frmt in format)
	{	# create the file
		if(!is.na(frmt))
		{	# set plot file name
			plot.filename <- filename
			if(toupper(substr(plot.filename, nchar(plot.filename)-2, nchar(plot.filename)))!=toupper(frmt))
				plot.filename <- paste0(plot.filename ,".",frmt)
			# handle format
			if(frmt=="PNG")
			{	
				png(filename=plot.filename,width=800,height=800,units="px",pointsize=20,bg="white")
			}
			if(frmt=="JPEG")
			{	
				jpeg(filename=plot.filename,width=800,height=800,units="px",pointsize=20,bg="white",quality="70")
			}
			else if(frmt=="PDF")
			{	pdf(file=plot.filename,bg="white")
			}
		}
		
		# create the plot
		# no need to use different colors, for now, blue is okay
		alluvial(x=partitions.df, freq=1, col="blue", border=border, axis_labels=x.labels, cex.axis=0.40)
#		print(plot.title)
		mtext(plot.title, 3, line=3, font=1.5)
#		title(plot.title, cex.main=0.5)
		
		
		# finalize plot file
		if(!is.na(frmt))
			dev.off()
	}

}