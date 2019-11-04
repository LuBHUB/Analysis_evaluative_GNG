## http://wiki.stdout.org/rcookbook/Graphs/Plotting%20means%20and%20error%20bars%20(ggplot2)/#Helper functions

## Summarizes data, handling within-subjects variables by removing inter-subject variability.
## It will still work if there are no within-S variables.
## Gives count, mean, standard deviation, standard error of the mean, and confidence interval (default 95%).
## If there are within-subject variables, calculate adjusted values using method from Morey (2008).
##   data: a data frame.
##   measurevar: the name of a column that contains the variable to be summarized
##   betweenvars: a vector containing names of columns that are between-subjects variables
##   withinvars: a vector containing names of columns that are within-subjects variables
##   idvar: the name of a column that identifies each subject (or matched subjects)
##   na.rm: a boolean that indicates whether to ignore NA's
##   conf.interval: the percent range of the confidence interval (default is 95%)
summarySEwithin <- function(data=NULL, measurevar, betweenvars=NULL, withinvars=NULL,
				    idvar=NULL, na.rm=FALSE, conf.interval=.95, .drop=TRUE) {
	
	# Ensure that the betweenvars and withinvars are factors
	factorvars <- sapply(data[, c(betweenvars, withinvars), drop=FALSE], FUN=is.factor)
	if (!all(factorvars)) {
		nonfactorvars <- names(factorvars)[!factorvars]
		message("Automatically converting the following non-factors to factors: ",
			  paste(nonfactorvars, collapse = ", "))
		data[nonfactorvars] <- lapply(data[nonfactorvars], factor)
	}
	
	# Norm each subject's data    
	data <- normDataWithin(data, idvar, measurevar, betweenvars, na.rm, .drop=.drop)
	
	# This is the name of the new column
	measureNormedVar <- paste(measurevar, "Normed", sep="")
	
	# Replace the original data column with the normed one
	data[,measurevar] <- data[,measureNormedVar]
	
	# Collapse the normed data - now we can treat between and within vars the same
	datac <- summarySE(data, measurevar, groupvars=c(betweenvars, withinvars), na.rm=na.rm,
				 conf.interval=conf.interval, .drop=.drop)
	
	# Apply correction from Morey (2008) to the standard error and confidence interval
	#  Get the product of the number of conditions of within-S variables
	nWithinGroups    <- prod(sapply(datac[,withinvars, drop=FALSE], FUN=nlevels))
	correctionFactor <- sqrt( nWithinGroups / (nWithinGroups-1) )
	
	# Apply the correction factor
	datac$sd <- datac$sd * correctionFactor
	datac$se <- datac$se * correctionFactor
	datac$ci <- datac$ci * correctionFactor
	
	return(datac)
}