###############################################################
# cNMTF
#	1. Preprocessing functions
# 1.1 Functions to create outcome matrix
#
# Corresponding author:
# Luis Leal, Imperial College London, email: lgl15@imperial.ac.uk
###############################################################


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#' @title Constructs the outcome matrix
#' @description Function to construct \code{Vo} (Low-rank matrix of outcome labels)
#'
# [Input]:
#' @param out Outcome vector for the patients
#' @param ki Number of clusters/columns for the outcome matrix
#'
# [Output]:
#' @return \code{Vo} outcome matrix of size [m x ki]
#'
#' @md
#' @author Luis G. Leal, \email{lgl15@@imperial.ac.uk}
#' @family Preprocessing functions
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


    construct.Vo <- function(out,ki){

      #Find number of levels in the outcome
    	nlevels.out = length(levels(as.factor(out)))
    	if( ki < nlevels.out ){print("Number of cluster must be as higher as the number of levels"); break}
    	if( ki < 2 ){print("Number of clusters lower than 2");  break}

    	Vo = matrix(0,nrow = length(out),ncol=ki)

    	#Assignation of outcome level to clusters
    	out.k = rep(1:nlevels.out,each = round(ki/nlevels.out))[1:ki]

    	#Fill the Vo matrix

    	for(i in 1:nlevels.out){
    		Vo[out == levels(as.factor(out))[i], out.k == i] <- 1
    	}
    	return(Vo)
    }
