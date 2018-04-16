###############################################################
# cNMTF
#	1. Preprocessing functions
# 1.2 Functions to work with ancestry/population origin
#
# Corresponding author:
# Luis Leal, Imperial College London, email: lgl15@imperial.ac.uk
###############################################################


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#' @title Calculates the kernel similarity matrix
#' @description Function to construct kernel matrices of ancestry/population origin
#' @references Li, Rakitisch, et al. "ccSVM: correcting SVMs for confounding factors in biological data classification" ISMB 2011
#'
# [Input]:
#' @param R Data matrix of size [n x m]
#' @param pop Population origin. \code{c("unknown", "known")}
#'
# [Output]:
#' @return Kernel matrices:
#' * \code{H} Centering matrix of size [m x m]
#' * \code{L} Side information kernel matrix (i.e. Kernel on the confounder random variable) [m x m]
#'
#' @md
#' @author Luis G. Leal, \email{lgl15@@imperial.ac.uk}
#' @family Preprocessing functions
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


    kernels.cnmtf = function(R, pop = c("unknown", "known")){

    	#Construct the Side information kernel matrix (L)
    	m = ncol(R)
    	n = nrow(R)
    	##If the population structure is unknown, then its obtained from the PCs
    	if (pop == "unknown"){
    		L = ( t(R) %*% R ) * 1/n

    	}else if (pop == "known"){ #Otherwise the vector of population membership is used to construct L
    		if(length(pop) == m){
    			L = matrix(0, nrow = m, ncol = m)
    			pop <- as.factor(pop)
    			for(i in levels(pop)){
    				L[ pop %in% i, pop %in% i] <- 1
    			}
    		}
    		else{
    			print("Vector of population memebership does not match the number of patients in R")
    			return(NULL)
    		}
    	}

    	#Construct the Centering matrix (H)
    	H = matrix((0 - 1/m), nrow = m, ncol = m)
    	diag(H) <- (1 - 1/m)

    	return(list(H,L))
    }



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#' @title Plot of PCA for patients
#' @description Function to plot axes of PCA and depict population structures
#[Input]
#' @param res.pca Results of PCA
#' @param var.col variable to colour the points (e.g., population)
#' @param axes.pca vector of axes to plot
#' @param plot.type Using own function or default function
#' @param outl Outliers vector to filter out PCA results
#[Output]
#' @return Plot of patients in the principal components space.
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


    plot.pca = function( res.pca,
                         var.col = NULL,
                         var.name = "Variable",
                         axes.pca,
                         show.legend = TRUE,
                         plot.type = c("ggplot", "default"),
                         col.manual = NULL,
                         outl = rep(TRUE, length(var.col)) )
      {

        #Calculate explained variance
        evar = round( res.pca$eig / sum(res.pca$eig) * 100, 1)

        #Labels and colours

        dpca = as.data.frame( cbind(res.pca$li[ outl ,axes.pca], var.col [ outl ]))
        print(dim(dpca))
        print(sum(outl))
        names(dpca) <- c("Axis1", "Axis2")
        xlab = as.character(paste("PC", axes.pca[1], " (EV = ", evar[1], "%)", sep = "" ))
        ylab = paste("PC", axes.pca[2], " (EV = ", evar[2], "%)", sep = "" )

        if(plot.type == "default"){

          cat("Plotting PCA with default function", "\n")

          s.class(res.pca$li[,c(1,2)], fac = as.factor(var.col), cellipse = 1, cstar = 1, col = rainbow(length(levels(as.factor(var.col)))), pch = 16)

        }else if(plot.type == "ggplot"){

          #Plot PCA

          cat("Plotting PCA with ggplot function", "\n")

          col.manual = if( is.null(col.manual) ) { rainbow( nlevels(var.col ) ) }else{ col.manual }

          p <- ggplot(dpca, aes(Axis1, Axis2)) +
            geom_point( aes(colour = var.col[ outl ] ), shape=16, show.legend = show.legend) +
            scale_colour_manual(values = col.manual, name=var.name) +
            geom_vline(xintercept = 0, lty = 2, col = "gray") +
            geom_hline(yintercept = 0, lty = 2, col = "gray") +
            theme_bw(base_size = 16) + labs(x = xlab, y = ylab)

          print(p)

        }

        return(p)

    }


