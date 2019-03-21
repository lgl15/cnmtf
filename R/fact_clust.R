###############################################################
# cNMTF
#	2. Factorisation functions
# 2.2 Function to perform factorisations
#
# Corresponding author:
# Luis Leal, Imperial College London, email: lgl15@imperial.ac.uk
###############################################################


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#' @title cNMTF algorithm
#' @description Function to perform one repetition of cNMTF
#' @references
#' * Shang, Fanha, et al. "Graph dual regularization non-negative matrix factorization for co-clustering" Pattern Recognition 45 (2012)
#' * Li, Rakitisch, et al. "ccSVM: correcting SVMs for confounding factors in biological data classification" ISMB 2011
#'
# [Input]
#' @param R Data matrix, \emph{n x m}
#' @param HLH Term in the correction for population structures. Dot product: H . L . H where: H is the Centering matrix \emph{(m x m)} and L is Ancestry matrix (i.e. Kernel on the confounder random variable) \emph{(m x m)}
#' @param Vo Phenotype/Outcome label matrix, \emph{m x k2}
#' @param Wu Adjacency matrix of the SNV-SNV network. \emph{n x n}
#' @param DWD Term in the normalised graph laplacian. Dot product:  1/sqrt(Du) . Wu . 1/sqrt(Du), where: Du is the Degree matrix and Wu is the adjacency matrix of the SNV-SNV network.
#' @param k Vector of rank parameters, \emph{k1 x k2}
#' @param iters Default number of itersatiors
#' @param calcObj Check convergency each X number of itersations
#' @param displ Logical. Print number of iterations
#' @param tof RelativeError
#' @param lparameters Vector of regularization parameters:  \eqn{ \gamma_{1}, \gamma_{2}, \gamma_{3} } for the SNV-SNV network, the Phenotype matrix and the Ancestry Kernel respectively.
#' @param init Initialize the matrices randomly (0) or by using PSVD (1)
#'
# [Output]
#' @return
#'   * \code{U} low-dimensional cluster indicator matrix for features observed data, \emph{n x k1}
#'   * \code{V} low-dimensional cluster indicator matrix for observed data, \emph{m x k2}
#'   * \code{S} cluster mapping matrix, \emph{k1 x k2}
#'   * \code{num.iters} number of itersation till convergence
#'   * \code{final.J} Final objective function value
#'   * \code{objF.vector} Objective function values across iterations
#'
#' @md
#' @author Luis G. Leal, \email{lgl15@@imperial.ac.uk}
#' @family Factorisation functions
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


    cnmtf = function(R, # Genotype data matrix [nxm]
    								 HLH, # HLH term in the correction for population structures. Product of H %*% L %*% H where: H is the Centering matrix [mxm] and L is the side information kernel matrix (i.e. Kernel on the confounder random variable) [mxm]
    								 Vo, # Phenotype/Outcome label matrix [mxk2]
    								 Wu, # Adjacency matrix of the SNV-SNV network.
    								 k, # [k1,k2] - Rank parameters
    								 iters = 1, # Default number of iteratiors
    								 calcObj, calcObj2, # Check convergency each X number of iterations after first Y iterations
    								 displ = TRUE, # Print number of iterations
    								 tof = 1e-5, # Minimum relativeError
    								 lparameters, # List of parameters
    								 init = 1, # Initialize the matrices randomly (0) or by using PSVD (1)
    								 V.init = NULL, U.init = NULL, # Initialisation of U and V using PSVD (1)
    								 init.S = FALSE, # Initialisation of S using my method
    								 DWD = NULL # DWD term in the normalised graph laplacian. Product of 1/sqrt(Du) %*% Wu %*% 1/sqrt(Du), where: Du is the Degree matrix and Wu is the adjacency matrix of the SNV-SNV network.
    )
    {


    	# Define dimmensiion
    	m = ncol(R)
    	n = nrow(R)
    	displ = F

    	# Replace zeros with small values in the relationship matrix to speed up process
    	R[ R == 0 ] = 2^(-16)


    	# Initializing low-rank matrices U,V and S
    	if (displ == TRUE){ cat("Initializing low-rank matrices U,V and S", "\n") }

    	if( init == 1 ){# PSDV initialization for U and V

    		S = matrix(runif(k[1] * k[2]), nrow=k[1], ncol=k[2])

    		if( is.null(V.init) | is.null(U.init) ){
    			V = psdv.init( t(R), k[2])
    			U = psdv.init( R, k[1])

    		}else{
    			V = V.init
    			U = U.init
    		}

    	} else { # Random initialization

    		U = matrix( runif( n * k[1] ), nrow = n, ncol = k[1])
    		V = matrix( runif( m * k[2] ), nrow = m, ncol = k[2])
    		S = matrix( runif( k[1] * k[2] ), nrow = k[1], ncol = k[2])

    	}


    	# Replace zeros with small values in the initialised matrices to speed up process
    	S[ S == 0 ] = 2^(-16)
    	V[ V == 0 ] = 2^(-16)
    	U[ U == 0 ] = 2^(-16)

    	# Extract parameters and asign them to variables in this environment
    	for(i in 1:3){

    		assign( names(lparameters)[i] , as.numeric( unlist(lparameters[i]) ) )

    	}

    	# Assess diagonal degree matrix and Laplacian matrix

    	if(gamma1 != 0){
    		if (displ == TRUE){ cat("Calculating Laplacian of Wu","\n") }
    		Du = diag(rowSums(Wu))

    		if ( !is.null(DWD) ){ # Normalised Laplacian

    			I = diag(rep(1, nrow(Wu))) # Identity matrix
    			Lu = I - DWD


    		}else{ # Basic Laplacian for the SNV-SNV network

    			Lu = Du - Wu

    		}

    	}


    	# Initialize objective function
    	objF.previous = Inf

    	# Display iteration message
    	if (displ == TRUE){ cat("Iterating","\n") }


    	# Number of times the function has increased in a row
    	# This is used when the function does not minimises to find an optimal minimum
    	count.increase = 0


    	# Population stratification is corrected till 'iters.gamma3' iterations. Then 'mode gamma2' is activated
    	iters.gamma3 = 50
    	if(gamma3 > 0){

    		mode.gamma2 = FALSE
    		if (displ == TRUE){ cat("Mode gamma2", mode.gamma2,  "\n"  ) }

    	}else{
    		mode.gamma2 = TRUE
    	}

    	# Run iterations
    	for (i in 1:iters){ #For each iteration


    		# Precalculate matrices to speed up process
    		#--------------------

    		US = U %*% S
    		tVV = t(V) %*% V
    		UStVV= US %*% tVV
    		RV = R %*% V

    		# Update rule for S
    		#--------------------

    		# Numerator of update rule
    		up_s = t(U) %*% RV

    		# Denominator of update rule
    		down_s = t(U) %*% UStVV

    		# Update S
    		S = S * (up_s/down_s)

    		# Update rule for U
    		#--------------------

    		# Numerator of update rule
    		up_u = RV %*% t(S)

    		# Denominator of update rule
    		down_u = UStVV %*% t(S)

    		# The denominator of this update rule also depends on the normalised/basic graph Laplacian option (DWD)
    		if ( gamma1 != 0 & ( i %% calcObj) == 0) {  # Updated every X iterations to speed up process

    			if ( !is.null(DWD) ){ # Normalised graph Laplacian

    				up_u = up_u + gamma1 * DWD %*% U
    				down_u = down_u + gamma1 * I %*% U


    			}else{ # Basic graph Laplacian

    				up_u = up_u + gamma1 * Wu %*% U
    				down_u = down_u + gamma1 * Du %*% U

    			}

    		}

    		# Update U
    		U = U * (up_u/down_u)


    		# Update rule for V
    		#--------------------

    		# Numerator of update rule
    		up_v = t(R) %*% US

    		# Denominator of update rule
    		down_v = V %*% t(S) %*% t(U) %*% US

    		# The denominator of this update rule also depends on the correction for population stratification (HLH)
    		# And the phenotypes Vo

    		if ( gamma3 != 0 & mode.gamma2 == FALSE){
    			down_v = down_v + gamma3 * HLH %*% V
    		}
    		if ( gamma2 != 0 & mode.gamma2 == TRUE){
    			up_v = up_v + gamma2 * Vo
    			down_v = down_v + gamma2 * V
    		}

    		# Update V
    		V = V * (up_v/down_v)




    		# Compute residual and print relative error and objective function

    		if( (( i %% calcObj) == 0 & (i >= calcObj2)) || i == iters ){ #If the iteration is a multiple of of calcObj or if last itersation was reached

    			# Compute terms of objective function

    			tmain = norm(R - U %*% S %*% t(V),"F")^2

    			if( gamma1 != 0 ){ tgamma1 = gamma1 * sum( diag(t(U) %*% Lu %*% U) ) } else { tgamma1 = 0 }
    			if( gamma3 != 0 & mode.gamma2 == FALSE){ tgamma3 = gamma3 * sum( diag(V %*% t(V) %*% HLH) ) } else { tgamma3 = 0 }
    			if( gamma2 != 0 & mode.gamma2 == TRUE){ tgamma2 = gamma2 * norm( V - Vo,"F")^2 } else { tgamma2 = 0 }


    			# Sum terms of objective function
    			objF.current = tmain + tgamma3 + tgamma2 + tgamma1

    			# Computing relative change of the objective function
    			relErr = (objF.previous - objF.current) / max(1,objF.previous) / calcObj
    			if(is.na(relErr)){relErr = Inf}


    			# Print relative error and objective function

    			if (displ == TRUE){ cat(paste("Iteration: ",i, "; Relative error: ", format(relErr,scientific = TRUE, digits = 3), "; Objective function: ", format(objF.current, scientific = TRUE, digits = 3), "\n", sep = "")) }


    			# If gamma3 or gamma2 are penalizing, check the number of iterations and change the mode from gamma3 to gamma2

    			if( (i >= iters.gamma3  | tof >= relErr) & mode.gamma2 == FALSE ){

    				mode.gamma2 = TRUE

    				if (displ == TRUE){ cat("Mode gamma2", mode.gamma2,  "\n"  ) }

    				objF.previous = Inf
    				next

    			}


    			# If the function have skipped convergency less than 2 times and its value is Inf
    			if( objF.current == Inf & count.increase < 2 ){

    				if (displ == TRUE){ cat("Objective function is Inf. Skipping convergency", count.increase, "time (s).", "\n") }
    				count.increase = count.increase + 1
    				next

    			}else if( objF.current == Inf & count.increase == 2 ){

    				if (displ == TRUE){ print(paste('Objective function is increasing from ', objF.previous, " to ", objF.current, sep = "")) }
    				break

    			}


    			# Check if objective function increases
    			if(objF.current > objF.previous){

    				# Number of times the function has increased in a row
    				count.increase = count.increase + 1

    				if (displ == TRUE){ print(paste('Objective function is increasing from ', objF.previous, " to ", objF.current, sep = "")) }


    				# If it has been increasing during the last two checks then stops
    				if( count.increase == 2){

    					if (displ == TRUE){ print(paste('Objective function is increasing from ', objF.previous, " to ", objF.current, sep = "")) }
    					break #stop (changed to break to return something, even unfinished optimizations)

    				}

    				if (displ == TRUE){ cat("Skiping the checking of convergency", count.increase, "time.", "\n") }

    				next #Skips the checking of convergency

    			}


    			# Update objective function
    			objF.previous <- objF.current

    			# Check convergency
    			if (tof >= relErr || i == iters){

    				if (displ == TRUE){ print(paste('Local minima found')) }

    				break
    			}
    		}


    	} # End iterations


    	# Conform cnmtf results
    	return(list(U, V, S, i, objF.current, relErr))
    }



