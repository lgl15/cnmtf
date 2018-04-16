###############################################################
# cNMTF
#	2. Factorisation functions
# 2.3 Functions to find consensus cNMTF results
#
# Corresponding author:
# Luis Leal, Imperial College London, email: lgl15@imperial.ac.uk
###############################################################


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#' @title Find consensus clustering and SNV score
#' @description Function for running the cNMTF algorithm a number of times and extract consensus results
#'
# [Input]
#' @param R Relationship matrix
#' @param Wv Prior knowledge networks for SNVs
#' @param HLH Matrices to correct for confounders
#' @param k Number of clusters
#' @param lparameters Penalization parameters
#' @param iters Number of iterations
#' @param labelsU Labels of SNVs
#' @param labelsV Labels of Patients
#' @param run.t Number of repetitions
#' @param calcObj Number of iterations to check convergency
#' @param calcObj2 Number of iterations to check convergency after \code{calcObj} iterations
#' @param ini Type of seeding/initialisation of matrices in the algorithm
#' @param parallel.opt Run the algorithm in parallel
#' @param n.cores Number of cores to use in the parallel processing
#' @param set.alg Set of algorithm results can be provided to find consensus results
#' @param set.clus.V Set of clustering of patients provided to find consensus results
#' @param export.res Whether the object with list of results is returned
#' @param V.init Initialisation of V using PSVD
#' @param U.init Initialisation of U using PSVD
#' @param do.U Logical. Perform consensus clustering on U matrix
#' @param do.V Logical. Perform consensus clustering on V matrix
#' @param do.O Logical. Perform consensus clustering on Omega matrix
#' @param find.parameters Logical. Run the NMTF algorithm to find estimations of parameters
#' @param estimate.penalisation Logical. Return the penalisation terms evaluated in the last iteration.
#' @param display.iters Logical. Display the iterations of function cnmtf
#
# [Output]
#' @return
#' * \code{rho} List of dispersion coefficients for the avg. consensus matrices of U and V
#' * \code{clus.U, clus.V} Clustering memberships
#' * \code{Os} Consensus SNV score
#' * \code{lniters, lJ, lrelErr} Variables to follow-up the algorithm performance
#' * \code{res.alg} Object with all the matrices at each repetition. Only IF parallel.opt = FALSE
#' * \code{max.parameters} sObject with the penalisation terms estimated in the last iteration of each repetition.
#' * \code{lparameters}  If \code{find.parameters = TRUE} this function only returns a list of estimations for the parameters
#'
#' @md
#' @author Luis G. Leal, \email{lgl15@@imperial.ac.uk}
#' @family Factorisation functions
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


	consensus.clust = function(R, #Relationship matrix
	                           Wv, Wu, #Prior knowledege networks for SNVs and Patients
	                           HLH, Vo, #Matrices to correct for confounders
	                           k, #Number of clusters
	                           lparameters, #Penalization parameters
	                           iters, #Number of iterations
	                           labelsU, labelsV, #Labels of matrices
	                           run.t, #Number of repetitions
	                           calcObj, calcObj2, #Check convergency each X number of iterations after first Y iterations
	                           init, #Type of seeding/initialisation of matrices in the algorithm
	                           parallel.opt, #Run the algorithm in parallel
	                           n.cores = NULL, #Number of cores to use in the parallel processing
	                           set.alg = NULL, #Set of algorithm results can be provided to find consensus results
	                           set.clus.V = NULL, #Set of clustering of patients provided to find consensus results
	                           export.res = TRUE, #Whether the object with list of results is returned
	                           V.init = NULL, U.init = NULL, #Initialisation of U and V using PSVD (1)
	                           do.U = TRUE, #Perform consensus clustering on U matrix
	                           do.V = TRUE, #Perform consensus clustering on V matrix
	                           do.O = TRUE, #Perform consensus clustering on Omega matrix
	                           find.parameters = FALSE, #Run the NMTF algorithm to find estimations of parameters
														 estimate.penalisation = FALSE, #Return the penalisation terms evaluated in the last iteration.
														 display.iters = FALSE #Display the iterations of function cnmtf

	){


	  #-------------------------------
	  # Run the algorithm: Repetitions in PARALLEL
	  #-------------------------------


				#If algorithm results are provided then extract them

				  if( !is.null(set.alg) ){

				    cat("Consensus V and U provided by user from object set.alg.", "\n")
				    rho = set.alg[[1]][[1]]
			      clus.U = set.alg[[1]][[2]]
			      clus.V = set.alg[[1]][[3]]

				  }



				#If results of algorithm are not provided AND the machine CAN store multiple runs of the algorithm in memory, then run the algorithm in parallel

					if( is.null(set.alg) & parallel.opt == TRUE ){

						res.alg = list();
						print(paste("Running cNMTF in PARALLEL with parameters: ", "k1:", k[1], "k2:", k[2], "lparameters:", paste(names(unlist(lparameters)),unlist(lparameters), collapse = " "), sep = " "))

						#Register the number of cores
			  			if( is.null(n.cores) ){
			  			  doMC::registerDoMC(cores = (detectCores() - 2 ))
			  			}else{
			  			  doMC::registerDoMC(cores = n.cores)
			  			}

						#Run each repetition of cNMTF in parallel
						res.alg <- foreach(runs = 1:run.t) %dopar% {
							res.alg = cnmtf(R=R, Wv=Wv, Wu=Wu,
															HLH= HLH, Vo=Vo,
															k=k, lparameters=lparameters,
															iters = iters,
															calcObj = calcObj, calcObj2 = calcObj2,
															init = init, V.init = V.init, U.init = U.init, displ = display.iters)
						}

					}else{
						print(paste("Running cNMTF in SERIES with parameters: ", "k1:", k[1], "k2:", k[2], "lparameters:", paste(names(unlist(lparameters)),unlist(lparameters), collapse = " "), sep = " "))
					}



			  #-------------------------------
			  # Calculate conectivity matrices for U and V and their dispersion coefficients
			  # Also calculate aggregated omega matricx
			  #-------------------------------


			  #Initialize connectivity and consensus matrices
			    n = nrow(R); m = ncol(R)
			    if(do.U == TRUE){ Cu = matrix(0, nrow = n, ncol = n); Ua = matrix(0, nrow = n, ncol = k[1]) }else{ Ua = NA }
			    if(do.V == TRUE){ Cv = matrix(0, nrow = m, ncol = m); Va = matrix(0, nrow = m, ncol = k[2]) }else{ Va = NA }
			    if(do.O == TRUE){ Or = matrix(0, nrow = n, ncol = m); Sa = matrix(0, nrow = k[1], ncol = k[2]) }else{ Or = NA; Sa = NA }


			  #Initialize vectors of convergency variables and penalisation terms
			    lniters = lJ = lrelErr <- rep(NA, run.t)

			  #Initialize vectors to estimate maximum penalisation terms
			    if(estimate.penalisation == TRUE) {
			    	malpha1 = malpha2 = mgamma1 = mgamma2 <- rep(0, run.t) #Vectors of maximum estimatedparameters at each run
			    }

			  #Extract results of each run
			  for (i in 1:run.t){#For each repetition of NMTF


			    #If the machine CANNOT store multiple runs of the algorithm in memory
			    #Run repetitions of the algorithm in SERIES and process the matrices at each run

			    if( is.null(set.alg) & parallel.opt == FALSE){

			      #Print progress
			    	cat("\n","Repetition: ", i, "\n")
			      cat("Calculating U,V,S matrices for repetition",i,"\n")

			      #Store/Overwrite the results in the first position of the list
			      res.alg = NULL
			      res.alg[[1]] = cnmtf(R=R, Wv=Wv, Wu=Wu,
			                           HLH, Vo=Vo,
			                           k=k, lparameters=lparameters,
			                           iters = iters,
			                           calcObj = calcObj, calcObj2 = calcObj2, init = init,
			                           V.init = V.init, U.init = U.init, displ = display.iters)


			      #Update position to be read in the results list
			      i.run <- 1

			    }else{
			      i.run <- i
			    } #End If


			    #Extract U,V,S matrices at this run
			      cat("Extracting U,V,S matrices for repetition",i,"\n", as.character(Sys.time()), "\n")
			      Ui = as.matrix(res.alg[[i.run]][[1]])
			      Vi = as.matrix(res.alg[[i.run]][[2]])
			      Si = as.matrix(res.alg[[i.run]][[3]])


			    #If tracking penalisation then re-calculate penalisation terms at each repetition

			      if(estimate.penalisation == TRUE) {

			      	#Assess main part of the objective function
			      	tmain = norm(R - Ui%*%Si%*%t(Vi),"F")^2

			      	#Compute parameters
			      	for(q in 1:length(unlist(lparameters) )){

			      			if( names(unlist(lparameters) )[q] == "alpha1" & !is.null(Wu)){
			      				cat("Finding maximum alpha1 at repetition",i,"\n")
			      				Du = diag(rowSums(Wu))
			      				Lu = Du - Wu
			      				malpha1[i] = tmain / sum( diag(t(Ui) %*% Lu %*% Ui) )
			      				cat("Maximum alpha1:", malpha1[i],"\n")
			      			}

			      			if( names(unlist(lparameters) )[q] == "alpha2" & !is.null(Wv)){
			      				cat("Finding maximum alpha2 at repetition",i,"\n")
			      				Dv = diag(rowSums(Wv))
			      				Lv = Dv - Wv
			      				malpha2[i] = tmain / sum( diag( t(Vi) %*% Lv %*% Vi) )
			      			}


			      			if( names(unlist(lparameters) )[q] == "gamma2" & !is.null(Vo)){
			      				cat("Finding maximum gamma2 at repetition",i,"\n")
			      				mgamma2[i] =  tmain / norm(Vi - Vo,"F")^2
			      				cat("Maximum gamma2:", mgamma2[i],"\n")

			      			}

				      		if( names(unlist(lparameters) )[q] == "gamma1" & !is.null(HLH) ){
				      			cat("Finding maximum gamma1 at repetition",i,"\n")
				      			mgamma1[i] =  tmain / sum( diag(Vi %*% t(Vi) %*% HLH) )
				      			cat("Maximum gamma1:", mgamma1[i],"\n")

				      		}



			      	}#End for 1:length(lparamenters)

			      }


			    #If not results of algorithm are provided then assess connectivity matrices for clustering
			    if(is.null(set.alg)) {


			      if(do.V == TRUE){

			        cat("Calculating connectivity matrix Cv for repetition",i,"\n")
			        #Calculate connectivity matrices (slow step!!)
			        Cvi = clus.membership(Vi)
			        #Add connectivity matrices to consensus matrices
			        Cv = Cv + Cvi
			        #Create matrix of accumulated V matrices
			        Va = Va + Vi


			      }

			      #Same for U
			      if(do.U == TRUE){

			        cat("Calculating connectivity matrix Cu for repetition",i,"\n")
			        #Calculate connectivity matrices (slow step!!)
			        Cui = clus.membership(Ui)
			        #Add connectivity matrices to consensus matrices
			        Cu = Cu + Cui
			        #Create matrix of accumulated U matrices
			        Ua = Ua + Ui

			      }
			    } #End If


			    if(do.O == TRUE){ #This procedure does not depend on "do.V"

			      cat("Calculating Omega matrix for repetition",i,"\n")

			      #Extract patient clustering membership at this run
			      clus.Vi = rep(NA, nrow(Vi))
			      for(a in 1:nrow(Vi)){
			        clus.Vi[a] = which.max(Vi[a,])
			      }
			      clus.Vi = cbind(labelsV, clus.Vi)

			      #Scale/normalize U, S matrices
			      #for(y in 1:ncol(Si)){
			       # Si[,y] = Si[,y] / sum(Si[,y])
			      #}

			      #for(x in 1:nrow(Ui)){
			       # Ui[x,] = Ui[x,] / sum(Ui[x,])
			      #}

			      #Calculate omega matrix
			      Oi = Ui %*% Si

			      #Create a matrix Or of SNPs by Patients with the omega values of each SNP at each patient in this run
			      for(x in 1:n){
			        for(z in 1:k[2]){
			          Or[x, clus.Vi[,2] == z] <- Or[x, clus.Vi[,2] == z] + Oi[x,z]
			        }
			      }


			    }else{
			      Or = NULL
			    }

			    lniters[i] = as.matrix(res.alg[[i.run]][[4]])
			    lJ[i] = as.matrix(res.alg[[i.run]][[5]])
			    lrelErr[i] = as.matrix(res.alg[[i.run]][[6]])

			    #Create matrix of accumulated S matrices
			    Sa = Sa + Si

			  }#End for through repetitions

			  #----------------------------------------------------
			  # Estimate maximum parameters
			  #----------------------------------------------------

				    if(estimate.penalisation == TRUE) {

				    	#Plot vectors of parameters
				    		#X11()
				    		par(mfrow = c(1,4))
				    		boxplot(malpha1,main = "alpha1", las = 2); boxplot(malpha2,main = "alpha2", las = 2)
				    		boxplot(mgamma1,main = "gamma1", las = 2); boxplot(mgamma2,main = "gamma2", las = 2)

				    	#Return median value of parameters as the maximum
				    	max.parameters = list(alpha1 = c(0,median(malpha1)), alpha2 = c(0,median(malpha2)),
				    												gamma1 = c(0,median(mgamma1)), gamma2 = c(0,median(mgamma2)))

				    }else{
				    	max.parameters = NULL
				    }


			  #----------------------------------------------------
			  # Calculate consensus matrices of U and V
			  # Also calculate the consensus omega
			  #----------------------------------------------------



			  #If not results of algorithm are provided then calculate hierarchical clustering from consensus matrices
			  if(is.null(set.alg)) {

			    cat("Calculating consensus U and V", "\n")

			    rho = c(NA,NA)

			    #Procedure for V

			        if(do.V == TRUE){

			          #Calculate average consensus matrices and dispersion coefficients
			            aCv = Cv/run.t
			            rho[2] = sum(sum(4*(aCv-0.5)^2))/(nrow(aCv)^2)
			          #Perform hierarchical clustering and print the clusters
			            clus.V = hierarchical.clust( aCv, k[2], labelsV )

			        }else if(!is.null(set.clus.V)){ #Only if clus.V is provided then assign it

			          cat("Consensus V provided by user from object set.clus.V", "\n")
			          clus.V <- set.clus.V

			        }else{

			          cat("Consensus V set to NULL", "\n")
			          clus.V = NA

			        }

			    #Procedure for U

			        if(do.U == TRUE){

			          #Calculate average consensus matrices and dispersion coefficients
			            aCu = Cu/run.t
			          #Perform hierarchical clustering and print the clusters
			            rho[1] = sum(sum(4*(aCu-0.5)^2))/(nrow(aCu)^2)
			            clus.U = hierarchical.clust( aCu, k[1], labelsU )

			        }else{
			          cat("Consensus U set to NULL", "\n")
			          clus.U = NA
			        }

			  }

			  if(do.O == TRUE){

			    cat("Calculating consensus Omega", "\n")

			    #Calculate consensus omega

			  	  #Create matrix with number of columns equals number of final patient clusters
			   	  n.clus.pat = max(as.numeric(clus.V[,2]))
			  	  On = matrix(0, nrow = n, ncol = n.clus.pat)

			  	  #Fill the matrix
			  	  for(x in 1:n){
			  	    for(z in 1:n.clus.pat){
			  	      On[x,z] <- median(Or[x, clus.V[,2] == z])
			  	    }
			  	  }

			  	  #Add rownames and colnames to consensus omega
			  	  rownames(On) <- labelsU
			  	  colnames(On) <- paste("p",1:ncol(On),sep="")

			  }else{
			    cat("Consensus Omega set to NULL", "\n")
			    On = NA
			    Or = NA
			  }

			  #Return results of consensus function
			  if( export.res == FALSE ){ #Decide whether the object with list of results is returned or not
			    res.alg = NA
			  }
			  #Special case for randomisations: Dont export Or either
				if( do.U == FALSE & do.V == FALSE){
					Or = NA
				}


			  return(list(rho = rho, clus.U = clus.U, clus.V = clus.V, On = On,
			  						lniters = lniters , lJ = lJ, lrelErr = lrelErr,
			  						res.alg = res.alg, max.parameters = max.parameters,
			  						Or = Or, Ua = Ua, Va = Va, Sa = Sa, lparameters = lparameters))

	}



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#' @title Calculate connectivity matrix
#' @description  Function to extract clustering assignments from a factorization matrix
#'
#	[Input]
#' @param H  Factorization matrix
#' @param indCluster  Clustering assignemnt vector
#'
# [Output]
#' @return \code{C} Binary connectivity matrix
#'
# [Example]
#' @examples
#' clus.membership(V)
#'
#' @md
#' @author Luis G. Leal, \email{lgl15@@imperial.ac.uk}
#' @family Factorisation functions
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


  	clus.membership = function(H){

  	  # Size of data matrix H
  	  a = nrow(H)
  	  b = ncol(H)

  	  # Declare connectivity matrix
  	  C = matrix(0, nrow = a, ncol = a)

  	  # Find largest value in each row of the matrix and return its position (index)
  	  (indCluster = as.factor(apply(H, 1, which.max) ))

  	  # Conform the connectivity matrix
  	  for(i in levels(indCluster)){
  	    pos <- indCluster == i
  	    C[pos,pos] <- 1
  	  }

  	  return(C)
  }



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#' @title Hierarchical clustering from consensus connectivity matrices
#' @description Function for extracting cluster labels based on hierarchical clustering
#'
# [Input]:
#' @param aC consensus matrix
#' @param k rank parameter (number of clusters)
#' @param linkF linkage function ('average','centroid','complete', etc.)
#' @param labels sample labels
#' @param out_filename exporting label assignment to clusters
# [Output]:
#' @return \code{res} cluster assignment
#'
#' @author Luis G. Leal, \email{lgl15@@imperial.ac.uk}
#' @family Factorisation functions
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


		hierarchical.clust = function(aC, ki, labels, linkF = c("average")){

				  #Transform consensus matrix into distance matrix
				  aC = aC - diag(diag(aC)) #Check this!
				  distC = as.dist(1.0 - aC)

				  #Compute linkage
				  tree = hclust(distC, method = linkF)

				  #Find cluster membership
				  j = 0; clust.id = 0
				  #If there are clusters of size < 5 then reduce the number of clusters
				  while( sum( table(clust.id) < 5 ) != 0 & (j < ki) ) {
				    clust.id = cutree(tree,ki - j)
				    j = j + 1
				  }

				  #Check labels length

				  if (length(clust.id) != length(labels)){
				    print('Wrong number of labels')
				    stop()
				  }
				  #Export clusters

				  res = cbind(labels,clust.id)
				  return(res)
		}


