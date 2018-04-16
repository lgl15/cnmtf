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
#' @param H Centering matrix, \emph{m x m}
#' @param L Side information kernel matrix (i.e. Kernel on the confounder random variable) [mxm]
#' @param Vo Outcome label matrix, \emph{m x k2}
#' @param k Vector of rank parameters, \emph{k1 x k2}
#' @param iters Default number of itersatiors
#' @param calcObj Check convergency each X number of itersations
#' @param displ Logical. Print number of iterations
#' @param tof RelativeError
#' @param lparameters Vector of regularization parameters
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


    cnmtf = function(R,HLH,Vo,Wu,Wv,k,
                     iters = 1, #Default number of iteratiors
                     calcObj, calcObj2, #Check convergency each X number of iterations after first Y iterations
                     displ = TRUE, #Print number of iterations
                     tof = 1e-5, #Minimum relativeError
                     lparameters, #List of parameters
                     init = 1, #Initialize the matrices randomly (0) or by using PSVD (1)
                     V.init = NULL, U.init = NULL, #Initialisation of U and V using PSVD (1)
                     init.S = FALSE #Initialisation of S using my method
    )
    {


      #Define dimmensiion
      m = ncol(R)
      n = nrow(R)

      #Replace zeros with small values in the relationship matrix to speed up process
      R[R==0] = 2^(-16)


      #Initializing low-rank matrices U,V and S
      if (displ == TRUE){ cat("Initializing low-rank matrices U,V and S", "\n") }

      if( init == 1 ){#PSDV initialization for U and V
        S = matrix(runif(k[1]*k[2]), nrow=k[1], ncol=k[2])

        if( is.null(V.init) | is.null(U.init) ){
          V = psdv.init(t(R),k[2])
          U = psdv.init(R,k[1])

        }else{
          V = V.init
          U = U.init
        }

        #Initialisation of S
        if( init.S == TRUE){
          S = matrix(NA, nrow=k[1], ncol=k[2])

          for(k1 in 1:k[1]){
            c.vec = rep(NA, m)
            for(j in 1:m){
              c.vec[j] = abs(cor( R[,j], U[,k1] ))
            }

            for(k2 in 1:k[2]){
              S[k1,k2] = abs(cor( c.vec, V[,k2] ))
            }
          }

        }

      } else { #Random initialization
        U = matrix(runif(n*k[1]), nrow=n, ncol=k[1])
        V = matrix(runif(m*k[2]), nrow=m, ncol=k[2])
        S = matrix(runif(k[1]*k[2]), nrow=k[1], ncol=k[2])
      }


      #Replace zeros with small values in the initialised matrices to speed up process
      S[S==0] = 2^(-16)
      V[V==0] = 2^(-16)
      U[U==0] = 2^(-16)

      #Extract parameters and asign them to variables in this environment
      for(i in 1:length(lparameters)){
        assign(names(lparameters)[i] , as.numeric(unlist(lparameters[i])) )
      }

      #Assess diagonal degree matrix and Laplacian matrix
      if(alpha1 != 0){
        if (displ == TRUE){ cat("Calculating Laplacian of Wu","\n") }
        Du = diag(rowSums(Wu))
        Lu = Du - Wu
      }
      if(alpha2 != 0){
        if (displ == TRUE){ cat("Calculating Laplacian of Wv","\n") }
        Dv = diag(rowSums(Wv))
        Lv = Dv - Wv
      }


      #Initialize objective function
      objF.previous = Inf

      if (displ == TRUE){ cat("Iterating","\n") }


      #Number of times the function has increased in a row
      count.increase = 0


      #Run iterations
      iters.gamma1 = 50
      if(gamma1 > 0){

        mode.gamma2 = FALSE
        if (displ == TRUE){ cat("Mode gamma2", mode.gamma2,  "\n"  ) }

      }else{
        mode.gamma2 = TRUE
      }


      for (i in 1:iters){ #For each iteration

        #cat("Iteration: ",i, "\n")
        #cat("Precalculate matrices to speed up process", "\n")

        #Precalculate matrices to speed up process
        US = U %*% S
        tVV = t(V) %*% V
        UStVV= US %*% tVV
        RV = R %*% V

        #Update rule for S
        #cat("Update rule for S: up", "\n")
        up_s = t(U) %*% RV
        #cat("Update rule for S: down", "\n")
        down_s = t(U) %*% UStVV
        S = S * (up_s/down_s)

        #Update rule for U
        #cat("Update rule for U: up", "\n")
        up_u = RV %*% t(S)
        #cat("Update rule for U: down", "\n")
        down_u = UStVV%*%t(S)
        if ( alpha1 != 0 & ( i %% calcObj) == 0 ){
          up_u = up_u + alpha1 * Wu %*% U
          down_u = down_u + alpha1 * Du %*% U
        }
        U = U*(up_u/down_u)

        #Update rule for V
        #cat("Update rule for V: up", "\n")
        up_v = t(R) %*% US
        #cat("Update rule for V: down", "\n")
        down_v = V %*% t(S) %*% t(U) %*% US
        if ( alpha2 != 0 ){
          up_v = up_v + alpha2 * Wv %*% V
          down_v = down_v + alpha2 * Dv %*% V
        }
        if ( gamma1 != 0 & mode.gamma2 == FALSE){ #& ( i %% calcObj) == 0
          down_v = down_v + gamma1 * HLH %*% V
        }
        if ( gamma2 != 0 & mode.gamma2 == TRUE){
          up_v = up_v + gamma2 * Vo
          down_v = down_v + gamma2 * V
        }
        V = V*(up_v/down_v)




        #Compute residual and print relative error and objective function

        if( (( i %% calcObj) == 0 & (i >= calcObj2)) || i == iters ){ #If the itersation is a multiple of of calcObj or if last itersation was reached

          #Compute objective function

          tmain = norm(R - U %*% S %*% t(V),"F")^2

          if( alpha1 != 0 ){ talpha1 = alpha1 * sum( diag(t(U) %*% Lu %*% U) ) } else { talpha1 = 0 }
          if( alpha2 != 0 ){ talpha2 = alpha2 * sum( diag(t(V) %*% Lv %*% V) ) } else { talpha2 = 0 }
          if( gamma1 != 0 & mode.gamma2 == FALSE){ tgamma1 = gamma1 * sum( diag(V %*% t(V) %*% HLH) ) } else { tgamma1 = 0 }
          if( gamma2 != 0 & mode.gamma2 == TRUE){ tgamma2 = gamma2 * norm( V - Vo,"F")^2 } else { tgamma2 = 0 }

          objF.current = tmain + tgamma1 + tgamma2 + talpha1 + talpha2

          #Computing relative change of the objective function
          relErr = (objF.previous - objF.current)/max(1,objF.previous)/calcObj
          if(is.na(relErr)){relErr = Inf}


          #Print relative error and objective function

          if (displ == TRUE){ print(paste("Iteration: ",i, "; Relative error: ", format(relErr,scientific = TRUE, digits = 3), "; Objective function: ", format(objF.current, scientific = TRUE, digits = 3),sep = "")) }


          #If gamma1 or gamma2 are penalizing, track their effect (Use out.cat and pop from global environment variables)
          if(gamma1 != 0 | gamma2 != 0){

            #Find cluster membership at this iteration
            #This section makes use of the environmental variables: pop and out.cat

            #Extract cluster membership
            clus.Vi = rep(NA, nrow(V))
            for(a in 1:nrow(V)){
              clus.Vi[a] = which.max(V[a,])
            }

            #Calculate similarities between partitions of individuals
            simVO.i = igraph::compare(as.numeric(as.factor(out.cat)), as.numeric(clus.Vi), method  =  c("nmi"))
            simVL.i = igraph::compare(as.numeric(as.factor(pop))[pop!="Unknown"], as.numeric(clus.Vi)[pop!="Unknown"], method  =  c("nmi"))
            simOL.i = igraph::compare(as.numeric(as.factor(pop))[pop!="Unknown"], as.numeric(as.factor(out.cat)[pop!="Unknown"]), method  =  c("nmi"))


            if (displ == TRUE){

              #Print similarities
              cat("SimVO", simVO.i, "SimVL", simVL.i, "SimOL", simOL.i, "\n"  )

              #Print number of significant SNVs
              cat("Number of associations retrieved :", round( length(sig.sc.i), digits = 2), "\n", sep=" ")
              cat("Number of known associations retrieved :", round( sum(sig.sc.i %in% snps.known), digits = 2), "\n", sep=" ")
            }

            #Find number of known associations retrieved at this iteration
            #This section makes use of the environmental variables: snps.known

            #Treshold
            tao.sc.i = 4

            #Calculate score matrix
            Oi = U %*% S
            Oi = scale(Oi) #The scaling of columns help to have a more normal alike distribution of delta scores

            #Delta score for a specific combination of patient clusters
            dOi = Oi[,2] - Oi[,1]

            #Calculate desviations from SD
            sd.sc.i = (dOi - mean(dOi))/sd(dOi)

            #Find significant SNVs at a given treshold
            sig.sc.i = rownames( R ) [ abs( sd.sc.i ) >  tao.sc.i ]



            #Check if SimVL increases

            if( simVL.i > 0.02 & 	mode.gamma2 == TRUE & gamma1!= 0){ #2% of similarity would reduce at maximum the number of false positives for PS
              if (displ == TRUE){ print(paste('SimVL is increasing above 0.01. Current objetive function ', objF.current, sep = "")) }
              break #stop (changed to break to return something, even unfinished optimizations)
            }


            #Check the Mode

            if( (i >= iters.gamma1  | tof >= relErr) & mode.gamma2 == FALSE ){
              mode.gamma2 = TRUE
              if (displ == TRUE){ cat("Mode gamma2", mode.gamma2,  "\n"  ) }
              objF.previous = Inf
              next
            }

          }

          #If the function have skipped convergency less than 2 times and its value is Inf
          if( objF.current == Inf & count.increase < 2 ){

            if (displ == TRUE){ cat("Objective function is Inf. Skipping convergency", count.increase, "time (s).", "\n") }
            count.increase = count.increase + 1
            next

          }else if( objF.current == Inf & count.increase == 2 ){

            if (displ == TRUE){ print(paste('Objective function is increasing from ', objF.previous, " to ", objF.current, sep = "")) }
            break
          }


          #Check if objective function increases
          if(objF.current > objF.previous){

            #Number of times the function has increased in a row
            count.increase = count.increase + 1

            if (displ == TRUE){ print(paste('Objective function is increasing from ', objF.previous, " to ", objF.current, sep = "")) }


            #If it has been increasing during the last two checks then stops
            if( count.increase == 2){

              if (displ == TRUE){ print(paste('Objective function is increasing from ', objF.previous, " to ", objF.current, sep = "")) }
              break #stop (changed to break to return something, even unfinished optimizations)

            }

            if (displ == TRUE){ cat("Skiping the checking of convergency", count.increase, "time.", "\n") }

            next #Skips the checking of convergency

          }




          #Update objective function
          objF.previous <- objF.current

          #Check convergency
          if (tof >= relErr || i == iters){

            if (displ == TRUE){ print(paste('Local minima found')) }

            break
          }
        }


      } #End iterations



      #Conform cnmtf results
      return(list(U, V, S, i, objF.current, relErr))
    }

