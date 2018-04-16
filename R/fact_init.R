###############################################################
# cNMTF
#	2. Factorisation functions
# 2.1 Functions to initialise matrices prior factorisation
#
# Corresponding author:
# Luis Leal, Imperial College London, email: lgl15@imperial.ac.uk
###############################################################


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#' @title SVD Initialisation for NMF
#' @description Function to initialize Non-negative Matrix Factorization algorithms
#' @references C. Boutsidis and E. Gallopoulos, SVD-based initialization: A head start for nonnegative matrix factorization, Pattern Recognition, Elsevier
#'
# [Input]:
#' @param A Nonnegative matrix A, \emph{n x n}
#' @param k Rank of the computed factors (number of clusters)
#'
# [Output]:
#' @return \code{W} nonnegative matrix, \emph{n x k}
#'
#' @md
#' @author Luis G. Leal, \email{lgl15@@imperial.ac.uk}
#' @family Factorisation functions
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


    psvd.init = function(A,k){

      #Check if matrix contains negative elements
      if (sum(A<0) > 0){
        print('Error: The input matrix contains negative elements.')
        stop()
      }
      #Size of data matrix A
      n = nrow(A)
      m = ncol(A)
      #Matrices of the factorization
      W = matrix(0,nrow = n, ncol = k)

      #1st SVD --> partial SVD rank-k to the input matrix A.
      svd.res = svd(A,k)
      U <- as.matrix(svd.res[[2]])
      S <- as.matrix(diag(svd.res[[1]]))
      V <- as.matrix(svd.res[[3]])


      #Choose the first singular triplet to be nonnegative
      W[,1] <- sqrt(S[1,1]) * abs(U[,1])

      #2nd SVD for the other factors
      for(i in 2:k){
        x = U[,i]
        y = V[,i]
        xp = pos.constrain(x); xn = neg.constrain(x); yp = pos.constrain(y); yn = neg.constrain(y)
        xpnrm = norm(xp, "2"); ypnrm = norm(yp, "2"); mp = xpnrm*ypnrm
        xnnrm = norm(xn, "2"); ynnrm = norm(yn, "2"); mn = xnnrm*ynnrm
        if (mp >= mn){
          u = xp/xpnrm; v = yp/ypnrm; sigma = mp
        }else{
          u = xn/xnnrm; v = yn/ynnrm; sigma = mn
        }
        W[,i] = sqrt(S[i,i]*sigma)*u
      }

      #Modify entries near to zero
      W[W<0.0000000001]<-0.1
      rownames(W) <- rownames(A)
      return(W)
    }



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#' @title Constraining negative entries in a matrix
#' @description This function sets to zero the negative elements of a matrix
#' @references C. Boutsidis and E. Gallopoulos, SVD-based initialization: A head start for nonnegative matrix factorization, Pattern Recognition, Elsevier
#'
# [Input]:
#' @param A Nonnegative matrix A, \emph{n x n}
#'
# [Output]:
#' @return \code{A} Constrained nonnegative A matrix, \emph{n x n}
#'
#' @md
#' @author Luis G. Leal, \email{lgl15@@imperial.ac.uk}
#' @family Factorisation functions
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


    pos.constrain = function(A){
      A [ A < 0 ] <- 0
      return(A)
    }



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#' @title Constraining positive entries in a matrix
#' @description This function sets to zero the positive elements of a matrix and takes the absolute value of the negative elements.
#' @references C. Boutsidis and E. Gallopoulos, SVD-based initialization: A head start for nonnegative matrix factorization, Pattern Recognition, Elsevier
#'
# [Input]:
#' @param A Nonnegative matrix A, \emph{n x n}
#'
# [Output]:
#' @return \code{A} Constrained nonnegative A matrix, \emph{n x n}
#'
#' @md
#' @author Luis G. Leal, \email{lgl15@@imperial.ac.uk}
#' @family Factorisation functions
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


    neg.constrain = function(A){
      A [ A > 0 ] <- 0
      A <- abs(A)
      return(A)
    }



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#' @title Initialising the low-rank matrices U and V
#' @description This function creates initialisations of U and V. It can also read stored initialisations to speed up the cNMTF algorithm.
#'
# [Input]:
#' @param R Relationship matrix, \emph{n x m}
#' @param k1,k2 Rank parameters (number of clusters for U and V, respectively)
#' @param name.init Name of the workspace with initialisations of U and V
#' @param work.dat Folder to read or write a workspace with initialisations
#' @param logfile Name of a file to write progess/status of this function
#'
# [Output]:
#' @return \code{U.init}, \code{V.init} intialisations of U and V.
#'
#' @md
#' @author Luis G. Leal, \email{lgl15@@imperial.ac.uk}
#' @family Factorisation functions
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


  initialise.UV = function(R, #Relationship matrix
                           name.init = NULL, #Name of workspace with initialisations of U and V
                           work.dat = NULL, #Folder to save and load workspaces
                           k2, #Number of clusters for U
                           k1, #Number of clusters for V
                           name.exp = NULL, #Name of experiment to save files
                           log.file = NULL
  )

  {

    if( !is.null(name.init)){

      if(file.exists(name.init) ){
        #Load initialisation from specific file
        msg <- paste("Loading initialisation of U and V from file.", "\n", as.character(Sys.time()), "\n", sep=" "); cat(msg); cat(msg, file = log.file, append = TRUE)
        load(file = name.init)
        if(exists("U.init") == FALSE){
          U.init = U.init.max[,1:k1]
        }
        U.init = U.init[,1:k1]
        init = 1
        if( !( all(rownames(U.init) == rownames(R)) | all(rownames(V.init) == colnames(R))) | is.null(rownames(U.init)) | is.null(rownames(V.init))){
          msg <- paste("Initialisations provided in file do not match current R.", "\n", as.character(Sys.time()), "\n", sep=" "); cat(msg); cat(msg, file = log.file, append = TRUE)
          stop()
        }

      }else{
        msg <- paste("Initialisation file not found.", "\n", as.character(Sys.time()), "\n", sep=" "); cat(msg); cat(msg, file = log.file, append = TRUE)
        stop()

      }


    }else{


      msg <- paste("Initialising U and V.", "\n", as.character(Sys.time()), "\n", sep=" "); cat(msg); cat(msg, file = log.file, append = TRUE)
      V.init = psdv.init(t(R),k2)
      U.init = psdv.init(R,k1)
      init = 1

      #Save initialisations
      save(list = c("V.init","U.init"),file = paste(work.dat,"init_", name.exp,".RData",sep=""))
      msg <- paste("Initialisation printed to workspace:" , paste(work.dat,"init_", name.exp,".RData",sep=""), "\n", as.character(Sys.time()), "\n", sep=" "); cat(msg); cat(msg, file = log.file, append = TRUE)

    }

    #Print dimensions of initialisations
    msg <- paste("Dimensions of U and V initialisations:" , paste(dim(U.init),collapse=" x "), "and" ,paste(dim(V.init),collapse=" x "), "\n", as.character(Sys.time()), "\n", sep=" "); cat(msg); cat(msg, file = log.file, append = TRUE)

    return( list(U.init, V.init))

  }


