###############################################################
# cNMTF
#	2. Factorisation functions
# 2.5 Main function to perform consensus cNMTF and randomisations
#
# Corresponding author:
# Luis Leal, Imperial College London, email: lgl15@imperial.ac.uk
###############################################################


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#' @title Calculate Omega matrix
#' @description Function to calculate the Omega matrix (SNV scores per cluster of patients)
#'
#[Input]:
#' @param R Relationship matrix
#' @param out.cat Categorical outcome variable
#' @param pop Population variable
#' @param log.file Log file to track progress of the function

#Variables to Save/Load data and workspaces
#' @param name.exp Name of experiment to save files
#' @param name.init Name of workspace with initialisations of U and V
#' @param work.dat Folder to save and load workspaces

#Number of clusters
#' @param define.k  Method to define \code{k1}. It can be: \code{c("user","NMTF","PCA")}
#' @param k1 Number of SNV clusters
#' @param k2 Number of patient clusters equals number of levels in the outcome


#Penalisation parameters
#' @param estimate.par Logical. Estimate penalisation parameters. If FALSE then provide the parameters in the argument "wparameters". If TRUE then provide preliminar weights to the parameters in the argument "wparameters".
#' @param wparameters Either the Penalization parameters (if \code{estimate.par == FALSE}) or Weights of penalisation terms to be computed (if estimate.par == TRUE)
#' @param save.parameters Logical. Save parameters to file
#' @param run.t.exp Number of repetitions for the experiment
#' @param run.t.par Number of repetitions for parameters fitting
#' @param range.parameters Range of parameters to be evaluated (if \code{estimate.par == FALSE})
#' @param sequential.estimation Set the parameters in a specific order and carry the optimals
#' @param max.try0 Maximum number of tries to fit the parameter
#' @param tao.sc  Trehsold of standard deviations for the SNV score
#' @param snps.known List of known associations

#Variables to control performance of the algorithm
#' @param calcObj,calcObj2 Number of iterations to check convergency. Check convergency each \code{calcObj} number of iterations after first \code{calcObj2} iterations
#' @param init Type of seeding/initialisation of matrices in the algorithm
#' @param parallel.opt Run some instances of the algorithm in parallel
#' @param n.cores Number of cores to use in the parallel processing
#' @param iters Number of iterations
#' @param do.U Perform clustering of SNPs
#' @param display.iters Display the iterations of function cnmtf

#Randomisations
#' @param score.pvalues Estimate p-values for the scores \code{c(TRUE, FALSE, "only")}.
#' @param randomisations Number of randomisations
#' @param random.parallel Logical. Run the radomisations in parallel

#Construction of penalisation terms
#' @param file.Gu Workspace with adjancency matrix for the SNV-SNV network.
#'
#[Output]:
#' @return The function internally generates:
#' * Estimation of number of SNP clusters (k1)
#' * Estimation of penalisation parameters (alpha1, alpha2, gamma1, gamma2)
#'
#' The function prints the following objects in workspaces:
#' * Initialisation of matrices U.init and V.init
#' * \code{res.cnmtf}: results of factorisations
#' * \code{lcnmtf.ran}: results of algorithm applied on randomisations of R
#'
#' @md
#' @author Luis G. Leal, \email{lgl15@@imperial.ac.uk}
#' @family Factorisation functions
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

	score.cnmtf = function( R, #Relationship matrix
	                        out.cat = NULL, #Categorical outcome variable
	                        pop = NULL, #Population variable
													log.file = NULL,

	                        #Variables to Save/Load data and workspaces
	                        name.exp = NULL, #Name of experiment to save files
	                        name.init = NULL, #Name of workspace with initialisations of U and V
	                        work.dat = NULL, #Folder to save and load workspaces

	                        #Number of clusters
	                        define.k = c("user","NMTF","PCA"), #Method to define k1
	                        k1 = NULL, #Number of clusters
	                        k2 = NULL, #Number of clusters equals number of levels in the outcome

	                        #Penalisation parameters
	                        estimate.par = TRUE, #Estimate penalisation parameters. If FALSE then provide the parameters in the argument "wparameters". If TRUE then provide preliminar weights to the parameters in the argument "wparameters".
	                        wparameters = NULL, #Either the Penalization parameters (if estimate.par == FALSE) or Weights of penalisation terms to be computed (if estimate.par == TRUE)
													save.parameters = TRUE, #Save parameters to file
													run.t.exp = NULL, #Number of repetitions for the experiment
	                        run.t.par = 4, #Number of repetitions for parameters fitting
	                        range.parameters = NULL, #Range of parameters to be evaluated (if estimate.par == FALSE)
													sequential.estimation = FALSE, #Set the parameters in a specific order and carry the optimals
													max.try0 = 4,  #Maximum number of tries to fit the parameter
													tao.sc  = 4, #Trehsold of standard deviations for the SNV score
													snps.known = NULL, #List of known associations (P)


	                        #Variables to control performance of the algorithm
	                        calcObj = NULL, calcObj2 = NULL, #Check convergency each X number of iterations after first Y iterations
	                        init = NULL, #Type of seeding/initialisation of matrices in the algorithm
	                        parallel.opt = FALSE, #Run some instances of the algorithm in parallel
	                        n.cores = 2, #Number of cores to use in the parallel processing
	                        iters = NULL, #Number of iterations
													do.U = FALSE, #Perform clustering of SNPs
													display.iters = TRUE, #Display the iterations of function cnmtf

	                        #Randomisations
	                        score.pvalues = c(TRUE, FALSE, "only"), #Estimate p-values for the scores
	                        randomisations = 50, #Number of randomisations
	                        random.parallel = FALSE, #Run the radomisations in parallel

	                        #Construction of penalisation terms
													file.Gu = file.Gu #Workspace with adjancency matrix
	)
	{

	  msg <- paste("\n","\n", "Running experiment", name.exp, "\n", as.character(Sys.time()), "\n", sep=" "); cat(msg); cat(msg, file = log.file, append = F)



		#---------------------------------------
		# 1. Select data and parameters
		#---------------------------------------


				  #Print dimensions
		  		labelsU <- rownames(R)
				  labelsV <- colnames(R)
				  (n = nrow(R)); (m = ncol(R))
				  msg <- paste("Dimensions of input R:", n, m,"\n", as.character(Sys.time()), "\n", sep=" "); cat(msg); cat(msg, file = log.file, append = TRUE)



			#-----------------------------------------------------------
			# 1.1 Define number of SNP clusters
			#-----------------------------------------------------------


				  msg <- paste("\n", "Setting numbers of clusters by ", define.k, "\n", as.character(Sys.time()), "\n", sep=" "); cat(msg); cat(msg, file = log.file, append = TRUE)

				  if(define.k == "NMTF"){ #Define number of SNP clusters based on dispersion coefficient of connectivity matrix

				  	#Define list of parameters to evaluate if not provided
				  	#If any of the parameters is fixed, its list will be of length 1

				  	#Parameters to evaluate for k1
					    if(!is.null(k1)){
					    	lk1 = k1
					    	if(length(k1) == 1){
					    		do.U.k1 = FALSE
					    	}else{
					    		do.U.k1 = TRUE
					    	}
					    }else{
					    	stop("Please provide the number of clusters k1")
					    }
					    nk1 = length(lk1)


				    #Parameters to evaluate for k2
					  	 if(!is.null(k2)){
					    	lk2 = k2
					    	if(length(k2) == 1){
					    		do.V.k2 = FALSE
					    	}else{
					    		do.V.k2 = TRUE
					    	}
					    }else{
					    	stop("Please provide the number of clusters k2")
					    }
					 		nk2 = length(lk2)


				    #Set number of iterations/id, repetitions of the algorithm and parameters = 0
				    lparameters.0 = list(alpha1 = 0, alpha2 = 0, gamma1 = 0, gamma2 = 0)

				    #Initialize U and V matrices for the highest number of clusters to evaluate
				    if(init == 1){
				      res.init = initialise.UV( R = R, k1 = max(lk1), k2 = max(lk2), name.exp = name.exp,
				                                name.init = name.init, work.dat = work.dat, log.file = log.file)
				      U.init.max = res.init[[1]]
				      V.init.max = res.init[[2]]

				    }else{
				      msg <- paste("Random initialisation of matrices U and V","\n", as.character(Sys.time()), "\n", sep=" "); cat(msg); cat(msg, file = log.file, append = TRUE)
				      U.init = NULL
				      V.init = NULL
				    }


				    #Create table of results for rho coefficients at each setting of parameters
				    t.res = matrix(NA, nrow = length(lk1) * length(lk2), ncol = 5) #res.clust = list()
				    z = 1 #Counter variable

				    #Run the algorithm for each setting of SNP-Patient clusters k1 and k2
				    for(m in nk1:1){
				      for(q in nk2:1){

				      	msg <- paste("Evaluating setting: k1 = ", lk1[m], ", k2 = ", lk2[q], "\n", as.character(Sys.time()), "\n", sep=" "); cat(msg); cat(msg, file = log.file, append = TRUE)

				      	#Set the parameters and aux. imput matrices
				      	k = c(lk1[m], lk2[q])

				      	#Filter columns of initialisation matrices
				      	if(init == 1){
				      		U.init = U.init.max[,1:lk1[m]]
				      		V.init = V.init.max[,1:lk2[q]]
				      	}

				      	#Run the algorithm and extract rho parameters
				      	res.clust = consensus.clust(R=R,
				      															Wv=0, Wu=0,
				      															HLH=0, Vo=0,
				      															k=k, lparameters = lparameters.0,
				      															iters = iters,
				      															labelsU = labelsU, labelsV = labelsV,
				      															run.t = run.t.par, calcObj = calcObj, calcObj2 = calcObj2,
				      															init = init,
				      															parallel.opt = parallel.opt, n.cores = n.cores,
				      															do.U = do.U.k1, do.V = do.V.k2, do.O = FALSE,
				      															V.init = V.init, U.init = U.init, export.res = FALSE, display.iters = display.iters)

				      	t.res[z,] = c(res.clust[[1]], mean(res.clust[[1]], na.rm = TRUE),k)
				      	z = z + 1 #Counter variable
				      }
				    }

				    #Explore results
				    #Add column names
				    colnames(t.res) <- c("rho1","rho2","avrho","k1","k2")
				    t.res<-as.data.frame(t.res)
				    print(t.res)

				    #Plot rho
				    pdf(paste(work.dat,"k1k2_", name.exp,".pdf",sep=""))
				    plot(t.res$k1, t.res$rho1, col = "blue", pch = 2, las=1,xlab="Number of SNP clusters (k1)", ylab="Dispersion coefficient (Rho1)", ylim = c(0,1))
				    dev.off()

				    #Define optimum k1
				    decision.k1k2 = as.numeric(readline("Check plot and input optimum k1:") )
				    if(decision.k1k2 > 0){
				    	k1 = decision.k1k2
				    }else{
				    	msg <- paste("Wrong k1 parameter. Processed interrumpted.")
				    	return()
				    }

				    #Define optimum k2
				    decision.k1k2 = as.numeric(readline("Check plot and input optimum k2:") )
				    if(decision.k1k2 > 0){
				    	k2 = decision.k1k2
				    }else{
				    	msg <- paste("Wrong k2 parameter. Processed interrumpted.")
				    	return()
				    }

				    #Save results
				    save(list = c("t.res"),file = paste(work.dat,"set_k1k2_", name.exp,".RData",sep=""))


				  }else if(define.k == "PCA"){

				    #Define number of SNP clusters based on number of significant principal components of SNPs
				    res.pca = dudi.pca(data.frame(R), center = TRUE,  scale = TRUE, scannf = F, nf = 3)
				    plot(cumsum(res.pca$eig)/sum(res.pca$eig))

				    #Define optimum k1
				    decision.k1k2 = as.numeric(readline("Check plot and input optimum k1:") )
				    if(decision.k1k2 > 0){
				    	k1 = decision.k1k2
				    }else{
				    	msg <- paste("Wrong k1 parameter. Processed interrumpted.")
				    	return()
				    }

				    #Define number of Patient clusters based on number of significant principal components of patients
				    res.pca = dudi.pca(data.frame(t(R)), center = TRUE,  scale = TRUE, scannf = F, nf = 3)
				    plot(cumsum(res.pca$eig)/sum(res.pca$eig))

				    #Define optimum k2
				    decision.k1k2 = as.numeric(readline("Check plot and input optimum k1:") )
				    if(decision.k1k2 > 0){
				    	k2 = decision.k1k2
				    }else{
				    	msg <- paste("Wrong k1 parameter. Processed interrumpted.")
				    	return()
				    }


				  }

				  #Set k1 and k2
				  k = c(k1,k2)
				  msg <- paste("Number of clusters set to:", k1, k2, "\n", as.character(Sys.time()), "\n", sep=" "); cat(msg); cat(msg, file = log.file, append = TRUE)



			#-----------------------------------------------------------
			# 1.2 Define penalisation parameters
			#-----------------------------------------------------------


				  msg <- paste("\n", "Setting penalisation parameters ","\n", as.character(Sys.time()), "\n", sep=" "); cat(msg); cat(msg, file = log.file, append = TRUE)

				  #Parameters to be set according to weights

				  (set.wparameters = names(unlist(wparameters))[unlist(wparameters) != 0])


				  #Construct the SNP-SNP network and calculate its laplacian matrix
				  if( "alpha1" %in% set.wparameters ){

				  	#Load workspace with Wu and Gu
				  	if( file.exists(file.Gu)){
				  		msg <- paste("Workspace with SNP-SNP network was found.","\n", as.character(Sys.time()), "\n", sep=" "); cat(msg); cat(msg, file = log.file, append = TRUE)
				  		load(file = file.Gu)
				  	}else{
				  		msg <- paste("Workspace with SNP-SNP network was not found.","\n", as.character(Sys.time()), "\n", sep=" "); cat(msg); cat(msg, file = log.file, append = TRUE)
				  		stop()
				  	}

			      #Match SNPs of R and Wu
			  		pos.snps.R = match(rownames(R), rownames(Wu))

			  		#Check that lists of SNPs match
			  		if( sum( is.na(pos.snps.R) ) > 0 ){
			  			stop("Some SNPs of input Wu and R do not match.")
			  		}

			  		#Filter adjacency matrix to map SNPs in R
			  		Wu = Wu[pos.snps.R ,pos.snps.R]
			  		Gu = graph.adjacency(Wu, mode="undirected")

			  		#Save graph object
			  		save(list = c("Wu","Gu"),file = paste(work.dat,"Gu_", name.exp,".RData",sep=""))


				  	#Print dimensions and node degree
				  	msg <- paste("Dimensions of input Wu:", nrow(Wu), ncol(Wu),"\n", as.character(Sys.time()), "\n", sep=" "); cat(msg); cat(msg, file = log.file, append = TRUE)
				    #msg <- paste("Summatory of node degrees in Wu:", sum(igraph::degree(Gu)), "\n", as.character(Sys.time()), "\n", sep=" "); cat(msg); cat(msg, file = log.file, append = TRUE)

				    if( nrow(Wu) != nrow(R) & ncol(Wu) != nrow(R)){
				      stop("Dimensions of input Wu and R do not match.")
				    }

				    #Calculate its laplacian matrix
				    Du = diag(rowSums(Wu))
				    Lu = Du - Wu

				  }else{
				    Wu = NULL
				  }

				  #Construct the Patient-Patient  network and calculate its laplacian matrix
				  if( "alpha2" %in% set.wparameters ){

				    msg <- paste("Constructing the Patient-Patient network","\n", as.character(Sys.time()), "\n", sep=" "); cat(msg); cat(msg, file = log.file, append = TRUE)
				    Wv = matrix(0, nrow = m, ncol = m)
				    Dv = diag(rowSums(Wv))
				    Lv = Dv - Wv


				  }else{
				    Wv = NULL
				  }

				  #Calculate the Side information kernel matrix (L) and the centering matrix
				  if( "gamma1" %in% set.wparameters ){

				    msg <- paste("Calculating the Side information kernel matrix (L) and the centering matrix (H)","\n", as.character(Sys.time()), "\n", sep=" "); cat(msg); cat(msg, file = log.file, append = TRUE)
				    res.kernels = kernels.cnmtf(R,"unknown")#kernels.cnmtf(R,pop) #
				    H = res.kernels[[1]]
				    L = res.kernels[[2]]
				    HLH = H %*% L %*% H

				  }else{
				    HLH = NULL
				  }

				  #Construct the outcome matrix
				  if( "gamma2" %in% set.wparameters ){

				    if( is.null(out.cat) ){
				    	stop("The outcome vector is NULL.")
				    }
				  	msg <- paste("Calculating the prior-knowledge matrix for patients (Vo)","\n", as.character(Sys.time()), "\n", sep=" "); cat(msg); cat(msg, file = log.file, append = TRUE)
				    Vo = construct.Vo( out.cat, k[2])
				    msg <- paste("Dimensions of input Vo:", nrow(Vo), ncol(Vo),"\n", as.character(Sys.time()), "\n", sep=" "); cat(msg); cat(msg, file = log.file, append = TRUE)


				  }else{
				    Vo = NULL
				  }


				  #Initialisation of U and V matrices

				    if(init == 1){
				      res.init = initialise.UV( R = R, k1 = k1, k2 = k2, name.exp = name.exp,
				                                name.init = name.init, work.dat = work.dat, log.file = log.file)
				      U.init = res.init[[1]]
				      V.init = res.init[[2]]

				    }else{

				      msg <- paste("Random initialisation of matrices U and V","\n", as.character(Sys.time()), "\n", sep=" "); cat(msg); cat(msg, file = log.file, append = TRUE)
				      U.init = NULL
				      V.init = NULL
				    }


				  #Run funtion to estimate parameters sequentially

				  if( estimate.par == TRUE){

				  	lparameters = parameters.cnmtf( R = R, #Relationship matrix
				  																out.cat = out.cat, #Categorical outcome variable
				  																pop = pop, #Population variable
				  																k = k, #Number of clusters
				  																log.file = log.file,

				  																#Penalisation terms
				  																Wv = Wv, Wu = Wu,
				  																HLH = HLH, Vo = Vo,

				  																#iGraph objects
				  																Gu = Gu, #iGraph for Wu (SNP-SNP network)
				  																Gv = Gv, #iGraph for Wv (Patient-Patient network)

				  																#Initialisations of matrices
				  																V.init = V.init, U.init = U.init,

				  																#Variables to Save/Load data and workspaces
				  																name.exp = name.exp, #Name of experiment to save files
				  																name.init = name.init, #Name of workspace with initialisations of U and V
				  																work.dat = work.dat, #Folder to save and load workspaces

				  																#Penalisation parameters
				  																sequential.estimation = sequential.estimation, #Set the parameters in a specific order and carry the optimals
				  																wparameters = wparameters, #Weights of penalisation terms to be computed
				  																run.t.par = run.t.par, #Number of repetitions for parameters fitting
				  																range.parameters = range.parameters, #Range of parameters to be evaluated (if estimate.par == FALSE)

				  																#Variables to control performance of the algorithm
				  																max.try0 = max.try0,  #Maximum number of tries to fit the parameter
				  																calcObj = calcObj, calcObj2 = calcObj2 , #Check convergency each X number of iterations after first Y iterations
				  																init = init, #Type of seeding/initialisation of matrices in the algorithm
				  																parallel.opt = parallel.opt, #Run some instances of the algorithm in parallel
				  																n.cores = n.cores, #Number of cores to use in the parallel processing
				  																iters = iters, #Number of iterations
				  																display.iters = display.iters, #Display the iterations of function cnmtf

				  																#Gamma2 auxiliar plot of known associations
				  																tao.sc  = 4, #Trehsold of standard deviations for the SNV score
				  																snps.known = snps.known #List of known associaitons (P)



				  	)

				    msg <- paste("Optimum parameters found:", "\n", as.character(Sys.time()), "\n", sep=" "); cat(msg); cat(msg, file = log.file, append = TRUE)
				  	print(unlist(lparameters))



				  }else{
				  	msg <- paste("Optimum parameters defined by user in object wparameters.", "\n", as.character(Sys.time()), "\n", sep=" "); cat(msg); cat(msg, file = log.file, append = TRUE)
				  	lparameters <- range.parameters
				  }

				  if( save.parameters == TRUE ){
				  	#Save parameters to file
				  	save(list = c("lparameters"), file = paste( substr( path.exp, start = 1, stop = nchar(path.exp) - nchar(name.exp) - 1 )  , trait.project, "/", trait.project, "_parameters.RData",sep=""))
				  }


		#-----------------------------------------------------------
		# 2. Run the algorithms
		#-----------------------------------------------------------

			#-----------------------------------------------------------
			# 2.1 Algorithm evaluated on real data
			#-----------------------------------------------------------

				  if(score.pvalues == "only"){

				    load(file = paste(work.dat,"cnmtf_",name.exp,".RData",sep=""))

				  }else{


				    #Create vector of results
				    res.cnmtf = list()

				    #Set parameters for the algorithm
				    run.t = run.t.exp

				    msg <- paste("\n", "Running", run.t, "repetitions of the algorithm ","\n", as.character(Sys.time()), "\n", sep=" "); cat(msg); cat(msg, file = log.file, append = TRUE)


				    #Run algorithm for experiment
				    res.cnmtf[[1]] = consensus.clust(R=R,
				                                     Wv=Wv, Wu=Wu,
				                                     HLH=HLH, Vo=Vo,
				                                     k=k, lparameters=lparameters,
				                                     iters = iters,
				                                     labelsU=labelsU, labelsV=labelsV,
				                                     run.t = run.t, calcObj = calcObj, calcObj2 = calcObj2,
				                                     init = init, parallel.opt = parallel.opt, n.cores = n.cores,
				                                     set.alg = NULL, export.res = FALSE,
				                                     do.U = do.U , V.init = V.init, U.init = U.init, display.iters = display.iters)

				    #Save results of algorithm
				    save(list = c("res.cnmtf"),file = paste(work.dat,"cnmtf_",name.exp,".RData",sep=""))
				    msg <- paste("Results of factorisation printed to file:" , paste(work.dat,"cnmtf_",name.exp,".RData",sep=""), "\n", as.character(Sys.time()), "\n", sep=" "); cat(msg); cat(msg, file = log.file, append = TRUE)


				    #Check clustering of patients
				    if( !is.null(out.cat) ){
				    	table(out.cat,res.cnmtf[[1]][[3]][,2] )
				    }

				  }



			#-----------------------------------------------------------
			# 2.2 Randomisations of imputed data
			#-----------------------------------------------------------


				  #Estimate p-values for the scores

				  if(score.pvalues %in% c(TRUE,"only")){

				    msg <- paste("\n", "Estimating p-values for the scores  ","\n", as.character(Sys.time()), "\n", sep=" "); cat(msg); cat(msg, file = log.file, append = TRUE)

				    #Initialise vector of randomizations
				    #There is not need of runing the algorihtm multiple times per randomisations (run.t = 1)
				    #Default initialization for randomizations is 0
				    run.t = 1
						iters = 120

				    #Run algorithm on randomisations

				    if(random.parallel == TRUE){

				    		 msg <- paste("Number of randomizations in PARALLEL: ", randomisations, "\n", as.character(Sys.time()), "\n", sep=" "); cat(msg); cat(msg, file = log.file, append = TRUE)
					       msg <- paste("Running randomizations in PARALLEL","\n", as.character(Sys.time()), "\n", sep=" "); cat(msg); cat(msg, file = log.file, append = TRUE)

					      #Register the number of cores
					      doMC::registerDoMC(cores = n.cores)

					      #Run each randomization of the experiment in parallel
					      l.cnmtf.ran <- foreach(ran = 1:randomisations) %dopar% {

					      			msg <- paste("\n", "Randomisation ", ran, "\n", as.character(Sys.time()), "\n", sep=" "); cat(msg); cat(msg, file = log.file, append = TRUE)

					      	  	#Randomize the R matrix
							        Rr <- R[, sample(1:m, size = m, replace = FALSE)]


							        #Perform consensus clustering
							        res.alg.ran =  consensus.clust(R=Rr,
							                                       Wv=Wv, Wu=Wu,
							                                       HLH=HLH, Vo=Vo,
							                                       k=k, lparameters=lparameters,
							                                       iters = iters,
							                                       labelsU=labelsU, labelsV=labelsV,
							                                       run.t = run.t, calcObj = calcObj, calcObj2 = calcObj2,
							                                       init = 0, parallel.opt = FALSE, #This option refers to performing repetititions of algorithm in parallel. As this is only 1 repetition then this is set to series.
							                                       set.alg = NULL, set.clus.V = res.cnmtf[[1]][[3]],
							                                       do.U = FALSE, do.V = FALSE, export.res = FALSE, display.iters = display.iters)

				      	} #End parallel

				    }else{

				    	msg <- paste("Number of randomizations in SERIES: ", randomisations, "\n", as.character(Sys.time()), "\n", sep=" "); cat(msg); cat(msg, file = log.file, append = TRUE)
				      msg <- paste("Running randomisations in SERIES","\n", as.character(Sys.time()), "\n", sep=" "); cat(msg); cat(msg, file = log.file, append = TRUE)

				      #Create object of results
				      l.cnmtf.ran = rep(list(NULL),length(randomisations))

				      for(ran in 1:randomisations){

				        msg <- paste("Randomisation ", ran, "\n", as.character(Sys.time()), "\n", sep=" "); cat(msg); cat(msg, file = log.file, append = TRUE)

				        #Randomize the R matrix
				        Rr <- R[, sample(1:m, size = m, replace = FALSE)]

				        #Perform consensus clustering
				        l.cnmtf.ran[[ran]] = consensus.clust(R=Rr,
				                                             Wv=Wv, Wu=Wu,
				                                             HLH=HLH, Vo=Vo,
				                                             k=k, lparameters=lparameters,
				                                             iters = iters,
				                                             labelsU=labelsU, labelsV=labelsV,
				                                             run.t = run.t, calcObj = calcObj, calcObj2 = calcObj2,
				                                             init = 0, parallel.opt = FALSE, #This option refers to performing repetititions of algorithm in parallel. As this is only 1 repetition then this is set to series.
				                                             set.alg = NULL, set.clus.V = res.cnmtf[[1]][[3]],
				                                             do.U = FALSE, do.V = FALSE, export.res = FALSE, display.iters = display.iters)
				      } #End for
				    }#End if random.parallel




				    #Save results of the algorithm
				    save(list = c("l.cnmtf.ran"),file = paste(work.dat,"randomizations_cnmtf_", name.exp,".RData",sep=""))
				    msg <- paste("Results of randomisations printed to file:" , paste(work.dat,"randomizations_cnmtf_", name.exp,".RData",sep=""), "\n", as.character(Sys.time()), "\n", sep=" "); cat(msg); cat(msg, file = log.file, append = TRUE)



				  }#End if score.pvalues



	}#End function

