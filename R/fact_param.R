###############################################################
# cNMTF
#	2. Factorisation functions
# 2.4 Functions to find optimal penalization parameters
#
# Corresponding author:
# Luis Leal, Imperial College London, email: lgl15@imperial.ac.uk
###############################################################


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#' @title Optimal penalization paramaters
#' @description Function to estimate the optimal parameters for the cNMTF algorithm.
#'
#[Input]
#' @param R Relationship matrix
#' @param out  Categorical outcome variable
#' @param pop Population variable
#' @param k Number of clusters
#' @param log.file Log file to track status of the function
#'
#Penalisation terms
#' @param Wv,Wu Adjacency matrices to guide clustering of U and V
#' @param HLH Term to correct for population stratification

#iGraph objects
#' @param Gu iGraph for Wu (SNP-SNP network)
#' @param Gv iGraph for Wv (Patient-Patient network)

#Initialisations
#' @param U.init,V.init Initialisations for matrices U and V

#Variables to Save/Load data and workspaces
#' @param name.exp Name of experiment to save files
#' @param name.init Name of workspace with initialisations of U and V
#' @param work.dat Folder to save and load workspaces

#Penalisation parameters
#' @param sequential.estimation Set the parameters in a specific order and carry the optimals
#' @param wparameters Weights of penalisation terms to be computed
#' @param run.t.par Number of repetitions for parameters fitting
#' @param range.parameters Range of parameters to be evaluated (if estimate.par == FALSE)

#Variables to control performance of the algorithm
#' @param max.try0 Maximum number of tries to fit the parameter
#' @param calcObj,calcObj2 Check convergency each X number of iterations after first Y iterations
#' @param init  Type of seeding/initialisation of matrices in the algorithm
#' @param parallel.opt Run some instances of the algorithm in parallel
#' @param n.cores Number of cores to use in the parallel processing
#' @param iters Number of iterations
#' @param display.iters Display the iterations of function cnmtf

#Auxiliar plots of known associations
#' @param tao.sc  Trehsold of standard deviations for the SNV score
#' @param snps.known List of known associations
#'
#[Output]
#' @return \code{opt.parameters} list of optimal parameters fitted: \code{gamma3, gamma2, gamma3}
#'
#' @md
#' @author Luis G. Leal, \email{lgl15@@imperial.ac.uk}
#' @family Factorisation functions
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

	parameters.cnmtf = function(  R, #Relationship matrix
	                              out = NULL, #Categorical outcome variable
	                              pop = NULL, #Population variable
	                              k, #Number of clusters
																log.file = NULL,

	                              #Penalisation terms
	                              Wv = NULL, Wu = NULL,
	                              HLH = NULL, Vo = NULL,

																#iGraph objects
																Gu = NULL, #iGraph for Wu (SNP-SNP network)
																Gv = NULL, #iGraph for Wv (Patient-Patient network)

																#Initialisations
																V.init = NULL, U.init = NULL,

	                              #Variables to Save/Load data and workspaces
	                              name.exp = NULL, #Name of experiment to save files
	                              name.init = NULL, #Name of workspace with initialisations of U and V
	                              work.dat = NULL, #Folder to save and load workspaces

	                              #Penalisation parameters
																sequential.estimation = FALSE, #Set the parameters in a specific order and carry the optimals
	                              wparameters = NULL, #Weights of penalisation terms to be computed
	                              run.t.par = NULL, #Number of repetitions for parameters fitting
	                              range.parameters = NULL, #Range of parameters to be evaluated (if estimate.par == FALSE)

	                              #Variables to control performance of the algorithm
	                              max.try0 = NULL,  #Maximum number of tries to fit the parameter
	                              calcObj = NULL, calcObj2 = NULL, #Check convergency each X number of iterations after first Y iterations
	                              init = NULL, #Type of seeding/initialisation of matrices in the algorithm
	                              parallel.opt = FALSE, #Run some instances of the algorithm in parallel
	                              n.cores = 2, #Number of cores to use in the parallel processing
	                              iters = NULL, #Number of iterations
																display.iters = TRUE, #Display the iterations of function cnmtf

																#Gamma2 auxiliar plot of known associations
																tao.sc  = 4, #Treshold of standard deviations for the SNV score
																snps.known = NULL #List of known associaitons (P)


	)

	{



		#Set the penalisation parameters to 0
		opt.parameters = list(gamma1 = 0, gamma3 = 0, gamma2 = 0)


			#Check if tracking clusters on PCs space
			track.pca = FALSE
			if(track.pca == TRUE){
				#Perform PCA
				res.pca = dudi.pca(data.frame(t(R)), center = TRUE,  scale = TRUE, scannf = F, nf = 3)

				#Calculate explained variance
				evar = round( res.pca$eig / sum(res.pca$eig) * 100, 0)

				#Labels and colours
				axes.pca = c(1,2)
				dpca = as.data.frame( cbind(res.pca$li[,axes.pca], pop))
				names(dpca) <- c("Axis1", "Axis2")
				xlab = as.character(paste("PC", axes.pca[1], " (EV = ", evar[1], "%)", sep = "" ))
				ylab = paste("PC", axes.pca[2], " (EV = ", evar[2], "%)", sep = "" )

			}


		  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	    #1. Fit parameter ranges (min,max)
	    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

	    #If there is not range of parameters to be evaluated, the algorithm will try to fit a range for each of them

		   if(is.null(range.parameters)){

	    				msg <- paste("\n","\n","Fitting range of parameters.","\n", as.character(Sys.time()), "\n", sep=" "); cat(msg); cat(msg, file = log.file, append = TRUE)

	            #Run the algorithm with parameters set to 0
	            res.param.0 = consensus.clust(R = R,
	                                          Wv = Wv, Wu = Wu,  #Providing all the penalisation terms to estimate maximum terms value
	                                          HLH = HLH, Vo = Vo,
	                                          k=k, lparameters = opt.parameters,
	                                          iters = iters,
	                                          labelsU = rownames(R), labelsV=colnames(R),
	                                          run.t = run.t.par, calcObj = calcObj, calcObj2 = calcObj2,
	                                          init = init, parallel.opt = parallel.opt, n.cores = n.cores,
	                                          V.init = V.init, U.init = U.init,
	                                          do.O = FALSE, do.U = FALSE, do.V = FALSE, export.res = FALSE,
	            															estimate.penalisation = TRUE, display.iters = display.iters)

	            range.parameters =  res.param.0[[9]]


	            msg <- paste("Range of parameters fitted to:", range.parameters, "\n", as.character(Sys.time()), "\n", sep=" "); cat(msg); cat(msg, file = log.file, append = TRUE)

	    }



		#Set the parameters in the following order # param.i = "gamma2"

			for(param.i in c("gamma2","gamma3","gamma1")){

							#Extract the weight (information contribution) for this parameter
							ww = as.numeric(wparameters[names(wparameters) %in% param.i])

							#Check if the weight is zero then goes to next parameter
							if( ww == 0){
								next
							}


							msg <- paste("\n", "Setting", param.i, "\n", as.character(Sys.time()), "\n", sep=" "); cat(msg); cat(msg, file = log.file, append = TRUE)





							#Check if optimum parameters were provided in the imput:
							range.i = unlist(range.parameters[names(range.parameters) %in% 	param.i ])

	    				if( length(range.i) == 1 ){

		    					msg <- paste("Optimum value for", param.i, "provided by user:", range.i, "\n", as.character(Sys.time()), "\n", sep=" "); cat(msg); cat(msg, file = log.file, append = TRUE)


		    					#If only one value is provided in the range then its taken as the optimum
		    					opt.parameters[ names(range.parameters) %in% 	param.i ] = range.i


		    					#Set gamma3 to be maximum 1 so it doesnt mess up the objective function
		    					if( param.i == "gamma3" ){

		    						msg <- paste("Optimum value for", param.i, "must be maximum 1. Setting this parameter to", min(1, range.i), "\n", as.character(Sys.time()), "\n", sep=" "); cat(msg); cat(msg, file = log.file, append = TRUE)
		    						opt.parameters[ names(range.parameters) %in% 	param.i ] <- min(1, range.i)

		    					}


		    					next

	    				}else if( length(range.i) > 2){
		    					msg <- paste("Wrong range of values for the parameter", "\n", as.character(Sys.time()), "\n", sep=" "); cat(msg); cat(msg, file = log.file, append = TRUE)
		    					stop()
	    				}


							#Check if the weight is one then takes the maximum value in the range as the optimum
							if( length(range.i) > 1 & ww == 1){

									#Select the maximum value
									opt.parameters[ names(range.parameters) %in% 	param.i ] = max(range.i)


									#Set gamma3 to be maximum 1 so it doesnt mess up the objective function
									if( param.i == "gamma3" ){

										msg <- paste("Optimum value for", param.i, "must be maximum 1. Setting this parameter to", min(1, max(range.i)), "\n", as.character(Sys.time()), "\n", sep=" "); cat(msg); cat(msg, file = log.file, append = TRUE)
										opt.parameters[ names(range.parameters) %in% 	param.i ] <- min(1, max(range.i))

										next
									}


									#Show optimum value in screen
									msg <- paste("Optimum", param.i, "was found to be", opt.parameters[ names(range.parameters) %in% 	param.i ] , "\n", as.character(Sys.time()), "\n", sep=" "); cat(msg); cat(msg, file = log.file, append = TRUE)
									next

							}

							#Number of parameter values to evaluate
							max.try <- max.try0

							#For gamma3 the range is predefined between 0 and 1
							if( param.i == "gamma3" ){
								range.i[1] <- 0
								range.i[2] <- 1
							}

							#Define order of magnitud of minimum and maximum value in the range
							mag.low = floor( log10( range.i[1] ))
							mag.up = floor( log10( range.i[2] ))

							#Set initial list of parameters to evaluate
							vparam.i = c(range.i[1], range.i[2],  range.i[1] + 1*10^seq( mag.up - 1, mag.up - max.try + 2 ) [ 1: (max.try - 2) ])
							vparam.i = unique( vparam.i[ vparam.i <= range.i[2]] )


							#Create vectors to save average clustering coefficients / average node degree and number of known associations retrieved
							y = y2 = rep(NA, max.try)
							y2.name = "Total known associations retrieved"


							#Set a marker variable to update the list of parameters to evaluate and a sub-range of parameter values to fit a linear model
							updated.vparam = FALSE
							range.fit = NULL

							#Check the weight for gamma2
							ww.gamma2 = as.numeric(wparameters[names(wparameters) %in% "gamma2"])

							#Set name of tracking variable and wheteher or not compute clustering of SNPs
							#Also set the maximum of the tracking variable

							if(param.i %in% c("gamma2")){

									y.name = "Similarity Cluster-Outcome"
									y.max = 1
									do.U.par = FALSE

							}else if(param.i %in% c("gamma1")){

										y.name = "Total node degree"
										y.max = sum(igraph::degree(Gu)) #Sum up node degrees in the graph
										do.U.par = TRUE


	            }else if(param.i %in% c("gamma3")){

	            		#If gamma2 was set, then the tracking variable is the simVO
	            		if( ww.gamma2  == 0) { #<---Check this

	            			y.name = "Similarity Cluster-Outcome"
	            			y.max = 1
	            			do.U.par = FALSE

	            		}else{ #If gamma2 was not set, then the tracking variable is the simVL

	            			y.name = "Similarity Cluster-Population"
	            			y.min = 0
	            			y.max = 1
	            			do.U.par = FALSE
	            		}


	            }

							#Print messages
							msg <- paste("Max", y.name, ":", y.max, "\n", as.character(Sys.time()), "\n", sep=" "); cat(msg); cat(msg, file = log.file, append = TRUE)
							msg <- paste("List of parameters", param.i, "to evaluate: ", paste( vparam.i, collapse = ", "),  "\n", as.character(Sys.time()), "\n", sep=" "); cat(msg); cat(msg, file = log.file, append = TRUE)

	            #If parameters are setted in a sequence or independently
	            if(sequential.estimation == TRUE){
	              #Declare list of penalisation parameters equal to the optimal parameters at this setting
	              lparameters.pp = 	opt.parameters
	            }else{
	              #All penalisation parameters are equal to 0 and only the parameter pp is updated in the loop
	              lparameters.pp = list(gamma1 = 0, gamma3 = 0, gamma2 = 0)
	            }


	            #Set counter variable of tries
	            pp = 1

	            #Set drawing plot
	            par( mfrow = c(1,1) )

	            #Run algorithm with each parameter:

	            while( pp <= max.try){

	                  msg <- paste("Setting", param.i, "to", vparam.i[pp], "\n", as.character(Sys.time()), "\n", sep=" "); cat(msg); cat(msg, file = log.file, append = TRUE)
	            			lparameters.pp[ names(range.parameters) %in% 	param.i ] <- vparam.i[pp]


	                  #Run algorithm for experiment
	                  res.param = consensus.clust(R = R,
	                                              Wv = Wv, Wu = Wu,
	                                              HLH = HLH, Vo = Vo,
	                                              k = k, lparameters = lparameters.pp,
	                                              iters = iters,
	                                              labelsU = rownames(R), labelsV=colnames(R),
	                                              run.t = run.t.par, calcObj = calcObj, calcObj2 = calcObj2,
	                                              init = init, parallel.opt = parallel.opt, n.cores = n.cores,
	                  														do.U = do.U.par, V.init = V.init, U.init = U.init,
	                                              do.O = TRUE, export.res = FALSE, display.iters = display.iters)



	                  #Tracking gamma3

	                  if(param.i == "gamma3"){

	                  	#If gamma2 was set previously the tracking variable is simVO
	                  	if(  ww.gamma2 == 0){ #<--- Check this

	                  		#Calculate similarity cluster-population
	                  		cat("SimVL", compare(as.numeric(as.factor(pop)) [pop!="Unknown"], as.numeric(res.param[[3]][,2]) [pop!="Unknown"], method  =  c("nmi")), "\n")

	                  		#Calculate similarity cluster-outcome
	                  		y[pp] = compare(as.numeric(as.factor(out)), as.numeric(res.param[[3]][,2]), method  =  c("nmi"))

	                  	}else{ #If only seting gamma3 the tracking variable is simVL

	                  		#Calculate similarity cluster-population
	                  		y[pp] = compare(as.numeric(as.factor(pop)) [pop!="Unknown"], as.numeric(res.param[[3]][,2]) [pop!="Unknown"], method  =  c("nmi"))


	                  	}


	                  }

	                  #Tracking gamma1

	                  #Calculate average clustering coefficient / average node degree at each SNP cluster
	                  if(param.i == "gamma1"){

												#Define list variable to save the node degree per cluster of SNPs
			                  lkd =  list()

			                  #Fill the list variable
			                  for(i in 1:k[1]){
			                  	Gu.i = induced_subgraph(graph = Gu, v = res.param[[2]][,2] == i)
			                  	lkd[i] = sum(igraph::degree(Gu.i))
			                  }

			                  #Sum aup the graph variable across SNP clusters
			                  y[pp] = sum(unlist(lkd))

	                  }


	                  #Tracking gamma2
	                  if(param.i == "gamma2"){

	                  	#Calculate similarity cluster-population
	                  	cat("SimVL", compare( as.numeric(as.factor(pop)) [pop!="Unknown"], as.numeric(res.param[[3]][,2]) [pop!="Unknown"], method  =  c("nmi")), "\n")

	                  	#Calculate similarity cluster-outcome
	                  	y[pp] = compare( as.numeric(as.factor(out)), as.numeric(res.param[[3]][,2]), method  =  c("nmi"))


	                  }


	                  #For all parameters an extra tracking variable is the number of known associations retrieved


		                  #Delta score for a specific combination of patient clusters
		                  Onp = scale( get("On",res.param) )
		                  dOnp = Onp[,2] - Onp[,1]

		                  #Calculate desviations from SD
		                  sd.sc = (dOnp - mean(dOnp))/sd(dOnp)

		                  #Find significant SNVs at a given treshold
		                  sig.sc = rownames( Onp ) [ abs( sd.sc ) > tao.sc ]
		                  y2[pp] = sum(sig.sc %in% snps.known)

		                  #Print number of significant SNVs
		                  msg <- paste("Number of associations retrieved", ":", paste( round(length(sig.sc), digits = 2), collapse = ", "), "\n", as.character(Sys.time()), "\n", sep=" "); cat(msg); cat(msg, file = log.file, append = TRUE)
		                  msg <- paste("Number of known associations retrieved", ":", paste( round(y2, digits = 2), collapse = ", "), "\n", as.character(Sys.time()), "\n", sep=" "); cat(msg); cat(msg, file = log.file, append = TRUE)


		                  #Plot parameter
		                  #plot.parameter(x = vparam.i, y = y2, logarithm.x = TRUE, y.max = length(snps.known), y.name = "Known associations retrieved (TP)",
		                  #							 x.name = param.i, fitmodel = F, print.file = NULL)


	                  #Update list of vparam.i to evaluate
	                  if(pp > 1){

		                    #Minimum value for the tracking variable
		                  		y.min = y[1]

		                  	#If gamma2 was not set previously the tracking variable is simVL
		                  		if(param.i %in% c("gamma3") & ww.gamma2 == 0){
		                  			y.max = y[1]
		                  			y.min = 0
		                  		}

		                		#Define optimal value of tracking variable
		                			y.opt = (y.max - y.min)*ww + y.min

	                	  	#If tracked variable is lower than the optimum, or if it is the last iteration, then update list of vparam.i and max.try
		                  	if( (y[pp] < y.opt | pp == max.try) & updated.vparam == FALSE ){

		                  		#Update list of param to evaluate in following tries:
		                  		max.try2 = pp + 3 #New max.try
		                  		vparam.i = 	c(vparam.i[1:pp], rep(NA,3)) #Add 3 more iterations to find the optimum

		                  		if( pp < max.try | y[pp] < y.opt ){ #If not reached first max.try: Makes a sequence between last two vparam.i

		                  			vparam.i[ (pp+1) : max.try2] = seq( vparam.i[pp-1], vparam.i[pp],length.out = max.try2 - pp + 2)[ 2 : (max.try2 - pp + 1)]
		                  			range.fit = c(vparam.i[pp-1], vparam.i[pp])

		                  		}else{ #If reached first max.try: Makes a sequence between last vparam.i and 0
		                  			vparam.i[ (pp+1) : max.try2] = seq( vparam.i[pp], 0,length.out = max.try2 - pp + 2)[ 2 : (max.try2 - pp + 1)]
		                  			range.fit = c(vparam.i[pp], 0)
		                  		}

		                  		max.try <- max.try2 #Update max.try
		                  		updated.vparam = TRUE #Update marker variable to dont do this procedure again


		                  		msg <- paste("Updated list of parameter", param.i,"to evaluate:",   paste( round(vparam.i, digits = 2), collapse = ", "), "\n", as.character(Sys.time()), "\n", sep=" "); cat(msg); cat(msg, file = log.file, append = TRUE)
		                  		msg <- paste("Sub-range to fit", param.i,":",  range.fit , "\n", as.character(Sys.time()), "\n", sep=" "); cat(msg); cat(msg, file = log.file, append = TRUE)

		                  	}

	                  }

	                  #Update number of tries
	                  msg <- paste("Tracked variable", y.name, ":", paste( round(y, digits = 2), collapse = ", "), "\n", as.character(Sys.time()), "\n", sep=" "); cat(msg); cat(msg, file = log.file, append = TRUE)
	                  pp = pp + 1





	            }



	            #Sort vectors
		            y = y[ order(vparam.i)]
		            y2 = y2[ order(vparam.i)]
		            vparam.i = vparam.i[ order(vparam.i) ]


	            #Save workspace for this parameter
		            save(list = c("vparam.i","param.i","y","snps.known" ,"y.max","y.name","range.fit","log.file", "ww", "name.exp","y2"),file = paste(work.dat, param.i,"_", name.exp,".RData",sep=""))


	            #Estimating optimal parameters and ploting changes in parameter versus tracking variable
		            print(vparam.i)
		            print(y)
		            (opt.param.i = plot.parameter(x = as.numeric(vparam.i), y = y, y.max = y.max, y.name = y.name, range.fit = if( param.i == "gamma3") { c(range.fit[2], range.fit[1]) }else{ range.fit }, log.file = log.file,
	  	                                            x.name = param.i, fitmodel = T, ww = ww, print.file = paste(work.dat,"parameter_", param.i ,"_", name.exp,".pdf",sep="")))

		          #For all the parameters creates an auxiliar plot with the number of known associations retrieved
	             	#plot.parameter(x = vparam.i, y = y2, logarithm.x = TRUE, y.max = length(snps.known), y.name = "Known associations retrieved (TP)",
	            	#							 x.name = param.i, fitmodel = F, print.file = paste(work.dat,"parameter_", param.i ,"_", name.exp,"_known.pdf",sep=""))


		          #Show the optimum param.i to user:
		            msg <- paste("Optimum", param.i, "was found to be", opt.param.i, "\n", as.character(Sys.time()), "\n", sep=" "); cat(msg); cat(msg, file = log.file, append = TRUE)

		            if( opt.param.i < 0 ){
		            	msg <- paste("Optimum", param.i, "can not be negative. Using minimum value tested in the optimisation (different than zero) :", min( vparam.i [vparam.i != 0] ), "\n", as.character(Sys.time()), "\n", sep=" "); cat(msg); cat(msg, file = log.file, append = TRUE)
		            	opt.parameters[ names(range.parameters) %in% 	param.i ] <- min( vparam.i [vparam.i != 0] )
		            }else{
		            	opt.parameters[ names(range.parameters) %in% 	param.i ] <- opt.param.i
		            }




					}

			return( opt.parameters )

	}



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#' @title Plot penalization parameters
#' @description Function to plot changes in parameter versus a tracking variable
#'
#[Input]
#' @param x Vector of parameter values
#' @param x.name Name of the parameter
#' @param y Tracking variable
#' @param y.max Maximum value of the tracking variable
#' @param y.name Name of tracking variable
#' @param	fitmodel Fit a linear model or a cubic spline on the data
#' @param	ww Information weigth to find optimal point between y.min and y.max
#' @param	print.file Print plots to file
#' @param	range.fit Range of values of x to fit the model
#' @param	log.file Log.file to print results
#' @param	logarithm.x Plot x axis in logarithmic scale

#[Output]
#' @return \code{opt.param.i} Vector of optimal parameters if \code{fitmodel == TRUE}, otherwise the function returns only a plot of both variables
#'
#' @md
#' @author Luis G. Leal, \email{lgl15@@imperial.ac.uk}
#' @family Factorisation functions
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


	plot.parameter = function(x, y, #Parameter and tracking variable
	                          y.max, #Maximum limit of tracking variable
	                          y.name, #Name of tracking variable
	                          x.name, #Name of parameter
	                          fitmodel = FALSE, #Fit a linear model or a cubic spline on the data
	                          ww, #Information weigth to find optimal point between y.min and y.max
														print.file = NULL, #Print plots to file
														range.fit = NULL, #Range of values of x to fit the model
														log.file = NULL, #Log.file to print restuls
														logarithm.x = FALSE #Plot x axe in logarithmic scale
	)
	  {


			#Define plot parameter for the x axes:
					if(logarithm.x == TRUE){
						log.x = "x"
						x.plot = x[ x!=0 ]; y.plot = y[ x!=0 ]
					}else{
						log.x = ""
						x.plot = x; y.plot = y
					}

	    #Plot the estimations of parameter
				#Open connection to print to file
				if(!is.null(print.file)){
					pdf(print.file, width = 6, height = 6)
					par(mar = c(6.5, 6.5, 0.5, 0.5), mgp = c(4, 1, 0))
				}

		    #plot(x.plot, y.plot, col="black", log = log.x , las = 1, ylab = y.name, xlab = x.name, ylim = c(min(y, na.rm = T),y.max))

		  #Fit a model or return NULL
	    if(fitmodel == FALSE){
	    	return(NULL)
	  	}

	    #Define x datapoints to be included in the model
		  if(!is.null(range.fit)){
		  	included = which(x <= range.fit[1] & x >= range.fit[2])

		  }else{
		  	included = which(x != 0 & y != y.max & y != 0)
		  }
		  included = included[ !is.na(included)]  #Remove NA due to second term: which( x <= range.fit[1] & y == y.max )[1]
	    msg <- paste("Including the following points in the linear model:",x.name, paste(x[included],collapse=", "),"and", y.name, paste(y[included],collapse=", "),"\n"); cat(msg)

	    #Fit a linear model or a cubic spline on the data
	    fit = lm(y [ included ] ~ x[ included ] ) # smooth.spline(x, y, df = 2 )

	    #Calculate the minimum tracking  variable and the optimal tracking variable (y.opt)
	    y.min = y[ x == 0 ] # y.min = predict(fit, x = 0)$y #y.min = y[ x == 0 ]
	    y.opt = (y.max - y.min)*ww + y.min

	    #Define estimations of y at different alphas using the model and plot the fitted
	    new.x = seq(min(x[ included ]),max(x[ included ]),length = 1000)
	    new.y =  new.x *  fit$coefficients[[2]]  + fit$coefficients[[1]]
	    plot(x.plot, y.plot, col="black", log = log.x, las = 1, ylab = y.name, xlab = x.name, ylim = c(min(y,y.min, na.rm = T),y.max))
	    lines(new.x, new.y, lty = 2, col = "red")  # lines(fit, lty = 2, col = "red") #w = c(1,rep(0.1,length(x)-1))

	    #Add lines showing maximum, minimum and medium node degree
	    abline(h = c(y.max, y.min, y.opt), lt = 2, col = "blue")

	    #Close file connection
	    if(!is.null(print.file)){
	    	dev.off()
	    }

	    #Find the optimum parameter at y.opt
	    (opt.param.i = (y.opt -  fit$coefficients[[1]] ) / fit$coefficients[[2]] )

	    #Return optimum parameter
	    return(opt.param.i)

  }


