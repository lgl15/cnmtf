###############################################################
# cNMTF
#	3. Delta score functions
# 3.2 Functions to calculate delta scores and their p-values
#
# Corresponding author:
# Luis Leal, Imperial College London, email: lgl15@imperial.ac.uk
###############################################################


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#' @title Delta score between clusters
#' @description Auxiliar function to calculate delta score between clusters of patients.
#'
#[Input]:
#' @param file.exp1 workspace of the experiment (created with function score.cnmtf)
#' @param file.exp2 Optional. Workspace of a second experiment.
#' @param snps.known List of known associations
#' @param snps.known1 List of known associations
#' @param snps.known2 Second list of known associations
#' @param print.file File to print plots
#' @param sig.snp.nmtf List of significant SNVs in this experiment
#' @param clus.a,clus.b List of clusters of patients to compute delta scores. Clusters in \code{clus.a} are compared with clusters in \code{clus.b}.

#'
#[Output]:
#' @return \code{dOn} Delta Standardised Omega matrix with the delta scores for each pair of clusters.
#'
#' @md
#' @author Luis G. Leal, \email{lgl15@@imperial.ac.uk}
#' @family Scoring functions
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


		delta.score = function(file.exp1, file.exp2 = NULL, #Workspace of the experiments 1 and 2 (created with function score.cnmtf)
													 snps.known1, #List of known associations
													 snps.known2, #Second list of known associations
												 	 print.file = NULL, #File to print plots
													 sig.snp.nmtf = NULL, #List of significant SNVs with cNMTF in this experiment
													 clus.a = list(c(1)), clus.b = list(c(2)) #List of patient clusters to use in the delta scores

		)
		{

			#-------------------------------------------------
			# 1. Load and correct results
			#-------------------------------------------------


					#Load experiment 1
					cat("Reading results of cNMTF experiment 1", "\n", as.character(Sys.time()), "\n")

					if(file.exists(file.exp1) ){
						load(file = file.exp1)
					}else{
						stop("Workspace of experiment not found")
					}


					#Load experiment 2
					if( !is.null(file.exp2) ){

						#Rename object with results
						res.cnmtf1 <- res.cnmtf

						cat("Reading results of cNMTF experiment 2", "\n", as.character(Sys.time()), "\n")
						if(file.exists(file.exp2) ){
							load(file = file.exp2)
						}else{
							stop("Workspace of experiment not found")
						}

						res.cnmtf2 <- res.cnmtf

						#Recalculate omega matrix on results of experiment 1 using patients clusters of experiment 2

						cat("Correcting results of cNMTF experiment 1 to match cluster conformation of experiment 2", "\n", as.character(Sys.time()), "\n")

						#Create matrix with number of columns equals number of final patient clusters
						clus.V2 = res.cnmtf2[[1]][[3]]
						Or1 = res.cnmtf1[[1]][[10]]
						n.clus.pat = max(as.numeric(clus.V2[,2]))
						n = nrow(Or1)
						On1 = matrix(0, nrow = n, ncol = n.clus.pat)

						#Fill the matrix
						for(x in 1:n){
							for(z in 1:n.clus.pat){
								On1[x,z] <- median(Or1[x, clus.V2[,2] == z])
							}
						}

						#Add rownames and colnames to consensus omega
						On1 = scale(On1)
						colnames(On1) <- paste("p",1:ncol(On1),sep="")

						#Extract Omega matrix from results
						On2 <- res.cnmtf2[[1]][[4]]
						rownames(On1) <- rownames(On2)

						#Scale Omega matrix
						On2 = scale(On2)
						cat("Dimensions of Omega matrix experiment 2:", dim(On2), "\n")

						#Dispersion plot of node degrees and SNV scores per cluster
						On.snps = rownames(On1)
						dOn = matrix(NA, nrow = nrow(On1), ncol = ncol(On1))
						color = rep("gray",nrow(On1))
						pos.known = which(On.snps  %in% snps.known)

						#Open connection to print PDF
						if(!is.null(print.file)){
							pdf(print.file, width = 8, height = 6)
						}

						#Plot delta omega
						par(mfrow = c(1,2))
						for(i in 1:ncol(On2)){
							dOn[,i] = On2[,i] - On1[,i]
							plot(dOn[,i], ylab = "Delta SNV score", xlab = "Index", col = color)
							points(x = pos.known, dOn[pos.known,i], pch = 16)
							#plot(dOn[,i] , kd.snps[match(On.snps, R.snps)], xlab = "Delta SNV score (SD)", ylab = "Node degree")
						}


					}else{ #Delta score between clusters of the same experiment

						#Extract Omega matrix from results
						On <- res.cnmtf[[1]][[4]]
						On = scale(On)

						#Open connection to print PDF
						if(!is.null(print.file)){
							pdf(print.file, width = 8, height = 6)
						}

						#Calculate delta omega
						for(i in 1:length(clus.a) ){

							cat("Plotting delta omega", "\n", as.character(Sys.time()), "\n")


							dOn = On[, clus.a[[i]] ] - On[, clus.b[[i]] ]


							#Define color and positions of known associations
							color = rep("gray",nrow(On))
							On.snps = rownames(On)
							pos.known1 = which(On.snps  %in% snps.known1)
							pos.known2 = which(On.snps  %in% snps.known2)
							pos.sig = which(On.snps  %in% sig.snp.nmtf [[ i ]] ) #List of significant SNVs with this method


							#Plot delta omega
							plot(dOn, las = 1, ylab = "Delta SNV score", xlab = "Index", col = color)
							points(x = pos.known1, dOn[pos.known1], pch = 16, col = "black")
							points(x = pos.known2, dOn[pos.known2], pch = 16, col = "green")
							points(x = pos.sig, dOn[pos.sig], pch = 16, col = "blue")

							abline(h = 0, lt = 2)


						}

						}

					#Close PDF file connection
					if(!is.null(print.file)){
						dev.off()
					}

					#rownames(dOn) <- rownames(On)
					return(dOn)


		}



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#' @title Calculate standardised delta score
#' @description Function to extract results from cNMTF and calculate delta scores
#'
#[Input]:
#' @param file.exp Workspace of the experiment (created with function score.cnmtf)
#' @param path.exp Path to files of experiment
#' @param name.exp Name of experiment (trait and id.exp)
#' @param snps.known1 List of known associations
#' @param snps.known2 Second List of known associations
#' @param out.cat Outcome vector
#' @param pop Population vector
#' @param clus.a List of patient clusters to use in the delta scores
#' @param use.randomisations Logical. Use randomisations to find the optimal cutoff points for the delta scores
#' @param file.ran Path to file of randomisations
#' @param alpha.cnmtf Level of significance
#'
#[Output]:
#' @return
#' * \code{ldOn} : list of delta scores per combination of clusters
#' * \code{lsd.sc} : list of SD per combination of clusters
#' * \code{tplot.sc} : table with frequencies of predicted associations
#'
#' @md
#' @author Luis G. Leal, \email{lgl15@@imperial.ac.uk}
#' @family Scoring functions
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


		sd.score = function(file.exp, #Workspace of the experiment (created with function score.cnmtf)
												 path.exp, #Path to files of experiment
												 name.exp, #Name of experiment (trait and id.exp)
												 snps.known1, #List of known associations
												 snps.known2, #Second List of known associations
												 out.cat = NULL, #Outcome vector
												 pop = NULL, #Population vector
												 clus.a = list(c(1)), clus.b = list(c(2)), #List of patient clusters to use in the delta scores
												 use.randomisations = TRUE, #Use randomisations to find the optimal cutoff points for the delta scores
												 file.ran = lfile.ran [[i]], #Path to file of randomisations
												 alpha.cnmtf = 0.005 #Level of significance
		)
		{



				#-------------------------------------------------
				# 1. Load and explore results
				#-------------------------------------------------


						#Load experiment
						cat("Reading results of cNMTF", "\n", as.character(Sys.time()), "\n")
						if(file.exists(file.exp) ){
							load(file = file.exp)
						}else{
							stop("Workspace of experiment not found")
						}


						#Extract Omega matrix from results
						Ono <- res.cnmtf[[1]][[4]]
						Ono = Ono[, !is.na( Ono[1,] )]
						On = scale(Ono) #The scaling of columns help to have a more normal alike distribution of delta scores
						cat("Dimensions of Omega matrix:", dim(On), "\n")




						#Extract patient clustering
						clus.V = res.cnmtf[[1]][[3]]
						toc = table(out.cat, clus.V[,2])

						#Print contigency table and similarity (simVO)
						cat("Contingency table of Patient clusters Vs Patient outcomes")
						print(toc)
						require("igraph")
						simVO.i =  igraph::compare(as.numeric(as.factor(out.cat)), as.numeric(clus.V[,2]), method  =  c("nmi"))
						cat("Similarity between Patient clusters Vs Patient outcomes", simVO.i, "\n")
						simVL.i = igraph::compare(as.numeric(as.factor(pop)), as.numeric(clus.V[,2]), method  =  c("nmi"))
						cat("Similarity between Patient clusters Vs Patient population", simVL.i, "\n")



				#-------------------------------------------------
				# 2. Find significant SNVs
				#-------------------------------------------------


						cat("Calculating delta SNVs scores", "\n")

						#Declare vectors to save delta scores, standard desviations from mean, sd tresholds and pvalues
						ldOn = lsd.sc = ltao.sd = lpval.sc = rec.list( length(clus.a) )

						#Load randomisations of the experiment

						if(file.exists(file.ran) ){
							cat("Reading randomisations of cNMTF from:", file.ran, "\n")
							load(file = file.ran)
						}else{
							stop("Workspace of randomisations not found")
						}

						#Number of randomizations
						(n.randoms = length(l.cnmtf.ran))


						#For each combination of cluster of patients find sd.scores

						for(i in 1:length(clus.a)){


								#-------------------------------------------------
								# 2.1 Explore dispersion of Delta scores
								#-------------------------------------------------


										cat("Calculating delta score between patient clusters:", clus.a[[i]], "and", clus.b[[i]], "\n", as.character(Sys.time()), "\n")

										#Delta score for this combination of patient clusters
										ldOn[[i]] <- dOn <- On[, clus.a[[i]] ] - On[, clus.b [[i]] ]

										#Calculate desviations from SD
										lsd.sc[[i]] =  (dOn - mean(dOn))/sd(dOn)
										ltao.q = seq(0.5,5,0.1)
										tplot.sc = matrix(0, nrow = length(ltao.q), ncol = 4)


										#Plot number of significant SNVs at different tresholds of SD
										for( k in 1:length(ltao.q)){
											sig.sc = rownames(On) [ abs( lsd.sc[[i]] ) > ltao.q[k] ]
											tplot.sc[k,2] = length(sig.sc)
											tplot.sc[k,3] = sum(sig.sc %in% snps.known1)
											tplot.sc[k,4] = sum(sig.sc %in% snps.known2)

										}


										#Plot table
										plot(ltao.q, tplot.sc[,2], col = "gray", t = "l", las = 1, xlab = "SD", ylab = "Number of predicted associations", ylim = c(0,30))
										lines(ltao.q, tplot.sc[,3], col = "black")
										lines(ltao.q, tplot.sc[,4], col = "green")

										#Plot ratio
										#plot(ltao.q, tplot.sc[,2]/(tplot.sc[,3] + tplot.sc[,4]), col = "black", t = "l", las = 1, xlab = "SD", ylab = "Ratio P/TP associations")
										#abline(h = 10, lty = 2)


										#Add colnames to table of results
										tplot.sc[,1] <-  ltao.q
										colnames(tplot.sc) <- c("sd", "total", "known1", "known2")



								#-------------------------------------------------
								# 2.2 Find the distribution of the score under the null hypothesis
								#-------------------------------------------------


										if(use.randomisations == TRUE){

												cat("Extracting randomised scores for this pair of clusters", "\n")


												#Declare list of randomised scores
												l.rdOn = list()

												#Counter variable for the number of effective randomisations without error
												j = 0

												#Extract and scale randomised scores across randomisations for this pair of clusters
												for(k in 1:n.randoms){

													#Extract randomised scores
													rOn = scale(l.cnmtf.ran[[k]][[4]])

													#Random Delta score for this combination of patient clusters
													rdOn <- rOn[, clus.a[[i]] ] - rOn[, clus.b [[i]] ]

													#Calculate desviations from SD. The Standardisation produces better results for all traits. Also the distribution of the non-standardised scores is always very symetric.
													rdOn <-  (rdOn - mean(rdOn))/sd(rdOn)

													#Check if the randomisation was effective (i.e., none NAs)
													if ( sum(is.na(rdOn)) > 0 ){
														next
													}

													#Save randomisations in a list
													l.rdOn[[k]] <- rdOn

													#Update counter
													j = j + 1

												}


												#Compares each score Vs all scores in all permutations
												u.rdOn = unlist(l.rdOn)
												g.rdOn = length(u.rdOn)

												#Declare vector of scores and probabilities
												x.score = seq(-10,10,by=0.01)
												prob = rep(0, length(x.score))


												cat("Finding the distribution of the score under the null hyphotesis", "\n")


												#Define the probability for each score
												for(m in 1:length(x.score)){

													prob[m] <- 	sum( u.rdOn < x.score[m])	/ g.rdOn

												}

												cat("Vector of probabilities ranging between :", min(prob), max(prob), "\n")
												cat(", for a score range of :", min(x.score), max(x.score), "\n")

												#Set thresholds for the delta score
												ltao.sd [[i]][1] = x.score[prob <= alpha.cnmtf/2] [sum(prob <= alpha.cnmtf/2)]
												ltao.sd [[i]][2] = x.score[prob >= (1 - alpha.cnmtf/2)][1]


												cat("Thresholds for the delta score set to", unlist( ltao.sd [[i]] ), "\n")


												#Plot the score versus the probability

												hist(u.rdOn, breaks = 300, freq = FALSE, xlab = "Delta SNV score", col = "gray", border = "darkgray", xlim = c(-6,6), las = 1, ylim = c(0,1), main = "")
												axis(side = 4, pretty(range(prob),5), las = 2)
												lines(x.score, prob, col = "red", lwd = 2, lt = 2)
												abline(h = c(alpha.cnmtf/2, (1 - alpha.cnmtf/2) ), lt = 2, col = "blue")
												abline(h = c(alpha.cnmtf/2, (1 - alpha.cnmtf/2) ), lt = 2, col = "blue")
												abline(v = c( ltao.sd [[i]][1] , ltao.sd [[i]][2] ), lt = 2, col = "black")



												#Check the threshold at probability = 0.5
												cat("Threshold at probability = 0.5 is", x.score[prob >= 0.5][1], "\n")


												#Those scores outside the interval [-10, 10] trimmed to the limits of the interval
												lsd.sc[[i]][ lsd.sc[[i]] < -10  ] <- (-10)
												lsd.sc[[i]][ lsd.sc[[i]] > 10  ] <- 10

												#Find p-values for each SNP (Rounds scores to 2 digits)
												pval.sc <- prob [ match( round(lsd.sc[[i]],2), round(x.score,2) )]


												#Unified pvalue (left and right tale of the distribution)
												pval.sc [pval.sc > 0.5] <- 1 - pval.sc[pval.sc > 0.5]

												#Add the pvalues to list
												lpval.sc [[i]] <- pval.sc

										} #End if use.randomisations

				} #End loop through each combination of clusters



				return( list( ldOn = ldOn, lsd.sc = lsd.sc, tplot.sc = tplot.sc, On = On, simVO.i = simVO.i, simVL.i = simVL.i, ltao.sd = ltao.sd, lpval.sc = lpval.sc) )

		}


