###############################################################
# cNMTF
#	3. Delta score functions
# 3.3 Functions for Manhattan plots, Venn diagrams and clusters
#
# Corresponding author:
# Luis Leal, Imperial College London, email: lgl15@imperial.ac.uk
###############################################################


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#' @title Plotting delta scores
#' @description  Function to plot dispersion plots, manhattan plots and venn diagrams of SNP sds
#'
#[Input]:
#' @param R.snps List of SNPs in the original seed files
#' @param object.sd Results from function sd.score
#' @param work.dat Working directory
#' @param pvals.lrm P-values obtained in a LRM
#' @param alpha.lrm Significance level
#' @param snps.known1 List of known associations
#' @param snps.known2 List of known associations 2
#' @param maf.snps Minor allele frequencies of SNPs in the original R matrix without filtering
#' @param kd.snps Node degress of SNPs in the original Wu matrix without filtering
#' @param snps.fp.lrm Potential False Positive association from LRM
#' @param print.file File to print plots
#' @param tao.sd Treshold of SD
#' @param clus.a List of patient clusters to use in the delta scores
#' @param trait.project Trait/outcome
#' @param tmap Mapping of SNPs to genes, chr and genomic position
#' @param ylim.sd Limits for \code{y} axis
#'
#[Output]:
#' @return
#' * Plots printed to \code{print.file}
#' * \code{sig.snp.nmtf, sig.snp.lrm} : Significant SNV ids from cNMTF and LRM
#'
#' @md
#' @author Luis G. Leal, \email{lgl15@@imperial.ac.uk}
#' @family Scoring functions
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


		sd.plot = function(  R.snps = rownames(R), #List of SNPs in the original seed files
												 object.sd = NULL, #Results from function sd.score
												 work.dat = NULL, #Working directory
												 pvals.lrm = NULL, #P-values obtained in a LRM
												 alpha.lrm = NULL, #Significance level
												 snps.known1 = NULL, #List of known associations
												 snps.known2 = NULL, #List of known associations 2
												 maf.snps = NULL, #Minor allele frequencies of SNPs in the original R matrix without filtering
												 kd.snps  = NULL, #Node degress of SNPs in the original Wu matrix without filtering
												 snps.fp.lrm = NULL, #Potential False Positive association from LRM
												 print.file = NULL, #File to print plots
												 tao.sd = NULL, #Treshold of SD
												 clus.a = list(c(1)), clus.b = list(c(2)), #List of patient clusters to use in the delta scores
												 trait.project = NULL, #Trait/outcome
												 tmap = NULL, #Mapping of SNPs to genes, chr and genomic position
												 ylim.sd = c(-10,10) #Limits for axe y

		)	{


					#-------------------------------------------------
					# 1. Load information
					#-------------------------------------------------

						#Extract information from object.sd
						On = get("On", object.sd)
						ldOn = get("ldOn", object.sd)
						lsd.sc = get("lsd.sc", object.sd)

						#Number of combinations of clusters
						z = length(ldOn)


						#Extract thresholds if not provided by user
						if( is.null(tao.sd) ){
							ltao.sd = get("ltao.sd", object.sd)
						}else{ #Or use provided taos for all combination of clusters
							ltao.sd = rep(list(tao.sd),z)
						}


						#Extract SNP names
						On.snps = rownames(On)


						#Define positions of known associations in matrix On
						pos.known1 = which(On.snps %in% snps.known1)
						pos.known2 = which(On.snps %in% snps.known2)

						#Open connection to PDF file
						if(!is.null(print.file)){
							pdf(print.file, width = 8, height = 6)
							#par(mar = c(6.5, 6.5, 0.5, 0.5), mgp = c(4, 1, 0))
						}



					#-------------------------------------------------
					# 2. Manhattan plot and dispersion plots of SNVs
					#-------------------------------------------------


							#Filter LRM p-values
							names(pvals.lrm) <- R.snps
							pvals.lrm = pvals.lrm[match(On.snps, R.snps)]


							cat("Preparing data for Manhattan plots", "\n")

							#Create table of chromosomal position used in manhatan plots
							obj.manhattan = manhattan.table(l.snps = On.snps, tmap = tmap)
							tmanhattan = obj.manhattan[[1]]
							med.point = obj.manhattan[[2]]

							#List of chromosomes
							lchr = sort(as.numeric(unique(tmanhattan$chr)))


							#Matching positions of On.snps in tmanhattan
							match.On = match( tmanhattan$snp, On.snps)

							#Position of all seeds and those found by Sabatti
							pos.snps.known1 = which( tmanhattan$snp %in% snps.known1 )
							pos.snps.known2 = which( tmanhattan$snp %in% snps.known2 )


							cat("Defining colour code of plots", "\n")


							#Vector of colours based on chromosomal position
							lcolor = llevels = color.levels = list()
							color.levels[[1]] = rep(gray.colors(2, start = 0.7, end = 0.8),12) [ as.numeric(levels( as.factor( tmanhattan$chr ) ))] 		#color.levels = rainbow( nlevels( as.factor(obj.manhattan[[1]][,1] ) ))
							llevels[[1]] = NULL
							lcolor[[1]] = color.levels[[1]][ as.factor( tmanhattan$chr ) ]

							#Repeat plots using color code for MAF (minor allele frequency)

							if(!is.null(maf.snps)){

								#Vector of colours based on MAF
								maf.snps = maf.snps[ match( On.snps, R.snps) ] #Filter vector
								maf.snps = maf.snps[ match.On ] #Order vector
								llevels[[2]] = levels(as.factor( maf.snps ))
								color.levels[[2]] = c("red","orange","yellow") #heat.colors( nlevels( as.factor( maf.snps ) ))
								lcolor[[2]] = color.levels[[2]][ as.factor( maf.snps ) ]


							}

							if(!is.null(snps.fp.lrm)){

								#Vector of colours based on potential false positives
								snps.fp.lrm = snps.fp.lrm[ match( On.snps, names(snps.fp.lrm)) ] #Filter vector
								snps.fp.lrm = snps.fp.lrm[ match.On ] #Order vector
								lcolor[[3]] <- lcolor[[1]]
								lcolor[[3]][ snps.fp.lrm == "FP" ] <- "brown"

							}


							#Loop over color codes
							for(lc in 1:length(lcolor)){


								#Plot pair of Manhatan plot and SNV score (SD)
								par(mfrow = c(1,2))


								#Create plot for each combination of clusters
								for(i in 1:z){


									cat("Plotting Manhattan plots of SNV p-values for combination of clusters", "\n")

									#Plot values from LRM
									manhattan.plot(color = lcolor[[lc]], lchr = lchr, pvalues = pvals.lrm[ match.On ], med.point = med.point, alpha = alpha.lrm,
																 pos.snps.known1 = pos.snps.known1, pos.snps.known2 = pos.snps.known2, print.file = NULL)


									#Extract SD of delta scores
									sd.sc = lsd.sc [[i]]


									#Extract thresholds of delta scores for this combination of clusters
									tao.sd = unlist( ltao.sd [[i]])
									cat("Tresholds for this combinations of clusters", tao.sd ,"\n")


									#Plot SD from cNMTF
									plot(sd.sc[ match.On ], xaxt = "n" ,ylab = "Delta SNV score", xlab = "Chromosomal position", col = lcolor[[lc]], las = 1, pch = 16, ylim = ylim.sd)
									axis(side = 1, at = med.point[as.numeric(lchr)], labels = paste("",lchr,sep=""),las = 1)
									abline(h = 0)
									abline(h = tao.sd[1],lty =3)
									abline(h = tao.sd[2],lty =3)


									#Highlight known associations
										if( !is.null(pos.snps.known1) ){
											points(y = sd.sc[ match.On ][pos.snps.known1], x = pos.snps.known1, pch =10, col = "black")
										}
										if( !is.null(pos.snps.known2) ){
											points(y = sd.sc[ match.On ][pos.snps.known2], x = pos.snps.known2, pch = 10, col = "green")
										}

										print(lc)
										if(lc %in% c(2)){
											legend(x = 1,y = 12, legend = llevels[[lc]], fill = color.levels[[lc]])
										}

									#Plot cNMTF vs LRM
										#plot( -log10(pvals.lrm[ match.On ]), sd.sc[ match.On ], ylab = "Delta SNV score (SD)", xlab = "-log(p-value)", pch = 16, col = unlist( lcolor[[lc]] ), ylim = ylim.sd, las = 1, xlim = c(0,10))
										#if( !is.null(pos.snps.known1) ){
										#	points(y = sd.sc[ match.On ][pos.snps.known1], x = -log10(pvals.lrm[ match.On ]) [ pos.snps.known1 ], pch =10, col = "black")
										#}
										#if( !is.null(pos.snps.known2) ){
										#	points(y = sd.sc[ match.On ][pos.snps.known2], x = -log10(pvals.lrm[ match.On ]) [ pos.snps.known2 ], pch = 10, col = "green")
										#}
										#abline(h = 0)
										#abline(h = tao.sd,lty =3)
										#abline(h = -tao.sd,lty =3)
										#abline(v = -log10(alpha.lrm),lty =3)


								} #End Loop over combination of clusters

							}#End Loop over color codes



					#-------------------------------------------------
					# 3. Plot venn diagrams
					#-------------------------------------------------


							cat("Plot venn diagrams for the union of clusters", "\n")

							#Declare vectors of results: significant associations
							sig.snp.nmtf = list()
							sig.snp.lrm = list()

							for(i in 1:z){

									#Extract SD of delta scores
									sd.sc = lsd.sc [[i]]


									#Extract list of significant SNPs
									sig.snp.nmtf[[i]] = sd.sc [ sd.sc  <= tao.sd[1] |  sd.sc  >= tao.sd[2] ]
									sig.snp.lrm[[i]] =  pvals.lrm[ pvals.lrm  < alpha.lrm ]

									#Delete NAs
									sig.snp.nmtf[[i]] <- sig.snp.nmtf[[i]][ !is.na( sig.snp.nmtf[[i]] )]
									sig.snp.lrm[[i]] <- sig.snp.lrm[[i]][ !is.na( sig.snp.lrm[[i]] )]

									#Create venn objects
									if( length(sig.snp.nmtf) != 0 & length(	sig.snp.lrm) != 0 ){

											#Plot venn diagrams for the union of clusters
											venn.plot <- venn.diagram(  x = list( "LRM" = names( sig.snp.lrm[[i]] ) , "cNMTF" = names( sig.snp.nmtf [[i]]) , "GWAS" = snps.known1 ),
																									filename = NULL,
																									output = TRUE,
																									imagetype="png" ,
																									resolution = 300,
																									compression = "lzw",
																									margin = 0.1)

											plot.new()
											par(mfrow = c(1,1))
											grid.draw( venn.plot )

									}

							}	#End Loop over combination of clusters

						#Close PDF file connection
						if(!is.null(print.file)){
							dev.off()
						}


					 return(list(sig.snp.nmtf = sig.snp.nmtf, sig.snp.lrm = sig.snp.lrm))

		}




#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#' @title Manhattan table
#' @description Function to create a table of chromosomal position used in manhatan plots
#'
#[Input]
#' @param l.snps List of SNPs to construct the table
#' @param tmap Mapping of SNPs to genes, chr and genomic position
#'
#[Output]
#' @return
#' * tmanhattan : table of chromosome position of SNPs.
#' * med.point : parameter to place the axes labels in a manhattan plot
#'
#' @md
#' @author Luis G. Leal, \email{lgl15@@imperial.ac.uk}
#' @family Scoring functions
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


		manhattan.table = function(l.snps, #List of SNPs to construct the table
																 tmap #Mapping of SNPs to genes, chr and genomic position
			)
			{

					#Define parameters to create the table

					tmanhattan = matrix(NA, nrow = length(l.snps), ncol = 3) # Table of chromosome position of SNPs. The SNPs positions are ordered by default from the reading files.
					med.point = rep(NA,22) # Parameter to place the axes labels in the plot
					k = 1 # Counter variable for the rows of the table

					#Read chromosome and position
					pos.match = match(as.character(l.snps), as.character(tmap$refsnp_id) )
					tmanhattan = data.frame(chr = as.numeric( as.character( tmap$chr[ pos.match ]) ),
																	position = as.numeric( as.character( tmap$position[ pos.match ]) ),
																	snp = l.snps, stringsAsFactors = FALSE)



					#Sort the table by genomic position and chromosome
					tmanhattan = tmanhattan[ order(tmanhattan$chr,tmanhattan$snp), ]

					#Find postion to place the axes labels in Manhattan plot
					for(i in 1:22){ # For each chromosome in the list of SNPs

							#Count number of SNPs ath this chromosome
							n.snps.chr = sum(tmanhattan$chr == i) #Do not add na.rm = T to discover mistakes

							#Calculate the parameter to place the axes labels in the plot
							med.point[i] <- ceiling( (k + n.snps.chr - 1) - k ) / 2 + k

							#Update counter variables and print progress
							k = k + n.snps.chr

					}


					return( list(tmanhattan, med.point) )
			}



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#' @title Manhattan plot
#' @description  Function to plot a Manhattan plot
#'
#[Input]
#' @param pvalues pvalues of the SNPs in tmanhattan
#' @param lchr List of chromosomes
#' @param color colour for each SNP
#' @param pos.snps.known1,pos.snps.known2 Positions of known associations
#' @param print.file PDF file to print plots
#' @param alpha Significance level
#' @param med.point parameter to place the axes labels in a manhattan plot
#'
#[Output]
#' @return Print manhattan plot either in console or in print.file
#'
#' @md
#' @author Luis G. Leal, \email{lgl15@@imperial.ac.uk}
#' @family Scoring functions
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


		manhattan.plot = function(pvalues,  #pvalues of the SNPs in tmanhattan
															lchr = 1:22, #List of chromosomes
															color, #colour for each SNP
															pos.snps.known1 = NULL, pos.snps.known2 = NULL, #Positions of known associations
															print.file = NULL, #PDF file to print plots
															alpha = 0.05, #Significance level
															med.point #parameter to place the axes labels in a manhattan plot

			)
			{

						if(!is.null(print.file)){
							pdf(print.file, width = 6, height = 6)
						}


						#Create plots
						plot( -log10(pvalues), xaxt = "n", ylim = c(0,10),ylab = "-log(p-value)", xlab = "Chromosomal position", col =color, las =1, pch = 16)
						axis( side = 1, at = med.point[as.numeric(lchr)], labels = paste("", lchr, sep=""),las = 1)
						abline( h = -log10(alpha), lty =3)

						#Highlight known associations

						if( !is.null(pos.snps.known1) ){
							points(y = -log10(pvalues)[pos.snps.known1], x = pos.snps.known1, pch =10, col = "black")
						}


						if( !is.null(pos.snps.known2) ){
							points(y = -log10(pvalues)[pos.snps.known2], x = pos.snps.known2, pch = 10, col = "green")
						}


						#abline(h=-log10(0.05/n),lty =3)
						#abline(h=-log10(1e-8),lty =3)
						#text(x = med.point[length(med.point)-4],y=-log10(0.05/n), labels = "adjusted")
						#text(x = med.point[length(med.point)-5],y=-log10(1e-8), labels = "treshold of significance")
						#identify(n = 13, y = -log10(pvalues.unc),x = seq(1:length(pvalues.unc)),labels = tmanhattan[,3])

						#Close file connection
						if(!is.null(print.file)){
							dev.off()
						}

			}



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#' @title Plot cNMTF clusters
#' @description Function to plot heatmaps and clusters of SNPs and patients
#'
#[Input]:
#' @param file.exp: Workspace of the experiment (created with function score.cnmtf)
#' @param snps.known: List of known associations
#' @param print.file File to print plots
#' @param maf.snps.cat Minor Allele Frequency
#' @param R.snps List of SNVs
#'
#[Output]:
#' @return Plots of clusters printed to \code{print.file}
#'
#' @md
#' @author Luis G. Leal, \email{lgl15@@imperial.ac.uk}
#' @family Scoring functions
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


    plot.clusters = function(file.exp, #Workspace of the experiment (created with function score.cnmtf)
    												 snps.known, #List of known associations
    												 print.file = NULL, #File to print plots
    												 maf.snps.cat, #Minor Allele Frequency
    												 R.snps #List of SNVs

    ){


    			#Load experiment
    			cat("Reading results of cNMTF", "\n", as.character(Sys.time()), "\n")
    			if(file.exists(file.exp) ){
    				load(file = file.exp)
    			}else{
    				stop("Workspace of experiment not found")
    			}

    			#Load matrices
    			Ua = res.cnmtf[[1]][[11]] #Average U matrix
    			Va = res.cnmtf[[1]][[12]] #Average V matrix
    			Or = res.cnmtf[[1]][[10]] #Expaned consensus On matrix
    			On = scale(res.cnmtf[[1]][[4]]) #Consensus On matrix
    			On.snps = rownames(On)	#SNPs in consensus On
    			clus.U = res.cnmtf[[1]][[2]] #Clustering membership of SNPs
    			clus.V = res.cnmtf[[1]][[3]]  #Clustering membership of patients

    			#Open connection to pdf file
    			if(!is.null(print.file)){
    				pdf(print.file, width = 8, height = 6)
    			}


    			#Plot Omega matrix

    			cat("Plotting heatmap of Omega matrix", "\n")

    			#Creat color index for labels
    			row.colors = c(rep("black",nrow(On)))
    			row.colors[ On.snps %in% snps.known] <- "blue"

    			#Print heatmap
    			Heatmap(On, col = colorRampPalette(c("cyan","black", "green"))( 8 ), cluster_rows = FALSE, cluster_columns = FALSE,
    							row_title = "SNVs", column_names_side = c("bottom"),  show_row_names = FALSE, row_names_side = c("left"),
    							row_names_gp = gpar(cex = 0.5, col = row.colors),
    							column_title = "Clusters of patients", heatmap_legend_param = list(title = "S.Score", color_bar = "discrete"),
    							width = 0.5)


    			cat("Plotting first and second factors of U and V", "\n")

    			#Plot first and second factors of U
    			par(mfrow = c(1,2))
    			plot( Ua[,1:2], las = 1, pch = 16, xlab = "First factor of U", ylab = "Second factor of U", col = topo.colors(nlevels(as.factor(clus.U[,2])))[as.factor(clus.U[,2])])
    			points( Ua[ rownames(Ua) %in% snps.known,1:2], pch = 10)
    			plot( Ua[,1:2], las = 1, pch = 16, xlab = "First factor of U", ylab = "Second factor of U",  col = c("red", "orange", "yellow")[maf.snps.cat[match(rownames(Ua),R.snps)]])
    			points( Ua[ rownames(Ua) %in% snps.known,1:2], pch = 10)

    			#Plot first and second factors of V
    			plot( Va[,1:2], las = 1, pch = 16, xlab = "First factor of V", ylab = "Second factor of V", col = rainbow(nlevels(as.factor(clus.V[,2])))[as.factor(clus.V[,2])])

    			#Plot first and second patient clusters of Omega
    			plot( On[,1:2],las = 1, pch = 16, xlab = "SNV score in cluster of patients 1", ylab = "SNV score in cluster of patients 2", col = c("red", "orange", "yellow")[maf.snps.cat[match(rownames(On),R.snps)]])
    			points( On[ rownames(On) %in% snps.known,1:2], pch = 10)

    			#Close PDF file connection
    			if(!is.null(print.file)){
    				dev.off()
    			}

    			return()

    }



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#' @title Plot Venn Diagrams
#' @description Function to plot customised the venn diagrams
#'
#[Input]
#' @param venn.list list of elements in each set
#' @param fill.venn colour for set
#' @param main title for the diagram
#'
#[Output]:
#' @return Printed venn diagram to console
#'
#' @md
#' @author Luis G. Leal, \email{lgl15@@imperial.ac.uk}
#' @family Scoring functions
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


		custom.venn = function( venn.list,  main = NULL, filename = NULL ){

				#If venn lists are of length 0, return NULL
			 	if ( sum(unlist(lapply( venn.list, FUN = length))) == 0 ) { return( NULL )}

				#Number of sets
				n.venn = length(venn.list)

				#Colours for diagrams
				fill.venn = c('blue', 'orange', 'green', 'yellow', 'brown')[1:(n.venn)]

				#Generate object with the geometry for the venn diagram
				venn.plot =	venn.diagram( x = venn.list ,
																	filename = filename,
																	main = main,
																	main.cex = 0.7,
																	height = 2500, width = 2500,
																	fill = fill.venn,
																	cex = rep( 1, c(1,3,7,15,31)[(n.venn)]),
																	cat.cex = 1, out = FALSE, cat.pos = c(0, 0, 180, 180)[1:n.venn],
																	euler.d = TRUE, cat.default.pos = "text", ext.text = TRUE)

		}



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#' @title Recursive lists
#' @description Function to create recursive lists
#'
#[Input]
#' @param len vector with the length of list levels.
#'
#[Output]
#' @return Return a nested list with different lengths at each level
#'
#' @md
#' @author Luis G. Leal, \email{lgl15@@imperial.ac.uk}
#' @family Scoring functions
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


		rec.list <- function(len){
			if(length(len) == 1){
				vector("list", len)
			} else {
				lapply(1:len[1], function(...) rec.list(len[-1]))
			}
		}



