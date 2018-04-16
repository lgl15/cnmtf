###############################################################
# cNMTF
#	1. Preprocessing functions
# 1.4 Functions to construct networks
#
# Corresponding author:
# Luis Leal, Imperial College London, email: lgl15@imperial.ac.uk
###############################################################


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#' @title Constructs the SNV-SNV network
#' @description Function to create a SNV-SNV network using a PPI reference network
#'
# [Input]
# Parameters for the reference network
#' @param net.type Type of reference network. Default "ppi"
#' @param file.dedges Filename with edges from reference network
#'
# Parameters for Linkage Disequilibrium
#' @param expand.highLD Expand SNV consequences to SNVs in high LD
#' @param ld.tao Treshold of LD. Default 0.8
#' @param res.ld Table of LD
#' @param keep.with.LD List of SNPs to include even if they are in LD
#'
# Parameters to counstruct Wu
#' @param R.snps List of SNPs in R
#' @param work.dat Working directory
#' @param trait.project Trait
#' @param n.cores Number of cores for parallel computing
#' @param tmap Mapping of SNPs to genes
#' @param alledges.snp Include edges between proteins for all SNPs or only for hubs
#' @param venn.diag Logical. Print Venn diagrams. Default = FALSE.
#' @param plot.file File to print Venn diagrams and node degree distribution
#'
#'
# [Output]
#' @return
#' * \code{Wu}, \code{G}: adjacency matrix and graph object of the network
#' * Table with node properties to be used in Cytoscape
#' * Venn diagrams of damaging variants and node degree distribution
#'
#' @md
#' @author Luis G. Leal, \email{lgl15@@imperial.ac.uk}
#' @family Preprocessing functions
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


	create.network = function(	#Parameters for the reference network
															net.type = "ppi", #Type of reference network
															file.dedges = NULL, #Define filename with edges from reference network paste(work.dat,"data/ref_networks/", net.type, "_edges_",trait.project,".txt",sep="")

															#Parameters for Linkage Disequilibrium
															remove.highLD = TRUE,
															ld.tao = 0.8, #Treshold of LD
															res.ld = NULL, #Table of LD
															keep.with.LD = NULL, #List of SNPs to include even if they are in LD

															#Parameters to counstruct Wu
															R.snps = NULL, #List of SNPs in R
															work.dat = NULL, #Working directory
															trait.project = NULL, #Trait
															n.cores = 4, #Number of cores
															tmap = NULL, #Mapping of SNPs to genes
															alledges.snp = FALSE, #Include edges between proteins for all SNPs or only for hubs
															plot.file = NULL, #File to print Venn diagrams and node degree distribution
															venn.diag = FALSE) #Print Venn diagrams?

	{

		#-------------------------------------------------
		# 1. Load/construct reference network ----
		#-------------------------------------------------


				#Read file with edges
				if( !file.exists(file.dedges) ){

						stop("The reference PPIN was not found.")

				}else{

					cat("Loading reference network: list of interactions between gene entrezids", "\n")
					dedges = read.table(paste(work.dat,"data/ref_networks/", net.type, "_edges_",trait.project,".txt",sep=""), head = TRUE)

				}

				cat("Dimmensions of table of edges (dedges):", "\n"); print(dim(dedges))




		#-------------------------------------------------
		# 2. Find SNP localization/consequences ----
		#-------------------------------------------------


				#Extract information of SNP consequences from ENSEMBL
				ensembl.snp <- useMart(biomart="ENSEMBL_MART_SNP", dataset = "hsapiens_snp", host = "grch37.ensembl.org" ) #listAttributes(ensembl.snp)
				snp.conseq = getBM(attributes=c("refsnp_id", "consequence_type_tv", "consequence_allele_string", "polyphen_prediction", "polyphen_score", "sift_prediction", "sift_prediction"), filters="snp_filter", values=R.snps, mart=ensembl.snp)

				#Check quality of data
				cat("Quality of the data", "\n")
				cat("SNV consequences ENSEMBL:", "\n"); print( table(snp.conseq$consequence_type_tv) )
				cat("SNV consequences Polyphen:", "\n"); print( table(snp.conseq$polyphen_prediction) )
				cat("SNV consequences Sift:", "\n"); print( table(snp.conseq$sift_prediction) )

				#Check how many hits per predictor
				cat("Number of SNVs with consequences ENSEMBL:", "\n");  print( length(unique(snp.conseq$refsnp_id[ snp.conseq$polyphen_prediction != "" ])))
				cat("Number of SNVs with consequences Polyphen:", "\n"); print( length(unique(snp.conseq$refsnp_id[ snp.conseq$sift_prediction != "" ])) )
				cat("Number of SNVs with consequences Sift:", "\n"); print( length(unique(snp.conseq$refsnp_id[ snp.conseq$consequence_type_tv != "" ])) )

				#Extract list of SNVs with damaging/deleterious effects
				snp.poly = unique(snp.conseq$refsnp_id[ snp.conseq$polyphen_prediction %in% c("possibly damaging", "probably damaging") ])
				snp.sift = unique(snp.conseq$refsnp_id[ snp.conseq$sift_prediction %in% c( "deleterious", "deleterious - low confidence" )  ])

				#Extract list of SNPs with HIGH impact according to ENSEMBL
				snp.impact = unique(snp.conseq$refsnp_id[ snp.conseq$consequence_type_tv %in% c("transcript_ablation", "splice_acceptor_variant", "splice_donor_variant", "stop_gained", "frameshift_variant", "start_lost", "transcript_amplification", "inframe_insertion", "inframe_deletion", "missense_variant", "protein_altering_variant") ] ) #Source: https://www.ensembl.org/info/genome/variation/predicted_data.html ])



		#-------------------------------------------------
		# 3. Share consequences to other SNPs in high LD ----
		#-------------------------------------------------


				if( remove.highLD == TRUE ){

						#Expand consequences to other SNPs in high LD
						snp.poly.exp = unique( c(snp.poly, as.character(res.ld$loc1) [ res.ld$loc2 %in% snp.poly & res.ld$r2 > ld.tao ], as.character(res.ld$loc2) [ res.ld$loc1 %in% snp.poly & res.ld$r2 > ld.tao ] ) )
						snp.sift.exp = unique( c(snp.sift, as.character(res.ld$loc1) [ res.ld$loc2 %in% snp.sift & res.ld$r2 > ld.tao ], as.character(res.ld$loc2) [ res.ld$loc1 %in% snp.sift & res.ld$r2 > ld.tao ] ) )

						#Define vector of high LD SNPs to remove
						set.out = unique(as.character(res.ld$loc2[ res.ld$r2 > ld.tao] ))

						#Include some SNPs specified by the user even if they are SNPs in LD
						set.out = setdiff( set.out, keep.with.LD )

						#Filter vectors
						snp.poly.exp = snp.poly.exp [ !(snp.poly.exp %in% set.out) ]
						snp.sift.exp = snp.sift.exp [ !(snp.sift.exp %in% set.out) ]
						snp.impact = snp.impact [ !(snp.impact %in% set.out) ]
						snps.known2 = snps.known2 [ !(snps.known2  %in% set.out) ]

				}



		#-------------------------------------------------
		# 4. Venn Diagram ----
		#-------------------------------------------------

				#Merge lists into one
				ldamaging = unique(c(snp.poly.exp, snp.sift.exp, snp.impact))

				#Check number of damaging variants and known associations
				length(ldamaging)
				cat( length(snps.known2), "SNVs known to be associated with the trait. source: GWAS catalogue.", "\n")
				cat( length(snp.impact), "SNVs with high or moderate impact (e.g, frameshift variant, stop gained) source: Ensembl.", "\n")
				cat( length(union(snp.poly.exp,snp.sift.exp)), "deleterious SNVs. source:  Sift and Polyphen.", "\n")
				cat( "Total:", length(union(ldamaging,snps.known2)), "damaging SNVs", "\n" )


				if( venn.diag == TRUE){

						#Create venn diagrams for each source
						venn.list = list(Sift.Poly = union( snp.poly.exp, snp.sift.exp), ENSEMBL = snp.impact, GWAS = snps.known2)
						fill.venn = c('blue', 'orange', 'green', 'yellow', 'brown')[1:(length(venn.list))]


						venn.plot =	venn.diagram( x = venn.list ,
																			filename = NULL,
																			output = TRUE ,
																			imagetype="png" ,
																			resolution = 300,
																			compression = "lzw",
																			lty = 'blank',
																			fill = fill.venn,
																			cex = rep( 2, c(1,3,7,15,31)[length(venn.list)]),
																			cat.cex = 2, margin = 0.1)

						#Plot diagram to file
						pdf(plot.file, width = 6, height = 6)
								plot.new()
								par(mfrow = c(1,1))
								grid.draw( venn.plot )


				}

		#-------------------------------------------------
		# 5. Construct the network guided by damaging SNVs ----
		#-------------------------------------------------


				#Hubs in the network
				lhubs = union( ldamaging, snps.known2 )

				#Run function to construct and print the network
				source("f_net.R")
				res.Wu = construct.Wu(R.snps = R.snps[ !(R.snps  %in% set.out) ], work.dat = work.dat,
															tmap = tmap, net.type = net.type,  dedges = dedges,  #Network as a list of edges
															trait.project = trait.project, n.cores = n.cores,
															alledges.snp = FALSE, lhubs = lhubs)


				#Check density and node degree distribution
				source("f_net.R")
				degree.distribution( res.Wu[[1]], gamma.kd = 2, weighted = TRUE)
				dev.off()

				#Number of damaging and known SNPs connected in the network
				cat("Number of damaging SNPs connected in the network", sum( R.snps[ rowSums( res.Wu[[1]] )>0] %in% ldamaging), "\n", as.character(Sys.time()), "\n")
				cat("Number of known SNPs connected in the network", sum( R.snps[ rowSums( res.Wu[[1]] )>0] %in% snps.known2), "\n", as.character(Sys.time()), "\n")


				#Create with node properties for Cytoscape
				snp.cat = rep("Candidate",length(R.snps))
				snp.cat[ R.snps %in% ldamaging ] <- "Damaging"
				snp.cat[ R.snps %in% snps.known2 ] <- "Known"

				#Add Gene
				snp.gene = as.character(tmap$entrezgene[ match(R.snps, tmap$refsnp_id)])
				R.snps[which(is.na(snp.gene))] #SNPs without linked gene

				#Write table
				cat("Writing table with node properties for Cytoscape", "\n");
				write.table( cbind(R.snps, snp.cat, snp.gene), file = paste(work.dat,"data/r_workspaces/", trait.project, "/Gu_",net.type, "_", trait.project, "_attributes.txt",sep=""), row.names = FALSE, quote = FALSE)



	}



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#' @title Creates edge list for the SNV-SNV network
#' @description Function to find pairs of SNVs linked in the network
#'
#[Input]
#' @param i counter variable over the list of genes
#' @param dedges network as a list of edges between gene entrezids
#' @param lgene list of genes mapping the SNPs in R
#' @param R.snps list of SNPs in the relationship matrix
#' @param tmap table for mapping SNPs to genes
#'
#[Output]
#' @return Nested list of SNP positions to be connected by edges
#'
#' @md
#' @author Luis G. Leal, \email{lgl15@@imperial.ac.uk}
#' @family Preprocessing functions
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


		rownet = function(i = NULL, dedges = NULL,lgene = NULL, R.snps = NULL, tmap = NULL){

		  print(cat("Finding SNP-SNP interactions for gene", i, "\n"))

		  #Declare lists of SNPs of gene i and SNPs of gene i's interactors
		  pos.a = NULL; pos.b = NULL

			#Interactors of protein i (Rows of protein i in the table dedges)
			a = c( which( dedges[,1] %in% lgene[i]), which( dedges[,2] %in% lgene[i] ))

			#Find the SNPs of gene i and map them to the list of snps total (i.e. rownames of R)
			pos.a = which( R.snps %in% tmap$refsnp_id[ tmap$entrezgene %in% lgene[i] ])


			for(j in i:length(lgene)){ #Including SNP-SNP interactions within the same gene

				##Interactors of protein j (Rows of protein j)
				b = c( which( dedges[,1] %in% lgene[j] ), which( dedges[,2] %in% lgene[j] ) )

				if( sum(a%in%b)>0 ){ # If they share any interaction

					##Find the SNPs of gene j and map them to the list of snps total (i.e. rownames of R)
					pos.b = c(pos.b, which( R.snps %in% tmap$refsnp_id[ tmap$entrezgene %in% lgene[j] ]))

				}
			}

			if(length(pos.b) > 0){
			  return(list(pos.a,pos.b))
			}else{
			  return()
			}


		}



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#' @title Constructs the SNV-SNV adjacency matrix
#' @description Function to construct the adjacency matrix Wu for the SNV-SNV network
#'
#[Input]
#' @param R.snps list of SNPs
#' @param trait.project project's name to be included in the file name
#' @param tmap mapping of SNPs to gene entrezids
#' @param ncores number of cores for parallel computing
#' @param lhubs list of damaging variants
#'
#[Output]
#' @return \code{Wu}, \code{G} adjacency matrix and graph object of the network
#'
#' @md
#' @author Luis G. Leal, \email{lgl15@@imperial.ac.uk}
#' @family Preprocessing functions
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


	construct.Wu = function(R.snps = NULL, work.dat = NULL, trait.project = NULL,
													dedges = NULL, n.cores = NULL, net.type,
													tmap = NULL, lhubs = NULL,
													alledges.snp = FALSE #Include edges between proteins for all SNPs or only for hubs
													)
	{

			cat("Filtering list of genes to match those mapped in R", "\n")


			#Filter tmap to include only snps in R
			tmap = tmap[tmap$refsnp_id %in% R.snps, ]

			#Define list of genes shared by the SNPs
			tmap = data.frame(tmap)
			lgene = unique(tmap$entrezgene)
			lgene = lgene[!is.na(lgene)]

			#Filter list of genes to match those mapped in R
			lgene.s = lgene[ lgene %in% tmap$entrezgene[ tmap$refsnp_id %in% R.snps] ]

			cat("Number of genes mapping R:", length(lgene.s), "\n", as.character(Sys.time()), "\n")


			#Register the number of cores for parallel computing
				if( is.null(n.cores) ){
					doMC::registerDoMC(cores = (detectCores() - 2 ))
				}else{
					doMC::registerDoMC(cores = n.cores)
				}

			cat("Finding edges in the PPI network for each gene", "\n", as.character(Sys.time()), "\n")

			#Create the adjacency matrix
			Wu = matrix(0, nrow = length(R.snps), ncol = length(R.snps))

			#Find edges

			if( alledges.snp == TRUE ){


						#Find the edges for each pair of genes
						cat("Creating SNP-SNP network", "\n")

						rows.Wu <- foreach(i = 1:(length(lgene.s))) %dopar% {
							rows.Wu = rownet(i = i, dedges = dedges, lgene = lgene.s, R.snps = R.snps, tmap = tmap)}

						cat("Filling the adjacency matrix", "\n", as.character(Sys.time()), "\n")


						#Crete edges between snps in the same protein and hubs

						for(i in 1:length(rows.Wu)){

							#Edges between snps in the same protein or between proteins
							if( length( rows.Wu[[i]][[1]] ) > 0 & length( rows.Wu[[i]][[2]] ) > 0 ){

								Wu[ rows.Wu[[i]][[2]] , rows.Wu[[i]][[1]]  ] <- 1.5
								Wu[ rows.Wu[[i]][[1]] , rows.Wu[[i]][[2]]  ] <- 1.5

								#Find edges between hubs
								hubs.1 = intersect( rows.Wu[[i]][[1]] , which(R.snps %in% lhubs))
								hubs.2 = intersect( rows.Wu[[i]][[2]] , which(R.snps %in% lhubs))

								#Find edges between non-hubs
								reg.1 = intersect( rows.Wu[[i]][[1]] , which(!(R.snps %in% lhubs)))
								reg.2 = intersect( rows.Wu[[i]][[2]] , which(!(R.snps %in% lhubs)))

								if( length( hubs.1 ) > 0 & length( hubs.2 ) > 0 ){
									Wu[ hubs.1, hubs.2 ] <- 2
									Wu[ hubs.2, hubs.1 ] <- 2
								}
								if( length( reg.1 ) > 0 & length( reg.2 ) > 0 ){
									Wu[ reg.1, reg.2 ] <- 1
									Wu[ reg.2, reg.1 ] <- 1
								}

							}

						} #End loop across pairs of genes




			}else{


						#Find the edges at each row
						cat("Creating SNP-SNP network with hubs", "\n")

						rows.Wu <- foreach(i = 1:(length(lgene.s))) %dopar% {
							rows.Wu = rownet.hub(i = i, dedges = dedges, lgene = lgene.s, R.snps = R.snps, tmap = tmap, lhubs = lhubs)}

						cat("Filling the adjacency matrix", "\n", as.character(Sys.time()), "\n")


						for(i in 1:length(rows.Wu)){

							#Edges between snps in the same protein
							if( length( rows.Wu[[i]][[1]] ) > 0 & length( rows.Wu[[i]][[2]] ) > 0 ){

								Wu[ rows.Wu[[i]][[2]] , rows.Wu[[i]][[1]] ] <- 1 / (length ( rows.Wu[[i]][[1]] ) - 1)
								Wu[ rows.Wu[[i]][[1]] , rows.Wu[[i]][[2]] ] <- 1 / (length ( rows.Wu[[i]][[1]] ) - 1)

							}

							#Edges between snps and hubs in the same protein
							if( length( rows.Wu[[i]][[1]] ) > 0 & length( rows.Wu[[i]][[3]] ) > 0 ){

								Wu[ rows.Wu[[i]][[1]] , rows.Wu[[i]][[3]] ] <- 2 / (length ( rows.Wu[[i]][[1]] ) - 1)
								Wu[ rows.Wu[[i]][[3]] , rows.Wu[[i]][[1]] ] <- 2 / (length ( rows.Wu[[i]][[1]] ) - 1)
							}

							#Edges between hubs (either in same protein or between proteins)
							if( length( rows.Wu[[i]][[3]] ) > 0  & length( rows.Wu[[i]][[4]] ) > 0 ){

								Wu[ rows.Wu[[i]][[3]] , rows.Wu[[i]][[4]] ] <- 2
								Wu[ rows.Wu[[i]][[4]] , rows.Wu[[i]][[3]] ] <- 2

							}

						} #End loop across pairs of genes



			}



			#Transform the matrix to an adjacency matrix
				cat("Changing adjacency matrix class", "\n", as.character(Sys.time()), "\n")
				Wu = as.matrix(Wu, matrix.type = c("adjacency"))
				diag(Wu) <- 0


			#Add rownames and column names
				rownames(Wu) <- colnames(Wu) <- R.snps

			#Print results
				Gu = graph.adjacency(Wu, mode="undirected", weighted = TRUE)

			#Write table with edges and weights
				cat("Printing list of edges to", paste(work.dat,"data/r_workspaces/", trait.project, "/Gu_",net.type, "_", trait.project, ".txt",sep=""), "\n", as.character(Sys.time()), "\n")
				write.graph(Gu,file = paste(work.dat,"data/r_workspaces/", trait.project, "/Gu_", net.type, "_", trait.project, ".txt",sep=""), format = "ncol")

			#Save graph object
				save(list = c("Wu","Gu"),file = paste(work.dat,"data/r_workspaces/", trait.project, "/Gu_", net.type, "_", trait.project,".RData",sep=""))

			return(list(Wu,Gu))

	}



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#' @title Plot node degree of a network
#' @description Function to plot node degree distribution and calculate some graph variables
#'
# [ Input ]
#' @param Gu graph object
#' @param gamma.kd gamma parameter for a law power distribution of node degree
#' @param weighted are the edges weighted?
#'
# [ Output ]
#' @return Print graph variables and plot node degree distribution
#'
#' @md
#' @author Luis G. Leal, \email{lgl15@@imperial.ac.uk}
#' @family Preprocessing functions
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


	degree.distribution = function(Gu, gamma.kd, weighted = FALSE){

			#Calculate network variables
			if( weighted == TRUE ){

					rsum = rowSums(Gu)
					nE = sum(Gu != 0)
					nV = sum( rsum > 0 )
					cat("Total edges", nE, "\n")
					cat("Total nodes", nV, "\n")
					cat("Average node degree", sum(rsum)/nV, "\n")

					hist(rsum, breaks = 50, xlab = "Weighted node degree", las = 1, main = "", col = "light blue")


			}else{

					nE = length(E(Gu))
					nV = length(V(Gu))
					cat("Total edges", nE, "\n")
					cat("Total nodes", nV, "\n")
					cat("Average node degree", nE/nV, "\n")

					#Calculate node degree distribution and its frequencies
					kd.snps = degree(Gu)
					kd = as.numeric(names(table(kd.snps)))
					p.kd = as.numeric(table(kd.snps))
					tkd = cbind(kd,p.kd)
					colnames(tkd) <- c("kd", "frequency")
					print(head(tkd))
					print(tail(tkd))

					#Plot degree distribution
					plot(kd, p.kd , xlab = "Node degree (k)", ylab = "Frequency", las = 1, pch = 16, col = "darkgrey")
					pred.p.kd = (kd ^ (-gamma.kd)) * nE
					lines(kd, pred.p.kd, lty = gamma.kd, col = "blue")


			}


	}

