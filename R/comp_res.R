###############################################################
# cNMTF
#	4. Comparing results
# 4.1 Functions to add metainformation on prioritised SNVs
#
# Corresponding author:
# Luis Leal, Imperial College London, email: lgl15@imperial.ac.uk
###############################################################


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#' @title Annotate prioritised SNVs
#' @description Function to map SNVs to genes and retrieve functional annotations from DAVID
#'
#[Input]:
#' @param lconsortiums List of data consortiums to merge
#' @param ltrait.project List of traits to merge
#' @param snps.known2 List of known associations
#' @param use.web.service Logical. Use DAVID web service or export/import manually
#' @param explore.gwas Logical. Explore mutations from GWAS cataloge
#' @param email.david Email account registered in DAVID.

#'
#[Output]:
#'
#' @return
#' * \code{cdav.all} Chart of annotations from DAVID. Printed to file "summary_david_chart.txt"
#' * \code{tsa.all} Table of annotations from DAVID. Printed to file "summary_david_table.txt"
#' * \code{tmut.all} Table merging annotations for a list of consortiums and traits. Printed to file "summary_david_table_genes.txt"
#' * \code{tres} Final table with metadata for all SNVs. Printed to file "explore_results.csv"
#'
#' @md
#' @author Luis G. Leal, \email{lgl15@@imperial.ac.uk}
#' @family Comparing functions
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


  annotate.results = function( lconsortiums, #List of data consortiums to analyse
                               ltrait.project,	#List of traits to analyse
                               snps.known2, #Second list of known associations
                               use.web.service = TRUE, #Use DAVID web service or export/import manually
                               explore.gwas = FALSE, #Logical. Explore mutations from GWAS cataloge
                               email.david #Email account registered in DAVID.
  )
  {


		#List of data consortiums to analyse
		lconsortiums = c("NFBC", "eMERGE", "NFBC2", "eMERGE2")[1:2]

		#List of traits to analyse
		ltrait.project = c("ldl","hdl","tg")[1:3]


		#Use web service or export/import manually
		use.web.service = TRUE
		explore.gwas = FALSE

		#Declare consolidate chart and table of functional annotations
		cdav.all = tsa.all = tmut.all =  NULL

		#Iterate across consortiums and traits
		for(data.consortium in lconsortiums){
			for(trait.project in ltrait.project){


				#------------------------------------------------
				# 4.1 Exporting/importing results to/from DAVID
				#------------------------------------------------


						#Define experiment id
						id.exp = "5002"

						cat("\n", "\n", rep("-",20), "\n", "\n")
						cat("Analysing results of ", data.consortium, "in trait" , trait.project, "(experiment", id.exp, ")", "\n", "\n")


						#Set working directory for the consortium
						work.dat = set.work.dat(data.consortium = data.consortium, trait.project = trait.project)


						#Load workspace with results
						fres = paste( work.dat, "data/r_workspaces/", trait.project, "/", id.exp,  "/score_results_", trait.project, "_", id.exp,".RData",sep="")

						#Load file
						load(file = fres)

						#Load SNP-gene mapping table
						load( file = paste(work.dat, "data/r_workspaces/", trait.project, "/",trait.project,"_tmap.RData",sep="") )



						cat("Loaded workspace with SNP-gene mapping table and preprocessing results.", "\n")


						#Extract list of significant SNVs for the first combination of clusters, and their p-values
						sig.snp.nmtf <- unlist(get( "sig.snp.nmtf" , res.sd.plot.i ) [[1]] )


						#Extract list of genes harbouring significant SNPS
						gene.cnmtf =  unique( as.character( tmap$entrezgene[ match( names( unlist( sig.snp.nmtf ))[ substr( names( unlist( sig.snp.nmtf )),1,1) == "r" ], tmap$refsnp_id ) ] ) )
						#gene.cnmtf = gene.cnmtf[ !is.na(gene.cnmtf)]

						#If there are no genes delete any old file and go next
						if(length(gene.cnmtf) == 0 ){

							file.annotations =  paste(work.dat, "data/r_workspaces/", trait.project,"/", id.exp, "/results_annotations_", id.exp, ".csv", sep = "" )

							if( file.exists( file.annotations ) ){
								file.remove( file.annotations )
							}

							file.david = paste(work.dat, "data/r_workspaces/", trait.project, "/", id.exp, "/chart_from_david_", trait.project, ".txt", sep = "" )

							if( file.exists( file.david ) ){
								file.remove( file.david )
							}

							next
						}

						if( use.web.service == TRUE ){


							cat("Setting connection to DAVID.", "\n")

							library("RDAVIDWebService")

							#Create object david
							david <- DAVIDWebService(email = email.david, url = "https://david.ncifcrf.gov/webservice/services/DAVIDWebService.DAVIDWebServiceHttpSoap12Endpoint/")

							#Add list of genes
							if( explore.gwas == TRUE){ #Sent all genes found mutated in GWAS catalogue

								addList( object = david, inputIds = lgenes.snp, idType = "ENTREZ_GENE_ID", listType = "Gene", listName = paste( data.consortium, trait.project, "GWAS", sep="-"))

							}else{ #Sent genes found mutated with cNTMF

								addList( object = david, inputIds = gene.cnmtf, idType = "ENTREZ_GENE_ID", listType = "Gene", listName = paste( data.consortium, trait.project, sep="-"))

							}

							cat("Quering table of functional annotations from DAVID.", "\n")

							#Extract functional annotations and print them to a file
							getFunctionalAnnotationTableFile( object = david, file = "table_from_david.tsv")

							#Load results from DAVID
							tdav = read.table( "table_from_david.tsv", header = T, sep = "\t", quote="\"",fill = TRUE)


							cat("Quering chart of functional annotations from DAVID.", "\n")

							#Extract chart of functional annotations
							cdav = getFunctionalAnnotationChart( object = david )
							cdav = as.data.frame(cdav)

							cat("Annotations from DAVID were loaded.", "\n")

						}else{ #If not using WEB service

							#Print list to file
							write.table(unique(gene.cnmtf), paste(work.dat, "data/r_workspaces/", trait.project,"/table_to_david_", trait.project, ".txt", sep = "" ),  row.names = FALSE, col.names = FALSE, quote = FALSE)

							cat("Manually export the list to DAVID, download the annotation chart and save it in the workspace of this trait. ", "\n")

							decision.file = readline("Results from DAVID ready? Press enter to continue.")

							#Load results from DAVID
							tdav = read.table( paste( "data/r_workspaces/", trait.project,"/table_from_david_", trait.project, ".txt", sep = "" ), header = T, sep = "\t", quote="\"", fill = TRUE)

							cat("Annotations from DAVID were loaded.", "\n")


						}



				#------------------------------------------------
				# 4.2 Working with results from DAVID
				#------------------------------------------------


						cat("Creating table of SNPs and gene annotations.", "\n")


						#Declare fields from DAVID to include in the results
						tdav.f = c("ID","Gene.Name", "OMIM_DISEASE", "KEGG_PATHWAY", "GOTERM_BP_DIRECT"  )
						tdav.f = tdav.f[ tdav.f %in% colnames(tdav)]

						#Load list of known SNPs
						load(paste(work.dat,"data/r_workspaces/", trait.project, "/preprocess_",trait.project,".RData", sep=""))


						#------------------------------------------------
						# 4.2.1 Table of mutated genes
						#------------------------------------------------


								#Create table of mutated genes
								tmut = tdav [ match(gene.cnmtf, tdav$ID),  tdav.f]
								tmut$ID <- gene.cnmtf
								tmut$known <- rep("",nrow(tmut))

								#Extract list of known mutated genes
								gene.known = tmap$entrezgene [ tmap$refsnp_id %in% snps.known  ]
								gene.known2 = as.character(tmap.o$entrezgene)

								#Depure list of known mutated genes
								gene.known = unique( gene.known[ !is.na(gene.known) ] )
								gene.known2 = unique( gene.known2[ !is.na(gene.known2) ] )

								#Create column of known mutated genes
								tmut$known [ tmut$ID %in% gene.known2 ] <- "known2"
								tmut$known [ tmut$ID %in% gene.known  ] <- "known1"

								#Add trait and data consortium
								tmut$trait = rep( trait.project, nrow(tmut))
								tmut$consortium = rep( data.consortium, nrow(tmut))

								#Transform table into dataframe
								tmut = as.data.frame(tmut)


								##Add concatenated list of significant SNVs
								for( g in 1:length(tmut$ID) ){

										#List of SNVs in gene g
										snvs.i = tmap$refsnp_id [ tmap$entrezgene == tmut$ID[g] ]

										#List of significant SNVs in gene g
										tmut$snvs[g]	<- paste( snvs.i[ snvs.i %in% unlist(names(sig.snp.nmtf)) ], collapse = ", ")

								}


								#Add table to consolidate table
								tmut.all = rbind.fill( tmut.all, tmut)

								cat("Created table of mutated genes.", "\n")



						#------------------------------------------------
						# 4.2.2 Table of mutated SNPs
						#------------------------------------------------



								#Paste fields to the table of SNPs and pvalues
								tsa = data.frame( dscore = sig.snp.nmtf ) #Table of snp and pvalue
								tsa$known <- rep("",nrow(tsa))
								tsa$known [ names(sig.snp.nmtf) %in% setdiff(snps.known2, snps.known) ] <- "known2"
								tsa$known [ names(sig.snp.nmtf) %in% snps.known ] <- "known1"

								#Add entrez gene
								tsa$entrezgene <- tmap$entrezgene [ match(names(sig.snp.nmtf), tmap$refsnp_id) ]

								#Add p-values for the first combination of clusters
								lpval.sc <- get( "lpval.sc" , res.sd.i )[[1]]
								lsd.sc <- get( "lsd.sc" , res.sd.i )[[1]]
								tsa$pvalue <- lpval.sc [ match( names(sig.snp.nmtf), names(lsd.sc) ) ]


								#Conform table of SNPs and gene annotations
								tsa = cbind( names(sig.snp.nmtf), tsa, tdav [ match( tsa$entrezgene, tdav$ID) ,  tdav.f]) #Add annotations from DAVID

								#Add trait and data consortium
								tsa$trait = rep( trait.project, nrow(tsa))
								tsa$consortium = rep( data.consortium, nrow(tsa))

								#Sort by dscore
								tsa = tsa[ order(tsa$dscore, decreasing = TRUE),]


								#Replace empty cells by "-"
								tsa = apply(tsa, MARGIN = 2, FUN = as.character )
								tsa[ tsa == "" ] <- "-"



								#Add colnames and print to file
								if( is.null(nrow(tsa)) ){ tsa <- as.data.frame(t(tsa)) }else{ tsa <- as.data.frame(tsa) }

								colnames(tsa) <- c("snp", "dscore", "known" , "gene", "pvalue", tolower(tdav.f), "trait", "consortium")
								write.csv(tsa, file = paste(work.dat, "data/r_workspaces/", trait.project,"/", id.exp, "/results_annotations_", id.exp, ".csv", sep = "" ),  row.names = F)

								cat("Printed table results_annotations_id.exp.csv", "\n")

								#Transform table into dataframe
								tsa = as.data.frame(tsa)


								#Add table to consolidate table
								tsa.all = rbind.fill( tsa.all, tsa)



						#------------------------------------------------
						# 4.2.3 Table of enrichment analysis
						#------------------------------------------------


								#Filter chart of anotations
								cdav = cdav[ cdav$Category %in% c("INTERPRO", "GOTERM_BP_DIRECT", "KEGG_PATHWAY", "OMIM_DISEASE"), ]
								cdav = cdav[ which(cdav$PValue <= 0.05) , ]

								#Add trait and data consortium
								cdav$trait = rep( trait.project, nrow(cdav))
								cdav$consortium = rep( data.consortium, nrow(cdav))

								#Print chart to a file
								write.table(cdav, paste(work.dat, "data/r_workspaces/", trait.project, "/", id.exp, "/chart_from_david_", trait.project, ".txt", sep = "" ),  row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")

								cat("Printed chart_from_david_trait.project.txt", "\n")


								#Add chart to consolidate chart
								cdav.all = rbind.fill( cdav.all, cdav)


			} #End loop across traits
		} #End loop across consortiums


		#Replace empty cells by "-"
		tmut.all = apply(tmut.all, MARGIN = 2, FUN = as.character )
		tmut.all[ tmut.all == "" ] <- "-"
		tmut.all <- as.data.frame(tmut.all)

		#Print consolidate chart and table of annotations
		write.table(cdav.all, "summary_david_chart.txt",  row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")
		write.table(tsa.all, "summary_david_table.txt",  row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")
		write.table(tmut.all, "summary_david_table_genes.txt",  row.names = F, col.names = TRUE,  sep = "\t", quote = TRUE)



    #------------------------------------------------
    # 6. Analyzing punctual mutations
    #------------------------------------------------


		#Id of experiment
		id.exp = "5002"

		#Table to consolidated total across consoritums and traits
		tres = NULL

		#------------------------------------------------
		# 6.1 Create information table of SNVs
		#------------------------------------------------


				#Iterate across consortiums and traits
				for( l in 1:length(lconsortiums)){
					for( k in 1:length(ltrait.project)){


						#Extract consortium, trait and experiment id
						data.consortium = lconsortiums[l]
						trait.project = ltrait.project[k]

						#Print script status
						cat("\n", "\n", rep("-",20), "\n", "\n")
						cat("Analysing SNVs raw data from ", data.consortium, "in trait" , trait.project, "\n", "\n")


						#Set working directory for the consortium
						work.dat = set.work.dat(data.consortium = data.consortium, trait.project = trait.project)


						#------------------------------------------------
						# 6.1.1 Create R matrix
						#------------------------------------------------


								#Load workspace with preprocessing variables
								load(paste(work.dat,"data/r_workspaces/", trait.project, "/preprocess_",trait.project,".RData", sep=""))
								load(paste(work.dat, "data/r_workspaces/", trait.project, "/",trait.project,"_tmap",".RData",sep=""))

								source("n_intro.R")

								#Load list of SNPs and Patients included after filtering R
								load( file = paste(work.dat,"data/r_workspaces/", trait.project, "/R.snps.s_", trait.project, ".RData", sep=""))
								R.snps <- R.snps.s
								load( file = paste(work.dat,"data/r_workspaces/", trait.project, "/R.pats.s_", trait.project, ".RData", sep=""))
								R.pats <- R.pats.s
								out.cat = out.cat[ colnames(R) %in% R.pats ]
								pop = pop[ colnames(R) %in% R.pats ]
								R = R[rownames(R) %in% R.snps, colnames(R) %in% R.pats]
								dim(R)


						#------------------------------------------------
						# 6.1.2 Fill information table of SNVs
						#------------------------------------------------

								#Check if base table of annotations exists. If not, then remove any old results file and go next
								file.annotations = paste(work.dat, "data/r_workspaces/", trait.project,"/", id.exp, "/results_annotations_", id.exp, ".csv", sep = "" )
								if( !file.exists(file.annotations) ) {

									file.results = paste(work.dat, "data/r_workspaces/", trait.project,"/explore_results", ".csv", sep = "" )

									if( file.exists(file.results )  ){
										file.remove( file.results )
									}

									next

								}

								#Read base table of annotations
								tfreq.snvs = read.csv(file = file.annotations)
								tfreq.snvs = as.data.frame(tfreq.snvs, stringsAsFactors = TRUE)

								#Fill table
								for(i in 1:nrow(tfreq.snvs)){

									pos.snv = which(tmap$refsnp_id %in% tfreq.snvs$snp[i])[1] #First position for SNV i
									maf.snv = sum( R[ rownames(R) %in% tfreq.snvs$snp[i], ] ) / ((ncol(R)*2))
									tfreq.snvs$alleleA [i] <- tmap$alleleA[ pos.snv ] #Allele A
									tfreq.snvs$alleleB [i] <- tmap$alleleB[ pos.snv ] #Allele B
									tfreq.snvs$freq.A [i]	<- 1 - maf.snv #Allele A MAF
									tfreq.snvs$freq.B [i]	<- maf.snv #Allele B MAF
									tcount = c( table(  factor( R[ rownames(R) %in% tfreq.snvs$snp[i],], levels = c("0","1","2") ), out.cat )  )
									tfreq.snvs$low.trait = paste(tcount[1:3], collapse = ",")
									tfreq.snvs$high.trait = paste(tcount[4:6], collapse = ",")

								}


								#Add GWAS catalog data for known associations
								tfreq.snvs = cbind(tfreq.snvs,
																	 dseed.s[ match( tfreq.snvs$snp, dseed.s$SNPS), c("FIRST.AUTHOR","INITIAL.SAMPLE.SIZE", "STRONGEST.SNP.RISK.ALLELE","RISK.ALLELE.FREQUENCY", "P.VALUE","OR.or.BETA", "X95..CI..TEXT.") ]
								)

								#Add other mapped traits to these SNPs
								dseed = read.table(paste(work.dat,"data/other_data/gwas_catalog_v1.0.1.tsv", sep=""), header=T,sep="\t",quote="\"",fill = TRUE,colClasses ="character")
								ctraits.snv = matrix(NA, nrow(tfreq.snvs), ncol = 2)

								#Specific traits to find match
								toMatch = c( "Hypertriglyceridemia","lipid measurement","lipoprotein measurement, blood metabolite measurement",
														 "response to statin","sphingolipid measurement","total cholesterol measurement", "high density lipoprotein cholesterol measurement",
														 "low density lipoprotein cholesterol measurement", "triglyceride measurement")


								for(i in 1:nrow(tfreq.snvs)){
									ltraits.snv = dseed$MAPPED_TRAIT [ dseed$SNPS %in% tfreq.snvs$snp[i] & dseed$MAPPED_TRAIT %in% toMatch ]
									ctraits.snv[i,] = c( length(ltraits.snv),
																			 paste( ltraits.snv, collapse = "; " )  )

								}
								tfreq.snvs = cbind(tfreq.snvs, ctraits.snv)


								#Find other SNVs in high LD with known associations within the same gene
								file.LD = paste(work.dat,"data/r_workspaces/", trait.project, "/linkage_table_",trait.project,"_", substr(id.exp,1,1),".RData", sep="")
								ld.tao = 0.8 #Treshold of LD
								load(file.LD)


								#Add SNVs in high LD
								for(i in 1:nrow(tfreq.snvs)){

									tfreq.snvs$LD.snps [ i ] =	paste( unique( c( as.character( res.ld$loc2[ res.ld$r2 > ld.tao & as.character(res.ld$loc1) %in% as.character(tfreq.snvs$snp[i]) ] ) ,
																															 as.character( res.ld$loc1[ res.ld$r2 > ld.tao & as.character(res.ld$loc2) %in% as.character(tfreq.snvs$snp[i]) ] )) ) ,
																										collapse = "; ")

								}


								#Add information of SNP consequences from ENSEMBL
								ensembl.snp <- useMart(biomart="ENSEMBL_MART_SNP", dataset = "hsapiens_snp", host = "grch37.ensembl.org" ) #listAttributes(ensembl.snp)
								snp.conseq = getBM(attributes=c("refsnp_id", "ensembl_transcript_stable_id","consequence_type_tv",  "polyphen_prediction", "sift_prediction"), filters="snp_filter", values=as.character(tfreq.snvs$snp), mart=ensembl.snp)

								#Load catalogue of consequences
								dconseq = read.table( paste(work.dat, "data/other_data/ensembl_consequences.txt",sep = "" ), header = TRUE, sep = "\t")

								#Find the most severe consequence
								for(i in 1:nrow(tfreq.snvs)){

									pos.consq = which(snp.conseq$refsnp_id %in% tfreq.snvs$snp[i])
									tfreq.snvs$consequence_type_tv[i] =  as.character( dconseq$consequence[ which( as.character(dconseq$consequence) %in% snp.conseq$consequence_type_tv[pos.consq] )] [1] )

								}


								#Add information of MAF in global population
								ncbi.query = ncbi_snp_query(SNPs = tfreq.snvs$snp)
								tfreq.snvs$global.maf = ncbi.query$MAF[ match(tfreq.snvs$snp, ncbi.query$Marker) ]

								#Add information of number of publications
								for(i in 1:nrow(tfreq.snvs)){

									annot.snvs = annotations( tfreq.snvs$snp [ i ], output = "all" )
									tfreq.snvs$npubs[ i ] =  length( !is.na( unique( annot.snvs$author) ) )

								}


							#Add table to consolidated total across consortiums and traits
							tres = rbind.fill(tfreq.snvs, tres)
							dim(tfreq.snvs)
							colnames(tres)

							#Clean empity cells
							tfreq.snvs = apply(tfreq.snvs, MARGIN = 2, FUN = as.character )
							tfreq.snvs[ is.na(tfreq.snvs) ] <- "-"
							tfreq.snvs[ tfreq.snvs == "" ] <- "-"

							#Print results to file
							write.csv(tfreq.snvs, file = paste(work.dat, "data/r_workspaces/", trait.project,"/explore_results", ".csv", sep = "" ),  row.names = F)



					} #End loop across traits
				} #End loop across consortiums



				#Find the most severe consequence among variants in high LD
				tres = as.data.frame(tres, stringsAsFactors = TRUE)

				for(i in 1:nrow(tres)){


					tres$consequence_type_tv2[i] <-  	tres$consequence_type_tv[i]
					tres$snp_cons_tv2[i] = ""

					if( tres$LD.snps [i] != ""){

						lsnps.ld =  unlist( strsplit( x = as.character(tres$LD.snps [ i ]), split =  "; "))
						lsnps.ld = lsnps.ld[ !(lsnps.ld %in% tres$snp [ tres$consortium == tres$consortium[i] &  tres$trait == tres$trait[i]  ] ) ]

						if( length(lsnps.ld ) > 0 ){


							snp.conseq = getBM(attributes = c("refsnp_id", "ensembl_transcript_stable_id","consequence_type_tv",  "polyphen_prediction", "sift_prediction"), filters="snp_filter",
																 values = lsnps.ld, mart=ensembl.snp)

							tres$consequence_type_tv2[i] =  as.character( dconseq$consequence[ which( as.character(dconseq$consequence) %in% c ( tres$consequence_type_tv[i] , snp.conseq$consequence_type_tv ) )] [1] )
							tres$snp_cons_tv2[i] =  paste( unique( snp.conseq$refsnp_id [ which(  snp.conseq$consequence_type_tv %in% tres$consequence_type_tv2[i] ) ] ), collapse = ", ")

						}

					}

				}


				#Print results to file
				write.csv(tres, file = "explore_results.csv",  row.names = F)

			return()

  } #End function



