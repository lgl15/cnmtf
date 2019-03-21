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
#' @param add.david.annotations Logical. Use DAVID web service or export/import manually
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


  annotate.results = function( name.exp = NULL, #Define experiment id
  														 snps.known = NULL, #Optional. List of SNVs known to be associated with the disease
  														 snps.known2 = NULL, #Optional. A second list of SNVs.
  														 add.david.annotations = TRUE, #Use DAVID web service or export/import manually
                               email.david = NULL, #Email account registered in DAVID.
  														 add.ensemble.conseq = FALSE, #Add SNV consequences from ENSEMBL
  														 work.dat = NULL, #Working directory
  														 tmap = NULL, #Mapping of SNPs to genes, chr and genomic position
  														 file.LD = NULL, #Workspace with a table of pairwise LD
  														 ld.tao = 0.8 #Treshold of LD

  )
  {




			#------------------------------------------------
			# 1. Reading list of prioritised variants
			#------------------------------------------------


						cat("\n", "\n", rep("-",20), "\n", "\n")
						cat("Annotating results of experiment", name.exp, "in trait", trait.project, "\n", "\n")


						#Load workspace with results
						fres = paste( work.dat, "score_results_", name.exp,".RData",sep="")

						#Load file
						load(file = fres)


						#Extract list of significant SNVs for the first combination of clusters, and their p-values
						sig.snp.nmtf <- unlist(get( "sig.snp.nmtf" , res.sd.plot.i ) [[1]] )


						#Extract list of genes harbouring significant SNPS
						gene.cnmtf =  unique( as.character( tmap$entrezgene[ match( names( unlist( sig.snp.nmtf ))[ substr( names( unlist( sig.snp.nmtf )),1,1) == "r" ], tmap$refsnp_id ) ] ) )


						#If there are no genes delete any old file and go next
						if(length(gene.cnmtf) == 0 ){

								file.annotations =  paste(work.dat, "results_annotations_", name.exp, ".csv", sep = "" )

								if( file.exists( file.annotations ) ){
									file.remove( file.annotations )
								}

								file.david = paste(work.dat, "chart_from_david_", trait.project, ".txt", sep = "" )

								if( file.exists( file.david ) ){
									file.remove( file.david )
								}

						next
						}



			#------------------------------------------------
			# 2. Importing annotations from DAVID
			#------------------------------------------------


						#If using DAVID web service to retrieve annotations of mutated genes

						if( add.david.annotations == TRUE ){

							cat("Setting connection to DAVID.", "\n")

							#Create object david
							david <- DAVIDWebService(email = email.david, url = "https://david.ncifcrf.gov/webservice/services/DAVIDWebService.DAVIDWebServiceHttpSoap12Endpoint/")

							#Add list of genes
							addList( object = david, inputIds = gene.cnmtf, idType = "ENTREZ_GENE_ID", listType = "Gene", listName = trait.project)


							cat("Quering table of functional annotations from DAVID.", "\n")

							#Extract functional annotations and print them to a file
							getFunctionalAnnotationTableFile( object = david, file = paste(work.dat, "table_from_david.tsv", sep = ""))

							#Load results from DAVID
							tdav = read.table( paste(work.dat, "table_from_david.tsv", sep = ""), header = T, sep = "\t", quote="\"",fill = TRUE)


							cat("Quering chart of functional annotations from DAVID.", "\n")

							#Extract chart of functional annotations
							cdav = getFunctionalAnnotationChart( object = david )
							cdav = as.data.frame(cdav)

							cat("Annotations from DAVID were loaded.", "\n")

							cat("Creating table of SNPs and gene annotations.", "\n")


							#Declare fields from DAVID to include in the results
							tdav.f = c("ID","Gene.Name", "OMIM_DISEASE", "KEGG_PATHWAY", "GOTERM_BP_DIRECT"  )
							tdav.f = tdav.f[ tdav.f %in% colnames(tdav)]

						}



			#------------------------------------------------
			# 3. Creating table of prioritised genes
			#------------------------------------------------


						#Reading annnotations
						if( add.david.annotations == TRUE ){

							tmut = tdav [ match(gene.cnmtf, tdav$ID),  tdav.f]
							tmut$ID <- gene.cnmtf


						}else{

									tmut =  as.data.frame(gene.cnmtf)
									colnames(tmut) <- c("ID")
						}

						#Add column for known associations
						tmut$known <- rep("",nrow(tmut))

						#Extract list of known mutated genes
						gene.known = tmap$entrezgene [ tmap$refsnp_id %in% snps.known  ]
						gene.known2 = tmap$entrezgene [ tmap$refsnp_id %in% snps.known2  ]

						#Depure list of known mutated genes
						gene.known = unique( gene.known[ !is.na(gene.known) ] )
						gene.known2 = unique( gene.known2[ !is.na(gene.known2) ] )

						#Create column of known mutated genes
						tmut$known [ tmut$ID %in% gene.known2 ] <- "known2"
						tmut$known [ tmut$ID %in% gene.known  ] <- "known1"

						#Add trait and data consortium
						tmut$trait = rep( trait.project, nrow(tmut))

						#Transform table into dataframe
						tmut = as.data.frame(tmut)


						##Add concatenated list of significant SNVs
						for( g in 1:length(tmut$ID) ){

								#List of SNVs in gene g
								snvs.i = tmap$refsnp_id [ as.character(tmap$entrezgene) == tmut$ID[g] ]

								#List of significant SNVs in gene g
								tmut$snvs[g]	<- paste( snvs.i[ snvs.i %in% unlist(names(sig.snp.nmtf)) ], collapse = ", ")

						}


						#Replace empty cells by "-"
						tmut = apply(tmut, MARGIN = 2, FUN = as.character )
						tmut[ tmut == "" ] <- "-"
						tmut <- as.data.frame(tmut)

						#Print consolidate chart and table of annotations
						write.table(tmut, paste(work.dat, "prioritised_genes_", name.exp, ".txt", sep = "" ),  row.names = F, col.names = TRUE,  sep = "\t", quote = TRUE)


						cat("Created table of mutated genes.", "\n")



			#------------------------------------------------
			# 4. Creating table of prioritised SNPs
			#------------------------------------------------


						#Paste fields to the table of SNPs and pvalues
						tsa = NULL
						tsa = data.frame( snp = names(sig.snp.nmtf), dscore = sig.snp.nmtf ) #Table of snp and pvalue
						tsa$known <- rep("",nrow(tsa))
						tsa$known [ names(sig.snp.nmtf) %in% setdiff(snps.known2, snps.known) ] <- "known2"
						tsa$known [ names(sig.snp.nmtf) %in% snps.known ] <- "known1"

						#Add entrez gene
						tsa$entrezgene <- tmap$entrezgene [ match(names(sig.snp.nmtf), tmap$refsnp_id) ]


						#Add p-values for the first combination of clusters
						lpval.sc <- get( "lpval.sc" , res.sd.i )[[1]]
						lsd.sc <- get( "lsd.sc" , res.sd.i )[[1]]
						tsa$pvalue <- lpval.sc [ match( names(sig.snp.nmtf), names(lsd.sc) ) ]


						#Add trait
						tsa$trait = rep( trait.project, nrow(tsa))

						#Add annotations from DAVID
						if( add.david.annotations == TRUE ){
							tsa = cbind(tsa, tdav [ match( tsa$entrezgene, tdav$ID) ,  tdav.f]) #Add annotations from DAVID
						}


						#Sort by dscore
						tsa = tsa[ order(tsa$dscore, decreasing = TRUE),]


						#Find other SNVs in high LD with known associations within the same gene
						load(file.LD)


						#Add SNVs in high LD
						for(i in 1:nrow(tsa)){

							tsa$LD.snps [ i ] =	paste( unique( c( as.character( res.ld$loc2[ res.ld$r2 > ld.tao & as.character(res.ld$loc1) %in% as.character(tsa$snp[i]) ] ) ,
																													 as.character( res.ld$loc1[ res.ld$r2 > ld.tao & as.character(res.ld$loc2) %in% as.character(tsa$snp[i]) ] )) ) ,
																								collapse = "; ")

						}

						#Replace empty cells by "-"
						tsa = apply(tsa, MARGIN = 2, FUN = as.character )
						tsa[ tsa == "" ] <- "-"



						#Add information of SNP consequences from ENSEMBL
						if( add.ensemble.conseq == TRUE){

							cat("Querying information ")

							#Create ensemble object
							ensembl.snp <- useMart(biomart="ENSEMBL_MART_SNP", dataset = "hsapiens_snp", host = "grch37.ensembl.org" ) #listAttributes(ensembl.snp)
							snp.conseq = getBM(attributes=c("refsnp_id", "ensembl_transcript_stable_id","consequence_type_tv",  "polyphen_prediction", "sift_prediction"), filters="snp_filter", values=as.character(tsa$snp), mart=ensembl.snp)

							#Load catalogue of consequences
							dconseq = read.table( paste(work.dat, "data/other_data/ensembl_consequences.txt",sep = "" ), header = TRUE, sep = "\t")

							#Find the most severe consequence
							for(i in 1:nrow(tsa)){

								pos.consq = which(snp.conseq$refsnp_id %in% tsa$snp[i])
								tsa$consequence_type_tv[i] =  as.character( dconseq$consequence[ which( as.character(dconseq$consequence) %in% snp.conseq$consequence_type_tv[pos.consq] )] [1] )

							}

						}


						#Transform, to dataframe (even if this a null table)
						if( is.null(nrow(tsa)) ){ tsa <- as.data.frame(t(tsa)) }else{ tsa <- as.data.frame(tsa) }


						#Print to file
						write.csv(tsa, file = paste(work.dat, "prioritised_snvs_", name.exp, ".csv", sep = "" ),  row.names = F)

						cat("Printed table results_annotations_name.exp.csv", "\n")



			#------------------------------------------------
			# 5. Table of enrichment analysis
			#------------------------------------------------

					if( add.david.annotations == TRUE ){

							#Filter chart of anotations
							cdav = cdav[ cdav$Category %in% c("INTERPRO", "GOTERM_BP_DIRECT", "KEGG_PATHWAY", "OMIM_DISEASE"), ]
							cdav = cdav[ which(cdav$PValue <= 0.05) , ]

							#Add trait and data consortium
							cdav$trait = rep( trait.project, nrow(cdav))

							#Print chart to a file
							write.table(cdav, paste(work.dat,"enrichement_analyisis_", trait.project, ".txt", sep = "" ),  row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")

							cat("Printed chart of enrichment analysis to chart_from_david_trait.project.txt", "\n")

					}

		return( list(tsa = tsa, tmut = tmut))

	}


