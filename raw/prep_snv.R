###############################################################
# cNMTF
#	1. Preprocessing functions
# 1.3 Functions to select SNVs and give format to genotypes
#
# Corresponding author:
# Luis Leal, Imperial College London, email: lgl15@imperial.ac.uk
###############################################################


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#' @title Select sample of SNVs
#' @description Function to retrieve list of seed SNVs and expanded list of genes
#'
#[Input]
#' @param seed.file.main  GWAS catalogue file
#' @param seed.file.sab  Table with SNPs found by Sabatti et al
#' @param seed.file.sab  Table with SNPs found by Talmut et al.
#' @param ppi.file  PPIN from BIOGRID
#' @param lauthors  List of authors from specific studiess
#' @param log.file File to save counts of SNPs
#'
#[Output]
#' @return A workspace named "\code{seed.file.results}" with the following objects:
#' * \code{lgenes.exp.attr} Expanded list of genes from PPIN
#' * \code{dseed.s} Depurated list of SNPs
#' * \code{tmap.o} Table to map seed SNPs to genes
#' * \code{lgenes.g.id} List of genes mapped from seed SNPs
#'
#' @md
#' @author Luis G. Leal, \email{lgl15@@imperial.ac.uk}
#' @family Preprocessing functions
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


seed.snps = function(  seed.file.main =  NULL,
											 seed.file.sab = NULL,
											 seed.file.tal = NULL,
											 ppi.file = NULL,
											 lauthors = NULL,
											 seed.file.results = NULL,
											 data.consortium = NULL,
											 add.random.snps = TRUE,
											 utrait = NULL,
											 trait = NULL)
	{

		#-------------------------------------------------
		# 1. Conform set of seed SNPs from GWAS catalogue
		#-------------------------------------------------


				#Load seeds of SNPs to be expanded
				dseed = read.table(seed.file.main, header=T,sep="\t",quote="\"",fill = TRUE,colClasses ="character")

				#Complement the seed list with all the SNPs found by Sabatti
				dseed.sab = read.table(seed.file.sab, header=T,sep="\t",quote="\"",fill = TRUE,colClasses ="character")
				dseed = rbind.fill( data.frame(dseed), data.frame(dseed.sab))

				#Complement the seed list with all the SNPs found by Talmud
				dseed.tal = read.table(seed.file.tal, header=T,sep="\t",quote="\"",fill = TRUE,colClasses ="character")
				dseed = rbind.fill( data.frame(dseed), data.frame(dseed.tal))


				#Filter out the seeds
				dseed.s = dseed[ dseed$MAPPED_TRAIT %in% trait, ] # Filter by trait
				dseed.s = dseed.s[ as.numeric(dseed.s$P.VALUE) <= 1e-8 | dseed.s$FIRST.AUTHOR %in% lauthors,] #Filter by p-value and keep Sabatti findings
				dseed.s = dseed.s[ !(dseed.s$CHR_ID%in%c("X","")),] # Remove chromosome X
				dseed.s = dseed.s[ match( unique(dseed.s$SNPS), dseed.s$SNPS ), ] # Remove repeated entries

				#Number of SNPs from the study of Sabatti (some SNPs are repeated between studies)
				cat( "Number of SNPs from specific publications:", sum(dseed.s$SNPS %in% dseed$SNPS [ dseed$FIRST.AUTHOR %in% lauthors & dseed$MAPPED_TRAIT %in% trait]), "\n")

				#Dimensions
				cat( "Number of SNP seeds: ", nrow(dseed.s), "\n")



		#-------------------------------------------------
		# 2. Map seed SNPs to genes
		#-------------------------------------------------


				#List of mapped and reported genes for the SNPs. Unlist the entries.
				lgenes.g = c(dseed.s$MAPPED_GENE, dseed.s$REPORTED.GENE.S.)
				lgenes.g = unlist(strsplit(lgenes.g,split = " - "))
				lgenes.g = unlist(strsplit(lgenes.g,split = ", "))
				lgenes.g = unique(lgenes.g)
				lgenes.g = lgenes.g[!is.na(lgenes.g)]

				cat( "Number of genes mapped from SNP seeds: ", length(lgenes.g), "\n")


				#Table that maps SNPs to gene entrez ids for SNP seeds
				tmap.o = NULL
				for(i in 1:nrow(dseed.s)){

					#Extract entrezid from upstream, downstrean and mapped genes
					lgenes.snp = c(dseed.s$UPSTREAM_GENE_ID[i], dseed.s$DOWNSTREAM_GENE_ID[i], dseed.s$SNP_GENE_IDS[i] )
					lgenes.snp = unlist(strsplit(lgenes.snp,split = " - "))
					lgenes.snp = unlist(strsplit(lgenes.snp,split = ", "))
					lgenes.snp = unique(lgenes.snp)
					lgenes.snp = lgenes.snp[!is.na(lgenes.snp)]
					nl = length(lgenes.snp)

					#If any gene was mapped for that SNP seed, then add the mapping to the table
					#Each SNP could map multiple genes because we are including closet upstream and downstream genes for SNPs in intronic regions
					if( nl > 0){
						tmap.o = rbind( tmap.o, cbind( rep(dseed.s$SNPS[i], nl), lgenes.snp, rep(dseed.s$CHR_ID[i], nl), rep(dseed.s$CHR_POS[i], nl) ) )
					}else{
						tmap.o = rbind( tmap.o, cbind( dseed.s$SNPS[i], NA, dseed.s$CHR_ID[i], dseed.s$position[i] ) )
					}

				}


				#Add colnames to the table
				tmap.o = as.data.frame(tmap.o)
				colnames(tmap.o) <- c("refsnp_id","entrezgene","chr","position")



		#-------------------------------------------------
		# 3. Expand set of genes from interactome
		#-------------------------------------------------


				#Number of chromosomes
				nchr = 22 #Chromosome 23 is not in the imputed data

				#Load PPI network data to expand the seed of SNPs
				dppi = read.table(ppi.file, header=T,sep="\t",quote="\"",fill = TRUE,colClasses ="character")
				dppi = dppi[dppi$Organism.Interactor.A == "9606", ] #This is the NCBI Taxonomy ID for Interactor A.

				#Retain columns of interactions: Entrez gene ID for interactor A and B respectively
				dppi = dppi[,2:3]


			  cat( "Number of interactions in the PPI: ", nrow(dppi), "\n")

				#Map the list of genes to entrezids
				ensembl <- useMart(biomart="ENSEMBL_MART_ENSEMBL", dataset="hsapiens_gene_ensembl", host = "grch37.ensembl.org" )
				lgenes.g.id <- unlist(getBM(attributes=c("entrezgene"), filters="external_gene_name", values=lgenes.g, mart=ensembl))

				#Expand the list of genes by adding upstream and downstream genes (entrezids from GWAS catalogue)
				lgenes.g.id.ud = unique(c(dseed.s$UPSTREAM_GENE_ID,dseed.s$DOWNSTREAM_GENE_ID))
				lgenes.g.id.ud = lgenes.g.id.ud[!is.na(lgenes.g.id.ud)]
				lgenes.g.id.ud = lgenes.g.id.ud[lgenes.g.id.ud != ""]
				lgenes.g.id = unique(c(lgenes.g.id,lgenes.g.id.ud))

				cat( "Number of genes mapped from SNP seeds (including upstream and downstream): ", length(lgenes.g.id), "\n")


				#Expand the list of genes by using the interactome
				lgenes.exp.id = unique( c(lgenes.g.id, dppi[ dppi[,1]%in%lgenes.g.id, 2], dppi[ dppi[,2]%in%lgenes.g.id, 1] ))


				#Retrieve the chromosome position of each gene candidate
				lgenes.exp.attr <- getBM(attributes=c("entrezgene","chromosome_name","start_position","end_position","external_gene_name"), filters="entrezgene", values=lgenes.exp.id, mart=ensembl)

				#Filter out by chromosomes
				lgenes.exp.attr = lgenes.exp.attr[lgenes.exp.attr$chromosome_name%in%c(1:nchr),]

				cat( "Number of genes expanded using the interactome: ", length(unique(lgenes.exp.attr$entrezgene)), "\n")


				#Obtain known associations
				if( data.consortium %in% c("NFBC","NFBC2")){
					snps.known = dseed.sab$SNPS[dseed.sab$MAPPED_TRAIT %in% trait]
				}else{
					snps.known = NULL
				}

				#List of known SNPS from GWAS catalog
				snps.known2 = dseed.s$SNPS[dseed.s$MAPPED_TRAIT %in% trait]


			#Save workspace with results
				save(list = c("lgenes.exp.attr","dseed.s", "tmap.o", "lgenes.g.id", "snps.known", "snps.known2"), file =  seed.file.results)



	}



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#' @title Select genetic model for the genotyping data
#' @description Function to change coding scheme in a relationship matrix
#'
#[Input]
#' @param R Relationship matrix
#' @param coding Desired coding scheme (genetic model) for the entries of \code{R}
#'
#[Output]
#' @return Relationship matrix \code{R} with a given coding scheme
#'
#' @md
#' @author Luis G. Leal, \email{lgl15@@imperial.ac.uk}
#' @family Preprocessing functions
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


coding.scheme <- function(R,
													coding = c("recessive","dominant"))
	{
			colnames(R) <- paste("snp",rep(1:ncol(R)),sep=".")

			if(coding == "recessive"){

				R[R == "1"] <- "0"
				R[R == "2"] <- "1"

			}else if(coding == "dominant"){

				R[R == "2"] <- "1"

			}
	return(R)
}




#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Function to calculate linkage disequilibrium between pairs of SNPs
#
#[References] Following guidelines on: https://snpinfo.niehs.nih.gov/cgi-bin/snpinfo/snptag.cgi
#[Input]
# - R : Relationship matrix
# - file.LD :	#File name for Linkage Disequilibrium table
# - min.maf : #Minimum Minor allele frequency
# - tmap : Table to map SNPs to genes
# - type.ld : c("gene","all") : Region to find SNPs in LD
#[Output]
# - Workspace with res.ld (table of pairwise LD)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



find.snps.ld = function( file.LD = NULL,	#File name for Linkage Disequilibrium table
												 type.ld = c("gene","all"), #Region to find SNPs in LD
												 parallel.opt = FALSE, #Run in Parallel
												 min.maf = 0.05, #Minimum Minor allele frequency
												 tmap = NULL, #Table to map SNPs to genes
												 R = NULL #Relationship matrix
){

	cat("Finding SNPs in high LD","\n", as.character(Sys.time()), "\n")

	#Define setting of algorithms
	maf.snps = rowSums(R)/((ncol(R)*2))

	#List of SNPs
	R.snps = rownames(R)

	#Extract list of genes
	lgene = unique(tmap$entrezgene[ tmap$refsnp_id %in% R.snps])
	length(lgene)

	if(type.ld == "gene"){

		#Declare table of results
		res.ld = NULL

		#Find LD among SNPs in the same gene region

		for(i in lgene){

			#Find SNPs in the gene
			pos.snps = which( ( R.snps %in% tmap$refsnp_id[tmap$entrezgene == i] ) & maf.snps > min.maf )

			if(length(pos.snps) > 1){

				#Calculate linkage disequilibrium between pairs of SNPs

				cat("Calculating LD for gene", which( lgene == i ), "out of", length(lgene), "\n")


				if(parallel.opt == TRUE){ #Parallel computing
					res.ld.chr = ld.table (R = R[pos.snps, ], #Relationship matrix
																 n.cores = 50, #Number of cores for parallel processing
																 step.width = 5e4, #Number of operations per cluster
																 option.parallel = 2)


				}else{

					res.ld.chr = LD.Measures(donnees = t(R[pos.snps, ]), V = NA, S = NA, data = "G", supinfo = FALSE, na.presence=TRUE)

				}


				#Merge with final table
				res.ld = rbind(res.ld, res.ld.chr)

			}


		}

	}else if(type.ld == "all"){


		#Calculate linkage disequilibrium between pairs of SNPs

		source("f_score.R")
		dim(R)
		res.ld = ld.table (R = R, #Relationship matrix
											 n.cores = 50, #Number of cores for parallel processing
											 step.width = 5e4, #Number of operations per cluster
											 print.file = print.file, #File to print results
											 option.parallel = 2)

	}


	#Save workspace with table
	save(list = c("res.ld"), file = file.LD)
	cat("Workspace with table of results saved to ", file.LD, "\n", as.character(Sys.time()), "\n")

}



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Function to calculate LD in PARALLEL
#[Input]
# - R : Relationship matrix
# - ncores : Number of cores for parallel processing
# - step.width : Number of operations per cluster
# - print.file : File to print results
# - option.parallel : Package for parallel processing
#[Output]
# - Workspace with res.ld (table of pairwise LD)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


ld.table = function(R = NULL, #Relationship matrix
										n.cores = 4, #Number of cores for parallel processing
										step.width = 1e3, #Number of operations per cluster
										print.file = NULL, #File to print results
										option.parallel = c(1,2) ) #Package for parallel processing: 1) mclapply 2)
{

	cat("Creating table of LD results", "\n", as.character(Sys.time()), "\n")

	#Create table z with indexes of rows
	n = nrow(R)
	z = matrix(NA, nrow = n*(n-1)/2, ncol = 2)
	z[,1] <- unlist( mclapply( 1:n, function(i) rep(i,n-i)) )
	z[,2] <- unlist( mclapply( 2:n, function(i) seq(i,n)) )

	#Create array x to save LD measures
	x = rep(NA, n*(n-1)/2)

	#Declare number of cores to work in option.parallel = 2 (the function for option.parallel = 1 has the argument mc.cores)
	if(option.parallel == 2){
		doMC::registerDoMC(cores = n.cores)
	}

	#Declare step: number of operations to distribute through cores
	step = min( step.width, nrow(z) )

	#Counter variable
	k = 1

	#Transpose of imput
	tR = t(R)

	#Fill table z
	while( k < nrow(z) ){

		#End position
		end = min( k + step, nrow(z) )

		#Calculate the LD measure for each pair of SNPs
		#Option parallel 1:
		if(option.parallel == 1){

			res.ld <- mclapply( k:end, function(i) { LD.Measures(donnees = tR[, z[i,] ] ) }, mc.cores = n.cores)

		}else if(option.parallel == 2){

			#Option parallel 2 (slower 1s per 1000):
			res.ld <- foreach(i = k:end) %dopar% ( LD.Measures(donnees = tR[, z[i,] ] ) )

		}

		#Extract results and save them in array x
		res.ld = unlist(lapply(res.ld, '[[', 3))
		x[ k:end ] <- res.ld


		cat("Calculated LD for ",end, "pairs of SNPs out of", nrow(z), "\n", as.character(Sys.time()), "\n")


		#Update counter variable
		k = k + step + 1

	}

	#Create table of results
	res.ld = cbind(rownames(R)[z[,1]], rownames(R)[z[,2]], x)
	colnames(res.ld) <- c("loc1","loc2","r2")
	res.ld = data.frame(res.ld, stringsAsFactors = FALSE)
	res.ld$r2 = as.numeric(res.ld$r2)

	return(res.ld)

}





