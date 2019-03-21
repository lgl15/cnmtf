###############################################################
# cNMTF
#	1. Preprocessing functions
# 1.3 Functions to select SNVs and give format to genotypes
#
# Corresponding author:
# Luis Leal, Imperial College London, email: lgl15@imperial.ac.uk
###############################################################




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





