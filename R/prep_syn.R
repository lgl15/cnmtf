###############################################################
# cNMTF
#	1. Preprocessing functions
# 1.5 Functions to generate synthetic genotyping data
#
# Corresponding author:
# Luis Leal, Imperial College London, email: lgl15@imperial.ac.uk
###############################################################


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#' @title Simulate synthetic genotypes
#' @description Function to create synthetic GWAS data sets
#' @references
#'  * Mingyao Li, \emph{et al.} (2010) "Correcting population stratification in GWAS using a phylogenetic approach" Bioinformatics
#'  * Alkes Price, \emph{et al.} (2006) "Principal components analysis corrects for stratification in genome-wide association studies" Nature genetics
#'
#[Input]
#' @param   n.snps number of SNPs in the synthetic dataset
#' @param   n.pats number of patients in the synthetic dataset
#' @param   n.pops number of patient populations in the synthetic dataset
#' @param   p.pops proportion of patients in each population
#' @param   c.pops proportion of cases in each patient population
#' @param   allele.freq.model  define the model to generate the allele frequency \code{ c("Balding-Nichols","fitted") }.
#' @param   fst Factor of differentiation between populations if \code{ allele.freq.model == "Balding-Nichols" }.
#' @param   pl.pop vector of fitted allele frequencies in each population if \code{ allele.freq.model == "fitted" }.
#' @param   r risk factor vector ( \emph{i.e.}, genotype relative risks) to generate the genotypes for cases/controls in each population
#' @param   pcases proportion of cases
#' @param   prefix.snp prefix to name the SNP columns
#'
# [Output]:
#' @return
#' * \code{db.g} synthetic GWAS data set
#' * \code{status} disease status of the patient
#' * \code{pop} population ancestry of the patient
#' * \code{snp.risk} genotype relative risk for each SNP
#'
#' @md
#' @author Luis G. Leal, \email{lgl15@@imperial.ac.uk}
#' @family Factorisation functions
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


		synthetic.gwas = function(n.snps = 10, #number of SNPs in the synthetic dataset
		                          n.pats = 100, #number of patients in the synthetic dataset
		                          n.pops = 3, #number of patient populations in the synthetic dataset
		                          p.pops = c(0.25,0.25,0.5),#proportion of patients in each population
		                          c.pops = c(0.2,0.2,0.2),#proportion of cases in each patient population
		                          r.pops = c(1.5,1,1.5), #risk factor vector (i.e., genotype relative risks for the SNPs causative of disease)
		                          fst = 0.01, #Factor of differentiation between populations
		                          allele.freq.model = c("Balding-Nichols","fitted"),
		                          pl.pop = c(0.1,0.1,0.9), #vector of fitted allele frequencies in each population
		                          prefix.snp = "snp"

		){

		  #Check inputs
		  if (n.pops!=length(pl.pop) & allele.freq.model == "fitted"){

		    print("Number of populations does not match the number of fitted allele frequencies")
		    return(NULL)

		  }

				  #Create output variables
				  pop = status = snp.risk = db.g = NULL #Population membership, patient status, genotype relative risk and genotypes

				  ##Percentage of low-frequency variants
				  lfv = 0
				  pa = runif(n.snps*(1-lfv), min = 0.05, max = 0.5) #Ancestral population frequency (used only if allele.freq.model == "Balding-Nichols")
				  pa = c(pa, runif(n.snps*lfv, min = 0.01, max = 0.04))
				  alpha = pa * (1 - fst)/fst #alpha parameter in the Beta distribution
				  beta = (1 - pa) * (1 - fst)/fst #beta parameter in the Beta distribution

				  #Generate population membership, patient status and genotypes for each (sub)population
				  for(i in 1:n.pops){ #For each (sub)population

				    #Simulate the allele frequencies according to the allele frequency model

				    if  (allele.freq.model == "Balding-Nichols"){

				      pl = rbeta(n.snps,alpha,beta) #allele frequencies for this population (generates a frequency for each SNP according with its alpha and betha parameters)

				    } else if (allele.freq.model == "fitted"){

				      #pl = rep(pl.pop[i],n.snps) #allele frequencies for this population
				      pl = runif(n.snps, min = pl.pop[i]-0.05, max = pl.pop[i]+0.05)

				    }

				    #Simulate the genotypes according to the risk model

				    n.pats.pop = round(p.pops[i] * n.pats, 0) #number of patients in this population
				    c.pats.pop = round(c.pops[i] * n.pats.pop, 0) #number of cases in this population
				    o.pats.pop = round((1 - c.pops[i]) * n.pats.pop, 0) #number of controls in this population
				    genotype.pop = NULL #genotypes for this population


				    #For each allele
				    for(j in 1:length(pl)){

				      #Built genotype probabilities matrix
				      gp = matrix(0,nrow = 3, ncol = 2)
				      rownames(gp) <- c("0","1","2"); colnames(gp) <- c("cases","controls")

				      #Generate genotypes for the causal allele conditionated on the disease status
				      gp["0","cases"] <-  (1 - pl[j])^2
				      gp["1","cases"] <-  2 * r.pops[i] * pl[j] * (1 - pl[j])
				      gp["2","cases"] <-  (r.pops[i])^2 * pl[j]^2
				      gp["0","controls"] <-  (1 - pl[j])^2
				      gp["1","controls"] <-  2  * pl[j] * (1 - pl[j])
				      gp["2","controls"] <-   pl[j]^2

				      #Scalate the probabilies in the cases group
				      gp[,"cases"] = gp[,"cases"] / ( (1 - pl[j])^2 +  2 * r.pops[i] * pl[j] * (1 - pl[j]) + (r.pops[i])^2 * pl[j]^2 )


				      #Generate genotypes for each group of patients
				      ##Cases
				      gf = round(gp[,"cases"] * c.pats.pop,0) #genotype frequencies
				      sgf = sum(gf)

				      #Print variables
				      #print(sgf)
				      #print(c.pats.pop)
				      #print(pl[j])

				      if(sgf != c.pats.pop){gf[1]<- (c.pats.pop - sum(gf[-1]))} #Correct difference between counting table and total patiens in this group (due to rounded digits)
				      genotype.l = rep(names(gf),gf)[sample(c.pats.pop,c.pats.pop)]

				      ##Controls
				      gf = round(gp[,"controls"] * o.pats.pop,0) #genotype frequencies
				      sgf = sum(gf)
				      if(sgf != o.pats.pop){gf[1]<- (o.pats.pop - sum(gf[-1]))}
				      genotype.l = c(genotype.l, rep(names(gf),gf)[sample(o.pats.pop,o.pats.pop)])


				      #Conform block of genotyping matrix for this population
				      genotype.pop = cbind(genotype.pop, genotype.l)

				    }

				    #Create vector of population membership and  patient status for this population of patients
				    status.pop = c(rep(1,c.pats.pop), rep(0,o.pats.pop))
				    pop.pop = rep(i,(c.pats.pop + o.pats.pop))
				    snp.risk.pop = rep(r.pops[i],(c.pats.pop + o.pats.pop))

				    #Add population membership, patient status and genotyping block to consolidate vector/table for all populations
				    status = c(status, status.pop)
				    pop = c(pop, pop.pop)
				    snp.risk = c(snp.risk, snp.risk.pop)
				    colnames(genotype.pop) <- paste(prefix.snp,rep(1:n.snps),sep=".")
				    db.g = rbind(db.g, genotype.pop)
				    #print(i)

				  }

		  return(list(db.g, status, pop, snp.risk))

		}



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#' @title Simulate synthetic genotypes
#' @description Function to create sythentic gwas dataset for a given snp type
#'
#[Input]
#' @param snpt SNP type
#' @inheritParams synthetic.gwas
#'
# [Output]:
#' @return
#' * \code{db.g} synthetic GWAS data set
#' * \code{status} disease status of the patient
#' * \code{pop} population ancestry of the patient
#' * \code{snp.risk} genotype relative risk for each SNP
#'
#' @md
#' @author Luis G. Leal, \email{lgl15@@imperial.ac.uk}
#' @family Preprocessing functions
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


		snpt.synthetic.gwas = function(snpt,n.snps,
                              		  n.pats = 2000, #number of patients in the synthetic dataset
                              		  n.pops = 2, #number of populations
                              		  p.pops = c(0.5,0.5), #proportion of patients in each population
                              		  c.pops = c(0.6,0.4), #proportion of cases in each patient population
                              		  r.pops1 = c(1,1), #Risk factor vector in each population for SNPs without association
                              		  r.pops2 = c(3,3), #Risk factor vector in each population for SNPs with association
                              		  pl.pop1 = c(0.1,0.1), #Fitted allele frequencies: It doesnot matter the percentages in this variable as the uniform distribution is used
                              		  pl.pop2 = c(0.15,0.45), #Fitted allele frequencies when SNP with populations structures
                              		  fst = 0.001 )
		{
		  if(snpt == 1){
		    #Generate data set of random SNP without association with the disease
		    ##N. Potential True Negatives
		    res.synt = synthetic.gwas(n.snps = n.snps, #number of SNPs in the synthetic dataset
		                              n.pats = n.pats, #number of patients in the synthetic dataset
		                              n.pops = n.pops, #number of patient populations in the synthetic dataset
		                              p.pops = p.pops, #proportion of patients in each population
		                              c.pops = c.pops, #proportion of cases in each patient population
		                              r.pops = r.pops1, #<-- risk factor vector in each population (i.e., genotype relative risks for the SNPs causative of disease in each population)
		                              fst = fst, #Factor of differentiation between populations
		                              allele.freq.model = c("Balding-Nichols"),
		                              pl.pop = pl.pop1, #<-- vector of fitted allele frequencies in each population
		                              prefix.snp = "snp.tn"
		    )
		  }

		  if(snpt == 2){
		    #Generate data set of population-differentiated SNP without association with the disease
		    ##N. Potential False Positives
		    res.synt = synthetic.gwas(n.snps = n.snps, #number of SNPs in the synthetic dataset
		                              n.pats = n.pats, #number of patients in the synthetic dataset
		                              n.pops = n.pops, #number of patient populations in the synthetic dataset
		                              p.pops = p.pops,#proportion of patients in each population
		                              c.pops = c.pops,#proportion of cases in each patient population
		                              r.pops = r.pops1, #<-- risk factor vector in each population (i.e., genotype relative risks for the SNPs causative of disease in each population)
		                              fst = fst, #Factor of differentiation between populations
		                              allele.freq.model = c("fitted"),
		                              pl.pop = pl.pop2, #<-- vector of fitted allele frequencies in each population
		                              prefix.snp = "snp.fp"
		    )
		  }


		  if(snpt == 3){
		    #Generate data set of random SNP with association with the disease
		    ##P. Potential True Positives
		    res.synt = synthetic.gwas(n.snps = n.snps, #number of SNPs in the synthetic dataset
		                              n.pats = n.pats, #number of patients in the synthetic dataset
		                              n.pops = n.pops, #number of patient populations in the synthetic dataset
		                              p.pops = p.pops,#proportion of patients in each population
		                              c.pops = c.pops,#proportion of cases in each patient population
		                              r.pops = r.pops2, #<-- risk factor vector in each population (i.e., genotype relative risks for the SNPs causative of disease in each population)
		                              fst = fst, #Factor of differentiation between populations
		                              allele.freq.model = c("Balding-Nichols"),
		                              pl.pop = pl.pop1, #<-- vector of fitted allele frequencies in each population
		                              prefix.snp = "snp.tp"
		    )
		  }

		  if(snpt == 4){
		    #Generate data set of population-differentiated SNP with association with the disease
		    ##P. Unknown mix. Mainly True Positives but might be also False Negatives if the population strata mask the association
		    res.synt = synthetic.gwas(n.snps = n.snps, #number of SNPs in the synthetic dataset
		                              n.pats = n.pats, #number of patients in the synthetic dataset
		                              n.pops = n.pops, #number of patient populations in the synthetic dataset
		                              p.pops = p.pops,#proportion of patients in each population
		                              c.pops = c.pops,#proportion of cases in each patient population
		                              r.pops = r.pops2, #<-- risk factor vector in each population (i.e., genotype relative risks for the SNPs causative of disease in each population)
		                              fst = fst, #Factor of differentiation between populations
		                              allele.freq.model = c("fitted"),
		                              pl.pop = pl.pop2, #<-- vector of fitted allele frequencies in each population
		                              prefix.snp = "snp.fn"
		    )
		  }
		  return(res.synt)
		}



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#' @title Simulate synthetic genotypes
#' @description Function to generate synthetic data for a combination of SNP types and number of SNPs
#'
#[Input]
#' @param comb.snpt Combination of SNP types in the synthetic dataset
#' @param comb.snpn Number of SNP types
#'
# [Output]:
#' @return
#' * \code{R} synthetic GWAS data set
#' * \code{out} disease status of the patient
#' * \code{pop} population ancestry of the patient
#' * \code{srisk} genotype relative risk for each SNP
#' * \code{snp.type} SNP type
#' * \code{snp.cat} SNP category (Association with the disease)
#'
#' @md
#' @author Luis G. Leal, \email{lgl15@@imperial.ac.uk}
#' @family Preprocessing functions
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


		#Declare
		mapply.snpt.synthetic.gwas <- function(comb.snpt, #Combination of SNP types in the synthetic dataset
		                                       comb.snpn# Number of SNP types
		                                       )
		{

    		  #Construct the relationship matrix (R)
    		  l.res.synt = mapply(FUN = snpt.synthetic.gwas, snpt = comb.snpt, n.snps = comb.snpn, SIMPLIFY = FALSE)
    		  Ri = lapply(l.res.synt, function(l) l[[1]])
    		  R <- t(do.call("cbind", Ri))
    		  R <- matrix(as.numeric(R), nrow = dim(R)[1], ncol = dim(R)[2], dimnames = list(row.names(R),paste("pat", seq(1:ncol(R)), sep = ".")))

    		  #Extract vectors of patient population and outcome membership
    		  pop <- l.res.synt[[1]][[3]]
    		  out <- l.res.synt[[1]][[2]]
    		  srisk <- l.res.synt[[1]][[4]]

    		  #Extract vectors of snp types
    		  snp.type = as.factor(substr(row.names(R),start = 1, stop = 6))

    		  #Create vector of snp categories (positive and negative association with disease)
    		  snp.cat = as.character(snp.type)
    		  snp.cat[snp.cat%in%c("snp.tp","snp.fn")]<-"p"
    		  snp.cat[snp.cat%in%c("snp.tn","snp.fp")]<-"n"
    		  snp.cat = as.factor(snp.cat)

    		  #Return results
    		  return(list(R,pop,out,snp.type,snp.cat,srisk))

		}
