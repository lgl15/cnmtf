
###############################################################
#	Full executable script for cNMTF
#
# Corresponding author:
# Luis Leal, Imperial College London, email: lgl15@imperial.ac.uk
###############################################################


#-------------------------------------------------
# 1. Installation
#-------------------------------------------------


			#Clean your current workspace

			rm(list=ls()) ; gc()

			#Load libraries

			library("devtools") #Development tools for packages in R
			library("biomaRt") #Query ENSEMBL for SNV functional consequences
			library("LDcorSV") #Linkage disequilibrium
			library("doParallel") #Processes in parallel
			library("igraph") #Network construction
			library("ade4") #Principal components analysis
			library("VennDiagram") #Venn diagrams


			#Installation of cnmtf package from github repository

			devtools::install_github("lgl15/cnmtf")


			#Create directory for your analysis

			dir.create("./test/")


#-------------------------------------------------
# 2. Preprocessing of data
#-------------------------------------------------


			#Find variants in high Linkage-disequilibrium

			find.snps.ld( file.LD = "./test/fileLD.RData",	type.ld = "gene", tmap = tmap, R = R)


			#Create SNV-SNV network

			create.network( #Parameters for the reference network
											net.type = "ppi", #Type of reference network
											dedges = dedges, #Define object with edges from reference network

											#Parameters for Linkage Disequilibrium
											remove.highLD = TRUE,
											ld.tao = 0.8, #Treshold of LD
											res.ld = "./test/fileLD.RData", #Table of LD
											keep.with.LD = NULL, #List of SNPs to include even if they are in LD

											#Parameters to construct Wu
											R.snps = rownames(R), #List of SNPs in R
											work.dat = ".", #Working directory
											trait.project = "test", #Trait
											n.cores = 3, #Number of cores
											tmap = tmap, #Mapping of SNPs to genes
											plot.file = "Gu_ppi_test_venn.pdf" ) #File to print Venn diagrams and node degree distribution



#-------------------------------------------------
# 3. Running cNMTF
#-------------------------------------------------


			score.cnmtf( R = R, #Relationship matrix
									 out = out, #Categorical outcome variable
									 pop = pop, #Population variable
									 log.file = "logfile_my_experiment", #Log-file to track the performance

									 #Variables to Save/Load data and workspaces
									 name.exp = "my_experiment", #Name of experiment to save files
									 work.dat = "./test/", #Folder to save and load workspaces

									 #Number of clusters
									 define.k = "method", #Method to define k1
									 k1 = c(2,3,5,8,10,15,20), #Number of clusters of SNVs
									 k2 = nlevels(out), #Number of clusters of patients

									 #Penalisation parameters
									 wparameters = list(gamma1 = 0.5, #Weight for the SNV network, Lu
									 									 gamma2 = 0.999, #Weight for the outcome matrix, Vo
									 									 gamma3 = 0.999), #Weight for the kernel matrix, A
									 save.parameters = TRUE, #Save parameters to file
									 run.t.par = 4, #Number of repetitions for parameters fitting
									 max.try0 = 4,  #Maximum number of tries to fit the parameter
									 snps.known = NULL, #List of known SNV associated with the trait

									 #Variables to control performance of the algorithm
									 parallel.opt = F, #Run some instances of the algorithm in parallel
									 n.cores = 3, #Number of cores to use in the parallel processing
									 init = 0, #Type of seeding/initialisation of matrices in the algorithm
									 do.U = TRUE, #Perform clustering of SNPs
									 calcObj = 20, #Check convergency every 20 iterations
									 calcObj2 = 40, #Start checking convergency after first 40 iterations
									 iters = 300, #Number of iterations
									 run.t.exp = 10, #Number of repetitions for the experiment
									 display.iters = FALSE, #Display the iterations of function cnmtf

									 #Randomisations
									 score.pvalues = T, #Estimate p-values for the scores
									 random.parallel = F, #Run each randomisation in parallel
									 randomisations = 100, #Number of randomisations

									 #Construction of penalisation terms
									 file.Gu = "./test/Gu_ppi_test.RData" #Workspace with SNV-SNV network

			)



#-------------------------------------------------
# 4. Prioritising SNVs
#-------------------------------------------------


		#Calculate the delta score and create Manhattan plots

			analyze.cnmtf (trait.project = "test", #Trait
										 name.exp = "my_experiment", #Name of experiment to analyse
										 work.dat = "./test/", #Working directory
										 alpha.cnmtf = 0.005, #Significance level for the delta SNV score
										 d.conf = NULL, #Optional. A dataframe of patients by confounder variables
										 snps.known = NULL, #Optional. List of known SNVs associated with the disease
										 snps.known2 = NULL #Optional. A second list of SNVs to depict on the Manhattan plots
			)


		#Map variants to genes and add functional annotations

			t.res = annotate.results( name.exp = "my_experiment", #Define experiment id
																snps.known = NULL, #Optional. List of SNVs known to be associated with the disease
																snps.known2 = NULL, #Optional. A second list of SNVs.
																add.david.annotations = TRUE, #Use DAVID web service or export/import manually
																email.david = email.david, #Email account registered in DAVID.
																add.ensemble.conseq = FALSE, #Add SNV consequences from ENSEMBL
																work.dat = "./test/", #Working directory
																tmap = tmap, #Mapping of SNPs to genes, chr and genomic position
																file.LD = "./test/fileLD.RData", #Workspace with a table of pairwise LD
																ld.tao = 0.8 #Treshold of LD
			)

		#Extract table of prioritised variants
		t.snvs = t.res[[1]]
		t.snvs [1:8]




