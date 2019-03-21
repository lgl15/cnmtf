###############################################################
# cNMTF
#	5. Examples
#
# Corresponding author:
# Luis Leal, Imperial College London, email: lgl15@imperial.ac.uk
###############################################################



#-------------------------------------------------
# 1. Create the sample data
#-------------------------------------------------


		rm(list=ls()) ; gc()

		#Load real data
		load("~/Desktop/eMERGE/data/r_workspaces/ldl/preprocess_ldl.RData")

		#Take a sample of SNVs and white American patients
		set.seed(22)
		pat.sample = sample(x = which(race == "White"), 50)
		snv.sample = sample(1:nrow(R), 100)
		R = R[ snv.sample , pat.sample ]

		#Anonymise SNVs and patients
		colnames(R) <- paste("patient.", 1:ncol(R), sep = "")
		rownames(R) <- sample(rownames(R), nrow(R))

		#Create tmap
		tmap = tmap[ match(rownames(R),tmap$refsnp_id), ]

		#Create outcome and population
		out = out.cat[ pat.sample]
		table(out)
		pop = pop[ pat.sample ]
		table(pop)

		#Save the data to the package
		devtools::use_data(out, pkg = ".")
		devtools::use_data(pop, pkg = ".")
		devtools::use_data(R, pkg = ".")
		devtools::use_data(tmap, pkg = ".")

		#Create small PPI edge list
		file.dedges = "../../lgl15/Desktop/eMERGE/data/ref_networks/ppi_edges_ldl.txt" #Define filename with edges from reference network
		dedges = read.table(file.dedges, head = TRUE)
		dedges = dedges[ (dedges$V2 %in% tmap$entrezgene) & (dedges$V3 %in% tmap$entrezgene), ]
		dedges = dedges[ !duplicated(paste(dedges$V2, dedges$V3, sep = "-")), ]
		dim(dedges)
		head(dedges)
		rownames(dedges) <- NULL

		#Save the data to the package
		devtools::use_data(dedges, pkg = ".",overwrite = TRUE )







#-------------------------------------------------
# 2. Examples of preprocessing
#-------------------------------------------------

		library("LDcorSV")
		cnmtf:::find.snps.ld( file.LD = "./test/fileLD.RData",	type.ld = "gene", tmap = tmap, R = R)


		library(doParallel)
		library('biomaRt')
		library("igraph")

		cnmtf:::create.network( #Parameters for the reference network
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
# 3. Examples of main clustering
#-------------------------------------------------

		rm(list=ls()) ; gc()

		#Load sampling data
		setwd("~/cnmtf/")

		load("data/R.rda")
		load("data/out.rda")
		load("data/pop.rda")
		load("data/tmap.rda")

		#Load libraries
		library('doParallel')
		library("igraph")

		#Load functions
		source("R/fact_clust.R")
		source("R/fact_cons.R")
		source("R/fact_init.R")
		source("R/fact_param.R")
		source("R/fact_main.R")
		source("R/prep_pop.R")
		source("R/prep_out.R")


		score.cnmtf( R = R, #Relationship matrix
								 out = out, #Categorical outcome variable
								 pop = pop, #Population variable
								 log.file = "logfile_my_experiment", #Log-file to track the performance

								 #Variables to Save/Load data and workspaces
								 name.exp = "my_experiment", #Name of experiment to save files
								 work.dat = "./test/", #Folder to save and load workspaces

								 #Number of clusters
								 define.k = "user", #Method to define k1
								 k1 = c(2,3,5,8,10,15,20)[5], #Number of clusters of SNVs
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
								 randomisations = 2, #Number of randomisations

								 #Construction of penalisation terms
								 file.Gu = "./test/Gu_ppi_test.RData" #Workspace with SNV-SNV network

		)


#-------------------------------------------------
# 4. Examples of prioritising SNVs
#-------------------------------------------------


		rm(list=ls()) ; gc()

		#Load sampling data
		setwd("~/cnmtf/")

		load("data/R.rda")
		load("data/out.rda")
		load("data/pop.rda")
		load("data/tmap.rda")

		#Load libraries
		library("igraph")
		library("ade4")
		library("ComplexHeatmap")
		library("VennDiagram")

		#Load functions
		source("R/fact_clust.R")
		source("R/fact_cons.R")
		source("R/fact_init.R")
		source("R/fact_param.R")
		source("R/fact_main.R")
		source("R/prep_pop.R")
		source("R/prep_out.R")
		source("R/comp_lrm.R")
		source("R/dscore_main.R")
		source("R/dscore_plot.R")
		source("R/dscore_aux.R")

		work.dat = "./test/"
		trait.project = "test"
		name.exp = "my_experiment"
		d.conf = NULL
		snps.known = NULL
		snps.known2 = NULL
		alpha.cnmtf = 0.005


		analyze.cnmtf (trait.project = "test", #Trait
									 name.exp = "my_experiment", #Name of experiment to analyse
									 work.dat = "./test/", #Working directory
									 alpha.cnmtf = 0.005, #Significance level for the delta SNV score
									 d.conf = NULL, #Optional. A dataframe of patients by confounder variables
									 snps.known = NULL, #Optional. List of known SNVs associated with the disease
									 snps.known2 = NULL #Optional. A second list of SNVs to depict on the Manhattan plots
)


		rm(list=ls()) ; gc()

		#Load sampling data
		setwd("~/cnmtf/")
		add.david.annotations = TRUE
		add.ensemble.conseq = FALSE
		email.david = "lgl15@imperial.ac.uk"
		work.dat = "./test/"
		trait.project = "test"
		name.exp = "my_experiment"
		snps.known = NULL
		snps.known2 = NULL
		file.LD = "./test/fileLD.RData"
		ld.tao = 0.8 #Treshold of LD

		load("data/R.rda")
		load("data/out.rda")
		load("data/pop.rda")
		load("data/tmap.rda")

		library("RDAVIDWebService")
		library('biomaRt')

		source("R/comp_res.R")

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
		head(t.snvs [1:8])

		#Extract table of mutated genes
		t.genes = t.results[[2]]
		head(t.genes)





