###############################################################
# cNMTF
#	4. Comparing functions
# 4.2 Functions to prioritise SNVs using regression models
#
# Corresponding author:
# Luis Leal, Imperial College London, email: lgl15@imperial.ac.uk
###############################################################


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#' @title Prioritise variants using regression models
#' @description Function to fit logistic or linear regression models
#'
#[Input]:
#' @param out clinical output
#' @param R relationship matrix
#' @param logistic.model Use logistic regression or linear regression models
#' @param coding Codding scheme for the R matrix
#' @param d.conf dataframe with confounder variables to correct for
#'
#[Output]:
#' @return \code{model.pvals} p-values for the association SNV-trait
#'
#' @md
#' @author Luis G. Leal, \email{lgl15@@imperial.ac.uk}
#' @family Factorisation functions
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


    regression.snps = function(out, #Outcome
                               R, #Relationship matrix
                               logistic.model = TRUE, #Use logistic regression
                               coding = "additive", #Codding scheme for the R matrix
                               d.conf = NULL # Dataframe with confounder variables to correct for (It might include the principal components to correct for population structures)
    ){

      #Filter out patients without outcome data
      pos.out = !is.na(out) #Because the rows of dcli and the columns of R match
      R.s = R[,pos.out]
      out = out[pos.out]


      #Save dimensions
      m = ncol(R.s)
      n = nrow(R.s)
      d.conf = d.conf[pos.out, ] #Remove patients witout trait data


      #Change the conding scheme to: Dominant effect coding
      if( coding == "dominant" ){
        R.s = coding.scheme(t(R.s),coding = c("dominant"))
      }


      #Define kind of regression model
      if( logistic.model ==  TRUE){ #Fit logit model
        family = binomial(link='logit')
      }else{
        family = gaussian
      }





      #Declare vector to save the p-values of the regression models
      model.pvals = rep(NA,n)

      #Fit models for each SNP
      for(i in 1:n){

        #If MAF = 0 then exclude that SNP from analysis
        if(sum(R.s[i,]) == 0 | sum(R.s[i,]) == m*2){
          model.pvals[i] = 1
          next
        }

        #Assing p-value = 1 to those alleles without variation
        if(sum(R.s[i,]) == 0 | sum(R.s[i,]) == m*2){#If MAF = 0
          model.pvals[i] = 1
          next
        }


        #If set of confounding variables is provided then use then directly in the model
        if( ! is.null(d.conf) ){

          covariates = data.matrix( cbind( as.numeric(R.s[i,]), d.conf) )


        }else{ #If not, then use only the genotyping data for that SNP

          covariates = data.matrix(as.numeric(R.s[i,]))

        }

        #Fit model
        model <- glm(out ~ 1 + covariates, family = family)

        #Extract results
        model.pvals[i] = coef(summary(model))[2,4]

      }#End loop across SNPs


      return(model.pvals)
    }
