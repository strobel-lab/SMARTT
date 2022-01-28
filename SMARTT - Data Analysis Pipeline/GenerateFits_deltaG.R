
# This script reads a csv file with %FL values for various ligand concentrations
# and fits the data to a modified binding equation for all variants. This is an
# alternative to analyzing the data using PRISM. This file also includes delta G
# values in the output.
#
#
# - Requires R (has been tested with version 3.4.4) and the minpack.lm package
# - Public release 1.0
# - Copyright 2018 Chad Torgerson

###################################################################################
# GPL statement:
#
# This file is part of SMART-DAP (SMARTT - Data Analysis Pipeline).
# 
# SMARTT-DAP is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# SMARTT-DAP is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with Foobar.  If not, see <http://www.gnu.org/licenses/>.
###################################################################################

library(minpack.lm)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~~~~~Variables to edit~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Name input file
inputfolder <- '/Users/Kathryn/Desktop/Toc/'
inputfileName <- 'Toc_singles_test.csv'
inFile <- paste(inputfolder, inputfileName, sep='')

# Name output file
outputfolder <- './'
outFileName <- 'Toc_singles_test_Fits.csv'
outFile <- paste(outputfolder, outFileName, sep='')


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Read input data
input_data <- read.csv(inFile)

R_val <- 0.001987216 #gas constant
Temp <- 310 #temp K

# Initialize variables
mat <- matrix(data = NA, nrow = I(NROW(input_data)-1), ncol = 14)
header <- c('K', 'K_SD', 'Ymin', 'Ymin_SD', 'Ymax', 'Ymax_SD', 'Amp', 'Amp_SD', 'dG_Amp', 'dG_Amp_SD', 
            'dG_K', 'dG_K_SD', 'dG_switch', 'dG_switch_SD')
colnames(mat) <- header
rownames(mat) <- input_data[-(1),1]
count <- 1

# Molar concentration values
X <- unname(unlist(input_data[1, seq(3, NCOL(input_data), 2)]))

# Constraints for K
K_lower <- min(X[X!=0])/100
K_upper <- max(X)*100

# Iteratively fit all variants in input file
for (row in 2:NROW(input_data)){
  count <- count + 1
  NormApprox = TRUE
  
  # Check that at least one read was found
  Reads <- unname(unlist(input_data[row, seq(2, NCOL(input_data), 2)]))   # number of reads
  if(max(Reads) == 0){
    next
  }
  
  #Set parameters
  Y <- unname(unlist(input_data[row, seq(3, NCOL(input_data), 2)]))  # % Full Length
  Weights <- rep(NA, length(Reads))   # inverse variance
  
  # Check whether the Normal Approximation or the Beta Binomial Distribution should be used
  if(min(min(Reads*Y/100, na.rm = TRUE),min(Reads*(1-(Y/100)), na.rm = TRUE)) <= 5){
    NormApprox = FALSE
  }
  
  for (i in 1:length(Reads)){
    n <- Reads[i]
    p <- (Y[i]/100)
    q <- (1-p)
    s <- round(n*p) # number of successes
    
    # Might need to deal with n = 0 cases somehow here
    if (n == 0){
      Weights[i] <- NA
      
      # Check if n*p>5 and n*q>5 at all positions
    } else if (NormApprox == TRUE){
      Weights[i] <- (1/(100*(p*q)/n))
      
      # If not, use Beta Binomial Distribution
    } else {
      # Set the %FL to the mean of the posterior distribution function
      Y[i] <- 100*(s+.5)/(n+1)
      
      # Determine the reciprocal of the variance
      Weights[i] <- 1/(100*((s+.5)*(n-s+.5))/(((n+1)^2)*(n+2)))
    }
    
  }
  
  # Set start values for nonlinear regression
  Ymin_est <- min(Y, na.rm = TRUE)
  Ymax_est <- max(Y, na.rm = TRUE)
  Amp_est <- (Ymax_est - Ymin_est)
  Ymid_est <- (Ymin_est + Amp_est/2)
  index = which.min(abs(Y-Ymid_est))
  
  # Don't allow the lowest or highest concentrations from being used as the 
  # start value for K
  if (index == 1){
    index <- 3
  } else if (index == length(X)){
    index <- index-2
  } 
  K_est = X[index]
  
  # Fit the data to a modified binding equation
  df = data.frame(X,Y)
  fit <- tryCatch({
    nlsLM(Y ~ I((Ymax-Ymin)*X/(K+X)+Ymin),
          data = df,
          start = list(K = K_est, Ymin = Ymin_est, Ymax = Ymax_est),
          control = nls.lm.control(maxiter = 500),
          weights = Weights,
          na.action = na.omit,
          lower = c(K_lower, 0, 0),
          upper = c(K_upper, 100, 100))
    
  }, warning = function(w) {
    print(paste("Line", count, "- MY_WARNING:  ", w))
    return(w)
  }, error = function(e) {
    print(paste("Line", count, "- MY_ERROR:  ", e))
    return(e)
  })
  
  if (inherits(fit, "warning") || inherits(fit, "error")){
    next
  }
  
  # Add fit values to output matrix
  fit_summary <- summary(fit)
  fit_params <- fit_summary$parameters
  
  K_vector <- fit_summary$parameters["K",c(1,2)]
  Ymin_vector <- fit_summary$parameters["Ymin",c(1,2)]
  Ymax_vector <- fit_summary$parameters["Ymax",c(1,2)]
  Amp_val <- (Ymax_vector[1] - Ymin_vector[1])
  Amp_SE <- sqrt(Ymin_vector[2]^2 + Ymax_vector[2]^2)
  
  #Change fit parameters into delta_G for each sequence
  dG_Amp_val <- Amp_val/(100-Amp_val)
  dG_Amp_SE <- abs(dG_Amp_val)*sqrt((Amp_SE^2/Amp_val^2)+(Amp_SE^4/(100-Amp_val)^2))
  dG_K_val <- R_val*Temp*log(K_vector[1]*10^-6)
  dG_K_SE <- R_val*Temp*(K_vector[2]/K_vector[1])
  dG_switch_val <- -1*R_val*Temp*log(dG_Amp_val)
  dG_switch_SE <- R_val*Temp*(dG_Amp_SE/dG_Amp_val)
  
  param_vector <- t(c(K_vector,Ymin_vector,Ymax_vector,Amp_val,Amp_SE,dG_Amp_val,dG_Amp_SE,dG_K_val,dG_K_SE
                      ,dG_switch_val,dG_switch_SE))
  mat[(row-1),] <- param_vector
}

#build a new matrix that compares delta_G values to the WT -- This is a delta_delta_G matrix

mat2 <- matrix(data = NA, nrow = I(NROW(mat)), ncol = 20)
header2 <- c('K', 'K_SD', 'Ymin', 'Ymin_SD', 'Ymax', 'Ymax_SD', 'Amp', 'Amp_SD', 'dG_Amp', 'dG_Amp_SD', 
            'dG_K', 'dG_K_SD', 'dG_switch', 'dG_switch_SD',  'ddG_K', 'ddG_K_SD','ddG_switch', 'ddG_switch_SD'
            ,'Function', 'Function_SD')
colnames(mat2) <- header2
rownames(mat2) <- rownames(mat)

for(x in 1:nrow(mat2)){
  mat2[x,1:14]=mat[x,1:14]
  mat2[x,15]=mat[x,'dG_K']-mat['WT','dG_K']
  mat2[x,16]=mat[x,'dG_K_SD']^2+mat['WT','dG_K_SD']^2
  mat2[x,17]=mat[x,'dG_switch']-mat['WT','dG_switch']
  mat2[x,18]=mat[x,'dG_switch_SD']^2+mat['WT','dG_switch_SD']^2
  mat2[x,19]=mat2[x,'ddG_K']+mat2[x,'ddG_switch']
  mat2[x,20]=mat2[x,'ddG_K_SD']^2+mat2[x,'ddG_switch_SD']^2
}


write.csv(mat2, file = outFile)
print('Complete!')
