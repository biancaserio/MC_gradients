
# ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------  

# Project: 28 Me 

# Supplementary Analyses

# Content: Estimate first- and second-order vector autoregressive (VAR) models of S-A axis loadings vs estradiol using parallel processing

# References: 
# https://github.com/tsantander/PritschetSantander2020_NI_Hormones/blob/master/analysis/jacobs28andMe_runEdgewiseVAR.R -> currently copied here
# https://github.com/tsantander/PritschetSantander2020_NI_Hormones/blob/master/analysis/jacobs28andMe_runGraphMetricVAR.R 

# ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------  




# Load in required packages.
#--------------------------------------------------------------------------

require(vars)
require(doMC)



# Specify number of cores for parallel estimation.
#--------------------------------------------------------------------------

registerDoMC(cores=24)



# Navigate to where the data live
#--------------------------------------------------------------------------


# Clear environment
rm(list = ls())

# set up directories
codedir = '/data/p_02667/28Me/scripts/'
datadir = '/data/p_02667/28Me/data/'
resdir = '/data/p_02667/28Me/results/'

# set directory to path of current script
setwd(codedir) 




# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# PREPARE DATA 
# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

# Aligned functional gradient loadings 
array_aligned_fc_G1_excl26 = read.csv(paste(resdir, 'array_aligned_fc_G1_excl26.csv', sep = ''), fileEncoding = 'UTF-8-BOM', header = FALSE)


# class(array_aligned_G1)
# typeof(array_aligned_G1)
dim(array_aligned_fc_G1_excl26)



# Descriptives 
hormones_28Me_excl26 = read.csv(paste(resdir, 'hormones_28Me_excl26.csv', sep = ''), fileEncoding = 'UTF-8-BOM')

str(hormones_28Me_excl26)



source('./permTestEdgeVAR.R')




#### original code 
# 
# zCoh <- read.csv('./edgewiseCoherence.csv', header = FALSE)
# dat  <- read.csv('./data.efficiency.csv', header = TRUE)
# 
# zCoh           <- data.frame(scale(zCoh))
# zDat           <- data.frame(scale(dat[,2:11]))
# colnames(zDat) <- colnames(dat)[2:11]
# 
# nEdge <- dim(zCoh)[2]
# 


### replacing with my variables

zCoh           <- data.frame(array_aligned_fc_G1_excl26)  # dim 29 400
zDat           <- data.frame(hormones_28Me_excl26)
#colnames(zDat) <- colnames(hormones_28Me_excl26)[2:11]

nEdge <- dim(zCoh)[2]  # i.e., 400 parcels




# Loop over edges in parallel, fit models.
#--------------------------------------------------------------------------

nptOut <- foreach(iEdge=1:nEdge, .combine=rbind, .packages='vars') %dopar% {
  
  modelDat           <- data.frame(zCoh[,iEdge], zDat$z_estradiol)
  colnames(modelDat) <- c('Edge', 'Estro')
  
  edgeVAR <- VAR(modelDat, p = 2)
  
  permResult <- permTestEdgeVAR(modelDat, edgeVAR, sample(1:5e5, 1))
  paramOut   <- permResult$Coefficients
  
  edgeCoef  <- paramOut$Estimate[paramOut$Model == 'Edge']
  estroCoef <- paramOut$Estimate[paramOut$Model == 'Estro']
  edgeNPT   <- paramOut$pValue[paramOut$Model == 'Edge']
  estroNPT  <- paramOut$pValue[paramOut$Model == 'Estro']
  
  npt <- c(edgeCoef, edgeNPT, estroCoef, estroNPT)
  
}

dim(nptOut)



# Save output.
#--------------------------------------------------------------------------

#save(nptOut, file = 'estroEdgeVAR.rda', compress = 'xz')
write.csv(nptOut, paste(resdir, 'estroEdgeVAR.csv', sep = ''), row.names = FALSE)
