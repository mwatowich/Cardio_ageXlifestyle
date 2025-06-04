### Generate joint-population DNAm clock using PC-Clock package 
# Code adopted from the PC-Clock package: Albert T., H., et al. (2022). A computational solution for bolstering reliability of epigenetic clocks: Implications for clinical trials and longitudinal tracking. Nature Aging
# Please see the PC-Clock package GitHub for more details and further instructions 

## Author: Marina Watowich
## Date last updated: June 3, 2025


library(glmnet)
library(dplyr)

# percMeth is a matrix of proportion methylation for Turkana samples, oa_meth is the same for Orang Asli data
# pheno is a data frame of age, sex, and estimated cell type frequency information per sample


set.seed(123)

# Subset methylation data to your set of CpGs
CpGs <- colnames(percMeth)

# Check order of data is correct
if(all(colnames(percMeth) == CpGs)){
  message("CpGs are all in order")
} else(message(paste("Only",sum(colnames(percMeth) == CpGs),"CpGs are in order")))

if(all(rownames(percMeth) == pheno$SampleID)){
  message("Samples are all in order")
} else(message(paste("Only",sum(rownames(percMeth) == pheno$SampleID),"Samples are in order")))


# Perform PCA and projections. Remove last PC. Columns should be CpGs, rows should be samples
PCA = prcomp(percMeth,scale.=F)
TrainPCData = PCA$x[,1:(dim(PCA$x)[2]-1)]
TrainAge <- pheno$Age

# Train PC clock
cv = cv.glmnet(TrainPCData, TrainAge, nfolds=10,alpha=0.5, family="gaussian")
fit = glmnet(TrainPCData, TrainAge, family="gaussian", alpha=0.5, nlambda=100)
plot(cv)

# Compress model
CalcPCAge_comb <- vector(mode = "list",length = 0)
temp = as.matrix(coef(cv,s = cv$lambda.min))
CalcPCAge_comb$model = temp[temp!=0,][-1]
CalcPCAge_comb$intercept = temp[1,1]
CalcPCAge_comb$center = PCA$center
CalcPCAge_comb$rotation = PCA$rotation[,names(CalcPCAge_comb$model)]

# Save model
save(CalcPCAge_comb,CpGs,file = paste0("my_file_", Sys.Date(), ".RData"))


# Calculate new clock
load(file = "my_file_date")
percMeth_2 = percMeth[,CpGs]
PCAge <- sweep(as.matrix(percMeth_2),2,CalcPCAge_comb$center) %*% CalcPCAge_comb$rotation %*% CalcPCAge_comb$model + CalcPCAge_comb$intercept
PCAge <- data.frame(SampleID = rownames(PCAge), TurkOA_Age = PCAge)

pheno_w_predictAge <- dplyr::left_join(pheno, PCAge, by = "SampleID")
cor.test(pheno_w_predictAge$Age, pheno_w_predictAge$TurkOA_Age)

# Calculate age residuals
calcPCClocks_Accel <- function(DNAmAge){
  clockColumns = c("TurkOA_Age")
  for (i in clockColumns){
    DNAmAge[,paste0(i,"Resid")] = resid(lm(DNAmAge[,i] ~ DNAmAge$Age))
  }
  return(DNAmAge)
}
pheno_w_predictAge <- calcPCClocks_Accel(pheno_w_predictAge)

