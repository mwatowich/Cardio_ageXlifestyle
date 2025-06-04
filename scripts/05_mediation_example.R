### Example mediation test

## Author: Marina Watowich
## Date last updated: June 3, 2025


df <- df[which(complete.cases(df$urbanicity_score) & 
                 complete.cases(df$Lymph) & 
                 complete.cases(df$Neutro) & 
                 complete.cases(df$tobacco) & 
                 complete.cases(df$TurkOA_ageResid)),]

out1 <- feols(scale(TurkOA_ageResid) ~ scale(urbanicity_score) + scale(age) + sex + scale(Lymph) + scale(Neutro), data = df)

out2 <- feols(scale(tobacco) ~ scale(urbanicity_score) + scale(age) + sex + scale(Lymph) + scale(Neutro), data = df) 

out3 <- feols(scale(TurkOA_ageResid) ~ scale(urbanicity_score) + scale(tobacco) + scale(age) + sex + scale(Lymph) + scale(Neutro), 
              data = df)

coef(out2)['scale(urbanicity_score)'] * coef(out3)['scale(tobacco)'] / coef(out1)['scale(urbanicity_score)'] * 100

library(boot)
set.seed(1)
mediation_fn <- function(data, i){
  df <- data[i,]
  
  a_path <- feols(scale(tobacco) ~ scale(urbanicity_score) + scale(age) + sex + scale(Lymph) + scale(Neutro), data = df)
  a <- coef(a_path)["scale(urbanicity_score)"]
  
  b_path <-  feols(scale(TurkOA_ageResid) ~ scale(tobacco) + scale(urbanicity_score) + scale(age) + sex + scale(Lymph) + scale(Neutro), data = df)
  b <- coef(b_path)["scale(tobacco)"]
  
  cp <- coef(b_path)["scale(urbanicity_score)"]
  
  ind_ef <- a*b
  total_ef <- a*b + cp
  return(c(ind_ef, total_ef))
}
boot_med <- boot(df, mediation_fn, R = 1000, parallel = "multicore", ncpus = 2)
boot.ci(boot_med, type = c("norm", "perc"))
colMeans(boot_med$t)
