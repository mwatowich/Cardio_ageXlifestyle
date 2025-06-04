### Model age x lifestyle effects for cardiometabolic biomarkers of interest
## Determine whether linear or quadratic age curve better fits the data - to look at age * lifestyle effects for the better age effect 

## Author: Marina Watowich
## Date last updated: June 3, 2025


### Generate 10 random subsample data downloaded from NHANES to match the the age distribution in our dataset
range(nhanes$age) # max is 80yo
full_data <- turkana[which(complete.cases(turkana$age) & complete.cases(turkana$sex) & turkana$age<=80),]

# Randomly subsample 10 times 
set.seed(2)
datasets_list <- vector("list", 10)
for (i in 1:10) {
  datasets_list[[i]] <- bind_rows(
    slice_sample(subset(nhanes, age>=18 & age<=20), n = nrow(subset(full_data, age>=18 & age<=20)), replace = F),
    slice_sample(subset(nhanes, age>20 & age<=25), n = nrow(subset(full_data, age>20 & age<=25)), replace = F),
    slice_sample(subset(nhanes, age>25 & age<=30), n = nrow(subset(full_data, age>25 & age<=30)), replace = F),
    slice_sample(subset(nhanes, age>30 & age<=35), n = nrow(subset(full_data, age>30 & age<=35)), replace = F),
    slice_sample(subset(nhanes, age>35 & age<=40), n = nrow(subset(full_data, age>35 & age<=40)), replace = F),
    slice_sample(subset(nhanes, age>40 & age<=45), n = nrow(subset(full_data, age>40 & age<=45)), replace = F),
    slice_sample(subset(nhanes, age>45 & age<=50), n = nrow(subset(full_data, age>45 & age<=50)), replace = F),
    slice_sample(subset(nhanes, age>50 & age<=55), n = nrow(subset(full_data, age>50 & age<=55)), replace = F),
    slice_sample(subset(nhanes, age>55 & age<=60), n = nrow(subset(full_data, age>55 & age<=60)), replace = F),
    slice_sample(subset(nhanes, age>60 & age<=65), n = nrow(subset(full_data, age>60 & age<=65)), replace = F),
    slice_sample(subset(nhanes, age>65 & age<=70), n = nrow(subset(full_data, age>65 & age<=70)), replace = F),
    slice_sample(subset(nhanes, age>70 & age<=75), n = nrow(subset(full_data, age>70 & age<=75)), replace = F),
    slice_sample(subset(nhanes, age>75 & age<=80), n = nrow(subset(full_data, age>75 & age<=80)), replace = F))
  
  write.table(datasets_list[[i]], 
              file = paste0("nhanes_subsample_", Sys.Date(), "_", i, ".txt"), 
              sep = "\t", row.names = F, col.names = T, quote = F)
}



### Model Turkana vs. US data
# Set lifestyle as "rural Turkana", "urban Turkana", or "US" (NHANES)
# "combined" = a data frame of Turkana (or OA) data and NHANES

# Model linear additive and linear interactive models for biomarkers best predicted by linear age terms, do for each NHANES permutation. Adapt for quadratic age models. 
linear_turkanaUS_list <- vector("list", length(combined))
for (i in 1:length(combined)) {
  combined_focal_iteration <- combined[[i]]
  # Run models
  linear_turkanaUS <- as.data.frame(do.call(rbind, lapply(linear_outcomes, function(x){
    do.call(rbind, lapply(c("rural","urban"), function(y){
      combined_focal_iteration <- combined_focal_iteration[which(complete.cases(combined_focal_iteration$age) & 
                                             complete.cases(combined_focal_iteration$sex) & 
                                             complete.cases(combined_focal_iteration$lifestyle)),]
      combined_focal_iteration$lifestyle <- relevel(as.factor(combined_focal_iteration$lifestyle), ref = y)
      
      # Run linear additive model
      out_lin <- lm(scale(combined_focal_iteration[[x]]) ~ scale(age) + lifestyle + sex, combined_focal_iteration)
      out_lin <- cbind(as.data.frame(summary(out_lin)$coefficients[-1,-3]),
                       AIC(out_lin))
      colnames(out_lin) <- c("beta", "se", "pval","AIC")
      out_lin$n <- nrow(combined_focal_iteration[complete.cases(combined_focal_iteration[[x]]),])
      out_lin <- out_lin %>% 
        tibble::rownames_to_column(var = "covariate") %>% 
        mutate(biometric = x) %>% 
        tidyr::pivot_wider(names_from = covariate, values_from = c(beta, se, pval))
      colnames(out_lin) <- gsub("scale(age)","age",colnames(out_lin), fixed = T)
      colnames(out_lin) <- gsub("rural","",gsub("urban","",colnames(out_lin)))
      out_lin$lifestyle_reference <- y
      out_lin$type <- "linear additive"
      
      # Run linear interactive model
      out_lin_inter <- lm(scale(combined_focal_iteration[[x]]) ~ scale(age)*lifestyle + sex, combined_focal_iteration)
      out_lin_inter <- cbind(as.data.frame(summary(out_lin_inter)$coefficients[-1,-3]),
                             AIC(out_lin_inter))
      colnames(out_lin_inter) <- c("beta", "se", "pval","AIC")
      out_lin_inter$n <- nrow(combined_focal_iteration[complete.cases(combined_focal_iteration[[x]]),])
      out_lin_inter <- out_lin_inter %>% 
        tibble::rownames_to_column(var = "covariate") %>% 
        mutate(biometric = x) %>% 
        tidyr::pivot_wider(names_from = covariate, values_from = c(beta, se, pval))
      colnames(out_lin_inter) <- gsub("scale(age)","age",colnames(out_lin_inter), fixed = T)
      colnames(out_lin_inter) <- gsub(":","_",colnames(out_lin_inter))
      colnames(out_lin_inter) <- gsub("rural","",gsub("urban","",colnames(out_lin_inter)))
      out_lin_inter$lifestyle_reference <- y
      out_lin_inter$type <- "linear interactive"
      
      # Combine 
      out <- dplyr::bind_rows(out_lin, out_lin_inter)
      out <- out %>% 
        dplyr::select(c("biometric","AIC","type",everything()))
      out$linear_additive_AIC <- out[out$type == "linear additive",]$AIC
      return(out)
    }))
  })))
  # Store the results
  linear_turkanaUS_list[[i]] <- linear_turkanaUS 
}

# Can extract age, lifestyle, and age*lifestyle for each comparison (rural v US, urban v US) 
