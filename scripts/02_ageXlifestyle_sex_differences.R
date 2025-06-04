### Sex differences in age*lifestyle effects 
## Model age*lifestyle effects in females and males separately, characterize whether interactive effects are shared or different between females and males

## Author: Marina Watowich
## Date last updated: June 3, 2025


# linear_outcomes = vector of biomarkers for which linear age term (ie. linear additive or linear interactive models) best fit the data 
# quadratic_outcomes = vector of biomarkers for which quadratic age term best fit the data


linear_interactive_models <- as.data.frame(do.call(rbind, lapply(linear_outcomes, function(x){
  do.call(rbind, lapply(c("F","M"), function(y){
    df <- df[which(complete.cases(df$age) & 
                          df$sex==y & 
                          complete.cases(df$urbanicity_score)),]
    
    # Run linear interactive model
    out <- lm(scale(df[[x]]) ~ scale(age)*scale(urbanicity_score), df)
    out <- cbind(as.data.frame(summary(out)$coefficients[-1,-3]),
                           AIC(out))
    colnames(out) <- c("beta", "se", "pval","AIC")
    out$n <- nrow(df[complete.cases(df[[x]]),])
    out <- out %>% 
      tibble::rownames_to_column(var = "covariate") %>% 
      mutate(biometric = x) %>% 
      tidyr::pivot_wider(names_from = covariate, values_from = c(beta, se, pval))
    colnames(out) <- gsub("scale(age)","age",colnames(out), fixed = T)
    colnames(out) <- gsub("scale(urbanicity_score)","urbanicity_score",colnames(out), fixed = T)
    colnames(out) <- gsub(":","_",colnames(out))
    out$sex <- y
    out$type = "linear interactive"
    return(out)
  }))
})))

quadratic_interactive_models <- as.data.frame(do.call(rbind, lapply(quadratic_outcomes, function(x){
  do.call(rbind, lapply(c("F","M"), function(y){
    df <- df[which(complete.cases(df$age) & 
                          df$sex==y & 
                          complete.cases(df$urbanicity_score)),]
    
    # Run quadratic interactive model
    out <- lm(scale(df[[x]]) ~ scale(poly(age,2))*scale(urbanicity_score), df)
    out <- cbind(as.data.frame(summary(out)$coefficients[-1,-3]),
                            AIC(out))
    colnames(out) <- c("beta", "se", "pval","AIC")
    out$n <- nrow(df[complete.cases(df[[x]]),])
    out <- out %>% 
      tibble::rownames_to_column(var = "covariate") %>% 
      mutate(biometric = x) %>% 
      tidyr::pivot_wider(names_from = covariate, values_from = c(beta, se, pval))
    colnames(out) <- gsub("scale(urbanicity_score)","urbanicity_score",colnames(out), fixed = T)
    colnames(out) <- gsub("scale(poly(age, 2))1","age",colnames(out), fixed = T)
    colnames(out) <- gsub("scale(poly(age, 2))2","age2",colnames(out), fixed = T)
    colnames(out) <- gsub(":","_",colnames(out))
    out$sex <- y
    out$type = "quadratic interactive"
    return(out)
  }))
})))

models_out <- dplyr::bind_rows(linear_interactive_models, quadratic_interactive_models)

# Adjust p-values for multiple columns
models_out <- models_out %>%
  mutate(across(starts_with("pval_"), ~ p.adjust(.x, method = "BH"), .names = "padj_{.col}"))
colnames(models_out) <- gsub("padj_pval","padj",colnames(models_out))

# Make df of interactive effects of interest
models_out_simple <- rbind(models_out[which(models_out$type=="linear interactive"),
                                            c("biometric","n","sex","type","beta_age_urbanicity_score",
                                              "se_age_urbanicity_score","pval_age_urbanicity_score","padj_age_urbanicity_score")] %>% 
                                dplyr::rename(ageLifestyle_effectSize = beta_age_urbanicity_score,
                                              ageLifestyle_se = se_age_urbanicity_score,
                                              ageLifestyle_pval = pval_age_urbanicity_score,
                                              ageLifestyle_padj = padj_age_urbanicity_score), 
                              models_out[which(models_out$type=="quadratic interactive"),
                                            c("biometric","n","sex","type","beta_age2_urbanicity_score",
                                              "se_age2_urbanicity_score","pval_age2_urbanicity_score","padj_age2_urbanicity_score")] %>% 
                                dplyr::rename(ageLifestyle_effectSize = beta_age2_urbanicity_score,
                                              ageLifestyle_se = se_age2_urbanicity_score,
                                              ageLifestyle_pval = pval_age2_urbanicity_score,
                                              ageLifestyle_padj = padj_age2_urbanicity_score))

biometrics_split <- split(models_out_simple, models_out_simple$biometric)

find_biometrics <- function(df) {
  # Ensure data has both F and M
  if(all(c("F", "M") %in% df$sex)) {
    padj_f <- df$ageLifestyle_padj[df$sex == "F"]
    padj_m <- df$ageLifestyle_padj[df$sex == "M"]
    pval_f <- df$ageLifestyle_pval[df$sex == "F"]
    pval_m <- df$ageLifestyle_pval[df$sex == "M"]
    beta_f <- df$ageLifestyle_effectSize[df$sex == "F"]
    beta_m <- df$ageLifestyle_effectSize[df$sex == "M"]
    
    # Check for the first condition
    condition1 <- (padj_f < 0.1 & pval_m >= 0.05) | (padj_m < 0.1 & pval_f >= 0.05) 
    
    # Check for the second condition
    condition2 <- (padj_f < 0.1 & padj_m < 0.1 & beta_f * beta_m < 0)
    
    # Return biometric if either condition is true
    if(condition1 | condition2) {
      return(df$biometric[1])
    }
  }
  return(NULL)
}
result <- lapply(biometrics_split, find_biometrics)
result <- Filter(Negate(is.null), result)
models_out_simple[which(models_out_simple$biometric %in% names(result)),]
