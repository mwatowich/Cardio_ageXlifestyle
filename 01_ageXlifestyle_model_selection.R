### Model age x lifestyle effects for cardiometabolic biomarkers of interest
## Determine whether linear or quadratic age curve better fits the data - to look at age * lifestyle effects for the better age effect 

## Author: Marina Watowich
## Date last updated: June 3, 2025


# outcomes = vector of biomarkers (eg. "ldl", "hdl")


# Define the priority order 
priority_order <- c("linear additive", "quadratic additive", "linear interactive", "quadratic interactive")

# Create a custom function to prioritize the model_type
priority_function <- function(x) {
  min_priority <- min(match(x, priority_order, nomatch = length(priority_order) + 1))
  return(priority_order[min_priority])
}

# Model all health variables 
all_models <- as.data.frame(do.call(rbind, lapply(outcomes, function(x){
  df <- df[which(complete.cases(df$age) & complete.cases(df$sex) & 
                        complete.cases(df$urbanicity_score)),]
  
  # Run linear additive model
  linear_additive <- lm(scale(df[[x]]) ~ scale(age) + scale(urbanicity_score) + sex, df)
  linear_additive <- cbind(as.data.frame(summary(linear_additive)$coefficients[-1,-3]),
                   AIC(linear_additive))
  colnames(linear_additive) <- c("beta", "se", "pval","AIC")
  linear_additive$n <- nrow(df[complete.cases(df[[x]]),])
  linear_additive <- linear_additive %>% 
    tibble::rownames_to_column(var = "covariate") %>% 
    mutate(biometric = x) %>% 
    tidyr::pivot_wider(names_from = covariate, values_from = c(beta, se, pval))
  colnames(linear_additive) <- gsub("scale(age)","age",colnames(linear_additive), fixed = T)
  colnames(linear_additive) <- gsub("scale(urbanicity_score)","urbanicity_score",colnames(linear_additive), fixed = T)
  linear_additive$type <- "linear additive"
  
  # Run quadratic additive model
  poly_additive <- lm(scale(df[[x]]) ~ scale(poly(age,2)) + scale(urbanicity_score) + sex, df)
  poly_additive <- cbind(as.data.frame(summary(poly_additive)$coefficients[-1,-3]),
                    AIC(poly_additive))
  colnames(poly_additive) <- c("beta", "se", "pval","AIC")
  poly_additive$n <- nrow(df[complete.cases(df[[x]]),])
  poly_additive <- poly_additive %>% 
    tibble::rownames_to_column(var = "covariate") %>% 
    mutate(biometric = x) %>% 
    tidyr::pivot_wider(names_from = covariate, values_from = c(beta, se, pval))
  colnames(poly_additive) <- gsub("scale(urbanicity_score)","urbanicity_score",colnames(poly_additive), fixed = T)
  colnames(poly_additive) <- gsub("scale(poly(age, 2))1","age",colnames(poly_additive), fixed = T)
  colnames(poly_additive) <- gsub("scale(poly(age, 2))2","age2",colnames(poly_additive), fixed = T)
  colnames(poly_additive) <- gsub(":","_",colnames(poly_additive))
  poly_additive$type <- "quadratic additive"
  
  # Run linear interactive model
  linear_interactive <- lm(scale(df[[x]]) ~ scale(age)*scale(urbanicity_score) + sex, df)
  linear_interactive <- cbind(as.data.frame(summary(linear_interactive)$coefficients[-1,-3]),
                         AIC(linear_interactive))
  colnames(linear_interactive) <- c("beta", "se", "pval","AIC")
  linear_interactive$n <- nrow(df[complete.cases(df[[x]]),])
  linear_interactive <- linear_interactive %>% 
    tibble::rownames_to_column(var = "covariate") %>% 
    mutate(biometric = x) %>% 
    tidyr::pivot_wider(names_from = covariate, values_from = c(beta, se, pval))
  colnames(linear_interactive) <- gsub("scale(age)","age",colnames(linear_interactive), fixed = T)
  colnames(linear_interactive) <- gsub("scale(urbanicity_score)","urbanicity_score",colnames(linear_interactive), fixed = T)
  colnames(linear_interactive) <- gsub(":","_",colnames(linear_interactive))
  linear_interactive$type <- "linear interactive"
  
  # Run quadratic interactive model
  poly_interactive <- lm(scale(df[[x]]) ~ scale(poly(age,2))*scale(urbanicity_score) + sex, df)
  poly_interactive <- cbind(as.data.frame(summary(poly_interactive)$coefficients[-1,-3]),
                          AIC(poly_interactive))
  colnames(poly_interactive) <- c("beta", "se", "pval","AIC")
  poly_interactive$n <- nrow(df[complete.cases(df[[x]]),])
  poly_interactive <- poly_interactive %>% 
    tibble::rownames_to_column(var = "covariate") %>% 
    mutate(biometric = x) %>% 
    tidyr::pivot_wider(names_from = covariate, values_from = c(beta, se, pval))
  colnames(poly_interactive) <- gsub("scale(urbanicity_score)","urbanicity_score",colnames(poly_interactive), fixed = T)
  colnames(poly_interactive) <- gsub("scale(poly(age, 2))1","age",colnames(poly_interactive), fixed = T)
  colnames(poly_interactive) <- gsub("scale(poly(age, 2))2","age2",colnames(poly_interactive), fixed = T)
  colnames(poly_interactive) <- gsub(":","_",colnames(poly_interactive))
  poly_interactive$type <- "quadratic interactive"
  
  # Combine 
  out <- dplyr::bind_rows(linear_additive, poly_additive, linear_interactive, poly_interactive)
  out <- out %>% 
    dplyr::select(c("biometric","AIC","type",everything()))
  out$linear_additive_AIC <- out[out$type == "linear additive",]$AIC
  return(out)
})))

# Adjust p-values for multiple columns
all_models <- all_models %>% 
  group_by(type) %>% 
  mutate(across(starts_with("pval_"), ~ p.adjust(.x, method = "BH"), .names = "padj_{.col}"))
colnames(all_models) <- gsub("padj_pval","padj",colnames(all_models))

# Select model with lower AIC when AIC difference > 2 and select the "priority" model when difference < 2.
best_all_model <- all_models %>%
  group_by(biometric) %>% 
  filter(!(type %in% c("quadratic additive") & pval_age2>0.05), # must at least pass pval threshold for the focal effect
         !(type %in% c("quadratic interactive") & pval_age2_urbanicity_score>0.05),
         !(type %in% c("linear interactive") & pval_age_urbanicity_score>0.05)) %>%
  filter(abs(AIC - min(AIC)) < 2) %>% 
  arrange(match(type, priority_order)) %>%
  filter(row_number() == 1) %>%
  ungroup() %>% 
  as.data.frame()

# From here, determine cardiometabolic biomarkers with age, lifestyle, and age*lifestyle effects. 
# For subsequent age*lifestyle analyses, we look at interactive effects from the model of the better-fit age term (eg. 'quadratic additive' model selected from as best model, age^2*lifestyle effects quantified for pertinent analyses)
