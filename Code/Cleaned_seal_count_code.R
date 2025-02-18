#### SETUP ####

## SET WORKING DIRECTORY TO WHEREVER THE DATA CSV IS LOCATED


## READ IN DATA
full_data=read.csv("../Data/full_data.csv", header=T)
full_data=full_data[1:134,]

## ADJUST DATASET FORMATTING

str(full_data)

full_data$date=as.Date(full_data$date)
full_data$month= as.factor(full_data$month) # FOR USING AS A RANDOM EFFECT
full_data$observ_id= as.numeric(full_data$observ_id) # FOR USING AS A SEQUENTIAL TIME VARIABLE IN ar()
full_data$time= as.numeric(full_data$time) # A DISCRETE NUMERIC VARIABLE

install.packages("renv")
library(renv)

## LOAD NEEDED PACKAGES

install.packages("ggplot2")
install.packages("writexl")
install.packages("openxlsx")
install.packages("posterior")
install.packages("brms")
install.packages("rstan")
install.packages("bayesplot")
install.packages("dplyr")

library(posterior) # FOR RANK-NORMALIZED EFFECTIVE SAMPLE SIZE AND RHAT CALCUALTIONS. ALSO ASSUMPTION CHECKING
library(writexl) # CREATING EXCEL SHEETS
library(openxlsx) # TURN DATA FRAME INTO XLSX FILE
library(ggplot2) # GRAPHING
library(dplyr) # FOR FORMATTING AND SUMMARIZE FUNCTIONS
library(brms) # FOR BRM MODEL AND LOO
library(rstan) # TO MAKE STAN RUN FASTER 
library(bayesplot) # PLOPT PARAMETERS IN MCMC_area

## HELP STAN RUN FASTER
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())


#### GENERAL DATA EXPLORATION ####

## QUICK CHECK FOR NORMALITY
shapiro.test(full_data$count) # NOT NORMAL

# CREATE A QQ PLOT 
qqnorm(full_data$count)
qqline(full_data$count, col = 2)  # SLIGHT S-SHAPE WITH UPWARD BEND


#### FINDING THE BEST DISTRIBUTION FOR FULL MODEL ####

### COUNT DATA SO NEGATIVE BINOMIAL OR POISSON
## EXAMINE DISTRIBUTION AND DETERMINE DEGREE OF OVERDISPERSION

## DISTRIBUTION HISTOGRAM 
ggplot(full_data, aes(x=count)) + 
  geom_histogram(binwidth=1, aes(y=..density..), colour="black", fill="white")+
  geom_density()+
  stat_density(alpha=.2,adjust = 1, fill="#FF6666")+
  xlab("Number of Seals Hauled-out")+ylab("Density")+
  theme(panel.background = element_blank())

## DOES THE MEAN EQUAL THE VARIANCE? Does the mean equal the variance

var(full_data$count)
mean(full_data$count)

## THE DISTRIBUTION RESEMBLES A POISSON DISTRIBUTION.
## THE VARIANCE IS ~2 TIMES LARGER THAN THE MEAN.
## WILL PROCEED WITH POISSON SINCE OVERDISPERSION IS ONLY INTERMEDIATE



### COLLINEARITY CHECK 
## CREATE A CORELATION MATRIX TO DETERMINE IF ANY VARIABLES ARE CORELATED WITH EACH OTHER

## FOR MATRIX MONTH NEED TO BE NUMERIC
full_data$month= as.numeric(full_data$month)

## RUN PAIRWISE cor BETWEEN ALL PREDICTORS WITH A CUT-OFF OF +/- 0.7
cor.matrix<-cor(full_data[,c("avgnoise", "month", "PC1", "time")])
cor.matrix[abs(cor.matrix)< 0.7]<-NA # KEEP ONLY LARGE CORRELATIONS IN SAME MODEL 
cor.matrix

## NO COLLINEARITY DETECTED
## RESET MONTH TO A FACTOR
full_data$month= as.factor(full_data$month)

#### TESTING FOR AUTOCORRELATION ####

### DETERMINE IF THERE IS ANY TEMPORAL AUTOCORRELATION IN OUR DATA
## SINCE WE IT IS MORE LIKELY WE MEASURED THE SAME INDIVIDUALS ON OBSERVATIONS TEMPORALLY CLOSER TO EACH OTHER THAN ONES FURTHER APART THERE IS INHERENT  AUTOCORRELATION
## CANNOT RUN A ACF PLOT ON THE COUNT DATA BECAUSE IT IS NON-CONTINUOUS (OBSERVATIONS UNEQUALLY SPACED APART)

## CHECK FOR AUTOCORREALTION BETWEEN MONTHS BY LOOKING AT THE RESIDUALS OF MONTH
fit1<- lm(count~ month, data = full_data)
plot(full_data$month, fit1$res, pch=20, col="blue")
abline(h=0) # Add the horizontal line at 0

## IF THE DATA POINTS WITHIN EACH MONTH STAY ABOVE OR BELOW 0 THEN THERE WOULD BE AUTOCORRELATION
## DOES NOT APPEAR THAT THEY DO
## WE ARE STILL GOING TO INCLUDE AN AUTOREGRESSIVE TERM IN THE MODEL BECAUSE OF THE INHERENT AUTOCORREALTION IN THE DATA


#### SET PRIORS FOR ALL MODELS ####

### DEFINE MY PRIORS FOR THE FULL MODEL AND CANDIDATE MODELS

## SET PRIOR FOR NULL MODEL
null_priors <- c(
  set_prior("cauchy(0, 2)", class = "sd", group = "month"),
  set_prior("uniform(-1, 1)", class = "ar")
)

## SET PRIOR FOR FULL MODEL AND CANDIDATE model1
full_priors <- c(
  set_prior("normal(0, 10)", class = "b", coef = "avgnoise"),
  set_prior("normal(0, 10)", class = "b", coef = "PC1"),
  set_prior("normal(0, 10)", class = "b", coef = "time"),
  set_prior("cauchy(0, 2)", class = "sd", group = "month"), # INCLUDE TO INDICATE THAT THERE IS LIKELY SOME VARIATION BETWEEN MONTHS
  set_prior("uniform(-1, 1)", class = "ar") # UNIFORM BECAUSE THERE IS EQUAL CHANCE OF THE EFFECT BEING ANYWHERE FROM -1 TO 1
)

## SET PRIOR FOR FULL MODEL AND CANDIDATE model2
model2_priors = c(
  set_prior("normal(0, 10)", class = "b", coef = "avgnoise"),
  set_prior("cauchy(0, 2)", class = "sd", group = "month"),
  set_prior("uniform(-1, 1)", class = "ar"))

## SET PRIOR FOR FULL MODEL AND CANDIDATE model3
model3_priors = c(
  set_prior("normal(0, 10)", class = "b", coef = "PC1"),
  set_prior("normal(0, 10)", class = "b", coef = "avgnoise"),
  set_prior("cauchy(0, 2)", class = "sd", group = "month"),
  set_prior("uniform(-1, 1)", class = "ar"))


#### FITTING A BAYESIAN FULL MODEL ####

### WE NEED TO FIT A FULL MODEL
## FOR EITHER MODEL WE ARE USING THE LOG LINK FUNCTION TO ENSURE THE PREDICTED MEAN REMAINS POSITIVE
## TO DO THIS WE NEED TO TEST WHETHER A ZERO-INFLATED OR REGULAR POISSON MODEL FITS THE DATA BEST

## FULL MODEL FITTED WITH A REGULAR POISSON DISTRIBUTION

fittedModel_bayes <- brm(count ~ avgnoise + PC1 + time + (1 | month) +
                           ar(time=observ_id, p=1), # INCLUDES AN AUTOREGRESSIVE TERM FOR observ_id
                         family = poisson(link = "log"), 
                         chains = 4,
                         iter = 4000,
                         data = full_data,
                         prior = full_priors,
                         control = list(adapt_delta = 0.95), # REDUCES DIVERGENT TRANSITIONS
                         save_pars = save_pars(all = TRUE)) # NEEDED FOR RUNNING LOO-CV

summary(fittedModel_bayes)

## CHECK FOR REGULAR POISSON MODEL CONVERGENCE
model <- fittedModel_bayes 
plot(model) # THIS SHOULD LOOK LIKE A NORMAL DIST AND FUZZY CATEPILLAR
pp_check(ZImodel, type = "dens_overlay") # SHOULD LINE UP

## FULL MODEL FITTED WITH A ZERO INFLATED POISSON DISTRIBUTION

ZImodel <- brm(count ~ avgnoise + PC1 + time + (1 | month) +
                 ar(time=observ_id, p=1),
               family = zero_inflated_poisson(link = "log"), 
               chains = 4,
               iter = 4000,
               data = full_data,
               prior = full_priors,
               control = list(adapt_delta = 0.95),
               save_pars = save_pars(all = TRUE)) 

summary(ZImodel)

## CHECK FOR ZI POISSON MODEL CONVERGENCE
model <- ZImodel 
plot(model) # THIS SHOULD LOOK LIKE A NORMAL DIST AND FUZZY CATEPILLAR
pp_check(ZImodel, type = "dens_overlay") # SHOULD LINE UP

## COMPARE REGULAR POISSON MODEL TO ZI POISSON MODEL
compare1 = loo_compare(loo(fittedModel_bayes), loo(ZImodel)) 
print(compare1, simplify=F)

## NO SIGNIFICANT DIFFERENCE BETWEEN THE TWO MODELS
## PROCEEDING WITH ZI POISSON FITTED MODEL BECAUSE IT IS LIKELY THERE IS SOME ZI IN OUR DATA EVEN IF IT IS MINOR

### ENSURE THAT THE ZI POISSON FITTED FULL MODEL CONVERGES WELL
## ASSESS MODEL CONVERGENCE USING ESS AND RHAT VALUES

## CHECK RANK-NORMALIZED EFFECTIVE SAMPLE SIZE AND RHAT VALUES 
stanfit_model <- ZImodel$fit # EXTRACT THE STANFIT OBJECT 
draws <- as_draws_array(stanfit_model) # EXTRACT THE SAMPLES INTO A FORMAT THE POSTERIOR CAN WORK WITH 
summary <- summarize_draws(draws) # SUMMARIZE DIAGNOSTICS FOR ALL PARAMETERS
summary[, c("variable", "ess_bulk", "ess_tail", "rhat")] # VIEW THE ESS DIAGNOSITICS
ess_threshold <- 500 # DEFINE ESS THRESHOLD
rhat_threshold <- 0.05  # DEFINE RHAT THRESHOLD
ess_bulk_criteria <- summary$ess_bulk > ess_threshold # CHECK ess_bulk CONDITIONS
ess_tail_criteria <- summary$ess_tail > ess_threshold  # CHECK ess_tail CONDITIONS
rhat_criteria <- abs(summary$rhat - 1) > rhat_threshold  # CHECK RHAT CONDITIONS
all_criteria <- ess_bulk_criteria & ess_tail_criteria & rhat_criteria # COMBINE ALL CRITERIA
n_parameters <- sum(all_criteria, na.rm = TRUE) # COUNT HOW MANY PARAMETERS SATIFY ALL CRITERIA
cat("Number of parameters meeting all criteria:", n_parameters, "\n") # PRINT THE RESULT

## BECAUSE ALL PARAMETERS HAVE FAVORABLE ESS AND RHAT VALUES THE ZI POISSON FITTED MODEL CONVERGES WELL


#### CHECK THE ASSUMPTIONS OF GLMMS WITH A ZI POISSON DISTRIBUTION ####

### BEFORE PROCEEDING FURTHER WITH OUR MODEL WE NEED TO MAKE SURE THE ASSUMPTIONS ARE NOT VIOLATED. 
## KEY ASSUMPTIONS OF THIS MODEL:
# 1. RESIDUALS FOLLOW THE DISTRIBUTION (ZIP MODELS OFTEN HAVE LOTS OF ZEROS AND A SKEWED DISTRIBUTION) 
# 2. ZERO-INFLATION COMPONENT IS PROPERLY ACCOUNTED FOR 
# 3. RANDOM EFFECTS ARE NORMAL
# 4. NO AUTOCORREALTION IN THE RESIDUALS

## 1. CHECK RESIDUAL DISTRIBUTION

# EXTRACT DEVIANCE RESIDUALS
residuals_ZI <- residuals(ZImodel, type = "pearson")[,1]  # Extract as numeric vector

# HISTOGRAM OF RESIDUALS
hist(residuals_ZI, breaks = 30, probability = TRUE, main = "Histogram of Pearson Residuals",
     xlab = "Residuals", col = "skyblue", border = "white")
lines(density(residuals_ZI), col = "darkblue", lwd = 2)

# QQ PLOT FOR RESIDUALS
qqnorm(residuals_ZI, main = "QQ Plot of Residuals")
qqline(residuals_ZI, col = "red", lwd = 2)

# CHECK FOR OVERDISPERSION (SUM OF SQUARED PEARSON RESIDUALS / DF)
overdispersion_stat <- sum(residuals_ZI^2) / (nrow(full_data) - length(fixef(ZImodel)))
print(paste("Overdispersion factor:", round(overdispersion_stat, 2)))


## 2. CHECK ZERO-INFLATION COMPONENT

# COMPARE OBSERVED VS. PREDICTED ZERO COUNTS
observed_zeros <- sum(full_data$count == 0)
predicted_zeros <- sum(predict(ZImodel, type = "response")[,1] < 1e-6)  # Close to zero
print(paste("Observed zeros:", observed_zeros))
print(paste("Predicted zeros:", predicted_zeros))


## 3. CHECK RANDOM EFFECTS DISTRIBUTION

# EXTRACT RANDOM EFFECTS FOR `month`
random_effects_ZI <- ranef(ZImodel)$month[, , "Intercept"]

# QQ PLOT FOR RANDOM EFFECTS
qqnorm(random_effects_ZI, main = "QQ Plot of Random Effects")
qqline(random_effects_ZI, col = "red", lwd = 2)

# HISTOGRAM OF RANDOM EFFECTS
hist(random_effects_ZI, breaks = 20, probability = TRUE, main = "Histogram of Random Effects",
     xlab = "Random Effect (Intercept)", col = "lightgreen", border = "white")
lines(density(random_effects_ZI), col = "darkgreen", lwd = 2)


## 4. CHECK AUTOCORRELATION IN RESIDUALS

# AUTOCORRELATION FUNCTION (ACF) PLOT
acf(residuals_ZI, main = "Autocorrelation of Residuals")



### OVERALL, ASSUMPTIONS ARE MET WELL ENOUGHT TO MOVE FORWARD WITH OUR ZI POISSON MODEL
## THERE MAY BE SOME UNDERDISPERSION BUT NOT VERY MUCH. SINCE OUR MODEL FITS WELL WE ARE NOT WORRIED
## THE ZERO INFLATION ASPECT OF THE MODEL DOES NOT SEEM TO PREDICT ANY ZEROS WHICH IS NOT VERY SUPRISING GIVEN THAT THE POISSON PROCESS ALONE ACCOUNTS FOR MOST IF NOT ALL ZEROS.
## GOING TO STICK WITH ZI POISSON TO ENSURE BECAUSE THE MODEL STILL DOES ESTIMATE AN EFFECT OF ZERO INFLATION



#### CREATE AND RUN CANDIDATE MODELS ####

### CREATE THE CANDIDATE MODELS AND RUN QUICK pp_check TO MAKE SURE THEY CONVERGE AND FIT

## NULL MODEL

null_model <- brm(count ~ 1 + (1 | month) + 
                  ar(time=observ_id, p=1), 
                  family = zero_inflated_poisson(link = "log"), 
                  chains = 4,
                  iter = 4000, 
                  data = full_data, 
                  prior = null_priors,
                  save_pars = save_pars(all = TRUE),
                  control = list(adapt_delta = 0.90))

summary(null_model)
plot(null_model) 
pp_check(null_model) 

## MODEL 1, SAME AS FITTED
model1 <- brm(count ~ avgnoise + PC1 + time + (1 | month) +
              ar(time=observ_id, p=1), 
              family = zero_inflated_poisson(link = "log"), 
              chains = 4, 
              iter = 4000, 
              data = full_data, 
              prior = full_priors,
              control = list(adapt_delta = 0.95), 
              save_pars = save_pars(all = TRUE))

summary(model1)
plot(model1) 
pp_check(model1)


## MODEL 2, JUST NOISE
model2 <- brm(count ~ avgnoise + (1 | month) +
              ar(time=observ_id, p=1), 
              family = zero_inflated_poisson(link = "log"), 
              chains = 4, 
              iter = 4000, 
              data = full_data, 
              prior = model2_priors, 
              save_pars = save_pars(all = TRUE)) # HELPS WITH LOO

summary(model2)
plot(model2)
pp_check(model2) 


## MODEL 3, PC1 AND AVERAGE NOISE, B/C I THINK IT WILL EXPLAIN MORE VARIANCE
model3 <- brm(count ~ avgnoise + PC1 + (1 | month) +
              ar(time=observ_id, p=1), 
              family = zero_inflated_poisson(link = "log"), 
              chains = 4, 
              iter = 4000, 
              data = full_data, 
              prior = model3_priors,
              control = list(adapt_delta = 0.90), # NO DIVERGENT TRANSITIONS AT JUST 0.90
              save_pars = save_pars(all = TRUE)) # HELPS WITH LOO

summary(model3)
plot(model3) 
pp_check(model3) 



#### MODEL SELECTION USING LOO-CV, ORIGINAL VERSION ####

### NOW WE NEED TO DETERMINE WHICH OF THE CANDIDATE MODELS IS THE "BEST"
## DOING THIS USING LOO-CV
## THE FOLLOWING VERSION OF THIS APPROACH HAS PARETO-K ISSUES
## FOLLOWING SECTIONS AVOID THESE

## RUN LOOS, WARNING MESSAGES FOR ALL
loo_null = loo(null_model)
loo1 = loo(model1)
loo2 = loo(model2)
loo3 = loo(model3)

## COMPARE THE FOUR LOOS
compare = loo_compare(loo_null, loo1, loo2, loo3) 
print(compare) # PRINT LOO COMPARISON RESULTS


#### NULL MODEL ROBUST LOO, SAMPLING FROM POSTERIOR DISTRIBUTION WITHOUT PROBLEMATIC OBSERVATIONS ####

### BECAUSE WE SAW WARNINGS ABOUT PARETO-K VALUES > 0.7 WE NEED TO TAKE A MORE ROBUST APPROACH
## WE WILL BE FOLLOWING THE APPROACH FROM Vehtari et al. 2017 TO DEAL WITH THE PROBLEMATIC OBSERVATIONS 

## 1. DEFINE THE FUNCTION TO REFIT THE MODEL EXLCUDING ONE PROBLEMATIC OBSERVATION AND EXTRACT POSTERIOR SAMPLES 
compute_posterior_without_i <- function(model, i) {
  # EXTRACT THE DATA FROM THE MODEL 
  model_data <- model$data
  
  # CREATE A COPY OF THE DATA EXCLUDING OBSERVATION i
  train_data <- model_data[-i, , drop = FALSE]
  
  # REFIT THE MODEL WITHOUT OBSERVATION i, KEEPING THE ar() STRUCTURE INTACT
  model_without_i <- brm(count ~ 1 + (1 | month) + ar(time = observ_id, p = 1),
                         family = zero_inflated_poisson(link = "log"),
                         chains = 4,
                         iter = 4000,
                         data = train_data,
                         prior = null_priors,
                         save_pars = save_pars(all = TRUE),
                         control = list(adapt_delta = 0.90)) # AVOID DIVERGENT TRANSITIONS
  
  # RETURN POSTERIOR SAMPLES FOR THE REFITTED MODEL
  posterior_samples <- as.data.frame(model_without_i)
  return(posterior_samples)
}

## 2. PERFORM INITIAL LOO-CV WITH PSIS
loo_null <- loo(null_model, save_psis = TRUE)

## 3. EXTRACT PARETO K DIAGNOSTICS 
pareto_k_values <- loo_null$diagnostics$pareto_k
print("Pareto k values:")
print(summary(pareto_k_values))

## 4. IDENTIFIY PROBLEMATIC OBSERVATIONS 
k_threshold <- 0.7
problematic_obs <- which(pareto_k_values > k_threshold)
print(paste("Number of problematic observations:", length(problematic_obs)))

## 5. REFIT MODELS EXCLUDING PROBLEMATIC OBSERVATIONS AND EXTRACT POSTERIOR SAMPLES
if(length(problematic_obs) > 0) {
  # STORE POSTERIOR SAMPLES FOR EACH PROBLEMATIC OBSERVATION
  posterior_samples_list <- lapply(problematic_obs, function(i) {
    compute_posterior_without_i(null_model, i)
  })
  
  # COMBINE ALL POSTERIOR SAMPLES FOR ANALYSIS 
  combined_posterior_samples <- do.call(rbind, posterior_samples_list)
  
  # NOW YOU CAN ANALYZE THE POSTERIOR SAMPLES OR COMPUTE OTHER METRICS
  print("Combined posterior samples from models excluding problematic observations:")
  print(head(combined_posterior_samples))
} else {
  print("No problematic observations found. Original LOO estimates are reliable.")
}

## 6. FUNCTION TO COMPUTE LOO FOR THE REFITTED MODELS 
compute_loo_for_posterior_samples <- function(posterior_samples, model_data, i) {
  # EXTRACT THE OBSERVED DATA POINT FOR THE PROBLEMATIC OBSERVATION
  test_data <- model_data[i, , drop = FALSE]
  
  # COMPUTE THE LOG-LIKELIHOOD FOR THE TEST OBSERVATION USING POSTERIOR SAMPLES
  log_lik_test_data <- apply(posterior_samples, 1, function(samples) {
    # EXTRACT THE RELEVANT PARAMETERS FOR THE MODEL
    # FOR SIMPLICITY, ASSUMING POSTERIOR SAMPLES CONTAIN THE RIGHT MODEL PARAMETERS (e.g., INTECEPT AND SLOPES) 
    # RECREATE THE LIKELIHOOD FOR THE i-th OBSERVATION USING THOSE SAMPLES 
    mu <- exp(samples[1])  # ASSUMING THE MODEL IS LOG-LINK
    log_lik <- dpois(test_data$count, lambda = mu, log = TRUE)  # ADJUST FOR YOUR MODEL FAMILY 
    return(log_lik)
  })
  
  # CALCULATE THE EXPECTED LOG PREDICTIVE DENSITY (ELPD)
  elpd <- mean(log_lik_test_data)
  return(elpd)
}

## 7. CALCULATE THE ELPD FOR EACH PROBLEMATIC OBSERVATION
elpd_values <- sapply(problematic_obs, function(i) {
  # GET POSTERIOR SAMPLES FOR THE REFITTED MODEL EXCLUDING THE i-th OBSERVATION
  posterior_samples_i <- posterior_samples_list[[which(problematic_obs == i)]]
  
  # COMPUTE ELPD FOR THIS PROBLEMATIC OBSERVATION
  compute_loo_for_posterior_samples(posterior_samples_i, null_model$data, i)
})

## 8. UPDATE THE LOO OBJECT WITH THE NEW ELPD VALUES FOR THE PROBLEMATIC OBSERVATIONS
loo_updated <- loo_null
loo_updated$pointwise[problematic_obs, "elpd_loo"] <- elpd_values

# RECOMPUTE THE TOTAL ELPD
loo_updated$estimates["elpd_loo", "Estimate"] <- sum(loo_updated$pointwise[, "elpd_loo"])

# UPDATE DIAGNOSTICS TO INDICATE EXACT COMPUTATION FOR PROBLEMATIC OBSERVATIONS
loo_updated$diagnostics$pareto_k[problematic_obs] <- NA  # MARK AS EXACTLY COMPUTED

# PRINT THE UPDATED LOO RESULTS 
print("Updated LOO results with exact ELPD contributions for problematic observations:")
print(loo_updated)

## 9. UPDATED THE LOO OBJECT WITH THE OTHER NEW QUANTITIES (e.g., LOOIC)

# a. RECALCULATE TOTAL ELPD AND STANDARD ERROR
loo_updated$estimates["elpd_loo", "Estimate"] <- sum(loo_updated$pointwise[, "elpd_loo"])
loo_updated$estimates["elpd_loo", "SE"] <- sqrt(sum((loo_updated$pointwise[, "elpd_loo"] - mean(loo_updated$pointwise[, "elpd_loo"]))^2) / (nrow(loo_updated$pointwise) - 1))

# b. RECALCULATE LOOIC AND ITS STANDARD ERROR
loo_updated$estimates["looic", "Estimate"] <- -2 * loo_updated$estimates["elpd_loo", "Estimate"]
loo_updated$estimates["looic", "SE"] <- 2 * loo_updated$estimates["elpd_loo", "SE"]

# c. RECALCULATE p_loo (EFFECTIVE NUMBER OF PARAMETERS) AND ITS STANDARD ERROR 
p_loo <- sum(loo_updated$pointwise[, "p_loo"])  # THIS IS HOW p_loo IS TYPICALLY COMPUTED
loo_updated$estimates["p_loo", "Estimate"] <- p_loo
loo_updated$estimates["p_loo", "SE"] <- sqrt(var(loo_updated$pointwise[, "p_loo"]))  # STANDARD ERROR FOR p_loo


## OPTIONAL COMPARE ORGINAL AND UPDATED ESTIMATES
cat("\nComparison of estimates:\n")
cat("Original total ELPD:", sum(loo_null$pointwise[, "elpd_loo"]), "\n")
cat("Updated total ELPD:", sum(loo_updated$pointwise[, "elpd_loo"]), "\n")

## OPTIONAL COMPARE ORGINAL AND UPDATED POINT WISE ESTIMATES
print(loo_null$pointwise[, "elpd_loo"])
print(sum(loo_updated$pointwise[, "elpd_loo"]))

test= loo_compare(loo_updated, loo_null)
print(test, simplify=F)
loo_null$diagnostics$pareto_k


#### MODEL 1 ROBUST LOO ####

### REPEATING THE SAME APPROACH AS TAKEN WITH THE NULL MODEL
## UPDATED APPROACH FOR PROBLEMATIC OBSERVATIONS WORKAROUND 

## 1. DEFINE THE FUNCTION TO REFIT THE MODEL EXCLUDING ONE OBSERVATION AND EXTRACT POSTERIOR SAMPLES 
compute_posterior_without_i <- function(model, i) {
   
  model_data <- model$data
 
  train_data <- model_data[-i, , drop = FALSE]
  
  model_without_i <- brm(count ~ avgnoise + PC1 + time + (1 | month) +
                           ar(time = observ_id, p = 1),
                         family = zero_inflated_poisson(link = "log"),
                         chains = 4,
                         iter = 4000,
                         data = train_data,
                         prior = full_priors,  
                         save_pars = save_pars(all = TRUE),
                         control = list(adapt_delta = 0.95))  
  
  posterior_samples <- as.data.frame(model_without_i)
  return(posterior_samples)
}

## 2. PERFORM INITIAL LOO-CV WITH PSIS USING THE ORIGINAL MODEL 
loo1 <- loo(model1, save_psis = TRUE)

## 3. EXTRACT PARETO K DIAGNOSITICS 
pareto_k_values <- loo1$diagnostics$pareto_k
print("Pareto k values:")
print(summary(pareto_k_values))

## 4. IDENTIFY PROBLEMATIC OBSERVATIONS
k_threshold <- 0.7
problematic_obs <- which(pareto_k_values > k_threshold)
print(paste("Number of problematic observations:", length(problematic_obs)))

## 5. REFIT MODELS EXCLUDING PROBLEMATIC OBSERVATIONS AND EXTRACT POSTERIOR SAMPLES 
if(length(problematic_obs) > 0) {

  posterior_samples_list <- lapply(problematic_obs, function(i) {
    compute_posterior_without_i(model1, i)
  })
  
  combined_posterior_samples <- do.call(rbind, posterior_samples_list)
  
  print("Combined posterior samples from models excluding problematic observations:")
  print(head(combined_posterior_samples))
} else {
  print("No problematic observations found. Original LOO estimates are reliable.")
}

## 6. FUNCTION TO COMPUTE LOO FOR THE REFITTED MODELS 
compute_loo_for_posterior_samples <- function(posterior_samples, model_data, i) {

  test_data <- model_data[i, , drop = FALSE]
  
  log_lik_test_data <- apply(posterior_samples, 1, function(samples) {
 
    lambda <- exp(samples[1])  
    log_lik <- dpois(test_data$count, lambda = lambda, log = TRUE) 
    return(log_lik)
  })
  
  elpd <- mean(log_lik_test_data)
  return(elpd)
}

## 7. CALCULATE THE ELPD FOR EACH PROBLEMATIC OBSERVATION
elpd_values <- sapply(problematic_obs, function(i) {
  
  posterior_samples_i <- posterior_samples_list[[which(problematic_obs == i)]]
  
  compute_loo_for_posterior_samples(posterior_samples_i, model1$data, i)
})

## 8. UPDATE THE LOO OBJECT WITH THE NEW ELPD VALUES FOR THE PROBLEMATIC OBSERVATIONS
loo_updated1 <- loo1
loo_updated1$pointwise[problematic_obs, "elpd_loo"] <- elpd_values

# RECOMPUTE THE TOTAL ELPD
loo_updated1$estimates["elpd_loo", "Estimate"] <- sum(loo_updated1$pointwise[, "elpd_loo"])

# UPDATE DIAGNOSITICS TO INDICATE EXACT COMPUTATION FOR PROBLEMATIC OBSERVATIONS
loo_updated1$diagnostics$pareto_k[problematic_obs] <- NA  

# PRINT THE UPDATED LOO RESULTS
print("Updated LOO results with exact ELPD contributions for problematic observations:")
print(loo_updated1)

## 9. UPDATE THE LOO OBJECT WITH THE OTHER NEW QUANTITIES (e.g., LOOIC)

# a. RECALCULATE TOTAL ELPD AND STANDARD ERROR 
loo_updated1$estimates["elpd_loo", "Estimate"] <- sum(loo_updated1$pointwise[, "elpd_loo"])
loo_updated1$estimates["elpd_loo", "SE"] <- sqrt(sum((loo_updated1$pointwise[, "elpd_loo"] - mean(loo_updated1$pointwise[, "elpd_loo"]))^2) / (nrow(loo_updated1$pointwise) - 1))

# b. RECALCULATE LOOIC AND ITS STANDARD ERROR
loo_updated1$estimates["looic", "Estimate"] <- -2 * loo_updated1$estimates["elpd_loo", "Estimate"]
loo_updated1$estimates["looic", "SE"] <- 2 * loo_updated1$estimates["elpd_loo", "SE"]

# c. RECALCULATE p_loo (EFFECTIVE NUMBER OF PARAMETERS) AND ITS STANDARD ERROR 
p_loo <- sum(loo_updated1$pointwise[, "p_loo"])  
loo_updated1$estimates["p_loo", "Estimate"] <- p_loo
loo_updated1$estimates["p_loo", "SE"] <- sqrt(var(loo_updated1$pointwise[, "p_loo"]))  

## OPTIONAL COMPARE ORIGINAL AND UPDATED ESTIMATES
cat("\nComparison of estimates:\n")
cat("Original total ELPD:", sum(loo1$pointwise[, "elpd_loo"]), "\n")
cat("Updated total ELPD:", sum(loo_updated1$pointwise[, "elpd_loo"]), "\n")

# OPTIONAL COMPARE ORGINAL AND UPDATED POINT WISE ESTIMATES
print(loo1$pointwise[, "elpd_loo"])
print(sum(loo_updated1$pointwise[, "elpd_loo"]))

# COMPARE LOO RESULTS 
test = loo_compare(loo_updated1, loo1)
print(test, simplify = FALSE)
loo1$diagnostics$pareto_k





#### MODEL 2 ROBUST LOO ####

## 1. DEFINE THE FUNCTION TO REFIT THE MODEL EXCLUDING ONE OBSERVATION AND EXTRACT POSTERIOR SAMPLES
compute_posterior_without_i <- function(model, i) {
  
  model_data <- model$data
  
  train_data <- model_data[-i, , drop = FALSE]
  
  model_without_i <- brm(count ~ avgnoise + (1 | month) + ar(time = observ_id, p = 1),
                         family = zero_inflated_poisson(link = "log"),
                         chains = 4,
                         iter = 4000,
                         data = train_data,
                         prior = model2_priors,
                         save_pars = save_pars(all = TRUE),
                         control = list(adapt_delta = 0.95))
  
  posterior_samples <- as.data.frame(model_without_i)
  return(posterior_samples)
}

## 2. PERFORM INITIAL LOO-CV WITH PSIS
loo2 <- loo(model2, save_psis = TRUE)

## 3. EXTRACT PARETO K DIAGNOSITICS 
pareto_k_values <- loo2$diagnostics$pareto_k
print("Pareto k values:")
print(summary(pareto_k_values))

## 4. IDENTIFY PROBLEMATIC OBSERVATIONS
k_threshold <- 0.7
problematic_obs <- which(pareto_k_values > k_threshold)
print(paste("Number of problematic observations:", length(problematic_obs)))

## 5. REFIT MODELS EXCLUDING PROBLEMATIC OBSERVATIONS AND EXTRACT POSTERIOR SAMPLES 
if(length(problematic_obs) > 0) {
  
  posterior_samples_list <- lapply(problematic_obs, function(i) {
    compute_posterior_without_i(model2, i)
  })
  
  combined_posterior_samples <- do.call(rbind, posterior_samples_list)
  
  print("Combined posterior samples from models excluding problematic observations:")
  print(head(combined_posterior_samples))
} else {
  print("No problematic observations found. Original LOO estimates are reliable.")
}

## 6. FUNCTION TO COMPUTE LOO FOR THE REFITTED MODELS
compute_loo_for_posterior_samples <- function(posterior_samples, model_data, i) {
  
  test_data <- model_data[i, , drop = FALSE]
  
  log_lik_test_data <- apply(posterior_samples, 1, function(samples) {
    mu <- exp(samples[1])
    log_lik <- dpois(test_data$count, lambda = mu, log = TRUE)  
    return(log_lik)
  })
  
  elpd <- mean(log_lik_test_data)
  return(elpd)
}

## 7. CALCULATE THE ELPD FOR EACH PROBLEMATIC OBSERVATION
elpd_values <- sapply(problematic_obs, function(i) {
  posterior_samples_i <- posterior_samples_list[[which(problematic_obs == i)]]
  compute_loo_for_posterior_samples(posterior_samples_i, model2$data, i)
})

## 8. UPDATE THE LOO OBJECT WITH THE NEW ELPD VALUES FOR THE PROBLEMATIC OBSERVATIONS
loo_updated2 <- loo2
loo_updated2$pointwise[problematic_obs, "elpd_loo"] <- elpd_values

# RECOMPUTE THE TOTAL ELPD
loo_updated2$estimates["elpd_loo", "Estimate"] <- sum(loo_updated2$pointwise[, "elpd_loo"])

# UPDATE DIAGNOSTICS TO INDICATE EXACT COMPUTATION FOR PROBLEMATIC OBSERVATIONS
loo_updated2$diagnostics$pareto_k[problematic_obs] <- NA

# PRINT THE UPDATED LOO RESULTS 
print("Updated LOO results with exact ELPD contributions for problematic observations:")
print(loo_updated2)

## 9. UPDATE THE LOO OBJECT WITH THE OTHER NEW QUATITIES (e.g., LOOIC)

# a. RECALCULATE TOTAL ELPD AND STANDARD ERROR
loo_updated2$estimates["elpd_loo", "Estimate"] <- sum(loo_updated2$pointwise[, "elpd_loo"])
loo_updated2$estimates["elpd_loo", "SE"] <- sqrt(sum((loo_updated2$pointwise[, "elpd_loo"] - mean(loo_updated2$pointwise[, "elpd_loo"]))^2) / (nrow(loo_updated2$pointwise) - 1))

# b. RECALCULATE LOOIC AND ITS STANDARD ERROR
loo_updated2$estimates["looic", "Estimate"] <- -2 * loo_updated2$estimates["elpd_loo", "Estimate"]
loo_updated2$estimates["looic", "SE"] <- 2 * loo_updated2$estimates["elpd_loo", "SE"]

# c. RECALCULATE p_loo AND ITS STANDARD ERROR
p_loo <- sum(loo_updated2$pointwise[, "p_loo"])
loo_updated2$estimates["p_loo", "Estimate"] <- p_loo
loo_updated2$estimates["p_loo", "SE"] <- sqrt(var(loo_updated2$pointwise[, "p_loo"]))

# OPTIONAL: COMPARE ORGINAL AND UPDATED ESTIMATES
cat("\nComparison of estimates:\n")
cat("Original total ELPD:", sum(loo2$pointwise[, "elpd_loo"]), "\n")
cat("Updated total ELPD:", sum(loo_updated2$pointwise[, "elpd_loo"]), "\n")

# OPTIONAL: COMPARE ORGINAL AND UPDATED POINTWISE ESTIMATES
print(loo2$pointwise[, "elpd_loo"])
print(sum(loo_updated2$pointwise[, "elpd_loo"]))

test <- loo_compare(loo_updated2, loo2)
print(test, simplify = FALSE)


#### MODEL 3 ROBUST LOO  ####

## 1. DEFINE THE FUNCTION TO REFIT THE MODEL EXCLUDING ONE OBSERVATION AND EXTRACT POSTERIOR SAMPLES
compute_posterior_without_i <- function(model, i) {
  
  model_data <- model$data
  
  train_data <- model_data[-i, , drop = FALSE]
  
  model_without_i <- brm(count ~ avgnoise + PC1 + (1 | month) + ar(time = observ_id, p = 1),
                         family = zero_inflated_poisson(link = "log"),
                         chains = 4,
                         iter = 4000,
                         data = train_data,
                         prior = model3_priors,
                         control = list(adapt_delta = 0.95),
                         save_pars = save_pars(all = TRUE))
  
  posterior_samples <- as.data.frame(model_without_i)
  return(posterior_samples)
}

## 2. PERFORM INITIAL LOO-CV WITH PSIS
loo3 <- loo(model3, save_psis = TRUE)

## 3. EXTRACT PARETO K DIAGNOSITICS
pareto_k_values <- loo3$diagnostics$pareto_k
print("Pareto k values:")
print(summary(pareto_k_values))

## 4. IDENTIFY PROBLEMATIC OBSERVATIONS
k_threshold <- 0.7
problematic_obs <- which(pareto_k_values > k_threshold)
print(paste("Number of problematic observations:", length(problematic_obs)))

## 5. REFIT MODELS EXCLUDING PROBLEMATIC OBSERVATIONS AND EXTRACT POSTERIOR SAMPLES
if(length(problematic_obs) > 0) {
  
  posterior_samples_list <- lapply(problematic_obs, function(i) {
    compute_posterior_without_i(model3, i)
  })
  
  combined_posterior_samples <- do.call(rbind, posterior_samples_list)
  
  print("Combined posterior samples from models excluding problematic observations:")
  print(head(combined_posterior_samples))
} else {
  print("No problematic observations found. Original LOO estimates are reliable.")
}

## 6. FUNCTION TO COMPUTE LOO FOR THE REFITTED MODELS
compute_loo_for_posterior_samples <- function(posterior_samples, model_data, i) {
  
  test_data <- model_data[i, , drop = FALSE]
  
  log_lik_test_data <- apply(posterior_samples, 1, function(samples) {
    mu <- exp(samples[1])
    log_lik <- dpois(test_data$count, lambda = mu, log = TRUE)  
    return(log_lik)
  })
  
  elpd <- mean(log_lik_test_data)
  return(elpd)
}

## 7. CALCULATE THE ELPD FOR EACH PROBLEMATIC OBSERVATION
elpd_values <- sapply(problematic_obs, function(i) {
  posterior_samples_i <- posterior_samples_list[[which(problematic_obs == i)]]
  compute_loo_for_posterior_samples(posterior_samples_i, model3$data, i)
})

## 8. UPDATE THE LOO OBJECT WITH THE NEW ELPD VALUES FOR THE PROBLEMATIC OBSERVATIONS
loo_updated3 <- loo3
loo_updated3$pointwise[problematic_obs, "elpd_loo"] <- elpd_values

# RECOMPUTE THE TOTAL ELPD
loo_updated3$estimates["elpd_loo", "Estimate"] <- sum(loo_updated3$pointwise[, "elpd_loo"])

# UPDATE DIAGNOSTICS TO INDICATE EXACT COMPUTATION FOR PROBLEMATIC OBSERVATIONS
loo_updated3$diagnostics$pareto_k[problematic_obs] <- NA

# PRINT THE UPDATED LOO RESULTS
print("Updated LOO results with exact ELPD contributions for problematic observations:")
print(loo_updated3)

## 9. UPDATE THE LOO OBKECT WITH THE OTHER NEW QUANTITIES (e.g., LOOIC)

# a. RECALCULATE TOTAL ELPD AND STANDARD ERROR
loo_updated3$estimates["elpd_loo", "Estimate"] <- sum(loo_updated3$pointwise[, "elpd_loo"])
loo_updated3$estimates["elpd_loo", "SE"] <- sqrt(sum((loo_updated3$pointwise[, "elpd_loo"] - mean(loo_updated3$pointwise[, "elpd_loo"]))^2) / (nrow(loo_updated3$pointwise) - 1))

# b. RECALCULATE LOOIC AND ITS STANDARD ERROR
loo_updated3$estimates["looic", "Estimate"] <- -2 * loo_updated3$estimates["elpd_loo", "Estimate"]
loo_updated3$estimates["looic", "SE"] <- 2 * loo_updated3$estimates["elpd_loo", "SE"]

# c. RECALCULATE p_loo AND ITS STANDARD ERROR
p_loo <- sum(loo_updated3$pointwise[, "p_loo"])
loo_updated3$estimates["p_loo", "Estimate"] <- p_loo
loo_updated3$estimates["p_loo", "SE"] <- sqrt(var(loo_updated3$pointwise[, "p_loo"]))

## OPTIONAL: COMPARE ORIGINAL AND UPDATED ESTIMATES
cat("\nComparison of estimates:\n")
cat("Original total ELPD:", sum(loo3$pointwise[, "elpd_loo"]), "\n")
cat("Updated total ELPD:", sum(loo_updated3$pointwise[, "elpd_loo"]), "\n")

## OPTIONAL: COMPARE ORIGINAL AND UPDATED POINTWISE ESTIMATES
print(loo3$pointwise[, "elpd_loo"])
print(sum(loo_updated3$pointwise[, "elpd_loo"]))

test <- loo_compare(loo_updated3, loo3)
print(test, simplify = FALSE)






#### ROBUST LOO COMPARISION RESULTS ####

### COMPARING THE LOO RESULTS OF MY CANDIDATE MODELS
## EACH TIME YOU RERUN THE CODE YOU WILL GET SLIGHTLY DIFFERENT RESULTS
## loo_updated IS FOR THE NULL MODEL
## loo_updated1 IS FOR MODEL 1 
## loo_updated2 IS FOR MODEL 2
## loo_updated3 IS FOR MODEL 3

compare = loo_compare(loo_updated, loo_updated1, loo_updated2, loo_updated3) 

print(compare, simplify = FALSE) # PRINT LOO COMPARISON RESULTS


#### FIGURES ####

### CREATE ALL FIGURES RELATING TO SEAL COUNT MODELS 

### FIGURE 2: DENSITY HISTOGRAM OF SEAL COUNT DATA

ggplot(full_data, aes(x = count)) + 
  geom_histogram(
    binwidth = 1, 
    aes(y = ..density..), 
    colour = "black", 
    fill = "white"
  ) +
  geom_density(
    colour = "red", 
    size = 1
  ) +
  stat_density(
    alpha = 0.2, 
    adjust = 1, 
    fill = "#FF6666", 
    colour = NA  # REMOVE BORDER AROUND THE FILLED AREA
  ) +
  xlab("Number of Seals") +
  ylab("Density") +
  theme_classic()

### FIGURE 3: POSTERIOR PREDICTIVE CHECK

ppc = pp_check(fittedModel_bayes) # SHOULD LINE UP

# CUSTOMIZE AXES
ppc + 
  xlab("Number of Seals") + 
  ylab("Density")+
  theme(text = element_text(family = "sans"))

### FIGURE 6: IN-AIR NOISE DISTRIBUTION

## ADD COLUMN TO full_data THAT HAS COMBINED MONTH AND YEAR NUMBERS

full_data$month <- as.numeric(as.character(full_data$month)) # NEED TO BE NUMERIC
full_data$year <- as.numeric(as.character(full_data$year)) # NEED TO BE NUMERIC
full_data$month <- as.factor(full_data$month) # RESET TO FACTOR FOR MDOELING
full_data$year <- as.factor(full_data$year) # RESET TO FACTOR FOR MDOELING

# MAP MONTHS > 12 BACK TO CORRECT MONTH NUMBERS
adjusted_month <- (full_data$month - 1) %% 12 + 1

# CREATE A NEW FORMATTED COLUMN FOR MONTH-YEAR
full_data$month_year <- paste0(month.abb[adjusted_month], " ", full_data$year)

# DEFINE THE DESIRED ORDER
desired_order <- c("Jun 2023", "Jul 2023", "Aug 2023", "Sep 2023", "Oct 2023", 
                   "Nov 2023", "Dec 2023", "Jan 2024", "Feb 2024", "Mar 2024", 
                   "Apr 2024", "May 2024", "Jun 2024") 

# REORDER THE FACTOR VARIABLE
full_data$month_year <- factor(full_data$month_year, levels = desired_order) 

## CALCULATE SAMPLE SIZE FOR EACH month_year
sample_sizes <- full_data %>%
  group_by(month_year) %>%
  summarize(sample_size = n(),
            mean_noise = mean(avgnoise, na.rm = TRUE), # CALCULATE MEAN
            sd_noise = sd(avgnoise, na.rm = TRUE)      # CALCULATE SD
  )

## FINAL GRAPH
ggplot(data = full_data, aes(x = month_year, y = avgnoise)) + 
  geom_violin(trim = FALSE, alpha = 0.3) + 
  geom_boxplot(outlier.colour = "black", outlier.shape = 20, outlier.size = 5, 
               notch = FALSE, alpha = 0.9, width = 0.4) + 
  stat_summary(fun = mean, geom = "point", shape = 17, size = 3, color = "black") + 
  geom_text(data = sample_sizes, aes(x = month_year, y = 71, 
                                     label = paste("n =", sample_size)), 
            inherit.aes = FALSE, size = 5, vjust = 0) + # ADD TEXT ABOVE PLOT
  theme_classic() + 
  theme(text = element_text(size = 20)) +
  labs(y = "In-air noise (dB)", x = "Month") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) # TITL LABELS


### FIGURE 7: SEAL COUNTS OVER THE MONTHS

## CALCULATE SAMPLE SIZE FOR EACH month_year
sample_sizes <- full_data %>%
  group_by(month_year) %>%
  summarize(sample_size = n(),
            mean_count = mean(count, na.rm = TRUE), 
            sd_count = sd(count, na.rm = TRUE)     
  )

## FINAL GRAPH
ggplot(data = full_data, aes(x = month_year, y = count)) + 
  geom_boxplot(outlier.colour = "black", outlier.shape = 20, outlier.size = 5, 
               notch = FALSE, alpha = 0.9, width = 0.4) + 
  stat_summary(fun = mean, geom = "point", shape = 17, size = 3, color = "black") + 
  geom_text(data = sample_sizes, aes(x = month_year, y = 15, 
                                     label = paste("n =", sample_size)), 
            inherit.aes = FALSE, size = 5, vjust = 0) + 
  theme_classic() + 
  theme(text = element_text(size = 20)) +
  labs(y = "Number of Seals", x = "Month") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

### FIGURE 8: POSTERIOR DISTRIBUTIONS

## GRAPH TO LOOK AT EFFECT SIZES OF COVARIATES
mcmc_plot <- mcmc_intervals(
  as.array(model1), 
  pars = c("b_avgnoise", "b_PC1", "b_time"),
  prob = 0.95, # 95% INTERVALS
  prob_outer = 0.99, # 99% INTERVALS
  point_est = "median" # USE MEDIAN AS THE POINT ESTIMATE
) +
  theme_minimal() 
  theme(
    text = element_text(family = "sans"), 
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(), 
    panel.background = element_blank(), 
    axis.line = element_line(color = "black") 
  ) +
  labs(x = "Log-Transformed Estimated Covariate Effect") 

mcmc_plot + scale_y_discrete(
  labels = c(
    "b_avgnoise" = "In-air Noise", 
    "b_PC1" = "Water Current",
    "b_time" = "Time of Day"
  )
)




#### TABLES ####

### CREATE ALL TABLES RELATING TO SEAL COUNT MODELS 

## TABLE 2

# CREATE A TABLE WITH THE LOO INFORMATION
results_table2 = loo_compare(loo_updated, loo_updated1, loo_updated2, loo_updated3) 
print(results_table2, simplify = FALSE) 

# CONVERT THE RESULTS TABLE TO A DATA FRAME
results_df2 <- as.data.frame(results_table2)
results_df2$Model <- rownames(results_table2) # ADD COLUMN WITH MODEL NAMES

# EXPORT RESULTS TABLE TO WORKING DIRECTORY
write_xlsx(results_df2, "Seal_count_results_robust.xlsx")



#### PRIORS SENSITIVITY ANALYSIS ####

### CHECK THAT OTHER VAGUE PRIORS PRODUCE SIMILAR POSTERIOR RESULTS

## DEFINE DIFFERENT COMBINATIONS OF PRIORS TO TEST

normal <- c(
  set_prior("normal(0, 10)", class = "b", coef = "avgnoise"),
  set_prior("normal(0, 10)", class = "b", coef = "PC1"),
  set_prior("normal(0, 10)", class = "b", coef = "time"),
  set_prior("cauchy(0, 2)", class = "sd", group = "month"),
  set_prior("uniform(-1, 1)", class = "ar") 
)

# DECREASE FIXED EFFECT PRIORS BY A SD 0F 5 
combination1 <- c(
  # Prior for avgnoise
  set_prior("normal(0, 5)", class = "b", coef = "avgnoise"),
  set_prior("normal(0, 5)", class = "b", coef = "PC1"),
  set_prior("normal(0, 5)", class = "b", coef = "time"),
  set_prior("cauchy(0, 2)", class = "sd", group = "month"), 
  set_prior("uniform(-1, 1)", class = "ar") 
)

# INCREASE FIXED EFFECT PRIORS BY A SD OF 5
combination2 <- c(
  # Prior for avgnoise
  set_prior("normal(0, 15)", class = "b", coef = "avgnoise"),
  set_prior("normal(0, 15)", class = "b", coef = "PC1"),
  set_prior("normal(0, 15)", class = "b", coef = "time"),
  set_prior("cauchy(0, 2)", class = "sd", group = "month"), 
  set_prior("uniform(-1, 1)", class = "ar") 
)

## MODELS TO COMPARE THE POSTERIOR RESULTS OF: 

normal_model <- brm(count ~ avgnoise + PC1 + time + (1 | month) +
                      ar(time=observ_id, p=1), 
                    family = zero_inflated_poisson(link = "log"), 
                    chains = 4,
                    iter = 4000,
                    data = full_data,
                    prior = normal,
                    control = list(adapt_delta = 0.95)) 

priortest1 = brm(count ~ avgnoise + PC1 + time + (1 | month) +
                   ar(time=observ_id, p=1), 
                 family = zero_inflated_poisson(link = "log"), 
                 chains = 4,
                 iter = 4000,
                 data = full_data,
                 prior = combination1,
                 control = list(adapt_delta = 0.95)) 

priortest2 = brm(count ~ avgnoise + PC1 + time + (1 | month) +
                   ar(time=observ_id, p=1), 
                 family = zero_inflated_poisson(link = "log"), 
                 chains = 4,
                 iter = 4000,
                 data = full_data,
                 prior = combination2,
                 control = list(adapt_delta = 0.95)) 

## POSTERIOR PREDICTIVE CHECKS

pp_check(normal_model)

pp_check(priortest1)

pp_check(priortest2)

## COMPARE LOOS

loonormal = loo(normal_model) # RUN A LOO FOR THE NULL MODEL, WARNING MESSAGES BUT NOT GOING TO DEAL WITH THIS FOR OUR CURRENT PURPOSE
looprior1 = loo(priortest1)
looprior2 = loo(priortest2)

loo_compare(loonormal, looprior1, looprior2) 

## THERE WERE NO LARGE DIFFERENCES BETWEEN THE MODELS WHEN WE TWEAKED OUR PRIORS
## THIS ENSURES THAT OUR PRIOR SPECIFICATION IS APPROPRIATE AND IS NOT GOING TO ALTER THE RESULTS UNDULY

