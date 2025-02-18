#### SETUP ####

## LOAD NEEDED PACKAGES

install.packages("ggplot2")
install.packages("writexl")
install.packages("openxlsx")
install.packages("posterior")
install.packages("brms")
install.packages("rstan")
install.packages("bayesplot")
install.packages("dplyr")
install.packages("gridExtra")
install.packages("tidybayes")
install.packages("tidyr")

library(posterior) # FOR RANK-NORMALIZED EFFECTIVE SAMPLE SIZE AND RHAT CALCUALTIONS
library(writexl) # CREATING EXCEL SHEETS
library(openxlsx) # TURN DATA FRAME INTO XLSX FILE
library(ggplot2) # GRAPHING
library(dplyr) # FOR FORMATTING AND SUMMARIZE FUNCTIONS
library(brms) # FOR BRM MODEL AND LOO
library(rstan) # TO MAKE STAN RUN FASTER 
library(bayesplot) # PLOPT PARAMETERS IN MCMC_area
library(tidybayes) # FOR ASSUMPTION CHECKING
library(gridExtra) # FOR ARRANGING GGPLOTS
library(tidyr) # TO USE drop_na()


## READ IN DATA
full_data=read.csv("../Data/full_data.csv", header=T)
full_data=full_data[1:134,]

## ADJUST DATASET FORMATTING
str(full_data) 

full_data$observ_id= as.numeric(full_data$observ_id) # FOR USING AS A SEQUANTIAL TIME VARIABLE IN AR
full_data$date=as.Date(full_data$date)
full_data$month= as.factor(full_data$month) # FOR USING AS A RANDOM EFFECT
full_data$time= as.numeric(full_data$time) # A DISCRETE NUMERIC VARIABLE
full_data$sealpresence= as.factor(full_data$sealpresence) # O = NO OTHER SEALS PRESENT, 1 = OTHERS PRESENT


## CREATE SURFACE DURATIONS DATASET
surfacetime=  full_data %>%
  drop_na(avgsurfacetime) # REMOVE NAs FROM avgsurfacetime

logsurfacetime= log(surfacetime$avgsurfacetime) # CREATE LOG VALUES FOR AVERAGE SURFACE TIME

surfacetime= cbind(surfacetime, logsurfacetime) # ADD LOG VALUES TO SURFACETIME DATASET

## HELP STAN RUN FASTER
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

#### GENERAL DATA EXPLORATION ####

## CHECK NORMALITY OF AVG SURFACE TIMINGS DISTRIBUTION
ggplot(surfacetime, aes(x=avgsurfacetime)) + 
  geom_histogram(binwidth=1, aes(y=..density..), colour="black", fill="white")+
  geom_density()+
  stat_density(alpha=.2,adjust = 1, fill="#FF6666")+
  xlab("avg time spent at surface")+ylab("Density")+
  theme(panel.background = element_blank())

# QQ PLOT OF AVG SURFACE TIMES
qqnorm(surfacetime$avgsurfacetime, main= "Avg surface time")
qqline(surfacetime$avgsurfacetime, col = 2)  

# VERY RIGHT SKEWED DISTRIBUTION BECAUSE OF TWO OUTLIERS. 
# THESE OUTLIERS ARE NOT ERRONEOUS.
# A LOG TRANSFORMATION WILL IMPROVE THE NORMALITY OF THE DATA.

## GRAPHICALLY CHECK DISTRIBUTION AFTER LOG TRANSFORM
ggplot(surfacetime, aes(x=logsurfacetime)) + 
  geom_histogram(binwidth=0.5, aes(y=..density..), colour="black", fill="white")+
  geom_density()+
  stat_density(alpha=.2,adjust = 1, fill="#FF6666")+
  xlab("log transformed avg time spent at surface")+ylab("Density")+
  theme(panel.background = element_blank())

# QQ PLOT OF LOG TRANSFORMED AVG SURFACE TIMES 
qqnorm(surfacetime$logsurfacetime, main="Logged avg surface time")
qqline(surfacetime$logsurfacetime, col = 2)  

# DISTRIBUTION APPEARS MORE NORMAL
# QQ PLOT IS NOT MUCH IMPROVED


#### IDENTIFY APPROPRIATE LINK FUNCTIONS ####

## PLANNING ON USING A GAUSSIAN DIST. 
## NEED TO CHECK IF A IDENTITY LINK FUNCTION IS APPROPRIATE

# PLOT FOR logsurfacetime vs avgnoise
lineartest1= ggplot(surfacetime, aes(x = avgnoise, y = logsurfacetime)) +
  geom_point() +  # Scatter plot
  labs(title = "Log Surface Time vs Avg Noise", 
       x = "Avg Noise", 
       y = "Log Surface Time")

# PLOT FOR logsurfacetime vs PC1
lineartest2= ggplot(surfacetime, aes(x = PC1, y = logsurfacetime)) +
  geom_point() +  # Scatter plot
  labs(title = "Log Surface Time vs PC1", 
       x = "PC1", 
       y = "Log Surface Time")

# PLOT FOR logsurfacetime vs time
lineartest3= ggplot(surfacetime, aes(x = time, y = logsurfacetime)) +
  geom_point() +  # Scatter plot
  labs(title = "Log Surface Time vs Time", 
       x = "Time", 
       y = "Log Surface Time")

# ARRANGE THE ABOVE 3 PLOTS
grid.arrange(lineartest1, lineartest2, lineartest3, ncol = 3, nrow = 1)

# NO DRAMATIC DEVIATION FROM A LINEAR RELATIONSHIP
# WILL USE AN IDENTITY LINK, AS LONG AS THE MODEL FITS WELL

#### FINDING THE BEST DISTRIBUTION FOR THE FULL MODEL ####

## PLAN TO USE gaussian() BUT WANT TO CHECK THAT THIS IS THE BEST FITTING MODEL DIST. 
## THE IDENTITY LINK FUNCTION IS THE DEFAULT FOR gaussian()


## CHECK MODEL WITH GAUSSIAN DIST. AND LOG TRANSFORMED DATA

fit_log_gaussian <- brm(logsurfacetime ~ avgnoise + PC1 + time + sealpresence + (1 | month) +
                        ar(time=observ_id, p=1), # INCLUDES AN AUTOREGRESSIVE STRUCTURE FOR observ_id
                        family = gaussian(), 
                        chains = 4,
                        iter = 4000,
                        data = surfacetime,
                        prior = STfull_priors, # MUST RUN THE CODE TO CREATE PRIORS BEFORE RUNNING THIS
                        control = list(adapt_delta = 0.95)) # HELPS WITH DIVERGENT TRANSITIONS

# CHECK MODEL WITH A STUDENT-T DIST. AND LOG TRANSFORMED DATA
fit_student_log <- brm(logsurfacetime ~ avgnoise + PC1 + time + sealpresence + (1 | month) +
                       ar(time=observ_id, p=1), 
                       family = student(link = "identity"), 
                       chains = 4,
                       iter = 4000,
                       data = surfacetime,
                       prior = STfull_priors,
                       control = list(adapt_delta = 0.90)) 

## COMPARE THE MODELS 

loo_log_gaussian = loo(fit_log_gaussian) # CALCULATE LOO INFO
loo_student_log = loo(fit_student_log)

compare = loo_compare(loo_log_gaussian, loo_student_log) 

print(compare, simplify = FALSE) # PRINT LOO COMPARISON RESULTS

## THE MODEL WITH A STUDENT-T DIST. AND THE LOG TRANSFORMED DATA OUTPERFORMED THE OTHERS.
## WILL BE PROCEEDING WITH THIS MODEL



#### CHECK THE ASSUMPTIONS OF GLMMS WITH A STUDENT-T DISTRIBUTION ####

## BEFORE PROCEEDING FURTHER WITH OUR MODEL WE NEED TO MAKE SURE THE ASSUMPTIONS ARE NOT VIOLATED. 
## KEY ASSUMPTIONS OF THIS MODEL:
# 1. MODEL RESIDUALS FOLLOW A STUDENT-T DIST.
# 2. DEGREES OF FREEDEM (nu) HAS A LOW MEDIAN VALUE
# 3. RANDOM EFFECTS FOLLOW A GAUSSIAN DIST.
# 4. NO AUTOCORRELATION IN THE RESIDUALS

## EXTRACT RESIDUALS
residuals <- residuals(fit_student_log)

## 1. RESIDUAL DIAGNOSTICS TO MAKE SURE THEY FOLLOW A THE DIST.

# HISTOGRAM OF RESIDUALS. HAS HEAVY TAILS
hist(residuals, breaks = 30, probability = TRUE, main = "Histogram of Residuals",
     xlab = "Residuals", col = "skyblue", border = "white")
lines(density(residuals), col = "darkblue", lwd = 2)

## 2. EVALUATE DEGREES OF FREEDOM

# POSTERIOR SUMMARY FOR DEGREES OF FREEDOM
summary(fit_student_log)$spec_pars  # EXTRACTS nu (df). LESS THAN 10 SUPPORTS t-distribution

# PLOT POSTERIOR FOR DEGREES OF FREEDOM
posterior_samples <- as_draws_df(fit_student_log)
hist(posterior_samples$nu, breaks = 30, probability = TRUE, 
     main = "Posterior Distribution of Degrees of Freedom",
     xlab = "Degrees of Freedom")
abline(v = median(posterior_samples$nu), col = "red", lwd = 2)

## 3. RANDOM EFFECTS ASSUMPTIONS
# EXTRACT RANDOM EFFECTS
random_effects <- ranef(fit_student_log)$month[, , "Intercept"]

# QQ PLOT FOR RANDOM EFFECTS
qqnorm(random_effects)
qqline(random_effects, col = "red", lwd = 2)

# HISTOGRAM OF RANDOM EFFECTS
hist(random_effects, breaks = 20, probability = TRUE, main = "Histogram of Random Effects",
     xlab = "Random Effect (Intercept)", col = "lightgreen", border = "white")
lines(density(random_effects), col = "darkgreen", lwd = 2)

## 4. AUTOREGRESSIVE STRUCTURE

# ACF
acf(residuals)



## ALSO NEED TO CHECK THAT THERE IS NO COLLINEARITY IN THE DATA

# MONTH AND SEAL PRESENCE NEEDS TO BE SET TO NUMERIC
surfacetime$sealpresence= as.numeric(surfacetime$sealpresence)
surfacetime$month= as.numeric(surfacetime$month)

# RUN PAIRWISE COR BETWEEN ALL INDEPENDENT VARIABLES
# CUT-OFF IS +/- 0.7
cor.matrix<-cor(surfacetime[,c("avgnoise", "month", "PC1", "time", "sealpresence")]) 
# KEEP ONLY LARGE CORRELATIONS IN THE SAME MODEL
cor.matrix[abs(cor.matrix)< 0.7]<-NA
cor.matrix

# MAKE SURE TO GO BACK AND RESET SEALPRESENCE AND MONTH TO FACTORS

## OVERALL, THERE ARE NO MAJOR ASSUMPTION VIOLATIONS SO WILL PROCEED WITH THIS FULL MODEL




#### TESTING FOR AUTOCORRELATION ####

### DETERMINE IF THERE IS ANY TEMPORAL AUTOCORRELATION IN OUR DATA
## SINCE WE IT IS MORE LIKELY WE MEASURED THE SAME INDIVIDUALS ON OBSERVATIONS TEMPORALLY CLOSER TO EACH OTHER THAN ONES FURTHER APART THERE IS INHERENT  AUTOCORRELATION
## CANNOT RUN A ACF PLOT ON THE COUNT DATA BECAUSE IT IS NON-CONTINUOUS (OBSERVATIONS UNEQUALLY SPACED APART)

## CHECK FOR AUTOCORREALTION BETWEEN MONTHS BY LOOKING AT THE RESIDUALS OF MONTH
fit1<- lm(count ~ month, data = surfacetime)
plot(surfacetime$month, fit1$res, pch=20, col="blue")
abline(h=0) # Add the horizontal line at 0

## IF THE DATA POINTS WITHIN EACH MONTH STAY ABOVE OR BELOW 0 THEN THERE WOULD BE AUTOCORRELATION
## DOES NOT APPEAR THAT THEY DO
## WE ARE STILL GOING TO INCLUDE AN AUTOREGRESSIVE TERM IN THE MODEL BECAUSE OF THE INHERENT AUTOCORREALTION IN THE DATA

#### SET PRIORS FOR ALL MODELS ####

### DEFINE MY PRIORS FOR THE FULL MODEL AND CANDIDATE MODELS

## SET PRIORS FOR NULL MODEL
STnull_priors <- c(
  set_prior("student_t(3, 0, 2.5)", class = "sd", group = "month"),
  set_prior("uniform(-1, 1)", class = "ar", lb = -1, ub = 1)
)

## SET PRIORS FOR FULL MODEL AND CANDIDATE model1
STfull_priors <- c(
  set_prior("normal(0, 10)", class = "b", coef = "avgnoise"),
  set_prior("normal(0, 10)", class = "b", coef = "PC1"),
  set_prior("normal(0, 10)", class = "b", coef = "time"),
  set_prior("normal(0, 10)", class = "b", coef = "sealpresence1"),  # USE "sealpresence1" FOR THE NON-REFERENCE LEVEL(SECOND VALUE WHEN YOU DO levels())
  set_prior("student_t(3, 0, 2.5)", class = "sd", group = "month"), # INCLUDE TO INDICATE THAT THERE IS LIKELY SOME VARIATION BETWEEN MONTHS
  set_prior("uniform(-1, 1)", class = "ar", lb = -1, ub = 1) # UNIFORM BECAUSE THERE IS EQUAL CHANCE OF THE EFFECT BEING ANYWHERE FROM -1 TO 1
)

## SET PRIORS FOR CANDIDATE model2
STmodel2_priors = c(
  set_prior("normal(0, 10)", class = "b", coef = "avgnoise"),
  set_prior("normal(0, 10)", class = "b", coef = "sealpresence1"),
  set_prior("student_t(3, 0, 2.5)", class = "sd", group = "month"),
  set_prior("uniform(-1, 1)", class = "ar", lb = -1, ub = 1))

## SET PRIORS FOR CANDIDATE model3
STmodel3_priors = c(
  set_prior("normal(0, 10)", class = "b", coef = "avgnoise"), 
  set_prior("student_t(3, 0, 2.5)", class = "sd", group = "month"),
  set_prior("uniform(-1, 1)", class = "ar", lb = -1, ub = 1))



#### FITTING A BAYESIAN FULL MODEL ####

### NOW THAT WE HAVE DEFINED PRIORS AND KNOW THAT OUR MODEL MEETS ASSUMPTIONS WE NEED TO FORMALLY CREATE A FULL MODEL
### WE ALSO NEED TO CHECK HOW WELL IT FITS AND CONVERGES


## CREATE FULL MODEL
STfittedModel_bayes <- brm(logsurfacetime ~ avgnoise + PC1 + time + sealpresence + (1 | month) +
                             ar(time=observ_id, p=1), 
                           family = student(link = "identity"), 
                           chains = 4,
                           iter = 4000,
                           data = surfacetime,
                           prior = STfull_priors,
                           control = list(adapt_delta = 0.90)) 

summary(STfittedModel_bayes)

## CHECK FOR MODEL CONVERGENCE
STmodel <- STfittedModel_bayes 
plot(STmodel) # THIS SHOULD LOOK LIKE A NORMAL DIST AND FUZZY CATEPILLAR
pp_check(STmodel) # SHOULD LINE UP




##  ASSESS MODEL CONVERGENCE USING ESS (RANK-NORMALIZED EFFECTIVE SAMPLE SIZE) AND RHAT VALUES
stanfit_model <- STfittedModel_bayes$fit # EXTRACT THE STANFIT OBJECT 
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


### THE FULL MODEL CONVERGES AND FITS WELL       


#### CREATE AND RUN CANDIDATE MODELS ####

### CREATE THE CANDIDATE MODELS AND RUN QUICK pp_check TO MAKE SURE THEY CONVERGE AND FIT

## NULL MODEL
STnull<- brm(logsurfacetime ~ 1 + (1 | month) +
             ar(time=observ_id, p=1), 
             family = student(link = "identity"), 
             chains = 4,
             iter = 4000,
             data = surfacetime,
             prior = STnull_priors,
             save_pars = save_pars(all = TRUE)) # HELPS WITH LOO

summary(STnull)
plot(STnull) 
pp_check(STnull) 

## MODEL 1, SAME AS FITTED MODEL
STmodel1<- brm(logsurfacetime ~ avgnoise + PC1 + time + sealpresence + (1 | month) +
               ar(time=observ_id, p=1), 
               family = student(link = "identity"), 
               chains = 4,
               iter = 4000,
               data = surfacetime,
               prior = STfull_priors,
               control = list(adapt_delta = 0.90),
               save_pars = save_pars(all = TRUE)) 


summary(STmodel1)
plot(STmodel1) 
pp_check(STmodel1)

## MODEL 2, WITHOUT TIME or PC1, SEE IF NOISE AND PRESENCE HAVE EFFECT
STmodel2<- brm(logsurfacetime ~ avgnoise + sealpresence + (1 | month) +
               ar(time=observ_id, p=1), 
               family = student(link = "identity"), 
               chains = 4,
               iter = 4000,
               data = surfacetime,
               prior = STmodel2_priors,
               control = list(adapt_delta = 0.90),
               save_pars = save_pars(all = TRUE)) 

summary(STmodel2)
plot(STmodel2) 
pp_check(STmodel2) 

## MODEL 3, ONLY AVGNOISE, IN CASE AL THE VARIANCE IS EXPLAINED BY NOISE
STmodel3<-brm(logsurfacetime ~ avgnoise + (1 | month) +
              ar(time=observ_id, p=1), 
              family = gaussian(), 
              chains = 4,
              iter = 4000,
              data = surfacetime,
              prior = STmodel3_priors,
              control = list(adapt_delta = 0.95),
              save_pars = save_pars(all = TRUE)) 

summary(STmodel3)
plot(STmodel3)
pp_check(STmodel3) 


#### MODEL SELECTION USING LOO-CV ####

## RUN A LOO FOR EACH MODEL
loo_STnull = loo(STnull) 
STloo1 = loo(STmodel1)
STloo2 = loo(STmodel2)
STloo3 = loo(STmodel3)

## COMPARE THE LOO RESULTS
compare = loo_compare(loo_STnull, STloo1, STloo2, STloo3) 
print(compare, simplify = FALSE) # PRINT LOO COMPARISON RESULTS


#### FIGURES ####

## FIGURE 4: DENSITY HISTOGRAM OF SURFACE TIMINGS DATA BEFORE AND AFTER LOG TRANSFORMATION 

# a: BEFORE TRANSFORM
a = ggplot(surfacetime, aes(x = avgsurfacetime)) + 
  geom_histogram(
    binwidth = 10, 
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
    colour = NA  # REMOVES BORDER AROUND THE FILLED AREA
  ) +
  xlab("Average Surface Duration") +
  ylab("Density") +
  theme_classic()

# b: AFTER TRANSFORM
b = ggplot(surfacetime, aes(x = logsurfacetime)) + 
  geom_histogram(
    binwidth = 0.5, 
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
    colour = NA 
  ) +
  xlab("Log Transformed Average Surface Duration") +
  ylab("") +
  theme_classic()

# ARRANGE
grid.arrange(a, b, ncol = 2)


## FIGURE 5: POSTERIOR PREDICTIVE CHECK

ppc = pp_check(STfittedModel_bayes) 

# CUSTOMIZE AXES
ppc + 
  xlab("Average Surface Duration") + 
  ylab("Density")+
  theme(text = element_text(family = "sans"))

## FIGURE 9: SURFACE DURATION OVER THE MONTHS

# ADD COLUMN TO surfacetime THAT HAS COMBINED MONTH AND YEAR NUMBERS

surfacetime$month <- as.numeric(as.character(surfacetime$month)) # NEED TO BE NUMERIC FOR THIS
surfacetime$year <- as.numeric(as.character(surfacetime$year)) # NEED TO BE NUMERIC FOR THIS
surfacetime$month <- as.factor(surfacetime$month) # RESET TO FACTOR FOR MODELING USE
surfacetime$year <- as.factor(surfacetime$year) # RESET TO FACTOR FOR MODELING USE

# MAP MONTHS > 12 BACK TO CORRECT MONTH NUMBERS 
adjusted_month <- (surfacetime$month - 1) %% 12 + 1

# CREATE A NEW FORMATTED COLUMN FOR MONTH-YEAR
surfacetime$month_year <- paste0(month.abb[adjusted_month], " ", surfacetime$year)

# DEFINE THE DESIRED ORDER
desired_order <- c("Jun 2023", "Jul 2023", "Aug 2023", "Sep 2023", "Oct 2023", 
                   "Nov 2023", "Dec 2023", "Jan 2024", "Feb 2024", "Mar 2024", 
                   "Apr 2024", "May 2024", "Jun 2024") 

# CONVERT TO FACTOR WITH THE CORRECT ORDER
surfacetime$month_year <- factor(surfacetime$month_year, levels = desired_order)

# CALCULATE SAMPLE SIZE FOR EACH month_year
sample_sizes <- surfacetime %>%
  group_by(month_year) %>%
  summarize(sample_size = n(),
            mean_surfacetime = mean(avgsurfacetime, na.rm = TRUE), # CALCULATE MEAN
            sd_surfacetime = sd(avgsurfacetime, na.rm = TRUE)      # CALCULATE SD
  )

# GRAPH 
ggplot(data = surfacetime, aes(x = month_year, y = avgsurfacetime)) + 
  geom_boxplot(outlier.colour = "black", outlier.shape = 20, outlier.size = 5, 
               notch = FALSE, alpha = 0.9, width = 0.4) + 
  stat_summary(fun = mean, geom = "point", shape = 17, size = 3, color = "black") + 
  geom_text(data = sample_sizes, aes(x = month_year, y = 160, 
                                     label = paste("n =", sample_size)), 
            inherit.aes = FALSE, size = 5, vjust = 0) + # ADD TEXT ABOVE THE PLOT
  theme_classic() + 
  theme(text = element_text(size = 20)) +
  labs(y = "Surfacing duration (s)", x = "Month") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) # TILT LABELS


## FIGURE 10: POSTERIOR DISTRIBUTION PLOT


mcmc_plot2 <- mcmc_intervals(
  as.array(STmodel1), 
  pars = c("b_avgnoise", "b_PC1", "b_time", "b_sealpresence1"),
  prob = 0.95, # 95% INTERVALS
  prob_outer = 0.99, # 99%
  point_est = "median" # USE MEDIAN AS THE POINT ESTIMATE
) +
  theme_minimal() + 
  theme(
    text = element_text(family = "sans"), 
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(), 
    panel.background = element_blank(), 
    axis.line = element_line(color = "black") 
  ) +
  labs(x = "Log-Transformed Estimated Covariate Effect") 

mcmc_plot2 + scale_y_discrete(
  labels = c(
    "b_avgnoise" = "In-air Noise", 
    "b_PC1" = "Water Current",
    "b_time" = "Time of Day",
    "b_sealpresence1" = "Seal Presence"
  )
)



#### TABLES ####

## TABLE 2

# CREATE A TABLE FOR THE LOO INFO
results_table <- loo_compare(loo_STnull, STloo1, STloo2, STloo3)

# PRINT THE TABLE
print(results_table, simplify=F)

# CONVERT THE RESULTS TABLE TO A DATA FRAME
results_df <- as.data.frame(results_table)
results_df$Model <- rownames(results_table) # ADD COLUMN WITH MODEL NAMES

# EXPORT RESULTS TABLE TO WORKING DIRECTORY
write_xlsx(results_df, "Surface_timings_results.xlsx")

#### PRIORS SENSITIVITY ANALYSIS ####

### CHECK THAT OTHER VAGUE PRIORS PRODUCE SIMILAR POSTERIOR RESULTS

## DEFINE DIFFERENT COMBINATIONS OF PRIORS TO TEST

STnormal <- c(
  set_prior("normal(0, 10)", class = "b", coef = "avgnoise"),
  set_prior("normal(0, 10)", class = "b", coef = "PC1"),
  set_prior("normal(0, 10)", class = "b", coef = "time"),
  set_prior("normal(0, 10)", class = "b", coef = "sealpresence1"),  
  set_prior("student_t(3, 0, 2.5)", class = "sd", group = "month"),
  set_prior("uniform(-1, 1)", class = "ar", lb = -1, ub = 1)
)

# DECREASE FIXED EFFECT PRIORS BY A SD 0F 5 
STcombination1 <- c(
  set_prior("normal(0, 5)", class = "b", coef = "avgnoise"),
  set_prior("normal(0, 5)", class = "b", coef = "PC1"),
  set_prior("normal(0, 5)", class = "b", coef = "time"),
  set_prior("normal(0, 5)", class = "b", coef = "sealpresence1"),  
  set_prior("student_t(3, 0, 2.5)", class = "sd", group = "month"),
  set_prior("uniform(-1, 1)", class = "ar", lb = -1, ub = 1) #
)

# INCREASE FIXED EFFECT PRIORS BY A SD OF 5
STcombination2 <- c(
  set_prior("normal(0, 15)", class = "b", coef = "avgnoise"),
  set_prior("normal(0, 15)", class = "b", coef = "PC1"),
  set_prior("normal(0, 15)", class = "b", coef = "time"),
  set_prior("normal(0, 15)", class = "b", coef = "sealpresence1"),
  set_prior("student_t(3, 0, 2.5)", class = "sd", group = "month"),
  set_prior("uniform(-1, 1)", class = "ar", lb = -1, ub = 1) 
)

## MODELS TO COMPARE THE POSTERIOR RESULTS OF: 

STnormal <- brm(logsurfacetime ~ avgnoise + PC1 + time + sealpresence + (1 | month) +
                  ar(time=observ_id, p=1), 
                family = student(link = "identity"), 
                chains = 4,
                iter = 4000,
                data = surfacetime,
                prior = STnormal,
                control = list(adapt_delta = 0.90))

STpriortest1 = brm(logsurfacetime ~ avgnoise + PC1 + time + sealpresence + (1 | month) +
                     ar(time=observ_id, p=1), 
                   family = student(link = "identity"), 
                   chains = 4,
                   iter = 4000,
                   data = surfacetime,
                   prior = STcombination1,
                   control = list(adapt_delta = 0.90))

STpriortest2 = brm(logsurfacetime ~ avgnoise + PC1 + time + sealpresence + (1 | month) +
                     ar(time=observ_id, p=1), 
                   family = student(link = "identity"), 
                   chains = 4,
                   iter = 4000,
                   data = surfacetime,
                   prior = STcombination2,
                   control = list(adapt_delta = 0.90))

## POSTERIOR PREDICTIVE CHECKS

pp_check(STnormal)

pp_check(STpriortest1)

pp_check(STpriortest2)

## COMPARE LOOS

STloonormal = loo(STnormal) # RUN A LOO FOR THE NULL MODEL
STlooprior1 = loo(STpriortest1)
STlooprior2 = loo(STpriortest2)

loo_compare(STloonormal, STlooprior1, STlooprior2) 

## THERE WERE NO LARGE DIFFERENCES BETWEEN THE MODELS WHEN WE TWEAKED OUR PRIORS
## THIS ENSURES THAT OUR PRIOR SPECIFICATION IS APPROPRIATE AND IS NOT GOING TO ALTER THE RESULTS UNDULY

