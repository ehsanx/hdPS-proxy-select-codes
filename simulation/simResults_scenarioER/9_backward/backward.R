library(autoCovariateSelection)
library(glmnet)
library(randomForest)
library(xgboost)
library(leaps)
library(GA)
library(MASS)
library(WeightIt)
library(cobalt)
library(ggplot2)
library(dplyr)

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
if (!dir.exists("result")) { dir.create("result") }
# ---------------------------------------------------------------

### global information
exposure <- "obese"
outcome <- "diabetes"
investigator.specified.covariates <- 
  c(# Demographic
    "age.cat", "sex", "education", "race", 
    "marital", "income", "born", "year",
    
    # health history related variables/access
    "diabetes.family.history", "medical.access",
    
    # behavioral
    "smoking", "diet.healthy", "physical.activity", "sleep",
    
    # Laboratory 
    "uric.acid", "protein.total", "bilirubin.total", "phosphorus",
    "sodium", "potassium", "globulin", "calcium.total", 
    "systolicBP", "diastolicBP", "high.cholesterol"
  )
covform <- paste0(investigator.specified.covariates, collapse = "+")
out.formula <- as.formula(paste0("outcome", "~", "exposure"))
path <- paste0("../../simData/scenarioER/data_1.rds")
data <- readRDS(path)
proxy.list <- names(data[, c(grep("^rec", names(data), value = TRUE))])
covarsTfull <- c(investigator.specified.covariates, proxy.list)
Y.form <- as.formula(paste0(c("outcome~ exposure", 
                              covarsTfull), collapse = "+") )
full.formula <- as.formula(paste0("outcome~exposure+",
                                  paste0(covarsTfull, collapse = "+"),
                                  collapse = "+"))
# ---------------------------------------------------------------

### initialization before for-loop
# 'a' ranges from 1 to 1000 to change the number of iterations
a <- 10

RD_backward <- data.frame(numProxy = integer(a), RD = numeric(a), SE = numeric(a))
num.proxy_backward <- c()
# ---------------------------------------------------------------

### scenarioER: for-loop generating RD & SE results
set.seed(42)

# Folder "scenarioER"
for (i in 1:a) {
  path <- paste0("../../simData/scenarioER/data_", i, ".rds")
  data <- readRDS(path)
  
  # found some id != idx
  data$idx <- data$id
  
  stepwise_backward <- regsubsets(full.formula, 
                                  data = data, 
                                  method = "backward", 
                                  nvmax = length(covarsTfull)+1, 
                                  nbest = 10)
  summary_stepwise <- summary(stepwise_backward)
  best_model <- which.max(summary_stepwise$adjr2)
  sel.variables <- names(summary_stepwise$which[best_model,])[summary_stepwise$which[best_model,]]
  proxy_backward <- proxy.list[proxy.list %in% sel.variables]
  proxyform <- paste0(proxy_backward, collapse = "+")
  rhsform <- paste0(c(covform, proxyform), collapse = "+")
  ps.formula <- as.formula(paste0("exposure", "~", rhsform))
  
  W.out_backward <- weightit(ps.formula,
                          data = data, 
                          estimand = "ATE",
                          method = "ps")
  fit.OR_backward <- glm(out.formula,
                      data = data,
                      weights = W.out_backward$weights,
                      family= binomial(link = "logit"))
  fit.RD_backward <- glm(out.formula,
                      data= data,
                      weights= W.out_backward$weights,
                      family=gaussian(link= "identity"))
  sum.RD_backward <- c(length(proxy_backward), 
                    summary(fit.RD_backward)$coef["exposure", c("Estimate")], 
                    sqrt(sandwich::sandwich(fit.RD_backward)[2,2]))
  
  RD_backward[i,] <- sum.RD_backward
  rownames(RD_backward)[i] <- paste0("data.", i, "_backward")
  num.proxy_backward <- c(num.proxy_backward, length(proxy_backward))
  names(num.proxy_backward)[i] <- paste0("data.", i, "_backward")
  
  results_backward <- c(i, sum.RD_backward)
  names(results_backward) <- c("iteration", "numProxy", "RD", "SE")
  results_backward.i <- paste0("results_backward.", i)
  assign(results_backward.i, results_backward)
  
  save_dir <- "result/"
  save_path <- paste0(save_dir, results_backward.i, ".RData")
  
  save(list = results_backward.i, file = save_path)
}
# ---------------------------------------------------------------

avg.num.proxy_backward <- mean(num.proxy_backward[complete.cases(num.proxy_backward)])

# ---------------------------------------------------------------

# save(RD_backward, num.proxy_backward, avg.num.proxy_backward,
#      file = "result/backward.RData")

