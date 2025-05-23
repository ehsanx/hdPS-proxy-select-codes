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
# ---------------------------------------------------------------

### initialization before for-loop
# 'a' ranges from 1 to 1000 to change the number of iterations for the for-loop
a <- 10

RD_rf <- data.frame(numProxy = integer(a), RD = numeric(a), SE = numeric(a))
num.proxy_rf <- c()
# ---------------------------------------------------------------

### scenarioER: for-loop generating RD & SE results
set.seed(42)

# Folder "scenarioER"
for (i in 1:a) {
  path <- paste0("../../simData/scenarioER/data_", i, ".rds")
  data <- readRDS(path)
  
  # found some id != idx
  data$idx <- data$id
  rf.fit <- randomForest(outcome ~ ., data = data, importance = TRUE)
  
  # Get important features
  rf.imp <- importance(rf.fit)
  ranked.variables <- rownames(rf.imp)[order(rf.imp[, 1], decreasing = TRUE)]
  proxy_rf <- ranked.variables[ranked.variables %in% proxy.list][1:100]
  proxyform <- paste0(proxy_rf, collapse = "+")
  rhsform <- paste0(c(covform, proxyform), collapse = "+")
  ps.formula <- as.formula(paste0("exposure", "~", rhsform))
  
  W.out_rf <- weightit(ps.formula,
                          data = data, 
                          estimand = "ATE",
                          method = "ps")
  fit.OR_rf <- glm(out.formula,
                      data = data,
                      weights = W.out_rf$weights,
                      family= binomial(link = "logit"))
  fit.RD_rf <- glm(out.formula,
                      data= data,
                      weights= W.out_rf$weights,
                      family=gaussian(link= "identity"))
  sum.RD_rf <- c(length(proxy_rf), 
                    summary(fit.RD_rf)$coef["exposure", c("Estimate")], 
                    sqrt(sandwich::sandwich(fit.RD_rf)[2,2]))
  
  RD_rf[i,] <- sum.RD_rf
  rownames(RD_rf)[i] <- paste0("data.", i, "_rf")
  num.proxy_rf <- c(num.proxy_rf, length(proxy_rf))
  names(num.proxy_rf)[i] <- paste0("data.", i, "_rf")
  
  results_rf <- c(i, sum.RD_rf)
  names(results_rf) <- c("iteration", "numProxy", "RD", "SE")
  results_rf.i <- paste0("results_rf.", i)
  assign(results_rf.i, results_rf)
  
  save_dir <- "result/"
  save_path <- paste0(save_dir, results_rf.i, ".RData")
  
  save(list = results_rf.i, file = save_path)
}
# ---------------------------------------------------------------

avg.num.proxy_rf <- mean(num.proxy_rf[complete.cases(num.proxy_rf)])

# ---------------------------------------------------------------

# save(RD_rf, num.proxy_rf, avg.num.proxy_rf,
#      file = "result/rf.RData")
