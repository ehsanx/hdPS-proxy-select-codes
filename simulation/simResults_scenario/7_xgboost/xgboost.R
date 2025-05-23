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
path <- paste0("../../simData/scenario/data_1.rds")
data <- readRDS(path)
proxy.list <- names(data[, c(grep("^rec", names(data), value = TRUE))])
covarsTfull <- c(investigator.specified.covariates, proxy.list)
Y.form <- as.formula(paste0(c("outcome~ exposure", 
                              covarsTfull), collapse = "+") )
# ---------------------------------------------------------------

### initialization before for-loop
# 'a' ranges from 1 to 1000 to change the number of iterations
a <- 10

RD_xgboost <- data.frame(numProxy = integer(a), RD = numeric(a), SE = numeric(a))
num.proxy_xgboost <- c()
# ---------------------------------------------------------------

### scenario: for-loop generating RD & SE results
set.seed(42)

# Folder "scenario"
for (i in 1:a) {
  path <- paste0("../../simData/scenario/data_", i, ".rds")
  data <- readRDS(path)
  
  # found some id != idx
  data$idx <- data$id
  covar.mat <- model.matrix(Y.form, data = data)[,-1]
  
  xgb.fit <- xgboost(data = covar.mat, label = data$outcome, 
                     max.depth = 30, eta = 1, nthread = 2, nrounds = 5, 
                     objective = "binary:logistic")
  xgb.imp <- xgb.importance(model = xgb.fit)
  sel.variables <- xgb.imp$Feature
  proxy_xgboost <- proxy.list[proxy.list %in% sel.variables]
  proxyform <- paste0(proxy_xgboost, collapse = "+")
  rhsform <- paste0(c(covform, proxyform), collapse = "+")
  ps.formula <- as.formula(paste0("exposure", "~", rhsform))
  
  W.out_xgboost <- weightit(ps.formula,
                          data = data, 
                          estimand = "ATE",
                          method = "ps")
  fit.OR_xgboost <- glm(out.formula,
                      data = data,
                      weights = W.out_xgboost$weights,
                      family= binomial(link = "logit"))
  fit.RD_xgboost <- glm(out.formula,
                      data= data,
                      weights= W.out_xgboost$weights,
                      family=gaussian(link= "identity"))
  sum.RD_xgboost <- c(length(proxy_xgboost), 
                    summary(fit.RD_xgboost)$coef["exposure", c("Estimate")], 
                    sqrt(sandwich::sandwich(fit.RD_xgboost)[2,2]))
  
  RD_xgboost[i,] <- sum.RD_xgboost
  rownames(RD_xgboost)[i] <- paste0("data.", i, "_xgboost")
  num.proxy_xgboost <- c(num.proxy_xgboost, length(proxy_xgboost))
  names(num.proxy_xgboost)[i] <- paste0("data.", i, "_xgboost")
  
  results_xgboost <- c(i, sum.RD_xgboost)
  names(results_xgboost) <- c("iteration", "numProxy", "RD", "SE")
  results_xgboost.i <- paste0("results_xgboost.", i)
  assign(results_xgboost.i, results_xgboost)
  
  save_dir <- "result/"
  save_path <- paste0(save_dir, results_xgboost.i, ".RData")
  
  save(list = results_xgboost.i, file = save_path)
}
# ---------------------------------------------------------------

avg.num.proxy_xgboost <- mean(num.proxy_xgboost[complete.cases(num.proxy_xgboost)])

# ---------------------------------------------------------------

# save(RD_xgboost, num.proxy_xgboost, avg.num.proxy_xgboost,
#      file = "result/xgboost.RData")

