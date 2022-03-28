##################-
### Author: Peter Kress
### Date: 2022/02/16
### Purpose: Functions to Estimate Intervention Effect
##################-

##################-
# Initialize Workspace ----
##################-

#+ include = F
## Paths ----
setwd("/Users/pkress/Documents/Personal Projects/ACIC/")

## Packages ----

if(!require("pacman")) install.packages("pacman")
library(pacman)
p_load(data.table, magrittr, stringr, ggplot2
       , fixest, MatchIt, BART, tmle, SuperLearner)

## Handy Functions ----
`%p%` = paste0

month_diff = function(d1, d2){
  12*(year(d2) - year(d1)) + (month(d2) - month(d1))
}

make_title = function(x){ str_to_title(gsub("_", " ", x))}

##################-
# Read in patient data ----
##################-
#+ include = F
read_files = function(folder, nos){
  lapply(list.files(folder, full.names = T, pattern = paste(nos%p%".csv",collapse = "|")), fread) %>% 
    rbindlist(idcol = "samp")
}

## Read in data
data = "data/track1b/"%p%c("patient", "patient_year", "practice", "practice_year") %>%
  setNames(.,c("patient", "patient_year", "practice", "practice_year")) %>%
  lapply(.,read_files, nos = 1201)

## Create patient data
patient_data = data$patient[
  data$patient_year
  , on = c("id.patient", "samp")
][
  data$practice
  , on = c("id.practice", "samp")
][
  data$practice_year[, .(id.practice, year, Z, post, samp)]
  , on = c("id.practice", "year", "samp")
][## For single samp value, drop samp indicator
  , samp:=NULL
]

## Create practice data
practice_data = patient_data[
  , c(lapply(.SD, mean)%>% 
        setNames(., names(.)%p%"_avg")
      , lapply(c("A", "B", "C"), \(x) sum(V5==x)/.N) %>% 
        setNames(., "V5_"%p%c("A", "B", "C")%p%"_avg")
      , .(n.patients = .N)
  )
  , .(id.practice, year, Z, post)
  , .SDcols = c("Y", "V"%p%c(1:4))
]

##################)
# BART + Conformal ----
##################)
#' # BART + Conformal inference
#' ## Outline

#' BART + conformal inference allows for  
#' inference.
#' 
#' I am using BART + TMLE because it flexibly models the treatment and response 
#' surfaces while giving valid inference. Additionally, this method 
#' was found to be extremely effective 
#' in terms of Bias, RMSE, and interval coverage rates
#' in a previous causal inference challenge (See Dorie
#' et al, 2017, https://arxiv.org/pdf/1707.02641.pdf). 
#' 
#' BART + TMLE has 4 main steps. 
#' 
#' * Initial Response Estimates
#' 
#'   * Use patient and practice covariates
#'   
#' * P-Score Estimates
#' 
#'   * Use only practice covariates, and (optionally) practice level patient 
#'   characteristic summaries from the first 2 years. 
#'   
#' * Response Estimate Updates
#' 
#' * Inference
#' 
#+ include = F

## Dummy Example ----
set.seed(1)
n <- 1000
d <- 5
X <- matrix(rnorm(n * d), nrow = n)
beta <- rep(1, 5)
Y <- X %*% beta + rnorm(n)

# Generate missing indicators
missing_prob <- pnorm(X[, 1])
if_missing <- missing_prob < runif(n)
# Y[if_missing] <- NA

data = data.table(X)
data[
  , `:=`(Y = Y
         , if_missing = if_missing
         , Y_mi = fifelse(if_missing, rep(NA_real_, .N), Y)
         )
]

obj <- conformalCf(X, data$Y, type = "mean",
                   outfun = "RF", useCV = FALSE)
predict(obj, Xtest, alpha = 0.1)
# Generate testing data
ntest <- 5
Xtest <- matrix(rnorm(ntest * d), nrow = ntest)

## Set up ----
set.seed(0)
samp_share = .01
init_sample = sample(nrow(patient_data), round(samp_share*nrow(patient_data)))

y_min = patient_data[, min(Y)]; y_range = patient_data[, max(Y) - min(Y)]

cleaned_data = copy(patient_data)[
  , `:=`(Y_transform=(Y - y_min)/y_range
         , z_post = Z * post)
][
  , lapply(.SD, \(x) {
    if(is.character(x)) x = as.factor(x)
    x
  })
]
x_init = cleaned_data[
  init_sample
  , !c("id.patient", "id.practice"
       , "Y", "Y_transform"
       , "Z", "post")
]
z_init = cleaned_data[
  init_sample
  , Z
]
y_init = cleaned_data[
  init_sample
  , Y_transform
]

aa = conformalCf(X = x_init, Y = y_init, estimand = "nonmissing", type = "mean"
                 , side = "two")
## Initial Response Estimates ----
## Run estimation
t1 = Sys.time()
init_out = mc.wbart(x.train = data.frame(x_init), y.train = y_init
                    , x.test = data.frame(x_init_all_untr)
                    , nskip = 100, ndpost = 250
                    , keepevery = 5
                    , mc.cores = 8)
t_elap = Sys.time() - t1
print(t_elap)