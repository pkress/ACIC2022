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
# BART + TMLE ----
##################)
#' # BART + TMLE Estimation method
#' ## Outline
#' 
#' BART + TMLE is a doubly robust ensemble method that models both the response
#' and treatment surfaces and exploits the efficient influence curve for 
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
## Set up ----
set.seed(0)
samp_share = 1
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
  # init_sample
  , !c("id.patient", "id.practice"
       , "Y", "Y_transform"
       , "Z", "post")
  ]
z_init = cleaned_data[
  # init_sample
  , Z
  ]
y_init = cleaned_data[
  #init_sample
  , Y_transform
  ]

x_init_all_untr = cleaned_data[
  # init_sample
][
  z_post==1
][
  , z_post:=0
][
  , !c("id.patient", "id.practice"
       , "Y", "Y_transform"
       , "Z", "post")
  ]

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
#' Note that estimation can take a long time. Even with sample of 
{{samp_share}}
#', the estimation with 200 post-burn chains took 
{{ t_elap}}
#' .

## Estimate Outcomes
init_ests = data.table(
  est_val = init_out$yhat.train.mean
  , id.practice = cleaned_data$id.practice
  , treat_flag = cleaned_data$Z
  , z_post = cleaned_data$z_post
  , year = cleaned_data$year
  , act_y = cleaned_data$Y_transform)
init_ests[## Add untreated coutnerfactural. 
  , est_untr := est_val
  ][
  z_post==1
  , est_untr:=init_out$yhat.test.mean
  ]


## Run convergence analysis
data.table(init_out$sigma) %>% 
  .[, id:=1:.N] %>% 
  melt.data.table(id.vars = "id") %>% 
  ggplot()+
  geom_line(aes(x = id, y = value, color = variable))


## Model Treatment ----

treatment_data = practice_data[
  year%in%1:2
  , -c("post")
  ] %>% 
  dcast(id.practice + Z~year, value.var = setdiff(names(.), c("id.practice", "Z", "year")))

treat_x = treatment_data[
  , -c('id.practice', "Z")
  ]
treat_y = treatment_data[, Z]

treat_out = mc.lbart(x.train = data.frame(treat_x), y.train = treat_y
                    , nskip = 50, ndpost = 200
                    , mc.cores = 8)

## Estimate Outcomes
treat_prob = data.table(treat = treat_y
           , prob = treat_out$prob.train.mean
           , id.practice = treatment_data$id.practice
           )
treat_prob[
  , `:=`(H_1 = 1/prob
         , H_0 = 1/(1-prob)
         , H_A = treat/prob - (1-treat)/(1-prob))
]

## Run convergence analysis 
treat_prob  %>% 
  ggplot()+
  geom_density(aes(x = prob, color = factor(true)))


## Fluctuation Param ----

fluct_data = init_ests %>% 
  merge(treat_prob, by.x = c("id.practice", "Z"), by.y = c("id.practice", "treat_y"))

glm_fit <- glm(act_y ~ -1 + offset(qlogis(est_val)) + H_A, data=fluct_data, family=binomial)

eps = coef(glm_fit)

fluct_data[, eps:=eps]

## Update Estimates ----

final_est = fluct_data[
  , `:=`(update_val = plogis(qlogis(est_val)) + eps*clev
         , update_untr = plogis(qlogis(est_untr)) + eps*H_0
  ]

satt = fluct_data[
  z_post==1
  , mean(update_)
]


##################)
# BART + Conformal ----
##################)



## BART mean, 5, 95 quantile Estimates on Training ----



## Q value from test data ----


## Confidence Interval from adjusting quantile with X weighted output ----


