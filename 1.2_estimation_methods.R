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
samp_share = 0.1
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
  ][
  init_sample
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
  # init_sample
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
save(init_out, file = "models/init_pred.RData", compress = F)
t_elap = Sys.time() - t1
print(t_elap)
#' Note that estimation can take a long time. Even with sample of 
{{samp_share}}
#', the estimation with 200 post-burn chains took 
{{ t_elap}}
#' .

## Estimate Outcomes
init_ests = cleaned_data[
  , .(id.practice, Z, z_post, year, Y_transform
      , X1, X2, X3, X4, X5
      , est_val = init_out$yhat.train.mean)
  ][
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
## Run convergence analysis 
treat_prob  %>% 
  ggplot()+
  geom_density(aes(x = prob, color = factor(treat)))


## Update Estimates ----
### Figure out how to implement this...
update_data = init_ests %>% 
  merge(treat_prob
        , by.x = c("id.practice", "Z")
        , by.y = c("id.practice", "treat")) %>% 
  .[
    , .(gn = prob, qn1 = est_val, qn0 = est_untr
        , Y_transform, Z, z_post)
    ]
update_data[
  , `:=`(hyn = 1/(sum(Z)/.N) - gn/((1-gn) * sum(Z)/.N)
         , hgn = (qn1 - qn0 - sum(gn * (qn1-qn0)/(sum(Z)/.N))/.N))
]

tol = 1/update_data[, .N]
hyn_score = update_data[, sum(hyn*(Y_transform - qn]
while()
gn = fluct_data$prob
qn = fluct_data$est_val
clev_y_n = fluct

## Fluctuation Param ----

treat_prob[
  , `:=`(H_1 = 1/(sum(treat_y)/.N)
         , H_0 = -1/(1-prob)
         , H_A = treat/prob - (1-treat)/(1-prob))
]


glm_fit <- glm(Y_transform ~ -1 + offset(qlogis(est_val)) + H_A, data=fluct_data, family=binomial)

eps = coef(glm_fit)

fluct_data[, eps:=eps]

## Final estimates ----

final_est = fluct_data[
  , `:=`(update_val = plogis(qlogis(est_val)) + eps*H_A
         , update_untr = plogis(qlogis(est_untr)) + eps*H_0)
  ][
  , `:=`(val_rescale = update_val*y_range + y_min
         , untr_rescale = update_untr*y_range + y_min)
  ]

satt = fluct_data[
  z_post==1
  , .(overall = mean(val_rescale - untr_rescale, na.rm = T)
      , year_3 = mean(val_rescale[year==3] - untr_rescale[year==3], na.rm = T)
      , year_4 = mean(val_rescale[year==4] - untr_rescale[year==4], na.rm = T))
]
for(v in "X"%p%1:5){
  vals = sort(unique(fluct_data[, get(v)]))
  for(val in vals){
    cat("\n", v%p%"_"%p%val)
    set(x = satt, i = 1L, j = v%p%"_"%p%val
        , value = (
          fluct_data[
            z_post==1 & get(v)==val
            , mean(val_rescale - untr_rescale, na.rm = T)
          ]))
  }
}
satt
