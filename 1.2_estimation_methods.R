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
setwd("/home/pkress/")

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

## Read in data
read_files = function(folder, nos){
  lapply(list.files(folder, full.names = T, pattern = paste(nos%p%".csv",collapse = "|")), fread) %>% 
    rbindlist(idcol = "samp")
}
dataset.nums = 1:2

data = "data/track1a_20220404/"%p%c("patient", "patient_year", "practice", "practice_year") %>%
  setNames(.,c("patient", "patient_year", "practice", "practice_year")) %>%
  lapply(.,read_files, nos = dataset.nums)

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
lapply(dataset.nums, \(dsetnum){
  
  cleaned_data = copy(patient_data)[
    samp==dsetnum
    ][
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
    , !c("samp", "id.patient", "id.practice"
         , "Y", "Y_transform"
         , "Z", "post")
    ]
  z_init = cleaned_data[
    , Z
    ]
  y_init = cleaned_data[
    , Y_transform
    ]
  
  x_init_all_untr = cleaned_data[
  ][
    z_post==1
  ][
    , z_post:=0
  ][
    , !c("samp", "id.patient", "id.practice"
         , "Y", "Y_transform"
         , "Z", "post")
    ]
  
  x_init_all_tr = cleaned_data[
  ][
    z_post==0
  ][
    , z_post:=1
  ][
    , !c("samp", "id.patient", "id.practice"
         , "Y", "Y_transform"
         , "Z", "post")
  ]
  ## Initial Response Estimates ----
  ## Run estimation
  nskip = 100; ndpost = 250
  t1 = Sys.time()
  init_out = mc.wbart(x.train = data.frame(x_init), y.train = y_init
                   , nskip = nskip, ndpost = ndpost
                   , keepevery = 5
                   , mc.cores = 8)
  save(init_out, file = "models/init_pred.RData", compress = F)
  pred_untr = predict(init_out, bartModelMatrix(data.frame(x_init_all_untr)))
  pred_tr = predict(init_out, bartModelMatrix(data.frame(x_init_all_tr)))
  t_elap = Sys.time() - t1
  print(t_elap)
  #' Note that estimation can take a long time. Even with sample of 
  {{samp_share}}
  #', the estimation with 
  {{ ndpost }}
  #' post-burn chains took 
  {{ t_elap}}
  #' .
  
  ## Estimate Outcomes
  init_ests = cleaned_data[
    , .(samp, id.practice, Z, z_post, year, Y_transform
        , X1, X2, X3, X4, X5
        , est_val = init_out$yhat.train.mean)
    ][
    , est_untr := est_val
    ][
    z_post==1
    , est_untr:=rowMeans(t(pred_untr))
    ][
    , est_tr:= est_val
    ][
    z_post==0
    , est_tr:=rowMeans(t(pred_tr))
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
      , .(samp, year, X1, X2, X3, X4, X5
          , gn = prob
          , qn = fifelse(est_val<=0.00001, 0.00001, est_val)
          , qn1 = fifelse(est_tr<=0.00001, 0.00001, est_tr)
          , qn0 = fifelse(est_untr<=0.00001,0.00001,est_untr)
          , Y_transform, Z, z_post)
      ]
  
  update_data[
    , `:=`(hyn = as.numeric(Z==1)/(sum(Z)/.N) - as.numeric(Z==0)*gn/((1-gn) * sum(Z)/.N)
           , hgn = (qn1 - qn0 - sum(gn * (qn1-qn0)/(sum(Z)/.N))/.N))
  ]
  
  tol = 1/update_data[, 100*.N]
  diff = Inf
  k = 0
  while(!all(abs(diff)<tol)){
    k = k +1
    cat("\nIter: ", k)
    g_fit = glm(Z ~ -1 + offset(qlogis(gn)) + hgn, data=update_data, family=binomial)
    update_data[, gn := plogis(qlogis(gn) + g_fit$coefficients * hgn)]
    update_data[
      , hyn:=as.numeric(Z==1)/(sum(Z)/.N)-as.numeric(Z==0)*gn/((1-gn) * sum(Z)/.N)
    ]
    q_fit = glm(Y_transform ~ -1 + offset(qlogis(qn)) + hyn, data=update_data, family=binomial)
    update_data[
      , `:=`(qn1 = plogis(qlogis(qn1) + q_fit$coefficients*hyn)
             , qn0 = plogis(qlogis(qn0) + q_fit$coefficients*hyn)
             , qn = plogis(qlogis(qn) + q_fit$coefficients*hyn))
      ]
    update_data[
      , hgn:=qn1 - qn0 - sum(gn * (qn1-qn0)/(sum(Z)/.N))/.N
    ]
    diff = c(hyn = update_data[, 1/.N * sum(hyn*(Y_transform - qn))]
             , hgn = update_data[, 1/.N * sum(hgn*(Z - gn))])
    print(max(abs(diff)))
  }
  cat("Solved for parameters in ", k, " iterations!")
  
  ## Generate final estimates
  ## Generate filters for each subgroup of interest
  subgroups = lapply("X"%p%1:5, \(v){
    vals = sort(unique(update_data[, get(v)]))
    sapply(vals, \(val,v){paste0(v,"==","'",val,"'")}, v = v) %>% 
      setNames(.,gsub("'", "", gsub("==", "_",.)))
  }) %>% unlist()
  
  scens = "z_post==1"%p%c(
    c("", "&year==3", "&year==4")
    , "&"%p%subgroups) %>% 
    setNames(.,c(c("all", "y3", "y4"), names(subgroups)))
  
  ## For each subgroup, calculate the satt and 90% confints
  outests = lapply(scens, \(x){
    out = update_data[
      eval(parse(text = x))
    ][
      , att := 1/.N * sum(gn/(sum(Z)/.N) * (qn1 - qn0))
    ][
      , .(dataset.num = first(samp)
          , att = first(att)
          , sd_dn = sd(hyn*(Y_transform - qn) + hgn*(Z - gn)
                       + gn/(sum(Z)/.N) * (qn1 - qn0 - att)))
    ][
      ,  `:=`(att = att*y_range + y_min
             , sd_dn = sd_dn*y_range)
    ][
      , .(dataset.num
          , satt = att
          , lower90 = att - 1.645*sd_dn
          , upper90 = att + 1.645*sd_dn)
    ]
  }) %>% 
    rbindlist(idcol = "scen")
  
  outcln = outests[
    , .(dataset.num
        , variable = fifelse(scen%in%c("overall", "y3", "y4"), "Overall"
                             , str_extract(scen, "X[1-5]"))
        , level = fifelse(scen%in%c("overall", "y3", "y4"), "NA"
                          , gsub("_", "", str_extract(scen, "_[ABC01]")))
        , year = fifelse(scen%in%c("y3", "y4"), str_extract(scen, "[34]$"), "NA")
        , satt, lower90, upper90)
  ]
  outcln
}) %>% 
  rbindlist()