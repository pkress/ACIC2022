##################-
### Author: Peter Kress
### Date: 2022/02/16
### Purpose: Functions to Estimate Intervention Effect
##################-

##################-
# Initialize Workspace ----
##################-

## Paths ----
setwd("/home/pkress/")

## Packages ----

if(!require("pacman")) install.packages("pacman")
library(pacman)
p_load(data.table, magrittr, stringr, ggplot2
       , fixest, BART, parallel)

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
  list.files(folder, full.names = T
             , pattern = paste(nos%p%".csv",collapse = "|")) %>% 
    setNames(.,nos) %>% 
    lapply(., fread) %>% 
    rbindlist(idcol = "samp")
}
dataset.nums = 2406:2407

data = "ACICdata/track1c_20220404/"%p%c("patient", "patient_year", "practice", "practice_year") %>%
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


all_out = lapply(dataset.nums, \(dsetnum){
  
  ## Set Up ----
  ## Rescale Y from 0 to 1, and limit to initial sample (if using a sample)
  y_min = patient_data[samp==dsetnum, min(Y)]
  y_range = patient_data[samp==dsetnum, max(Y) - min(Y)]
  
  set.seed(0)
  # samp_share = 1
  # init_sample = patient_data[samp==dsetnum, sample(.N, round(samp_share*.N))]
  init_sample = 1:patient_data[samp==dsetnum, .N]
  
  ## limimt to dataset of current iteration
  cleaned_data = copy(patient_data)[
    samp==dsetnum
    ][## Rescale Y
    , `:=`(Y_transform=(Y - y_min)/y_range
           , z_post = Z * post)
  ][## convert strings to factors
    , lapply(.SD, \(x) {
      if(is.character(x)) x = as.factor(x)
      x
    })
    ][## Limit to sample
    init_sample
    ]
  
  ## Covariates for estimation
  x_init = cleaned_data[
    , !c("samp", "id.patient", "id.practice"
         , "Y", "Y_transform"
         , "Z", "post")
    ]
  
  ## Treatment
  z_init = cleaned_data[
    , Z
    ]
  
  ## Outcome
  y_init = cleaned_data[
    , Y_transform
    ]
  
  ## Create version where treated are untreated in post period
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
  
  ## Create version where untreated are treated in post period
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
  ## Run bart estimation for response
  ncores = parallel::detectCores()
  nskip = 200; ndpost = 200
  t1 = Sys.time()
  init_out = mc.wbart(x.train = data.frame(x_init), y.train = y_init
                   , nskip = nskip, ndpost = ndpost
                   , keepevery = 5
                   , mc.cores = ncores, seed = 0)
  # save(init_out, file = "models/init_pred.RData", compress = F)
  pred_untr = predict(init_out, bartModelMatrix(data.frame(x_init_all_untr)))
  pred_tr = predict(init_out, bartModelMatrix(data.frame(x_init_all_tr)))
  t_elap = Sys.time() - t1
  print(t_elap)
  
  ## Estimate Outcomes
  init_ests = cleaned_data[## Take variables of importance and estimated outputs
    , .(samp, id.practice, Z, z_post, year, Y_transform
        , X1, X2, X3, X4, X5
        , est_val = init_out$yhat.train.mean)
    ][
    , est_untr := est_val
    ][## Create untreated counterfactual estimates
    z_post==1
    , est_untr:=rowMeans(t(pred_untr))
    ][
    , est_tr:= est_val
    ][## Create treated counterfactual estimates
    z_post==0
    , est_tr:=rowMeans(t(pred_tr))
    ]
  
  ## Run convergence analysis
  # data.table(init_out$sigma) %>% 
  #   .[, id:=1:.N] %>% 
  #   melt.data.table(id.vars = "id") %>% 
  #   ggplot()+
  #   geom_line(aes(x = id, y = value, color = variable))
  
  
  ## Model Treatment Probability ----
  
  ## Only use practive level summary statistics to determine treatment probability
  treatment_data = practice_data[
    year%in%1:2
    , -c("post")
    ] %>% 
    dcast(id.practice + Z~year, value.var = setdiff(names(.), c("id.practice", "Z", "year")))
  
  ## covariates for estimation
  treat_x = treatment_data[
    , -c('id.practice', "Z")
    ]
  ## Outcome is treatment
  treat_y = treatment_data[, Z]
  
  ## Run estimation
  ncores = parallel::detectCores()
  nskip = 100; ndpost = 200
  treat_out = mc.lbart(x.train = data.frame(treat_x), y.train = treat_y
                      , keepevery = 5
                       , nskip = nskip, ndpost = ndpost
                      , mc.cores = ncores)
  
  ## Save outcomes of treatment estimate
  treat_prob = data.table(treat = treat_y
             , prob = treat_out$prob.train.mean
             , id.practice = treatment_data$id.practice
             )
  ## Run convergence analysis 
  # treat_prob  %>% 
  #   ggplot()+
  #   geom_density(aes(x = prob, color = factor(treat)))
  
  ## TMLE Targeting: Update Estimates ----
  ## First, take in response and treatment estimates
  ## Then, create initial values for 
  ### gn: treatment probabity
  ### qn: estimated outcome
  ### qn1: estimated treated counterfactual outcome
  ### qn0: estimated untreated counterfactual outcome
  
  ## We treat bottom code estimates to 0.00001 to avoid NaNs
  update_data = init_ests %>% 
    merge(treat_prob
          , by.x = c("id.practice", "Z")
          , by.y = c("id.practice", "treat")) %>% 
    .[
      , .(samp, id.practice, year, X1, X2, X3, X4, X5
          , gn = prob
          , qn = fifelse(est_val<=0.00001, 0.00001, est_val)
          , qn1 = fifelse(est_tr<=0.00001, 0.00001, est_tr)
          , qn0 = fifelse(est_untr<=0.00001,0.00001,est_untr)
          , Y_transform, Z, z_post)
      ]
  
  ## Create initial values for "clever covariates"
  ### hyn: clever covariate for Y
  ### hgn: clever covariate for Z
  update_data[
    , `:=`(hyn = as.numeric(Z==1)/(sum(Z)/.N) - as.numeric(Z==0)*gn/((1-gn) * sum(Z)/.N)
           , hgn = (qn1 - qn0 - sum(gn * (qn1-qn0)/(sum(Z)/.N))/.N)/(sum(Z)/.N))
  ]
  
  ## Iteratively target the ATT (heuristic is 1/N, we choose 1/(100N))
  tol = 1/update_data[, 100*.N]
  diff = Inf
  k = 0
  while(!all(abs(diff)<tol)){
    ## For each iteration, update the estimate for gn using the hgn clever covar
    k = k +1
    cat("\nIter: ", k)
    g_fit = glm(Z ~ -1 + offset(qlogis(gn)) + hgn, data=update_data, family=binomial)
    update_data[, gn := plogis(qlogis(gn) + g_fit$coefficients * hgn)]
    ## Update the estimate for hyn clever covar, and update qn estimates using 
    ## the updated value
    update_data[
      , hyn:=as.numeric(Z==1)/(sum(Z)/.N)-as.numeric(Z==0)*gn/((1-gn) * sum(Z)/.N)
    ]
    q_fit = glm(Y_transform ~ -1 + offset(qlogis(qn)) + hyn, data=update_data, family=binomial)
    update_data[
      , `:=`(qn1 = plogis(qlogis(qn1) + q_fit$coefficients*hyn)
             , qn0 = plogis(qlogis(qn0) + q_fit$coefficients*hyn)
             , qn = plogis(qlogis(qn) + q_fit$coefficients*hyn))
      ]
    ## Update hgn clever covar with the updated qn estimates
    update_data[
      , hgn:=qn1 - qn0 - sum(gn * (qn1-qn0)/(sum(Z)/.N))/.N
    ]
    ## Repeat until covar adjusted total differences for Y and Z are small
    diff = c(hyn = update_data[, 1/.N * sum(hyn*(Y_transform - qn))]
             , hgn = update_data[, 1/.N * sum(hgn*(Z - gn))])
    print(max(abs(diff)))
  }
  cat("Solved for parameters in ", k, " iterations!")
  
  ## Generate final estimates
  calc_out = function(filter_text){
    out = update_data[## limit to the subsample of interest
      eval(parse(text = filter_text))
    ][## Calculate the satt estimate
      , att := 1/.N * sum(gn/(sum(Z)/.N) * (qn1 - qn0))
    ][## Calculate the SD of the influence function for inference
      , .(dataset.num = first(samp)
          , att = first(att)
          , sd_dn = sd(hyn*(Y_transform - qn) + hgn*(Z - gn)
                       + gn/(sum(Z)/.N) * (qn1 - qn0 - att)))
    ][## Rescale to original bounds
      ,  `:=`(att = att*y_range + y_min
              , sd_dn = sd_dn*y_range)
    ][## Return final outputs
      , .(dataset.num
          , satt = att
          , lower90 = att - 1.645*sd_dn
          , upper90 = att + 1.645*sd_dn)
    ]
  }
  
  ## Generate filters for each subgroup of interest
  ### Start with full sample and specified subgroups based on years and covars
  subgroups = lapply("X"%p%1:5, \(v){
    vals = sort(unique(update_data[, get(v)]))
    sapply(vals, \(val,v){paste0(v,"==","'",val,"'")}, v = v) %>% 
      setNames(.,gsub("'", "", gsub("==", "_",.)))
  }) %>% unlist()
  
  scens = "z_post==1"%p%c(
    c("", "&year==3", "&year==4")
    , "&"%p%subgroups) %>% 
    setNames(.,c(c("all", "y3", "y4"), names(subgroups)))
  
  ### Also create scenarios for each practice
  prac_scens = "z_post==1&id.practice=="%p%update_data[, unique(id.practice)] %>% 
    setNames(.,update_data[, unique(id.practice)])
  
  ## Estimate overall outcomes
  outests = lapply(scens, calc_out) %>% 
    rbindlist(idcol = "scen")
  ### Clean output
  outcln = outests[
    , .(dataset.num
        , variable = fifelse(scen%in%c("all", "y3", "y4"), "Overall"
                             , str_extract(scen, "X[1-5]"))
        , level = fifelse(scen%in%c("all", "y3", "y4"), "NA"
                          , gsub("_", "", str_extract(scen, "_[ABC01]")))
        , year = fifelse(scen%in%c("y3", "y4"), str_extract(scen, "[34]$"), "NA")
        , satt, lower90, upper90)
  ]
  
  ## Estimate practice specific outputs
  outprac = lapply(prac_scens, calc_out) %>% 
    rbindlist(idcol = "id.practice") %>% 
    setcolorder(c("dataset.num", "id.practice", "satt", "lower90", "upper90"))
  
  ## Return overall and practice specific estimates for this draw
  return(list(overall = outcln, practice = outprac))
})

## Combine outputs across draws
overalls = lapply(all_out, `[[`, "overall") %>% 
  rbindlist()
practices = lapply(all_out, `[[`, "practice") %>% 
  rbindlist()

## Save outputs for draws
fwrite(overalls, "ACICdata/output_overall_"%p%paste(range(dataset.nums), collapse = "_")%p%".csv")
fwrite(practices, "ACICdata/output_practices_"%p%paste(range(dataset.nums), collapse = "_")%p%".csv")
