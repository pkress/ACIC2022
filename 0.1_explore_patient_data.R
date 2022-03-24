##################-
### Author: Peter Kress
### Date: 2022/02/16
### Purpose: ACIC Data Processing
##################-

##################-
# Initialize Workspace ----
##################-
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
# Read in Practice level data ----
##################-
run_nos = c(1201)
read_files = function(folder, nos){
  lapply(list.files(folder, full.names = T, pattern = paste(nos%p%".csv",collapse = "|")), fread) %>% 
    rbindlist(idcol = "samp")
}
data = "data/track1b/"%p%c("patient", "patient_year", "practice", "practice_year") %>% 
  setNames(.,c("patient", "patient_year", "practice", "practice_year")) %>% 
  lapply(.,read_files, nos = run_nos)

merged_data = data$patient[
  data$patient_year
  , on = c("id.patient", "samp")
][
  data$practice
  , on = c("id.practice", "samp")
][
  data$practice_year[, .(id.practice, year, Z, post, samp)]
  , on = c("id.practice", "year", "samp")
]

##################-
# Basic Exploration ----
##################-

## Patient level vars ----
merged_data[
  , .N
  , c("samp", "id.patient", "v"%p%1:5)
][, .N, c("samp", "id.patient")][,.N,N]
### Nonvarying by time

### V1 and V4 are continuous ----
merged_data %>% 
  ggplot()+
  geom_histogram(aes(x = V1))+
  facet_wrap(~samp)+
  scale_x_log10()
merged_data %>% 
  ggplot()+
  geom_density(aes(x = V1, fill = factor(samp)), alpha = 0.2)+
  scale_x_log10()
merged_data[, .(sd(log10(V1)), mean(log10(V1))), .(samp)]
### Lognormal around 10, sd (log10) of 0.16

merged_data %>% 
  ggplot()+
  geom_histogram(aes(x = V4))+
  facet_wrap(~samp)
merged_data %>% 
  ggplot()+
  geom_density(aes(x = V4, fill = factor(samp)), alpha = 0.2)
merged_data[, .(sd(V4), mean(V4)), .(samp)]
### Left skewed mean 0 variance 1. 


### V2, and V5 are categorical ----

v2_sum = merged_data[
  , .N
  ,, .(samp, V2)
  ][
  , share:=N/sum(N)
  , samp
  ] 
v2_wide = v2_sum %>% 
  dcast(V2~samp, value.var = "share") %>% 
  print()
v2_sum %>% 
  ggplot()+
  geom_col(aes(x = factor(V2), y = share, fill = factor(samp)), position = "dodge")
### Very similar distributions. rounded lognormal it seems

v5_sum = merged_data[
  , .N
  ,, .(samp, V5)
][
  , share:=N/sum(N)
  , samp
] %>% print()
v5_wide = v5_sum %>% 
  dcast(V5~samp, value.var = "share") %>% 
  print()
v5_sum %>% 
  ggplot()+
  geom_col(aes(x = factor(V5), y = share, fill = factor(samp)), position = "dodge")
### Very similar distributions. Approx 69, 20.5, 10.5. 

### V3 is binary ----
merged_data[
  , .N
  ,,.(samp, V3)
][
  , share:=N/sum(N)
  , samp
] %>% print()
### 47% = 0 and 53% = 1

## Practice covariates ----
### X1, 3, 5 ----
merged_data[
  , .(N = uniqueN(id.practice))
  ,,.(samp, X1)
][
  , share:=N/sum(N)
  , samp
] %>% print()
### about 50/50

merged_data[
  , .(N = uniqueN(id.practice))
  ,,.(samp, X3)
][
  , share:=N/sum(N)
  , samp
] %>% print()
### differs by samp, 44-46% = 0

merged_data[
  , .(N = uniqueN(id.practice))
  ,,.(samp, X5)
][
  , share:=N/sum(N)
  , samp
] %>% print()
### differs by samp, 50-55% = 0

### X2, X4 ----
outsum = lapply("X"%p%c(2,4), \(x){
  sumdata =merged_data[
    , .(N = uniqueN(id.practice))
    ,, .(samp, get(x))
  ][
    , share:=N/sum(N)
    , samp
  ] %>% 
    setnames(c("samp", x, "N", "share")) %>% 
    print()
  plot = sumdata %>% 
    ggplot()+
    geom_col(aes(x = factor(get(x)), y = share, fill = factor(samp)), position = "dodge")+
    labs(x = x)
  return(list(sumdata = sumdata, plot = plot))
})
outsum
### Differs by samp. 
### X2 falling in shares from A to C. A between 30 and 50, C between 23 and 25
### X4 mostly in A (~77%), B around 17% and C around 6%

### X6:9 ----
merged_data[
  ,.N
  , .(samp, X6, id.practice)
] %>% 
  ggplot()+
  geom_density(aes(x = X6, fill = factor(samp)), alpha = 0.2)
### right skewed and inconsistent across practices.
### Generally 0 to 60, mean near 30 with sd ~ 8

merged_data[
  ,.N
  , .(samp, X7, id.practice)
] %>% 
  ggplot()+
  geom_density(aes(x = X7, fill = factor(samp)), alpha = 0.2)
### right skewed and inconsistent across practices.
### Generally 0 to 60, mean near 8.5 with sd ~ 5


merged_data[
  ,.N
  , .(samp, X8, id.practice)
]  %>% 
  ggplot()+
  geom_density(aes(x = X8, fill = factor(samp)), alpha = 0.2)
merged_data[, .(sd(X8), mean(X8)), .(samp)]
### very slightly right skewed skewed and inconsistent across practices.
### Generally 0.1 to 0.9, mean near 0.4 with sd ~ 1.4

merged_data[
  ,.N
  , .(samp, X9, id.practice)
]  %>% 
  ggplot()+
  geom_density(aes(x = X9, fill = factor(samp)), alpha = 0.2)
merged_data[, .(sd(X9), mean(X9)), .(samp)]
### very slightly left skewed skewed and inconsistent across practices.
### Generally -100 to 150, mean near ~20 with sd ~ 40

## outcome and treatment data  ----

merged_data %>% 
  ggplot()+
  geom_density(aes(x = Y, fill = factor(samp)), alpha = 0.2)
merged_data %>% 
  ggplot()+
  geom_density(aes(x = Y, fill = factor(samp)), alpha = 0.2)+
  scale_x_log10()
merged_data[, .(sd(Y), mean(Y)), .(samp)]
merged_data[, .(sd(log10(Y), na.rm = T), mean(log10(Y), na.rm =T)), .(samp)]
### Basically lognormal distributed with mean ~10^2.7 and sd(log10) of 0.49. 
### Slightly skewed right (faster ascent left)

merged_data %>% 
  ggplot()+
  geom_boxplot(aes(x = factor(year), y = Y, fill = factor(Z)))+
  facet_wrap(~samp)+
  scale_y_log10()

### Every group has growth in Y overall. 

##################-
# Clean Data ----
##################-
## Encode types ----
catvars = c("V2", "V5", "X2", "X4")
contvars = c("V1", "V4", "X"%p%6:9)
binvars = c("V3", "X1", "X3", "X5")
merged_data[
  , c(binvars, "Z", "post"):=lapply(.SD, as.logical)
  , .SDcols = c(binvars, "Z", "post")
  ][
  , c(catvars):=lapply(.SD, as.factor)
  , .SDcols = c(catvars)
  ]
## Rename Vars ----
new_names = c("V1" = "pat_cont_lognorm"
              , "V4" = "pat_cont_left_skew_0"
              , "V2" = "pat_cat_lognorm_pois"
              , "V5" = "pat_cat_abc_exp"
              , "V3" = "pat_bin_4753" #0/1 ratios
              , "X1" = "prac_bin_5050"
              , "X3" = "prac_bin_4555"
              , "X5" = "prac_bin_5248"
              , "X2" = "prac_cat_abc_453525"
              , "X4" = "prac_cat_abc_751510"
              , "X6" = "prac_cont_sl_right_skew_30"
              , "X7" = "prac_cont_right_skew_10"
              , "X8" = "prac_cont_norm_narrow_0.5"
              , "X9" = "prac_cont_norm_wide_20"
              )
setnames(merged_data, old = names(new_names), new = new_names, skip_absent = T)

##################-
# Quick regressions ----
##################-
## Linear all ----
preds = lapply(merged_data[, unique(samp)][1], \(x){
  feols(as.formula("Y~Z+"%p%paste(new_names[c(catvars, binvars, contvars)], collapse = "+")%p%"+i(year)")
        , data = merged_data[samp==x], cluster = "id.patient")
})

sigsigns = lapply(preds, \(x){
  out = data.table(sign(x$coeftable$Estimate) * fifelse(x$coeftable$`Pr(>|t|)`<=0.01, 1, 0)) %>% 
    transpose()
  names(out) = rownames(x$coeftable)
  out
}) %>% 
  rbindlist(idcol = "samp", fill = T)
etable(preds)
coefplot(preds, keep = c("year"))
### Year seems very linear with about 125 per year. 
coefplot(preds, keep = c("prac"))
### Signs and sig vary across model, esp for the continuous vars (not subgroup) and the X4 (cat) and X5 (bin) vars
### Somewhat more consistant for some cat/bin vars 
### "X4/prac_cat_abc_751510" (neg for both), "X1prac_bin_5050" (pos)
### "X3" = "prac_bin_4555" (pos) are farily consistent (except for one or two. )
coefplot(`preds`, keep = c("pat"))
### Mostly imprecise and around zero, except pat cont left skew 


## Log linear ----
log_preds = lapply(merged_data[, unique(samp)], \(x){
  feols(as.formula("log(Y)~"%p%paste(new_names[c(catvars, binvars, contvars)], collapse = "+")%p%"+i(year)")
        , data = merged_data[samp==x], cluster = "id.patient")
})

log_sigsigns = lapply(log_preds, \(x){
  out = data.table(sign(x$coeftable$Estimate) * fifelse(x$coeftable$`Pr(>|t|)`<=0.01, 1, 0)) %>% 
    transpose()
  names(out) = rownames(x$coeftable)
  out
}) %>% 
  rbindlist(idcol = "samp", fill = T)
etable(log_preds)

##################-
# Quick ests ----
##################-

## DiD ----
diff_in_diff = lapply(merged_data[, unique(samp)][1], \(x){
    feols(as.formula("Y~i(year, ref = 2):Z+"%p%paste(grep("^prac_", new_names, value = T), collapse = "+")%p%"|id.patient + year"), data = merged_data[samp==x]
          , cluster = "id.practice")
  })
diff_in_diff = lapply(merged_data[, unique(samp)][1], \(x){
  feols(as.formula("log(Y)~i(year, ref = 2):Z|id.patient + year")
        , data = merged_data[samp==x]
        , cluster = "id.practice")
})
coefplot(diff_in_diff)


## linear ----

eff = lapply(merged_data[, unique(samp)][2], \(x){
  feols(as.formula("Y~i(Z, post)+"%p%paste(new_names[c(catvars, binvars, contvars)], collapse = "+")%p%"+i(year, ref = 1)")
        , data = merged_data[samp==x], cluster = "id.patient")
})
etable(eff)


##################-
# Matching Estimates ----
##################-

## CEM ----
### Prolly not best
# Function to perform cem estimations
perform_cem = function(match_vars, in_data, out_var){
  cem_data = in_data[
    , c(.(out_var = get(out_var)
          , eoc_id = eoc_id)
        , lapply(.SD, as.factor)
    )
    , .SDcols = c("treatment", all_vars)
  ]
  
  imb = imbalance(group = cem_data[, treatment], cem_data[, ..match_vars], )
  cem_out = cem(treatment = "treatment",  data = cem_data, drop = c(setdiff(all_vars, match_vars), "out_var", "eoc_id")
                , eval.imbalance = T)
  cem_data = cem_data[
    , `:=`(weight=cem_out$w
           , strata = cem_out$mstrata)
  ][
    weight>0
  ]
  reg_out = lm(out_var ~ treatment, data = cem_data, weights = weight)
  reg_cln = coeftest(reg_out, vcov. = vcovHC(reg_out, type = "HC1"))
  out = list(imb = imb, cem = cem_out, reg = reg_out, reg_cln = reg_cln, data = cem_data)
  return(out)
}

cem_data = merged_data[samp==1][sample(1:.N, round(.N/10))]
cem_data[
  , names(cem_data)[sapply(cem_data, class)=="logical"]:=
    lapply(.SD, as.factor)
  , .SDcols = names(cem_data)[sapply(cem_data, class)=="logical"]
]
imb = imbalance( group = cem_data[, Z]
                 , merged_data[samp==1, .SD, .SDcols = c(paste(new_names[c(catvars, binvars, contvars)]))], )
cem_out = cem(treatment = "Z"
              ,  data = data.frame(cem_data[, .SD, .SDcols = c("Z", new_names[c(catvars, binvars)])])
              , drop = c("Y", "id.patient", "id.practice", "id.samp", "year", "post")
              , eval.imbalance = T)
cem_data = cem_data[
  , `:=`(weight=cem_out$w
         , strata = cem_out$mstrata)
][
  weight>0
]
reg_out = lm(out_var ~ treatment, data = cem_data, weights = weight)
reg_cln = coeftest(reg_out, vcov. = vcovHC(reg_out, type = "HC1"))

m_out = matchit(as.formula("Z~"%p%paste(new_names[c(catvars, binvars, contvars)], collapse = "+")%p%"+factor(year)")
                , data = merged_data[samp==1], method = "cem")


##################-
# BART (Joint) + TMLE ----
##################-
### Use BART to jointly model treatment and response. 
### For treatment, use cross validation to pick optimal hyperparameter (struggles with binary)
### then apply TMLE treatment. 
## Step 1: Initial Estimates ----
## Use wbart for initial, then lbart or pbart for PScore
?BART::abart()


##################-
# SL + TMLE ----
##################-



