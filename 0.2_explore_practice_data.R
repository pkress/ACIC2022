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
       , fixest)

## Handy Functions ----
`%p%` = paste0

month_diff = function(d1, d2){
  12*(year(d2) - year(d1)) + (month(d2) - month(d1))
}

make_title = function(x){ str_to_title(gsub("_", " ", x))}

##################-
# Read in Practice level data ----
##################-
run_nos = c(1201:1218)
read_files = function(folder, nos){
  lapply(list.files(folder, full.names = T, pattern = paste(nos%p%".csv",collapse = "|")), fread) %>% 
    rbindlist(idcol = "samp")
}
data = "data/track1b/"%p%c("practice", "practice_year") %>% 
  setNames(.,c("practice", "practice_year")) %>% 
  lapply(.,read_files, nos = run_nos)

merged_data = data$practice[
  data$practice_year
  , on = c("id.practice", "samp")
]

##################-
# Clean Data ----
##################-

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

##################-
# Quick regressions ----
##################-
### TRy CEM To match similar facilities, and then check how treatment affects spending controlling for patient chars
## CEM ----
# Function to perform cem estimations
perform_cem = function(match_vars, in_data, out_var){
  
  cat_names = names(new_names)[!grepl("pat", new_names)&grepl("cat|bin", new_names)]
  cont_names = "X"%p%6:9
  pat_names = names(new_names)[grepl("pat", new_names)]
  
  cem_data = copy(in_data)[
    , c("Z", cat_names):=lapply(.SD, as.factor)
    , .SDcols = c("Z", cat_names)
  ]
  
  
  imb = imbalance(group = cem_data[year==1, Z]
                  , cem_data[year==1, .SD, .SDcols = c(cat_names, cont_names)])
  cem_out = cem(treatment = "Z"
                ,  data = cem_data[year==1, .SD, .SDcols = c("Z", cont_names)]
                , eval.imbalance = T, verbose = T, keep.all = T)
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
