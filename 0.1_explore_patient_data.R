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
       , )

## Handy Functions ----
`%p%` = paste0

month_diff = function(d1, d2){
  12*(year(d2) - year(d1)) + (month(d2) - month(d1))
}

make_title = function(x){ str_to_title(gsub("_", " ", x))}

##################-
# Read in Practice level data ----
##################-
run_nos = c(1201:1203)
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
]

##################-
# Basic Exploration ----
##################-

## Patient level vars ----
### V1 and V4 are continuous
### V2, and V5 are categorical
### V3 is binary
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
### Right skewed mean 0 variance 1. 
