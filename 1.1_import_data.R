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
p_load(data.table, magrittr, stringr, ggplot2)

## Handy Functions ----
`%p%` = paste0

month_diff = function(d1, d2){
  12*(year(d2) - year(d1)) + (month(d2) - month(d1))
}

make_title = function(x){ str_to_title(gsub("_", " ", x))}

##################-
# Read in Practice level data ----
##################-

practice_year_data = lapply(list.files("data/track2/practice_year", full.names = T), fread) %>% 
  rbindlist(idcol = "samp")
