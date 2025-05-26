if (!require(optparse,quietly = TRUE)) {
  install.packages("optparse")
  library(optparse)
} else {
  library(optparse)}
if (!require(magrittr)) {
  install.packages("magrittr")
  library(magrittr)
} else {
  library(magrittr)}
if (!require(dplyr)) {
  install.packages("dplyr")
  library(dplyr)
} else {
  library(dplyr)}
if (!require(tidyr)) {
  install.packages("tidyr")
  library(tidyr)
} else {
  library(tidyr)}
if (!require(data.table)) {
  install.packages("data.table")
  library(data.table)
} else {
  library(data.table)}
if(!require(ggplot2)){
  install.packages("ggplot2")
  library(ggplt2)
}else{library(ggplot2)}
