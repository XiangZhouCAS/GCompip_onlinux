if(!require(optparse,quietly = TRUE)){
  install.packages("optparse")
  library(optparse)
}else{
  library(optparse)}
option_list <- list(
  make_option(c("--help","-h"),type = "character",default = F,
              help = "Show this help message and exit"))
opt_parser <- OptionParser(
  usage = "usage: gcompip custom diy    To calculate gene abundance in microbial community (GAM,%) by customed database.\n       gcompip custom ter    To calculate terminal enzyme gene abandance in microbial community (GAM,%).\n       gcompip custom hyd    To calculate single copy Hydrogenase abandance in microbial community.",
  add_help_option = TRUE,
  prog=NULL ,
  description = "This page is to show how to run the program.")
opt <- parse_args(opt_parser)
