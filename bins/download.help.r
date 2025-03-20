if(!require(optparse,quietly = TRUE)){
  install.packages("optparse")
  library(optparse)
}else{
  library(optparse)}
option_list <- list(
  make_option(c("--help","-h"),type = "character",default = F,
              help = "Show this help message and exit"))
opt_parser = OptionParser(
  usage = "gcompip download",
  add_help_option = TRUE,
  prog=NULL ,
  description = "This page is to show how to download the databases if you need, database types include 'USCG', 'ter', 'hyd' and 'all'.")
opt <- parse_args(opt_parser)
