if(!require(optparse,quietly = TRUE)){
  install.packages("optparse")
  library(optparse)
}else{
  library(optparse)}
option_list <- list(
  make_option(c("--help","-h"),type = "character",default = F,
              help = "Show this help message and exit"))
opt_parser = OptionParser(
  usage = "gcompip geneset uscg     To calculate RPKM abandance of 14 universal single copy ribosomal genes (USCGs).\n       gcompip geneset res      To convert RPKM to GAM of single copy gene through GeneSet.",
  add_help_option = TRUE,
  prog=NULL ,
  description = "This page is to show how to run the comts.")
opt <- parse_args(opt_parser)
