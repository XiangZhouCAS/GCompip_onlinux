if (!require(optparse,quietly = TRUE)) {
  install.packages("optparse")
  suppressMessages(library(optparse))
} else {
  suppressMessages(library(optparse))}
option_list <- list(
   make_option(c("--db_filepath","-d"),type = "character",default = F,
               help = "Set the download filepath."),
   make_option(c("--db_type","-t"),type = "character",default = "USCG",
               help = "Select the type of database to download."))
opt_parser <- OptionParser(
  usage = "usage: gcompip download",
  option_list = option_list,
  add_help_option = TRUE,
  prog=NULL ,
  description = "This page is to show how to download the databases if you need, database types include 'USCG', 'ter', 'hyd' and 'all'.")
opt <- parse_args(opt_parser)
db_filepath <- opt$db_filepath
db_type <- opt$db_type
setwd(db_filepath)
if(db_type == "USCG"){
    system(sprintf("wget -c https://github.com/XiangZhouCAS/GCompip/raw/refs/heads/main/database/Ribo_14.dmnd"))
    print("USCG database has been downloaded.")}
if(db_type == "hyd"){
    system(sprintf("wget -c https://github.com/XiangZhouCAS/GCompip/raw/refs/heads/main/database/hyd_id-name.script.txt"))
    system(sprintf("wget -c https://github.com/XiangZhouCAS/GCompip/raw/refs/heads/main/database/hyddb.all.dmnd"))
    print("Hydrogense database and hyd_id-name.script.txt have been downloaded.")}
if(db_type == "ter"){
    system(sprintf("wget -c https://github.com/XiangZhouCAS/GCompip/raw/refs/heads/main/database/ter.dmnd.gz"))
    print("Terminal database has been downloaded.")}
if(db_type == "all"){
    system(sprintf("wget -c https://github.com/XiangZhouCAS/GCompip/raw/refs/heads/main/database/Ribo_14.dmnd"))
    system(sprintf("wget -c https://github.com/XiangZhouCAS/GCompip/raw/refs/heads/main/database/hyd_id-name.script.txt"))
    system(sprintf("wget -c https://github.com/XiangZhouCAS/GCompip/raw/refs/heads/main/database/hyddb.all.dmnd"))
    system(sprintf("wget -c https://github.com/XiangZhouCAS/GCompip/raw/refs/heads/main/database/ter.dmnd.gz"))
    print("All databases have been downloaded.")}
