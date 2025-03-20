if(!require(optparse,quietly = TRUE)){
  install.packages("optparse")
  suppressMessages(library(optparse))
}else{
  suppressMessages(library(optparse))}
option_list <- list(
  make_option(c("--input_geneset","-i"),type = "character",default = F,
              help = "Please set the directory of geneset RPKM abundance file"),
  make_option(c("--input_ribo_rpkm","-r"),type = "character",default = F,
	      help = "Please set the directory of USCGs' RPKM file e.g.,ribo_rpkm.txt."),
  make_option(c("--output","-o"),type = "character",default = F,
	      help = "Please set the directory of output"))
opt_parser = OptionParser(
  usage = "usage: comts rpkm com [options]",
  option_list = option_list,
  add_help_option = TRUE,
  prog=NULL ,
  description = "This Script is to convert rpkm to gene abundance in microbial community (GAM,%).")
opt <- parse_args(opt_parser)
geneset <- opt$input_geneset
rpkm <- opt$input_ribo_rpkm
output <- opt$output
all_data <- read.table(rpkm,header = T,sep = "\t")
geneset <- read.table(geneset,header = T,
                  sep = "\t",check.names = F)
if(!require(dplyr,quietly = TRUE)){
   install.packages("dplyr")
   suppressMessages(library(dplyr))}else{
   suppressMessages(library(dplyr))}
if(!require(tidyr,quietly = TRUE)){
   install.packages("tidyr")
   suppressMessages(library(tidyr))}else{
   suppressMessages(library(tidyr))}
all_data <- all_data%>%
  group_by(sample_name)%>%
  mutate(geomean = exp(mean(log(T_RPKM))))
all_data <- unique(all_data[,c(1,4)])
all_data <- all_data[match(colnames(geneset[,-1]),
		           all_data$sample_name),]
for(i in 1:nrow(all_data)){
  ribo_rpkm <- as.numeric(all_data[i,2])
  geneset[,i+1] <- geneset[,i+1]*100/ribo_rpkm}
df <- data.frame(lapply(geneset[,-1],
                        function(x)ifelse(x > 100, 100, x)),check.names = F)
res <- basename(output)
colnames(geneset)[1] <- "GeneID"
df$GeneID <- geneset$GeneID
df$GeneID <- geneset$GeneID
df <- df[,c(ncol(df),1:ncol(df)-1)]
write.table(df,res,sep = "\t",quote = F,row.names = F)
