if (!require(optparse,quietly = TRUE)) {
  install.packages("optparse")
  suppressMessages(library(optparse))
} else {
  suppressMessages(library(optparse))}
option_list <- list(
  make_option(c("--input_reads","-i"),type = "character",default = F,
                                help = "Please set the directory of reads"),
  make_option(c("--result","-o"),type = "character",default = F,
              help = "Please set the result file name"),
  make_option(c("--threads","-t"),type = "numeric",default = 1,
              help = "Setting the threads of CPU,default is 1"),
  make_option(c("--USCG_db","-u"),type = "character",default = F,
              help = "Please set the directory of universal single copy genes (USCGs) database"),
  make_option(c("--skip_fastp","-s"),action = "store_true", default = F,
              help = "If you have already filtered the reads, you can set this parameter to skip running fastp. The default is to run fastp."),
  make_option(c("--min_length","-m"),type = "numeric",default = 100,
	      help = "Set the minimum length required for filtering reads, the default is 100, but it is recommended to set this parameter to 140 if hydrogenases or hydrogen metabolism terminal enzymes are to be calculated."),
  make_option(c("--run_seqkit","-k"),type = "character",default = "run",
	      help = "If you have already counted the total number of reads using seqkit, you can specify the directory of the seqkit results (e.g., 'sample_name.all.reads.txt') to skip running seqkit. By default, seqkit will be executed."),
  make_option(c("--keep_samples","-e"), action = "store_true", default = F,
              help = "By default, the temporary results will be deleted unless this parameter is setted."))
opt_parser = OptionParser(
  usage = "usage: comts rpkm ribo [options]",
  option_list = option_list,
  add_help_option = TRUE,
  prog=NULL ,
  description = "This Script is to calculate rpkm abandance of singleM's universal single copy ribosomal genes.")
opt <- parse_args(opt_parser)
input_reads <- opt$input_reads
skip_fastp <- opt$skip_fastp
min_length <- opt$min_length
input_geneset <- opt$input_geneset
res <- opt$result
singleM <- opt$USCG_db
threads <- opt$threads
run_seqkit <- opt$run_seqkit
keep_samples <- opt$keep_samples
singleM_out <- paste0(basename(res),".USCG.txt")
fastp_output <- paste0(basename(res),".filtered.fq.gz")
seqkit_out <- paste0(res,".all.reads.txt")
seqkit <- sprintf("seqkit stat %s > %s",
                  fastp_output,seqkit_out)
seqkit2 <- sprintf("seqkit stat %s > %s",
                  input_reads,seqkit_out)
fastp <- sprintf("fastp -i %s -o %s --length_required %d -w %d",
                 input_reads,fastp_output,min_length,threads)
diamond_singleM <- sprintf("diamond blastx --db %s --query %s --out %s --threads %d --outfmt 6 slen stitle qcovhsp bitscore --max-target-seqs 1 --max-hsps 1 > /dev/null 2>&1",
                           singleM, fastp_output, singleM_out, threads)
diamond_singleM2 <- sprintf("diamond blastx --db %s --query %s --out %s --threads %d --outfmt 6 slen stitle qcovhsp bitscore --max-target-seqs 1 --max-hsps 1 > /dev/null 2>&1",
                            singleM, input_reads, singleM_out, threads)
if(skip_fastp == F){
  print("fastp is Running.")
  system(fastp)
  print("fastp is completed.")}else{
    print("Not run fastp because you set the directory of filtered reads to skip it.")}
if(skip_fastp == F){
  print("diamond is Running (USCGs).")
  system(diamond_singleM)
  print("diamond is completed (USCGs).")}else{
    print("diamond is Running (USCGs).")
    system(diamond_singleM2)
    print("diamond is completed (USCGs).")
  }
if(skip_fastp == T){
   if(run_seqkit != "run"){
   print("Not run seqkit because you set the directory of seqkit result to skip it.")}else{
	   print("seqkit is Running.")
	   system(seqkit2)
           print("seqkit is completed.")}}else{
	   print("seqkit is Running.")
	   system(seqkit)
           print("seqkit is completed.")}
if(!require(magrittr)){
  install.packages("magrittr")
  suppressMessages(library(magrittr))
}else{
  suppressMessages(library(magrittr))}

if(!require(dplyr)){
  install.packages("dplyr")
  suppressMessages(library(dplyr))
}else{
  suppressMessages(library(dplyr))}

d2 <- read.table(singleM_out,sep = "\t")
colnames(d2)[c(1,2,3,4)] <- c("slen","sseq_id","qcovhsp","bitscore")
d2 <- filter(d2,qcovhsp > 80 & bitscore > 40)
d2 <- d2[,-c(3,4)]
d2 <- d2%>%
  group_by_all()%>%
  count()
s_out <- paste0(basename(res),".RUSCG.txt")
if(run_seqkit == "run"){
  p <- read.table(seqkit_out,header = T)
  }else{p <- read.table(run_seqkit,header = T)}
d2$all_reads_nums <- p$num_seqs
d2$all_reads_nums <- gsub(",","",d2$all_reads_nums)
d2$all_reads_nums <- as.numeric(d2$all_reads_nums)
d2$RPKM <- d2$n*10^9/(d2$slen*d2$all_reads_nums*3)
d2$sample_name <- res
d2 <- d2[,c(6,2,1,3,4,5)]
d2$sseq_id <- gsub("-.*","",d2$sseq_id)
d2 <- d2%>%
  group_by(sseq_id)%>%
  mutate(T_RPKM = sum(RPKM))
d2 <- unique(d2[,c(1,2,7)])
write.table(d2,s_out,quote = F,sep = "\t",row.names = F)
system(sprintf("cat %s >> RUSCG.txt",s_out))
total <- read.table("RUSCG.txt",header = T,sep = "\t")
total <- filter(total,T_RPKM != "T_RPKM")
write.table(total,"RUSCG.txt",quote = F,sep = "\t",row.names = F)

system(sprintf("mkdir %s",res))
if(run_seqkit == "run"){
  system(sprintf("mv %s %s",seqkit_out,res))
  system(sprintf("mv %s %s",s_out,res))
  system(sprintf("mv %s %s",singleM_out,res))}else{
    system(sprintf("mv %s %s",s_out,res))
    system(sprintf("mv %s %s",singleM_out,res))}
if(skip_fastp == F){
  system(sprintf("mv %s %s",fastp_output,res))}else{
    print("There is no fastp result.")}
wd1 <- paste0(getwd(),"/",res)
wd2 <- paste0(getwd())
if(keep_samples == T){
  print("samples are keeped")}else{if(file.exists("samples") == T){
    setwd(wd1)
    system(sprintf("mv %s ../samples",s_out))
    setwd(wd2)
    system(sprintf("rm %s -rf",res))}else{
      system(sprintf("mkdir samples"))
      setwd(wd1)
      system(sprintf("mv %s ../samples",s_out))
      setwd(wd2)
      system(sprintf("rm %s -rf",res))}}
print("All Completed!")
