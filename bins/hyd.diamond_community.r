if (!require(optparse,quietly = TRUE)) {
  install.packages("optparse")
  suppressMessages(library(optparse))
} else {
  suppressMessages(library(optparse))}
option_list <- list(
  make_option(c("--input_reads","-i"),type = "character",default = F,
              help = "Please set the reads to include only the forward reads if the data is paired-end (PE)"),
  make_option(c("--result","-o"),type = "character",default = F,
              help = "Please specify the name for the result file."),
  make_option(c("--threads","-t"),type = "numeric",default = 1,
              help = "Set the number of CPU threads, with the default value being 1."),
  make_option(c("--hyd_db","-d"),type = "character",default = F,
              help = "Please set the directory of diamond database (e.g.'hyddb.all.dmnd')."),
  make_option(c("--script","-c"),type = "character",default = F,
	      help = "Please set the directory of hyd_id-name.script"),
  make_option(c("--USCG_db","-u"),type = "character",default = F,
              help = "Please specify the directory for universal single-copy genes database (e.g., 'Ribo_14.dmnd'). If you have already calculated the RPKM for these genes, you may instead specify the directory for the results (e.g., 'sample_name.UCSG.hits.txt') to skip this step."),
  make_option(c("--skip_fastp","-s"),action = "store_true", default = F,
              help = "If you have already filtered the reads, you can set this parameter to skip running fastp. The default is to run fastp."),
  make_option(c("--min_length","-m"),type = "numeric",default = 50,
	      help = "Set the minimum length required for filtering reads, the default is 50."),
  make_option(c("--run_seqkit","-k"),type = "character",default = "run",
	      help = "If you have already counted the total number of reads using seqkit, you can specify the directory of the seqkit results (e.g., 'sample_name.all.reads.txt') to skip running seqkit. By default, seqkit will be executed."),
  make_option(c("--keep_samples","-e"), action = "store_true", default = F,
              help = "By default, the temporary results will be deleted unless this parameter is setted."))
opt_parser = OptionParser(
  usage = "usage: gcompip custom hyd [options]",
  option_list = option_list,
  add_help_option = TRUE,
  prog=NULL ,
  description = "This Script is to calculate Hydrogenase community abandance.")
opt <- parse_args(opt_parser)
diamond_db <- normalizePath(opt$hyd_db)
singleM <- normalizePath(opt$USCG_db)
input_reads <- normalizePath(opt$input_reads)
if(opt$run_seqkit != "run"){
  run_seqkit <- normalizePath(opt$run_seqkit)
}else{run_seqkit <- opt$run_seqkit}
setwd(dirname(opt$result))
outpath <- basename(opt$result)
threads <- opt$threads
min_length <- opt$min_length
fastp_output <- paste0(outpath,".filtered.fq.gz")
script <- opt$script
system(sprintf("file %s > tmp.txt",singleM))
system(sprintf("sed 's/.*://' -i tmp.txt"))
system(sprintf("sed 's/ ASCII //g' -i tmp.txt"))
tmp <- read.table("tmp.txt",sep = "\t")
skip_fastp <- opt$skip_fastp
keep_samples <- opt$keep_samples
diamond_out <- paste0(outpath,".hits.txt")
singleM_out <- paste0(outpath,".USCG.hits.txt")
seqkit_out <- paste0(outpath,".all.reads.txt")
fastp <- sprintf("fastp -i %s -o %s --length_required %d -w %d",
                 input_reads,fastp_output,min_length,threads)
diamond <- sprintf("diamond blastx --db %s --query %s --out %s --threads %d --outfmt 6 slen stitle pident qcovhsp --max-target-seqs 1 --max-hsps 1 > /dev/null 2>&1",
                   diamond_db, fastp_output, diamond_out, threads)
diamond_singleM <- sprintf("diamond blastx --db %s --query %s --out %s --threads %d --outfmt 6 slen stitle qcovhsp bitscore --max-target-seqs 1 --max-hsps 1 > /dev/null 2>&1",
                           singleM, fastp_output, singleM_out, threads)
seqkit <- sprintf("seqkit stat %s > %s",
                  fastp_output,seqkit_out)
diamond2 <- sprintf("diamond blastx --db %s --query %s --out %s --threads %d --outfmt 6 slen stitle pident qcovhsp --max-target-seqs 1 --max-hsps 1 > /dev/null 2>&1",
		   diamond_db, input_reads, diamond_out, threads)
diamond_singleM2 <- sprintf("diamond blastx --db %s --query %s --out %s --threads %d --outfmt 6 slen stitle qcovhsp bitscore --max-target-seqs 1 --max-hsps 1 > /dev/null 2>&1",
			    singleM, input_reads, singleM_out, threads)
seqkit2 <- sprintf("seqkit stat %s > %s",
		   input_reads,seqkit_out)
if(skip_fastp == F){
	print("fastp is running.")
	system(fastp)
        print("fastp is completed.")}else{
		print("Not run the fastp because you set the directory of filtered reads to skip it.")}
if(skip_fastp == F){
	print("diamond is Running (target genes).")
	system(diamond)
        print("diamond is completed (target genes).")}else{
		print("diamond is Running (target genes).")
		system(diamond2)
		print("diamond is completed (target genes).")}
if(skip_fastp == F){
  if(tmp$V1 != "text"){
    print("diamond is Running (USCGs).")
    system(diamond_singleM)
    print("diamond is completed (USCGs).")
  }else{print("Not count the RPKM of singleM marker genes.")}
  }else{if(tmp$V1 != "text" ){
  print("diamond is Running (USCGs).")
  system(diamond_singleM2)
  print("diamond is completed (USCGs).")
    }else{print("Not count the RPKM of singleM marker genes.")}}
if(skip_fastp == F){
  if(run_seqkit == "run"){
    print("seqkit is Running.")
    system(seqkit)
    print("seqkit is completed.")}else{
      print("Not run the seqkit because you set the directory of seqkit result")}
  }else{if(run_seqkit == "run"){
  print("seqkit is Running.")
  system(seqkit2)
  print("seqkit is completed.")
}else{print("Not run the seqkit because you set the directory of seqkit result")}}
if (!require(magrittr)) {
  install.packages("magrittr")
  suppressMessages(library(magrittr))
} else {
  suppressMessages(library(magrittr))}
if(!require(dplyr)){
  install.packages("dplyr")
  suppressMessages(library(dplyr))}else{
  suppressMessages(library(dplyr))}
if(!require(tidyr)){
  install.packages("tidyr")
  suppressMessages(library(tidyr))}else{
  suppressMessages(library(tidyr))}
if(!require(data.table)){
  install.packages("data.table")
  suppressMessages(library(data.table))}else{
  suppressMessages(library(data.table))}
if(run_seqkit == "run"){
  p <- read.table(seqkit_out,header = T)
  }else{p <- read.table(run_seqkit,header = T)}
d1 <- fread(diamond_out,sep = "\t")
colnames(d1)[c(1:4)] <- c("slen","id","pident","qcovhsp")
d1 <- d1%>%
  group_by_all()%>%
  count()
d1$all_reads_nums <- p$num_seqs
d1$all_reads_nums <- gsub(",","",d1$all_reads_nums)
d1$all_reads_nums <- as.numeric(d1$all_reads_nums)
d1$tmp_RPKM <- d1$n*10^9/(d1$slen*d1$all_reads_nums*3)
d1$sample_name <- outpath
name <- read.table(script,sep = "\t",header = T)
d1$id <- gsub(" .*","",d1$id)
d1 <- left_join(d1,name,by = "id")
d12 <- as.data.frame(d1[grepl("NiFe4", d1$gene),])
d13 <- as.data.frame(d1[grepl("FeFe", d1$gene),])
d12 <- rbind(d12,d13)
d11 <- filter(d1,!id %in% d12$id)%>%
	as.data.frame()
d12 <- filter(d12,pident >= 60 & qcovhsp >= 80)
d11 <- filter(d11,pident >= 50 & qcovhsp >= 80)
d1 <- rbind(d11,d12)
d1 <- d1[,c(8,9,7)]
d1 <- d1%>%
  group_by(gene)%>%
  mutate(RPKM = sum(tmp_RPKM))
d1 <- unique(d1[,-3])
if(tmp$V1 != "text"){
  d2 <- fread(singleM_out,sep = "\t")
  }else{
  d2 <- fread(singleM,sep = "\t")
  }

colnames(d2)[c(1,2,3,4)] <- c("slen","sseq_id","qcovhsp","bitscore")
d2 <- filter(d2,qcovhsp >= 80 & bitscore >= 40)
d2 <- d2[,-c(3,4)]
d2 <- d2%>%
  group_by_all()%>%
  count()
d2$all_reads_nums <- p$num_seqs
d2$all_reads_nums <- gsub(",","",d2$all_reads_nums)
d2$all_reads_nums <- as.numeric(d2$all_reads_nums)
d2$RPKM <- d2$n*10^9/(d2$slen*d2$all_reads_nums*3)
d2$sample_name <- outpath
d2 <- d2[,c(6,2,1,3,4,5)]
d2$sseq_id <- gsub("-.*","",d2$sseq_id)
d2 <- d2%>%
  group_by(sseq_id)%>%
  mutate(T_RPKM = sum(RPKM))
d2 <- unique(d2[,c(1,2,7)])
s_out <- paste0(outpath,".USCG.hits.comb.txt")
write.table(d2,s_out,quote = FALSE,sep = "\t",row.names = F)
geo <- exp(mean(log(d2$T_RPKM)))
d1$GAM <- d1$RPKM/geo*100
d1$GAM <- ifelse(d1$GAM > 100,100,d1$GAM)
r_out <- paste0(outpath,".all.abd.txt")
colnames(d1)[4] <- "GAM"
write.table(d1,r_out,quote = FALSE,sep = "\t",row.names = F)
system(sprintf("mkdir %s",outpath))
if(skip_fastp == F){
	system(sprintf("mv %s %s",fastp_output,outpath))}else{
		print("……")}
system(sprintf("mv %s %s",diamond_out,outpath))
if(tmp$V1 != "text"){
  system(sprintf("mv %s %s",singleM_out,outpath))}else{
	  print("There is no need to move singleM result to outpath")}
if(run_seqkit == "run"){
  system(sprintf("mv %s %s",seqkit_out,outpath))
  }else{
    print("There is no need to move seqkit result to outpath")}
system(sprintf("cat %s >> all.abd.txt",r_out))
tmp <- read.table("all.abd.txt",sep = "\t",header = T)
tmp <- filter(tmp,sample_name != "sample_name")
write.table(tmp,"all.abd.txt",quote = F,sep = "\t",row.names = F)
system(sprintf("mv %s %s",s_out,outpath))
system(sprintf("mv %s %s",r_out,outpath))
if(skip_fastp == F){
	system("rm \"fastp.html\" \"fastp.json\" \"tmp.txt\"")}else{
		system("rm  \"tmp.txt\"")}
wd1 <- paste0(getwd(),"/",outpath)
wd2 <- paste0(getwd())
if(keep_samples == T){
  print("samples are keeped.")}else{if(file.exists("samples") == T){
    setwd(wd1)
    system(sprintf("mv %s ../samples",r_out))
    setwd(wd2)
    system(sprintf("rm %s -rf",outpath))}else{
      system(sprintf("mkdir samples"))
      setwd(wd1)
      system(sprintf("mv %s ../samples",r_out))
      setwd(wd2)
      system(sprintf("rm %s -rf",outpath))}}
tmp <- read.table("all.abd.txt",header = T,sep = "\t")
rpkm <- pivot_wider(tmp,id_cols = gene,
              values_from = "RPKM",
              names_from = "sample_name")
rpkm <- as.data.frame(rpkm)
rpkm[is.na(rpkm)] <- 0
write.table(rpkm,"rpkm.abd.txt",
                    sep = "\t",
                    quote = F,row.names = F)
com <- pivot_wider(tmp,id_cols = gene,
              values_from = "GAM",
              names_from = "sample_name")
com <- as.data.frame(com)
com[is.na(com)] <- 0
com <- com[!duplicated(com[,c(1,2)]),]
write.table(com,"GAM.abd.txt",
                    sep = "\t",
                    quote = F,row.names = F)
if(!require(ggplot2)){
  install.packages("ggplot2")
  suppressMessages(library(ggplt2))
}else{suppressMessages(library(ggplot2))}
lv <- colnames(com)[-1]
lv <- as.factor(lv)
samplept <- com%>%
        pivot_longer(cols = !gene,values_to = "val",names_to = "sample")
samplept$sample <- factor(samplept$sample,levels = lv)
samplept$val <- as.numeric(samplept$val)
hp <- ggplot(samplept,aes(sample,gene))+
  geom_tile(aes(fill = val),color = "grey50")+
  scale_fill_gradientn(colours = c("white","#FF471A",
                                   "#EF1A22","#924EA4", "#262479"),
                       name = "GAM (%)")+
  theme(panel.background = element_blank(),
        axis.title = element_blank(),
        axis.text.x = element_text(size = 15,color = "black",angle = 30,hjust = 1),
        legend.position = "bottom",
        legend.title = element_text(vjust = 0.75),
        axis.text.y = element_text(size = 12,color = "black"),
        axis.ticks.y = element_line(color = "black"))
w <- length(unique(samplept$sample))
h <- length(unique(samplept$gene))
ggsave(plot = hp,"Heatmap_GAM.pdf",width = w,height = h/4)
print("All Completed!")
