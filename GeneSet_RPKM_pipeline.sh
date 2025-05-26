#Step1: Remove host’s sequence.
bowtie2 -p threads -x hostgenome -1 sample_name.raw.1.fq.gz -2 sample_name.raw.2.fq.gz -S sample_name.sam --un-conc sample_name.raw.clean.fq

#Step2: Assemble contigs from sequences. (rename_fasta.sh at the end of the document)
megahit -1 sample_name.clean.1.fq.gz -2 sample_name.clean.2.fq.gz --min-contig-len length -t threads -o sample_name
sh rename_fasta.sh final.contigs.fa sample_name sample_name.contig.ok.fa


#Step3: Predicting genes from contigs.
prodigal -i sample_name.contig.ok.fa -f gff -o sample_name.gff -p meta -q -d sample_name.temp.orf.ffn -a sample_name.temp.orf.faa
perl gene-filter.pl sample_name.temp.orf.ffn  sample_name.temp.orf.faa sample_name.orf.ffn sample_name.orf.faa

#Step4: Generate non-redundant representative sequences.
cd-hit-est -i sample_name.orf.ffn -o sample_name.geneSet.ffn -n 9 -c 0.95 -G 0 -M 0 -d 0 -aS 0.9 -r 0 -T 40
Cat *geneSet.ffn >> DNA.geneSet.tmp.ffn
cd-hit-est -i DNA.geneSet.tmp.ffn -o DNA.geneSet.ffn -n 9 -c 0.95 -G 0 -M 0 -d 0 -aS 0.9 -r 0 -T 40

#Step5: Clean sequences aligning with gene set.
bwa index DNA.geneSet.ffn
bwa mem -t #threads -M -R '@RG\tID:$i\tSM:$i\tLB:$i\tPL:Illumina\tPI:$i' DNA.geneSet.ffn sample_name.clean.1.fq.gz sample_name.clean.2.fq.gz > sample_name.sam

Step6: Calculate genes’ RPKM in gene set. (rpkm.r at the end of the document)
# Counting total read counts (TRC) and length of genes.
perl read_counts.pl DNA.geneSet.ffn sample_name.sam sample_name.sam.abd
sh ./list.sh
perl abd.merge.pl list TPM.txt reads_counts.txt
# Calculate genes’ RPKM in gene set.
cut -f 2 sample.sam.abd (only need one sample to generate it) > genelength.txt
sed -i '1d' genelength.txt
Rscript rpkm.r -i reads_counts.txt -g genelenth.txt -o RPKM.txt

rename_fasta.sh:
#!/bin/bash
if [ "$#" -ne 3 ]; then
    echo "Usage: $0 <input_fasta_file> <prefix> <output_file>"
    exit 1
fi
input_file="$1"
prefix="$2"
output_file="$3"
awk -v prefix="$prefix" '
BEGIN { counter = 1 }
/^>/ {print ">" prefix "_" counter;counter++}
!/^>/ {print}
' "$input_file" > "$output_file"
echo "Renamed sequences have been written to $output_file"

gene-filter.pl
echo "Renamed sequences have been written to $output_file"
die "perl input.ffn input.faa  output.ffn output.faa\n" if(@ARGV!=4);
open IA, "$ARGV[0]" or die "can not open file: $ARGV[0]\n";
open IB, "$ARGV[1]" or die "can not open file: $ARGV[1]\n";
open OA, ">$ARGV[2]" or die "can not open file: $ARGV[2]\n";
open OB, ">$ARGV[3]" or die "can not open file: $ARGV[3]\n";
$/=">";<IA>;<IB>;
while($line=<IA>){
        chomp $line;
        my @ele=split /\n/,$line;
        my $cds="";
        my @inf=split /\s+/,$ele[0];
        for(my $i=1;$i<=$#ele;$i++){
                $cds.=$ele[$i];
        }
        $line=<IB>;
        chomp $line;
        @ele=split /\n/,$line;
        my $aa="";
        for(my $i=1;$i<=$#ele;$i++){
                $aa.=$ele[$i];
        }
        if(length($cds)>100){
                $aa=~s/\*//g;
                print OA ">$inf[0]\n$cds\n";
                print OB ">$inf[0]\n$aa\n";
        }
}
close IA;
close OA;
close IB;
close OB;

read_counts.pl
use FindBin qw($Bin);
die "perl uniqGene-ProfilebySoap.pl gene.fa sam out\n" if(@ARGV!=3);
my ($line,@inf,%genelen,%abu,$sum);
open IN, "$ARGV[0]" or die "can not open $ARGV[0]\n";
$/=">";
<IN>;
while($line=<IN>){
        chomp $line;
        @inf=split /\n/,$line;
        my @id=split /\s/,$inf[0];
        my $seq="";
        for(my $i=1;$i<=$#inf;$i++){
                $seq.=$inf[$i];}
        $genelen{$id[0]}=length($seq);
        $abu{$id[0]}=0;}
close IN;
$/="\n";
#sam results
if($ARGV[1]=~/\.bam$/){
        open IN,"$Bin/samtools view $ARGV[1]|" or die "can not open $ARGV[1]\n";
}
else{
        open IN,"$ARGV[1]" or die "can not open $ARGV[1]\n";
}
while($line=<IN>){
        chomp $line;
        @inf=split /\t/,$line;
        next if($line=~/^@/ || $inf[2] eq "*" ||$inf[5] eq "*");
        my ($m,$l)=(0,0);
        my @ele=split /[A-Z]/,$inf[5];
        for(my $i=0;$i<=$#ele;$i++){
                $l+=$ele[$i];
        }
        my $mn=$inf[5]=~s/M/M/g;
        if($mn==1){
                $inf[5]=~/(\d+)M/;
                $m=$1;
        }
        if($mn==2){
                $inf[5]=~/(\d+)M\d+\w(\d+)M/;
                $m+=$1+$2;
        }
        if($mn==3){
                $inf[5]=~/(\d+)M\d+\w(\d+)M\d+\w(\d+)M/;
                $m+=$1+$2+$3;
        }
#       print "$inf[5]\t$m\t$l\n";
        if($m>=50 && ($m/$l)>0.95){
                $sum++;
                $abu{$inf[2]}++;
        }
}
close IN;

open OA, ">$ARGV[2]" or die "can not open $ARGV[2]\n";
print OA "GeneID\tGeneLen\tReadsNum\tAbundance\tRelativeAbundance\tTotalAbundance\n";
my $temp=0;
foreach my $i (sort keys %genelen){
        $temp+=$abu{$i}/$genelen{$i}/$sum;
}
foreach my $i (sort keys %genelen){
        printf OA "$i\t%d\t%ld\t%.8e\t%.8e\t%.8e\n",$genelen{$i},$abu{$i},$abu{$i}/$genelen{$i}/$sum,$abu{$i}/$genelen{$i}/$sum/$temp,$temp;
}
close OA;

list.sh
close OA;
for i in sample_name1 sample_name2 sample_name3 ...
do
echo "$i        $i.sam.abd" >> list
done

abd.merge.pl
die "perl Profile-Merge.pl sam.list out.profile out.read \n" if(@ARGV!=3);
my ($line,@inf,%ff,@sample);
open IN, "$ARGV[0]" or die "can not open $ARGV[0]\n";
open OA, ">$ARGV[1]" or die "can not open $ARGV[1]\n";
open OB, ">$ARGV[2]" or die "can not open $ARGV[2]\n";
print OA "GeneID";      print OB "GeneID";
while($line=<IN>){
        chomp $line;    @inf=split /\s+/,$line;
        print OA "\t$inf[0]";   print OB "\t$inf[0]";
        push @sample,$inf[0];
        open $ff{$inf[0]}, "$inf[1]" or die "can not open $inf[1]\n";
        my $temp=$ff{$inf[0]};<$temp>;}
print OA "\n";  print OB "\n";  close IN;
open IN, "$inf[1]" or die "can not open $inf[1]\n";
<IN>;
while($line=<IN>){
        chomp $line;@inf=split /\t/,$line;
        print OA "$inf[0]";     print OB "$inf[0]";
        for(my $j=0;$j<=$#sample;$j++){
                my $temp=$ff{$sample[$j]};
                $line=<$temp>;chomp $line;@inf=split /\t/,$line;
                printf OA "\t%.6f",$inf[4]*1000000;
                print  OB "\t$inf[2]";}
        print OA "\n";  print OB "\n";}
close OA;       close IN;       close OB;

rpkm.r:
if(!require(optparse,quietly = TRUE)){
  install.packages("optparse")
  library(optparse)
}else{
  library(optparse)}
if(!require(data.table,quietly = TRUE)){
  install.packages("data.table")
  library(data.table)
}else{
  library(data.table)}
option_list <- list(
  make_option(c("--reads_count","-i"),type = "character",
              default = F,
              help = "Please set the directory of reads count file."),
  make_option(c("--gene_length","-l"),type = "character",
              default = F,
              help = "Please set the directory of gene_length file"),
  make_option(c("--output","-o"),type = "character",
              default = F,
              help = "Please set the result name's prefix."))
opt_parser <- OptionParser(
  usage = "usage: comts rpkm [option]",
  option_list = option_list,
  add_help_option = TRUE,
  prog = NULL,
  description = "To calculate the RPKM abundance of geneset.")
opt <- parse_args(opt_parser)
reads_count <- opt$reads_count
gene_length <- opt$gene_length
output <- opt$output
exprSet <- read.table(reads_count,header=T,row.names=1)
lengths <- fread(gene_length)
lengths <- as.vector(unlist(lengths))
total_count<- colSums(exprSet)
rpkm <- t(do.call(rbind,
                  lapply(1:length(total_count),
                         function(i){
                           10^9*exprSet[,i]/lengths/total_count[i]
                         }) ))
rownames(rpkm) <- rownames(exprSet)
colnames(rpkm) <- colnames(exprSet)
rpkm<-as.data.frame(rpkm)
write.table(rpkm,output,quote=F,sep = "\t")
