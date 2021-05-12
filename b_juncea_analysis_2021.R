#This is the code that is to be uploaded to Git - Open date May 12, 2021
#Brassica project 
#April 21st, 2021
#This project is based off of data collected by the Kentville research station AAFC

install.packages('ggplot2')
install.packages('ggfortify')
library(ggplot2)
library(ggfortify)

############################################################
#!# DATA EXPLORATION #!#
############################################################
setwd("~/Desktop/DALMSc/brassica_project")
metadata = read.csv('brassica_metadata.csv')
sensorydata = read.csv('brassica_sensory_data.csv')

#METADATA
#First, I want to do a general exploration of the data to see what I am dealing with 
total_sample_no = length(metadata$E_UNIT)
#132 samples in total
table(metadata$GENUS)
# All samples are from the genus Brassica
table(metadata$SPECIES)
#hirta     juncea     JUNCEA   narinosa perviridis       rapa 
#1         24         51          1          2         14 
#Total = 93
# hirta no longer belogs to the genus Brassica

change_juncea = which(metadata$SPECIES == "JUNCEA")
metadata$SPECIES[change_juncea] = "juncea"
table(metadata$SPECIES)
slices <- c(1, 75,1,39, 2, 14)
lbls <- c("hirta (1)", "juncea (75)", "narinosa (1)","unlabeled (39)", "perviridis (2)", "rapa (14)")
pie(slices, labels = lbls, main="Pie Chart of Species", init.angle = 80)


total_unlabelled = (sum(1*(is.na(metadata$SPECIES))))
#39 samples do not have species labels.

#A total of 75 Juncea samples. (56.8% of the samples are B. juncea)

table(metadata$COUNTRY)
# a total of 10 different labeled countries over 51 samples

table(metadata$E_UNIT)
# All E-units are unique

name_table = table(metadata$NAME)
# There are some replicates. Eg, there are 3 'Tatsoi' samples. 
#There are 95 unique accession names

length(table(metadata$SOURCE))


#SENSORY DATA
table(sensorydata$ASSESSOR)
#B_Amyotte      B_Lees S_MacKinnon 
#83          51         130 
table(sensorydata$NAME)
# Each name was sampled at least 2 times. 


#SENSORY DATA
table(sensorydata$LEAF_COLOUR)
#one entry of leaf colour is registering as a different G. This may have a space in it or some other invisible character. 
which(sensorydata$LEAF_COLOUR == "G ")
sensorydata$LEAF_COLOUR[14] = "G"
table(sensorydata$LEAF_COLOUR)

table(sensorydata$LEAF_SHAPE)
#There is a 6. This must be a typo. 
which(sensorydata$LEAF_SHAPE == '6')
sensorydata$LEAF_SHAPE[260] = '5'
table(sensorydata$LEAF_SHAPE)

table(sensorydata$LEAF_SIZE)
#Acceptable.

table(sensorydata$SWEET)
#There is a 12? This must be excluded
which(sensorydata$SWEET == '12')
sensorydata = sensorydata[-c(234),]
table(sensorydata$SWEET)

table(sensorydata$SOUR)
#There is a 23. This must be excluded
which(sensorydata$SOUR == '23')
sensorydata = sensorydata[-c(205),]
table(sensorydata$SOUR)

table(sensorydata$BITTER)
# three 0s, must be a typo
which(sensorydata$BITTER == '0')
sensorydata$BITTER[191] = '1'
sensorydata$BITTER[210] = '1'
sensorydata$BITTER[211] = '1'
table(sensorydata$BITTER)
#accpetable

table(sensorydata$AROMA)
#There are two 6s. These must be typos
which(sensorydata$AROMA == '6')
sensorydata$AROMA[94] = '5'
sensorydata$AROMA[95] = '5'
table(sensorydata$AROMA)

table(sensorydata$GRASSY)
#There are three 6s. These must be typos
which(sensorydata$GRASSY == '6')
sensorydata$GRASSY[5] = '5'
sensorydata$GRASSY[8] = '5'
sensorydata$GRASSY[9] = '5'
table(sensorydata$GRASSY)

table(sensorydata$HOT)
#There are two 6s, and 9 0s. These must be typos
which(sensorydata$HOT == '6')
sensorydata$HOT[94] = '5'
sensorydata$HOT[95] = '5'
zeros = which(sensorydata$HOT == '0')
sensorydata$HOT[zeros] = '1'
table(sensorydata$HOT)

table(sensorydata$CRISP)
#acceptable

table(sensorydata$DISSOLVING)
#acceptable (no 5s)

table(sensorydata$ASTRINGENT)
#acceptable

#There are a few rows that contain rows of NAs 



sensorydata_leaf_data = sensorydata[,6:17]
keep_index = which(complete.cases(sensorydata_leaf_data))
# only these rows should be kept
sensorydata = sensorydata[keep_index,]

#write.csv(sensorydata,file = '~/Desktop/DALMSc/brassica_project/final_brassica_sensory_data_table.csv',row.names = FALSE)


# some traits should also be encoded as binaries, so that they can be used to performa GWAS. 

sensorydata = read.csv('final_brassica_sensory_data_table.csv')
# The phenotype that should be changed from ordinal (which they arent) to binary, are: 
# Colour
# Original data has the phenotypes encoded as Light Green, Medium Green, Dark Green, Green-purple, and Purple. 
# This should be changed to: Level of green-ness (1,2,3) *NOTE: there will be some zeros because samples that score on the pruple scale will not have scored on the green scale and will therefor be 0s. 
# Purple: Binary (Individuals rated as green-purple AND purple are 1, those not are 0)
#Therefor this trait needs to be delt with in 2 columns, one that ranks 1-3 and one that is binary
#First, I will deal with the level of green-ness
table(sensorydata$LEAF_COLOUR)
# The following changes to encoding should be made: G = 1, M = 2, D = 3
#Get index of each: 
light_green_index = which(sensorydata$LEAF_COLOUR == 'G')
med_green_index = which(sensorydata$LEAF_COLOUR == 'M')
dark_green_index = which(sensorydata$LEAF_COLOUR == 'D')
purple_index1 = which(sensorydata$LEAF_COLOUR == 'GP')
purple_index2 = which(sensorydata$LEAF_COLOUR == 'P')
purple_index_total = c(purple_index1, purple_index2)

zero_vec = c(rep(0,93))
sensorydata$GREEN=zero_vec
#now fill in column with 1,2 or 3 for appropriate samples
sensorydata$GREEN[light_green_index] = 1
sensorydata$GREEN[med_green_index] = 2
sensorydata$GREEN[dark_green_index] = 3

#create purple column that will be binary
sensorydata$PURPLE=zero_vec
sensorydata$PURPLE[purple_index_total] = 1

#remove the leaf colour column as it is redundant. 
sensorydata = subset(sensorydata, select = -LEAF_COLOUR)
#re-write as new phenotype table

write.csv(sensorydata,file = '~/Desktop/DALMSc/brassica_project/final_brassica_sensory_data_table_with binary_traits.csv',row.names = FALSE)



############################################################
#!# DATA CURATION #!#
############################################################
#It would be useful to do some basic visualizations with the phenotype data. 
# I have decided that it would make the most sense to go with a dataset that was from a single sampler, because using 2 or 3 doesnt give that much power, and also narrows the sample size because only Mackinnon sampled all plants. 

sensorydata = read.csv('final_brassica_sensory_data_table.csv')

index_amyotte = which(sensorydata$ASSESSOR == 'B_Amyotte')
amyotte_samples = as.data.frame(sensorydata[index_amyotte,])
index_mackinnon = which(sensorydata$ASSESSOR == 'S_MacKinnon')
mackinnon_samples = as.data.frame(sensorydata[index_mackinnon,])
index_lees = which(sensorydata$ASSESSOR == 'B_Lees')
lees_samples = as.data.frame(sensorydata[index_lees,])

break_no = c(0,1,2,3,4,5)
mackinnon_samples_plotting = mackinnon_samples[,7:17]
which(is.na(mackinnon_samples_plotting))
# There are no NAs

#Basic Hists: 
layout(matrix(1:4, 2, 2, byrow = T))
for (i in 1:ncol(mackinnon_samples_plotting)) {
  hist(mackinnon_samples_plotting[,i], xlab = colnames(mackinnon_samples_plotting[i]), breaks = break_no, main = NULL)
}
dev.off()

#ggplot bar plots 
for (i in 1:ncol(mackinnon_samples_plotting)) {
  gg = ggplot( mackinnon_samples_plotting, aes(x=mackinnon_samples_plotting[,i])) +geom_histogram(binwidth = 1)+xlab(colnames(mackinnon_samples_plotting[i]))
  print(gg)
}

#Leaf colour barplot
colour_table = table(sensorydata$LEAF_COLOUR)
colour_count = c(rep(1,9),rep(2,69),rep(3,47),rep(4,126),rep(5,6))
colour_frame = data.frame(colour_count,row.names=NULL)
hh = ggplot(colour_frame, aes(x=colour_count)) +geom_histogram(binwidth = 1)+xlab('LEAF_COLOUR')
hh




install.packages("ggcorrplot")
library(ggcorrplot)
corr <- round(cor(mackinnon_samples_plotting), 3)
p.mat <- p.adjust(cor_pmat(mackinnon_samples_plotting))
pstuff <- (cor_pmat(mackinnon_samples_plotting))
ggcorrplot(corr, p.mat = pstuff,type = 'lower',insig = 'blank')
ggcorrplot(corr, p.mat = p.mat,insig = 'blank')




#############################################################
#!# KEYFILE #!#
#############################################################
setwd("~/Desktop/DALMSc/brassica_project")
GBS_data = read.delim('mustardGreenKeyfile_edited.txt')
#There is GBS data for 96 samples. Therefor not every one of the 134 samples that were taste-tested were genotyped. I will only be presenting the phenotype data for those samples that I also have genotype data for. 

#I need to find out which samples have GBS data, and then use that to trim down the phenotype data. 

#What is the species breakdown for the GBS data?

gbs_species_index_table = mackinnon_samples[mackinnon_samples$E_UNIT %in% GBS_data$FullSampleName,]

meta_index = metadata$E_UNIT %in% GBS_data$FullSampleName

meta_table_trimmed = metadata[metadata$E_UNIT %in% GBS_data$FullSampleName,]
table(meta_table_trimmed$SPECIES)

#This will produce a pie chart showing the species breakdown for GBS data 
slices2 <- c(1, 64,1,19, 2, 9)
lbls2 <- c("hirta (1)", "juncea (64)", "narinosa (1)","unlabeled (19)", "perviridis (2)", "rapa (9)")
pie(slices2, labels = lbls2, main="Pie Chart of Species (GBS)", init.angle = 80)


gbs_species_index = gbs_species_index_table$E_UNIT

for (i in gbs_species_index){
  metadata[i,]
}

#############################################################
#!# EDIT KEYFILE #!#
#############################################################
setwd("~/Desktop/DALMSc/brassica_project")
GBS_data = read.delim('mustardGreenKeyfile_edited.txt')

# I need a table that has the the following layout:
#Sample names, barcodes, RE, .....everything else 
GBS_for_GBX = GBS_data
GBS_for_GBX["RE"]= 'PstI'
GBS_for_GBX= GBS_for_GBX[,c(4,3,9,1,2,5,6,7,8)]
true_gbx = GBS_for_GBX[,c(1,2,3)]

write.table(true_gbx, file ="~/Desktop/DALMSc/brassica_project/barcodes_for_GBSX.txt", row.names = F, col.names = F, sep = '\t', quote = F)

#############################################################
#!# DEMULTIPLEXING #!#
#############################################################
#The fastq files supplied from Platform Genetics were demultiplexed using GBSX via the following code: 

#FILE NAME : run_demilti_v2
#!/usr/bin/bash

#SBATCH --ntasks=1
#SBATCH --mem-per-cpu=24G
#SBATCH --account=def-smyles
#SBATCH --job-name=6-demux-mustard
#SBATCH --output=demulti_mustard1.out
#SBATCH --time=1-00:00 


module load java

java \
-jar /project/6003429/myles_lab/bin/GBSX-master/releases/GBSX_v1.3/GBSX_v1.3.jar \
--Demultiplexer \
-f1 /project/def-smyles/myles_lab/mustard_GBS/tom-snp-calling-test/data/HC3M3CCX2_1.fastq.gz \
-f2 /project/def-smyles/myles_lab/mustard_GBS/tom-snp-calling-test/data/HC3M3CCX2_2.fastq.gz \
-i /home/tdavies/projects/def-smyles/myles_lab/mustard_GBS/tom-snp-calling-test/barcodes_for_GBSX.txt \
-gzip true \
-o /home/tdavies/projects/def-smyles/myles_lab/mustard_GBS/tom-snp-calling-test/data\
-mb 0 \
-me 1 \

#The demultiplexer finished with a total of 229,308,764 reads demultiplexed

#This produced 2 fastq files for each sample that are filed in the /data folder 
#This also produced a stats file, which I will load in and take a look at below

correction_stats = read.delim('correction.stats', header = F)
correction_stats = correction_stats[-c(1),]
#turns out this file has two distinct sections that need to be looked at separately, I cant really make sense of what each of these mean

gbs_stats = read.delim('gbsDemultiplex.stats')
read_total_sum = (gbs_stats$total.count)
sum(read_total_sum)
#229308764, this is the total number of reads 
#the final row in the dataframe shows how many were undetermined 
gbs_stats$total.count[97]
#7611234 , 7.6 million reads were undetermined 
total_perc_sum = gbs_stats$total.perc
sum(total_perc_sum)
#1, this sums to 1. I think this column is telling me what % of the reads each sampled had 
gbs_stats$total.perc[97]
# 3.319% of all reads were undetermined

#I want to take a look at the histograms from each relevant column from the gbsDemultiplex.stats file. Columns 5-10. 

for (i in colnames(gbs_stats[,5:10])){
  
  hist(gbs_stats[,i], main = i)
}
#Interpretations of histograms
#TOTAL PERC: The proportion of the reads that this sample accounts for. 
#makes a normal looking distribution, with a n outlier at each end. The low outlier is sample 148. This sample has very few reads, and thus the total percent of reads that in accounts for in the data is low. The high outlier is 'undetermined' reads. This makes sense, as undetermined reads had the highest read count. (samples 148,176,502,525)
#MISMATCH.0.COUNT: The number of reads with 0 barcode mismatches. 
# This produces a normal looking distribution. THis number will be the same as the total.count. Thus, the samples with the lowest read counts, and highest read counts are on the tails. 
#MISMATCH.0.PERC: The reads with 0% barcode mismatch. % of reads with no mismatch aka proportion of reads that have no mismatches
#This is 1 for all samples, which checks out with previous histogram. The hist made by the loop by default makes a massive square, but this is just because of default binning, and all the values are 1.
# I will attempt to create the proper histogram for this now: 
hist(gbs_stats$mismatch.0.perc, main = 'mismatch.0.perc', xlim = c(-2,3), ylim = c(0,100), freq = TRUE)
gbs_mismatch = (table((gbs_stats$mismatch.0.perc)*100))
barplot(gbs_mismatch, ylim = c(0,100), xlim = c(0,4), ylab = "Number of samples (N=96)", xlab = "?", main ="mismatch.0.perc" )

#BASECALL.COUNT: the total number of bases called for each sample
#this creates a normal looking distribution. As expected, the same samples end up at the tail of the distribution (samples 148,176,502,525). 
#BASECALL.ABOVE.30.PERC: the proportion of the reads for that sample that are above 30 percent ????
#There is a single outlier, sample 148. All other samples are between 0.86-0.89. 
#BASECALL.QUAL.AVG: The average phred score for the reads for each sample. 
#There is a single outlier, sample 148. All other samples are between 36.9-37.9. This is just shy of 99.99% accuracy, which is pretty solid. 

#CONCLUSIONS: 
#Samples 148,176,502,525 should probably be removed, but these will will be eliminated in downstream pipeline steps that are built for filtering out low quality samples. 


#below is a bar chart that will display the read counts for each sample. 
bar_x = as.integer(gbs_stats$sampleID[1:96])
bar_y = (gbs_stats$total.count[1:96])
y_axis_buckets = c(seq(0,8000000,1000000))
barplot(bar_y, xlab='Samples', ylab = 'Read Count', ylim = c(0,8000000))
#version with bar for undetermined
bar_y = (gbs_stats$total.count)
y_axis_buckets = c(seq(0,8000000,1000000))
barplot(bar_y, xlab='Samples', ylab = 'Read Count', ylim = c(0,8000000))

#############################################################
#!# TRIMMING #!#
#############################################################
# the following shell script was used to trim reads. It is saved as run_trimming2.sh

#!/usr/bin/bash


#SBATCH --ntasks=1
#SBATCH --mem-per-cpu=24G
#SBATCH --account=def-smyles
#SBATCH --job-name=3-trimmomatic
#SBATCH --output=mustard-trim.out
#SBATCH --time=1-00:00


module load java
module load trimmomatic


NUM_CPU=20
OUTPUT_DIR="/project/def-smyles/myles_lab/mustard_GBS/tom-snp-calling-test"
for file in /project/def-smyles/myles_lab/mustard_GBS/tom-snp-calling-test/data/*.R1.fastq.gz
do
fname=$(basename $file)
base=$(echo "$fname" | cut -d '.' -f1)

java -jar $EBROOTTRIMMOMATIC/trimmomatic-0.39.jar PE -threads $NUM_CPU \
-trimlog $OUTPUT_DIR/trimlog.txt \
-summary $OUTPUT_DIR/summaryStats.txt \
$OUTPUT_DIR/data/$base.R1.fastq.gz \
$OUTPUT_DIR/data/$base.R2.fastq.gz \
-baseout $OUTPUT_DIR/data/$base \
SLIDINGWINDOW:4:20 MINLEN:25

#The output file indicated that trimming was done successfully. 
#The files produced from this command were 2 output files per input file. EX: 110.R1.fastq.gz was output as 110_1P 110_1U. P = paired reads, U = Unpaired reads 

#These files were combined into 2 fastqfiles, named all_1P and all_2P. 1P were the forward paired reads, and 2P were the reverse paired reads. 



#############################################################
#!# ALIGNMENT #!#
#############################################################
#BWA was used to align reads to the B. juncea reference genome, producing .amb, .ann, .bwt, .pac, .sa files. 

#!/usr/bin/bash

#SBATCH --ntasks=1
#SBATCH --mem=30G
#SBATCH --time=1-00:00:00
#SBATCH --cpus-per-task=30
#SBATCH --account=def-smyles
#SBATCH --job-name=align
#SBATCH --output=align.out

module load bwa/0.7.17
module load samtools/1.12
module load sambamba/0.8.0

REF_FILE="/home/tdavies/projects/def-smyles/myles_lab/mustard_GBS/tom-snp-calling-test/4-alignment/GCA_015484525.1_UDSC_Var_1.1_genomic.fna.gz"
FASTQ_FILE1="/home/tdavies/projects/def-smyles/myles_lab/mustard_GBS/tom-snp-calling-test/data/all_1P.fastq"
FASTQ_FILE2="/home/tdavies/projects/def-smyles/myles_lab/mustard_GBS/tom-snp-calling-test/data/all_2P.fastq"
NUM_CPU=20
OUTPUT_BAM="/home/tdavies/projects/def-smyles/myles_lab/mustard_GBS/tom-snp-calling-test/4-alignment/mustard_GBS.bam"
COMPRESSION_LEVEL=5


# ALIGNMENT
bwa index $REF_FILE
bwa mem $REF_FILE $FASTQ_FILE1 $FASTQ_FILE2 \
-R '@RG\tID:myreadgroup\tSM:sid1\tLB:null\tPL:illumina\tCN:null' \
-t $NUM_CPU -v 3 -Y -M | samtools view -1 > $OUTPUT_BAM



#############################################################
#!# SORT, INDEX, MARK DUPLICATES #!#
#############################################################
#SORTING 

#SBATCH --ntasks=1
#SBATCH --mem=30G
#SBATCH --time=1-00:00:00
#SBATCH --cpus-per-task=30
#SBATCH --account=def-smyles
#SBATCH --job-name=sort
#SBATCH --output=sort.out

module load bwa/0.7.17
module load samtools/1.12
module load sambamba/0.8.0

REF_FILE="/home/tsoomro/projects/def-smyles/myles_lab/tayab-snp-calling-test/data/ref/GCA_015484525.1_UDSC_Var_1.1_genomic.fna"
NUM_CPU=20
OUTPUT_BAM="/home/tdavies/projects/def-smyles/myles_lab/mustard_GBS/tom-snp-calling-test/4-alignment/mustard_GBS.bam"
OUTPUT_SORTED_BAM="/home/tdavies/projects/def-smyles/myles_lab/mustard_GBS/tom-snp-calling-test/5-sort-index-markdups/mustard.sort.bam"
OUTPUT_SORTED_NODUP_BAM="/home/tdavies/projects/def-smyles/myles_lab/mustard_GBS/tom-snp-calling-test/5-sort-index-markdups/mustard.sort.nodup.bam"
COMPRESSION_LEVEL=5

# SORTING
echo "Starting Sorting..."
#time sambamba index -m 500M --tempdir ./ -l $COMPRESSION_LEVEL --show-progress --nthreads $NUM_CPU
time samtools sort -@ $NUM_CPU -m 500M $OUTPUT_BAM -o $OUTPUT_SORTED_BAM -T ./ -O BAM -l $COMPRESSION_LEVEL
echo "Finished Sorting"

#This produced the "mustard.sort.bam" file 

#INDEX 

#!/usr/bin/bash

#SBATCH --ntasks=1
#SBATCH --mem=30G
#SBATCH --time=1-00:00:00
#SBATCH --cpus-per-task=30
#SBATCH --account=def-smyles
#SBATCH --job-name=index
#SBATCH --output=index.out

module load bwa/0.7.17
module load samtools/1.12
module load sambamba/0.8.0

REF_FILE="/home/tsoomro/projects/def-smyles/myles_lab/tayab-snp-calling-test/data/ref/GCA_015484525.1_UDSC_Var_1.1_genomic.fna"
NUM_CPU=20
OUTPUT_BAM="/home/tdavies/projects/def-smyles/myles_lab/mustard_GBS/tom-snp-calling-test/4-alignment/mustard_GBS.bam"
OUTPUT_SORTED_BAM="/home/tdavies/projects/def-smyles/myles_lab/mustard_GBS/tom-snp-calling-test/5-sort-index-markdups/mustard.sort.bam"
OUTPUT_SORTED_NODUP_BAM="/home/tdavies/projects/def-smyles/myles_lab/mustard_GBS/tom-snp-calling-test/5-sort-index-markdups/mustard.sort.nodup.bam"
COMPRESSION_LEVEL=5


# INDEXING
echo "Starting Indexing..."
time sambamba index --show-progress -t $NUM_CPU $OUTPUT_SORTED_BAM
echo "Finished Indexing"

#This produced the mustard.sort.bam.bai file. 

#PLACE DUPLICATE MARKING COMMAND HERE 

#After marking duplicates, Picard spits out a duplicate statistics file, which I will look at below: 

setwd("~/Desktop/DALMSc/brassica_project/Prelim_data")
duplicate_stats = read.delim("marked_dups_metrics.txt", row.names = NULL)


