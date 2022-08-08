# LASKY-MORRIS LAB: Sorghum project

## Table of contents    
* [Page 1: 2020-12-05](#id-section1). Chapter 1: Getting resequencing data (by Aayudh)

* [Page 2: 2020-12-05](#id-section2). Chapter 2: SnpEff (by Luke)

* [Page 3 2020-12-06](#id-section3). Chapter 3: Running Beagle (by Luke) 

* [Page 4 2020-12-07](#id-section4). Chapter 4: Environmental GWAS (by Aayudh)

* [Page 5 2020-12-08](#id-section5). Chapter 5: GWAS using vcftools

* [Page 6 2020-12-08](#id-section6). Chapeter 6: 

* [Page 7 2020-12-14](#id-section7). Chapeter 7: 

* [Page 8 2020-12-16](#id-section8). Chapeter 8: 

* [Page 9 2020-12-22](#id-section9). Chapter 9: 

* [Page 10 2020-12-22](#id-section10). Chapter10: 

* [Page 11 2021-04-06](#id-section11). Chapter11: 
* [Page 12 2021-04-06](#id-section12). Chapter12: 

------
<div id='id-section1'/>

## Chapter 1: Getting resequencing data (by Aayudh)

#Getting resequencing data
Login info: ssh -l azd6024 submit.aci.ics.psu.edu

### Data details
Here is a vcf for variants (SNP, MNP, indel independently) for all the resequencing libraries of Sorghum. This includes
375 BAP-TERRA
90 CASP-TERRA-PNNL
914 Gates
21 GatesRef
14 JGI_Ref
109 Mullet
230 TERRA-Cornell
Download LInks:
http://hagsc.org/restricted_access/sorg_wga/Sorghum_1757g_AllChr_1757g.tar.gz
http://hagsc.org/restricted_access/sorg_wga/Sorghum_all_1757g_metadata.txt


1. Code to wget
```
wget --user xyz --password xyz https://hagsc.org/restricted_access/sorg_wga/Sorghum_1757g_AllChr_1757g.tar.gz
```
2. Unzip data

```
tar -zxvf Sorghum_1757g_AllChr_1757g.tar.gz
```
**Getting GBS SNP data 10.3835/plantgenome2018.06.0044**

ssh azd6024@128.118.42.11

```
cd /data/Sorghum/Hu_GBS_2019/
```
**DATA STORAGE IN LASKY LAB**

```
/gpfs/group/jrl35/default/aayudh
```

-----
<div id='id-section2'/>

## Chapter 2: Running Beagle (by Luke) 

**Rationale**
Beagle performs genomic imputation on variant calls, which will be useful for downstream analyses. The goal is for everyone to use the same set of 
imputations. The problem with Beagle in this analysis is that the input vcf file for the full resequencing dataset is huge (>430GB), which necessitates breaking it up into smaller pieces to avoid memory issues. The separate imputed chunks then need to be concatenated back together again to obtain a single imputed file. 

```
# Run Beagle on the relatively small indels vcf file 
java -jar beagle.05May22.33a.jar gt=Sorghum_1757g_AllChr.polymorphic.indel.noRepeats.5pctMasked.vcf.gz out=Sorghum_1757g_AllChr.polymorphic.indel.noRepeats.5pctMasked.imputed

# Split the very large (432GB unzipped) vcf file for SNPs into ~10,000-line chunks
java -jar SnpSift.jar split -l 10000 Sorghum_1757g_AllChr.polymorphic.snp.noRepeats.5pctMasked.vcf.gz

# Run Beagle on each of the subset vcf files output from SnpSift
for f in *.vcf
do
	java -jar /projects/luwh7529/software_builds/Beagle/beagle.05May22.33a.jar gt=$f out=$f".imputed"
done

# Run a series of steps to bgzip, index, and bcftools concatenate the separate imputed vcf files

#gunzip the existing gz files
gunzip *.imputed.vcf.gz

#rezip the files using bgzip (otherwise bcftools won't parse the files)
for f in *.vcf
do
	bgzip $f
done


#index the files using bcftools (this way they can be concatenated together)
for f in *.imputed.vcf.gz
do
	bcftools index $f
done


#concatenate the resulting bgzipped and index vcf files into a single combined file
ls *imputed.vcf.gz > vcfoutfiles.txt

bcftools concat --file-list vcfoutfiles.txt > Sorghum_1757g_AllChr.polymorphic.snp.noRepeats.5pctMasked.imputed.combined.vcf

bgzip Sorghum_1757g_AllChr.polymorphic.snp.noRepeats.5pctMasked.imputed.combined.vcf

```

------


<div id='id-section3'/>

## Chapter 3: SnpEff (by Luke)

**Introduction**
SnpEff is an open source tool that annotates variants and predicts their effects on genes by using an interval forest approach. This program takes pre-determined variants listed in a data file that contains the nucleotide change and its position and predicts if the variants are deleterious.

**Rationale**
Running SnpEff provides predicted effect sizes for the observed SNPs in the variant dataset

**Codes**

```
# Run SnpEff on the full resequencing variant calls vcf SNPs
java -jar snpEff.jar -c snpEff.config -v Sorghum_bicolor Sorghum_1757g_AllChr.polymorphic.snp.noRepeats.5pctMasked.vcf > Sorghum_1757g_AllChr.polymorphic.snp.noRepeats.5pctMasked.ann.vcf

# Run SnpEff on the full resequencing variant calls vcf indels
java -jar snpEff.jar -c snpEff.config -v Sorghum_bicolor Sorghum_1757g_AllChr.polymorphic.indel.noRepeats.5pctMasked.vcf.gz > Sorghum_1757g_AllChr.polymorphic.indel.noRepeats.5pctMasked.ann.vcf

# Run snpEff on the *imputed* full resequencing vcf variant calls SNPs
java -jar snpEff.jar -c snpEff.config -v Sorghum_bicolor Sorghum_1757g_AllChr.polymorphic.snp.noRepeats.5pctMasked.imputed.combined.vcf > Sorghum_1757g_AllChr.polymorphic.snp.noRepeats.5pctMasked.imputed.combined.ann.vcf

# Run snpEff on the *imputed* full resequencing vcf variant calls indels
java -jar snpEff.jar -c snpEff.config -v Sorghum_bicolor Sorghum_1757g_AllChr.polymorphic.indel.noRepeats.5pctMasked.imputed.vcf.gz > Sorghum_1757g_AllChr.polymorphic.indel.noRepeats.5pctMasked.imputed.ann.vcf

```


------


<div id='id-section4'/>

## Chapter 4: Environmental GWAS (by Aayudh)

Let’s take a peek inside the vcf file first. 
Note: zcat lets us open a .gz (gzipped) file; we then “pipe” | this output from zcat to the head command and print as many lines as we want -n #

**Details of the vcf.gz file**

```
zcat Sorghum_1757g_AllChr.polymorphic.snp.noRepeats.5pctMasked.imputed.combined.vcf.gz | head -n 9
```

```
##fileformat=VCFv4.2
##FILTER=<ID=PASS,Description="All filters passed">
##filedate=20220706
##source="beagle.05May22.33a.jar"
##INFO=<ID=AF,Number=A,Type=Float,Description="Estimated ALT Allele Frequencies">
##INFO=<ID=DR2,Number=A,Type=Float,Description="Dosage R-Squared: estimated squared correlation between estimated REF dose [P(RA) + 2*P(RR)] and true REF dose">
##INFO=<ID=IMP,Number=0,Type=Flag,Description="Imputed marker">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=DS,Number=A,Type=Float,Description="estimated ALT dose [P(RA) + 2*P(AA)]">
```



### STEPS for GWAS

**1. First create gds file from vcf or vcf.gz file.**

```
#!/usr/bin/env Rscript

#PBS -l nodes=1:ppn=8
#PBS -l walltime=24:00:00
#PBS -l pmem=48gb
#PBS -M azd6024@psu.edu
#PBS -A open
#PBS -j oe

setwd("~/scratch/test_reseq_gunzip")

library(gdsfmt)
library(SNPRelate)

vcf.fn <- "/storage/home/azd6024/scratch/test_reseq_gunzip/Sorghum_1757g_AllChr.polymorphic.snp.noRepeats.5pctMasked.imputed.combined.vcf"

snpgdsVCF2GDS(vcf.fn, "Sorghum_1757_Reseq.gds", method="biallelic.only")
```

**Create the id matching files with trait value** 

```
setwd("~/Library/CloudStorage/OneDrive-UniversityofVermont/PENN STATE/resequencing data_GWAS")
list.files()

lasky_2015 <- read.csv("Lasky_2015_Accessions_cords.csv")
head(lasky_2015)

Reseq_matadata <- read.csv("Sorghum_all_1757g_metadata.csv")
head(Reseq_matadata)
final_data_1 <- gsub("_","",Reseq_matadata$sample)
head(final_data_1)
Reseq_matadata_1 <- cbind (Reseq_matadata,final_data_1)
head(Reseq_matadata_1)
names(Reseq_matadata_1)[names(Reseq_matadata_1) == 'final_data_1'] <- 'Accession'

Reseq_lasky<- merge(Reseq_matadata_1,lasky_2015, by="Accession")
dim(Reseq_lasky)
head(Reseq_lasky)
library(data.table)
MASTER_data_final <- unique(setDT(Reseq_lasky)[order(Accession, -Accession)], by = "Accession")
dim(MASTER_data_final)
head(MASTER_data_final)

library(raster)
library(sp)
library(rgdal)
library(tidyverse)

setwd("~/OneDrive - University of Vermont/PENN STATE/RAstor data")
list.files()

preds.all <- raster(paste0("~/OneDrive - University of Vermont/PENN STATE/RAstor data/preds.all.tif"))
preds.all

Suitability_data <- raster::extract(preds.all, MASTER_data_final[,c("Lon","Lat")])
ReseqGWAS_traits <- cbind(MASTER_data_final,Suitability_data)
setwd("~/Library/CloudStorage/OneDrive-UniversityofVermont/PENN STATE/resequencing data_GWAS")
write.table(ReseqGWAS_traits, "1. ReseqGWAS_traits.csv", sep=",")
```

**3. From genofile to ibs.matrix,mySNPbed.bed and HS_score_BimBam.txt file.** 

```
#!/usr/bin/env Rscript

#PBS -l nodes=1:ppn=8
#PBS -l walltime=24:00:00
#PBS -l pmem=48gb
#PBS -M azd6024@psu.edu
#PBS -A open
#PBS -j oe

setwd("~/scratch/gwas_reseq")

library(gdsfmt)
library(SNPRelate)

genofile <- snpgdsOpen('Sorghum_1757g_2ndtry.gds')
sample.id <- read.gdsn(index.gdsn(genofile, 'sample.id'))
snp.id <- read.gdsn(index.gdsn(genofile, 'snp.id'))
snp.position <- read.gdsn(index.gdsn(genofile, 'snp.position'))
chromosome.id <- read.gdsn(index.gdsn(genofile, 'snp.chromosome'))
snp.allele <- read.gdsn(index.gdsn(genofile, 'snp.allele'))
snp.rs.id <- read.gdsn(index.gdsn(genofile, 'snp.rs.id'))

#you can use your traits file to get SNPs only for genotypes with phenotypic data
#the trait file should have a column 'sample.id' that corresponds to the id in the genofile
traits <- data.frame(read.csv('1. ReseqGWAS_traits.csv', header=TRUE))
head(traits)

#no need to account for LD when getting SNPs for GWAS. 
mySNPs <- snpgdsSelectSNP(genofile, traits$LIB, maf = 0.05, missing.rate = 0)
length(mySNPs)

mySNPmatrix <- snpgdsGetGeno(genofile, sample.id = traits$LIB, snp.id = mySNPs, with.id = TRUE)
dim(mySNPmatrix$genotype)

#these lines create the bed file for gemma
rs <- data.frame(chr = chromosome.id[mySNPs], positions = snp.position[mySNPs])
SNPs_bim <- data.frame(paste0('rs_', rs$chr, '_', rs$positions), 'A', 'T', t(mySNPmatrix$genotype)) 
write.table(SNPs_bim, 'mySNPbed.bed', row.names = FALSE, quote = FALSE, col.names = FALSE, sep = ',') 

#kinship matrix for gemma obtained with 'identity by state'
ibs <- snpgdsIBS(genofile, sample.id=traits$LIB, snp.id=mySNPs, maf = 0.05, missing.rate = 0)
ibs.matrix <- ibs$ibs
dim(ibs.matrix)

write.table(ibs.matrix, 'ibs_matrix.txt', row.names = FALSE, col.names = FALSE, quote = FALSE, sep = ', ')

#prepare the trait data for gemma, accessions should be in the same order as in the genetic matrices
order <- as.data.frame(ibs$sample.id)
dim(order)
colnames(order) <- c("sample.id")

order.idx <- match(order$sample.id, traits$LIB)
order.idx

ordered <- traits[order.idx,]

write.table(ordered$Suitability_data, 'HS_score_BimBam.txt', row.names = FALSE, quote = FALSE, col.names = FALSE)
```




<div id='id-section5'/>

## Chapter 5: GWAS using vcftools



**Intalling vcftools**

```
cd ~/work
git clone https://github.com/vcftools/vcftools.git
cd vcftools/
./autogen.sh 
./configure --prefix=/storage/work/azd6024/vcftools
make
make install
   
echo "export PATH=/storage/work/azd6024/vcftools/bin:$PATH" >> ~/.bashrc
source ~/.bashrc
vcftools --help
 ```
**1. Getting info with vcf.gz file**

Create a .sh file (vi test.sh) with a code embedded to run it in server

```
#!/bin/bash

#PBS -l nodes=1:ppn=8
#PBS -l walltime=10:00:00
#PBS -l pmem=48gb
#PBS -A open

WORKINGDIR=/storage/home/azd6024/work/test_vcffiles
cd $WORKINGDIR

vcftools --gzvcf Sorghum_1757g_AllChr.polymorphic.snp.noRepeats.5pctMasked.imputed.combined.vcf.gz
```
Paste this code in the test.sh and save by *wq* then *qsub test.sh*

*WHAT THIS VCF.GZ FILE CONTAINS?*

After filtering, kept 1757 out of **1757 Individuals**
After filtering, kept 64618398 out of a possible **64618398 Sites**
 
