# LASKY-MORRIS LAB: Sorghum project

## Table of contents    
* [Page 1: 2020-12-05](#id-section1). Chapter 1: Getting resequencing data (by Aayudh)

* [Page 2: 2020-12-05](#id-section2). Chapter 2: SnpEff (by Luke)

* [Page 3 2020-12-06](#id-section3). Chapter 3: Running Beagle (by Luke) 

* [Page 4 2020-12-07](#id-section4). Chapter 4: Environmental GWAS (by Aayudh)

* [Page 5 2020-12-08](#id-section5). Chapter 5: Mapping SNPs and finding SNPs nearby (by Aayudh)

* [Page 6 2020-12-08](#id-section6). Chapeter 6: Finding nearby SNPs 

* [Page 7 2020-12-14](#id-section7). Chapeter 7: Vcftools related functions 

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

**2. Create the id matching files with trait value** 

We are using accessions from Lasky et al (2015). NOTE: You need to match the genofile$sample.id format with the metadata column.

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
#PBS -l walltime=12:00:00
#PBS -l pmem=24gb
#PBS -M azd6024@psu.edu
#PBS -A jrl35_c_g_sc_default
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

Excluding 3,238,815 SNPs on non-autosomes
Excluding 53,629,157 SNPs (monomorphic: TRUE, MAF: 0.05, missing rate: 0)
Genotype matrix: 339 samples X 7750426 SNPs
Identity-By-State (IBS) analysis on genotypes:
Excluding 56,867,972 SNPs (non-autosomes or non-selection)
Excluding 0 SNP (monomorphic: TRUE, MAF: 0.05, missing rate: 0)
Working space: 339 samples, 7,750,426 SNPs
IBS:    the sum of all selected genotypes (0,1,2) = 3811307892


**4. Running gemma** 

gemma: https://github.com/genetics-statistics/gemma-wrapper/blob/master/README.md 

Install gemma
```
wget https://github.com/genetics-statistics/GEMMA/releases/download/v0.98.5/gemma-0.98.5-linux-static-AMD64.gz
gzip -d gemma-0.98.5-linux-static-AMD64.gz 
mv gemma-0.98.5-linux-static-AMD64 gemma-0.98.5-linux-static
chmod 700 gemma-0.98.5-linux-static
```
*Run gemma*

```
#!/bin/bash

#PBS -l nodes=1:ppn=8
#PBS -l walltime=12:00:00
#PBS -l pmem=24gb
#PBS -M azd6024@psu.edu
#PBS -A jrl35_c_g_sc_default
#PBS -j oe

WORKINGDIR=/storage/home/azd6024/work/gwas_reseq
cd $WORKINGDIR

./gemma-0.98.5-linux-static -g mySNPbed.bed -k ibs_matrix.txt -lmm 4 -miss 0.1 -p HS_score_BimBam.txt -o Reseq_gwas_HS_score_out
```
**preds.all**

**number of total individuals = 339**

number of analyzed individuals = 339
number of covariates = 1
number of phenotypes = 1

**number of total SNPs/var        =  7750426**

**number of analyzed SNPs         =  7750411**

Start Eigen-Decomposition...
pve estimate =0.902493
se(pve) =0.109123

REMLE log-likelihood in the null model = 98.4581

MLE log-likelihood in the null model = 99.5659

pve estimate in the null model = 0.58953

se(pve) in the null model =  0.103705

vg estimate in the null model = 0.0933044

ve estimate in the null model = 0.0175318

beta estimate in the null model =   0.270463

se(beta) =   0.00719141


**5. Clean gemma txt file** 

```
setwd("~/Library/CloudStorage/OneDrive-UniversityofVermont/PENN STATE/resequencing data_GWAS")
list.files()


gemma_output <- read.table("Reseq_gwas_HS_score_out.assoc.txt", header = TRUE, sep = "", dec = ".")
head(gemma_output)
tail(gemma_output)
dim(gemma_output)

rs <- gemma_output[,2]
tmp <- sub('rs_', '', rs)
chr <- sub('_.*$','',tmp)
ps <- sub('._', '', tmp)
gemma_output.clean <- cbind(chr,rs,ps,gemma_output[,c(4:15)])
head(gemma_output.clean)
write.table(gemma_output.clean, 'Reseq_preds_all_HS_score_out.assoc.clean.txt', row.names = FALSE, quote = FALSE, col.names = TRUE, sep = '\t')
```

**6. FDR corrections** 

Same script can be run for p_lrt and p_score.

```
#!/usr/bin/env Rscript

#PBS -l nodes=1:ppn=8
#PBS -l walltime=12:00:00
#PBS -l pmem=24gb
#PBS -M azd6024@psu.edu
#PBS -A open
#PBS -j oe

setwd("~/work/preds_all_gwas/output")

glm_stats <- read.table("Reseq_preds_all_HS_score_out.assoc.clean.txt", header = T, sep = "\t")
head(glm_stats)

library(dplyr)

# Calculate Bonferroni Correction and False Discovery Rate 

adj_glm <- glm_stats %>%
  transmute(rs, chr, ps,af,p_wald,
            p_Bonferroni =  p.adjust(glm_stats$p_wald,"bonferroni"),
            p_FDR = p.adjust(glm_stats$p_wald,"fdr")
  )

head(adj_glm)

setwd("~/work/preds_all_gwas/fdr_correction")

write.csv(adj_glm, file="Reseq_preds_all_HS_score_adj_p_GLM.csv", quote = T, eol = "\n", na= "NA")

library(qqman)

adj_glm_KRN_4 <- read.csv("Reseq_preds_all_HS_score_adj_p_GLM.csv", header = T)
head(adj_glm_KRN_4)

tiff("Reseq_preds_all_HS_score_FDR.tiff", width = 11, height = 7, units = 'in', res = 300)
par(mfrow=c(1,3))
qq(adj_glm_KRN_4$p_wald, main = "non-adjusted P-value")
qq(adj_glm_KRN_4$p_Bonferroni, main = "Bonferroni")
qq(adj_glm_KRN_4$p_FDR, main = "FDR")
par(mfrow=c(1,1))
dev.off()

png(file = 'Reseq_preds_all_HS_score_FDR.png', width=1400, height=960, res=300)
par(mfrow=c(1,3))
qq(adj_glm_KRN_4$p_wald, main = "non-adjusted P-value")
qq(adj_glm_KRN_4$p_Bonferroni, main = "Bonferroni")
qq(adj_glm_KRN_4$p_FDR, main = "FDR")
par(mfrow=c(1,1))
dev.off()
```
![alt text](https://github.com/aayudh454/Lasky-Morris-Lab-Sorghum-project/blob/main/fdr%20correction.png)

**7. Making manhattan plot** 

```
setwd("~/Library/CloudStorage/OneDrive-UniversityofVermont/PENN STATE/resequencing data_GWAS/manhattan")
list.files()

library(tidyverse)
library(ggtext)
library(normentR)

gwas_data_load <- read.csv("Reseq_gwas_HS_score_adj_p_GLM.csv", header=T)
head(gwas_data_load)

sig_data <- gwas_data_load %>% 
  subset(p_wald < 0.05)
notsig_data <- gwas_data_load %>% 
  subset(p_wald >= 0.05) %>%
  group_by(chr) %>% 
  sample_frac(0.1)
gwas_data <- bind_rows(sig_data, notsig_data)

data_cum <- gwas_data %>% 
  group_by(chr) %>% 
  summarise(max_bp = max(ps)) %>% 
  mutate(bp_add = lag(cumsum(max_bp), default = 0)) %>% 
  select(chr, bp_add)

gwas_data <- gwas_data %>% 
  inner_join(data_cum, by = "chr") %>% 
  mutate(bp_cum = ps + bp_add)

axis_set <- gwas_data %>% 
  group_by(chr) %>% 
  summarize(center = mean(bp_cum))

ylim <- gwas_data %>% 
  filter(p_wald == min(p_wald)) %>% 
  mutate(ylim = abs(floor(log10(p_wald))) + 2) %>% 
  pull(ylim)

sig <- 5e-8

setwd("~/Library/CloudStorage/OneDrive-UniversityofVermont/PENN STATE/resequencing data_GWAS/manhattan")
tiff("Reseq_gwas_HS_score_ggplot.tiff", width = 11, height = 7, units = 'in', res = 300)
ggplot(gwas_data, aes(x = bp_cum, y = -log10(p_wald), 
                                  color = as_factor(chr), size = -log10(p_wald))) +
  geom_hline(yintercept = -log10(sig), color = "grey40", linetype = "dashed") + 
  geom_point(alpha = 0.75) +
  scale_x_continuous(label = axis_set$chr, breaks = axis_set$center) +
  scale_y_continuous(expand = c(0,0), limits = c(0, ylim)) +
  scale_color_manual(values = rep(c("#276FBF", "#183059"), unique(length(axis_set$chr)))) +
  scale_size_continuous(range = c(0.5,3)) +
  labs(x = NULL, 
       y = "-log<sub>10</sub>(p)") + 
  theme_minimal() +
  theme( 
    legend.position = "none",
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    axis.title.y = element_markdown(),
    axis.text.x = element_text(angle = 60, size = 8, vjust = 0.5))
dev.off()
```
![alt text](https://github.com/aayudh454/Lasky-Morris-Lab-Sorghum-project/blob/main/manhattan_plot.png)

**8. Finding top SNPs**

```
setwd("~/Library/CloudStorage/OneDrive-UniversityofVermont/PENN STATE/eGWAS_revised list/manhattan plot")
list.files()

library(tidyverse)
library(ggtext)
library(normentR)
library(ggplot2)

Reseq_preds.all <- read.csv("1. Reseq_preds_all_chr10_corrected.csv", header=T)
head(Reseq_preds.all)


top_SNPs <- Reseq_preds.all[with(Reseq_preds.all,order(p_wald)),]
head(top_SNPs)

Reseq_topSNPS <- top_SNPs[1:1000,]
head(Reseq_topSNPS)
setwd("~/Library/CloudStorage/OneDrive-UniversityofVermont/PENN STATE/eGWAS_revised list/TOP SNPs")
write.csv(Reseq_topSNPS, file="Reseq_preds_all_top1000SNPS.csv", quote = T, eol = "\n", na= "NA")

Reseq_top100SNPS <- top_SNPs[1:100,]
head(Reseq_top100SNPS)

write.csv(Reseq_top100SNPS, file="Reseq_preds_all_top100SNPS.csv", quote = T, eol = "\n", na= "NA")
```

**9. Annotating top SNPs**

First find gff file with annotations from phytazome. For reseq data we used Sbicolor_454_v3.1.1.gene_exons.gff3 and Sbicolor_454_v3.1.1.annotation_info files https://data.jgi.doe.gov/refine-download/phytozome?organism=Sbicolor&expanded=454 

```
#!/usr/bin/env Rscript

#PBS -l nodes=1:ppn=8
#PBS -l walltime=12:00:00
#PBS -l pmem=24gb
#PBS -M azd6024@psu.edu
#PBS -A open
#PBS -j oe

setwd("~/work/preds_all_gwas/annotation")

#setwd("~/Library/CloudStorage/OneDrive-UniversityofVermont/PENN STATE/eGWAS_revised list/Annotation")
#list.files()
library(tidyverse)
#library(ggtext)
library(normentR)
library(dplyr)
library(fuzzyjoin)

Sorghum_annot <- read.delim("Sbicolor_454_v3.1.1.gene.gff3", header=T, comment.char="#")
head(Sorghum_annot)
names(Sorghum_annot)[9]<-paste("Gene_name")
names(Sorghum_annot)[4]<-paste("start")
names(Sorghum_annot)[5]<-paste("stop")
names(Sorghum_annot)[7]<-paste("direction")
names(Sorghum_annot)[1]<-paste("chr")
Sorghum_annot = within(Sorghum_annot, rm(.,..1,phytozomev12))
head(Sorghum_annot)

chr_new = substring(Sorghum_annot$chr, 4)
#replace column
Sorghum_annot[, "chr"] <- chr_new
head(Sorghum_annot)
dim(Sorghum_annot)

#Sorghum_annot[389067,1]
#Sorghum_annot[389068,1]

chr10_Sorghum_annot <- Sorghum_annot[Sorghum_annot$chr == '10',]
dim(chr10_Sorghum_annot)
chr1_9_Sorghum_annot <- Sorghum_annot[1:389067,]
chr1_9_Sorghum_annot_1 = substring(chr1_9_Sorghum_annot$chr, 2)
head(chr1_9_Sorghum_annot_1)

chr1_9_Sorghum_annot[, "chr"] <- chr1_9_Sorghum_annot_1
head(chr1_9_Sorghum_annot)
dim(chr1_9_Sorghum_annot)



Sbicolor_annot <- rbind(chr1_9_Sorghum_annot,chr10_Sorghum_annot)
head(Sbicolor_annot)
dim(Sbicolor_annot)

#adding these 2 colums to the dataframe will give you the +/- 5 kb window
Sbicolor_annot$start_5kb <- Sbicolor_annot$start - 5000
Sbicolor_annot$stop_5kb <- Sbicolor_annot$stop + 5000
head(Sbicolor_annot)


#setwd("~/Library/CloudStorage/OneDrive-UniversityofVermont/PENN STATE/eGWAS_revised list/Annotation")
#list.files()
top100SNPs <- read.csv("Reseq_preds_all_top100SNPS.csv", header=T)
top100SNPs = within(top100SNPs, rm(X))
head(top100SNPs)

library(dplyr)
library(fuzzyjoin)
gwas_annot <- fuzzy_left_join(top100SNPs, Sbicolor_annot, by=c("ps"="start_5kb", "ps"="stop_5kb", "chr"="chr"),
                            match_fun=list(`>=`, `<=`, `==`)) %>% select(Gene_name,ps,rs,p_wald,start,stop,direction)


gwas_annot_split <- data.frame(do.call("rbind", strsplit(as.character(gwas_annot$Gene_name), ";", fixed = TRUE)))
pacId <- sub('pacid=', '', gwas_annot_split$X3)
head(gwas_annot_split)

gwas_annot_1 <- cbind (pacId,gwas_annot)
gwas_annot_1 = within(gwas_annot_1, rm(Gene_name))
head(gwas_annot_1)
gwas_annot_2 <- gwas_annot_1[- grep("ancestorIdentifier", gwas_annot_1$pacId),]
head(gwas_annot_2)

gwas_annot_3 <- gwas_annot_2[- grep("ID=Sobic", gwas_annot_2$pacId),]
head(gwas_annot_3)

library(data.table)
gwas_annot_4 <- unique(setDT(gwas_annot_3)[order(ps, -ps)], by = "ps")
dim(gwas_annot_4)
head(gwas_annot_4)

write.table(gwas_annot_4, "Annotation_reseq_top100_5kb_preds_all.csv", sep=",")

```
Now download this and run it in R to get gene info

```
list.files()
library(tidyverse)
#library(ggtext)
library(normentR)
library(dplyr)
library(fuzzyjoin)

gwas_annot <- read.csv("Annotation_reseq_top100_5kb_preds_all.csv", header = T)
head(gwas_annot)

library(tidyr)
gwas_annot <- gwas_annot %>% drop_na()
dim(gwas_annot)
head(gwas_annot)


Sorghum_mdata <- read.csv("Sbicolor_454_v3.1.1.annotation_info.csv", header = T)
head(Sorghum_mdata)

annotation_Loc <- merge(Sorghum_mdata,gwas_annot, by = "pacId")
dim(annotation_Loc)
head(annotation_Loc)
annotation_Loc_1 <- annotation_Loc[order(annotation_Loc$p_wald),]

write.table(annotation_Loc_1, "1. Annotation_reseq_top100_preds_all.csv", sep=",")

```

<div id='id-section5'/>

## Chapter 5: Mapping SNPs and finding SNPs nearby

Based on your top hit based on annotation you can map this top SNPs

```
#!/usr/bin/env Rscript

#PBS -l nodes=1:ppn=8
#PBS -l walltime=12:00:00
#PBS -l pmem=24gb
#PBS -M azd6024@psu.edu
#PBS -A open
#PBS -j oe

setwd("/gpfs/group/jrl35/default/aayudh/gwas_reseq")

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
setwd("~/work/preds_all_gwas")
traits <- data.frame(read.csv('1_ReseqGWAS_traits_revissed_predsALL.csv', header=TRUE))
head(traits)

#no need to account for LD when getting SNPs for GWAS. 
mySNPs <- snpgdsSelectSNP(genofile, traits$LIB, maf = 0.05, missing.rate = 0)
length(mySNPs)

mySNPmatrix <- snpgdsGetGeno(genofile, sample.id = traits$LIB, snp.id = mySNPs, with.id = TRUE)
dim(mySNPmatrix$genotype)

#these lines create the bed file for gemma
rs <- data.frame(chr = chromosome.id[mySNPs], positions = snp.position[mySNPs], allele = snp.allele[mySNPs])
SNPs_bim <- data.frame(paste0('rs_', rs$chr, '_', rs$positions, ';', rs$allele), t(mySNPmatrix$genotype)) 
dim(SNPs_bim)
length(SNPs_bim)
SNPs_bim[1,1]

data_2 <- SNPs_bim[,2:340]
colnames(data_2) <- mySNPmatrix$sample.id
data_3 <- cbind(SNPs_bim[,1],data_2)
head(data_3)
dim(data_3)
data_3[1,1]
names(data_3)[1]<-paste("SNPs")

data_4 <- data.frame(do.call("rbind", strsplit(as.character(data_3$SNPs), ";", fixed = TRUE)))

data_5 <- cbind (data_4,data_3)
names(data_5)[1]<-paste("position")
names(data_5)[2]<-paste("allele")
data_5 = within(data_5, rm(SNPs))


data_49813813 <- subset(data_5, position=="rs_04_49813813")
data_new <- data_49813813[ , colSums(is.na(data_49813813)) < nrow(data_49813813)] 
data_snp_49813813 <- t(data_new)

SNP_49813813 <- as.data.frame(data_snp_49813813)
setwd("~/work/preds_all_gwas/map_snps")
write.table(SNP_49813813, "SNP_rs_04_49813813.csv", sep=",")
```
Now next part is in your R

```
#-------------------ggplot-part---------------------------------
setwd("~/Library/CloudStorage/OneDrive-UniversityofVermont/PENN STATE/eGWAS_revised list/map SNPs")
list.files()
SNP_49813813 <- read.csv("SNP_rs_04_49813813.csv")
head(SNP_49813813)

metadata <- read.csv("1. ReseqGWAS_traits_revissed_predsALL.csv")
head(metadata)

SNP_49813813_full<- merge(SNP_49813813,metadata, by="LIB")

SNP_49813813_full$ref_alt <- ifelse(SNP_49813813_full$rs_04_49813813_C.G=="2",'REF','ALT')
head(SNP_49813813_full)

##MAP
library(tidyverse)
library(sf)
library(mapview)
library(ggplot2)

SNP_49813813_full %>% 
  select(Accession, Lon, 
         Lat,ref_alt,population) %>%
  head()

world <- map_data("world")

setwd("~/Library/CloudStorage/OneDrive-UniversityofVermont/PENN STATE/eGWAS_revised list/map SNPs")
tiff("SNP_49813813.tiff", width = 6, height = 6, units = 'in', res = 300)
ggplot() + 
  geom_map(data = world, map = world,aes(long, lat, map_id = region),
           color = "black", fill = "white", size = 0.1) +
  geom_point(data = SNP_49813813_full,aes(Lon, Lat, color = ref_alt),alpha = 1) +
  coord_sf(xlim = c(-20, 50), ylim = c(-35, 35), expand = FALSE)+
  scale_color_manual(values=c('red', '#56B4E9'))+
  labs(x = NULL, y = NULL) 
dev.off()
```
**Density plot**
```
setwd("~/Library/CloudStorage/OneDrive-UniversityofVermont/PENN STATE/eGWAS_revised list/map SNPs")
list.files()
SNP_45529877 <- read.csv("SNP_45529877.csv")
head(SNP_45529877)

metadata <- read.csv("1. ReseqGWAS_traits_revissed_predsALL.csv")
head(metadata)

SNP_45529877_full<- merge(SNP_45529877,metadata, by="LIB")

SNP_45529877_full$ref_alt <- ifelse(SNP_45529877_full$rs_10_45529877_G.A=="2",'REF','ALT')
head(SNP_45529877_full)


library(raster)
library(rgdal)
library(classInt)
library(RColorBrewer)

##Here, lgs1 is a table, with the suitability variable ‘enm’ and the lgs1 alleles in ‘V4’. 
##I’m collapsing the two different deletion alleles into the purple bars.
### make allele freq plot ####
delz <- hist(SNP_45529877_full$HS_score[SNP_45529877_full$ref_alt %in% 'ALT'], breaks = seq(0, 1, by = 0.05), plot = F)$counts 
intz <- hist(SNP_45529877_full$HS_score[!SNP_45529877_full$ref_alt %in% 'ALT'], breaks = seq(0, 1, by = 0.05), plot = F)$counts 
sez <- sqrt(((delz / colSums(rbind(delz, intz))) * (intz / colSums(rbind(delz, intz)))) / colSums(rbind(delz, intz)))

setwd("~/Library/CloudStorage/OneDrive-UniversityofVermont/PENN STATE/eGWAS_revised list/allele_freq")
tiff("SNP_45529877.tiff", width = 6, height = 6, units = 'in', res = 300)
plot(seq(0, 0.95, by = 0.05) + 0.025, delz / colSums(rbind(delz, intz)))
plot.new()
plot.window(xlim = c(0, 0.9), ylim = c(0, 1))
for(i in 1:length(delz)){
  polygon(seq(0, 1, by = 0.05)[c(i,i+1, i+1, i)], c(0, 0, rep(intz[i]/(intz[i] + delz[i]), 2)), col = gray(0.7))
  polygon(seq(0, 1, by = 0.05)[c(i,i+1, i+1, i)], c(1, 1, rep(intz[i]/(intz[i] + delz[i]), 2)), col = brewer.pal(4, 'Set1')[4])
  lines(seq(0, 1, by = 0.05)[rep(i,2)] + 0.025, rep(intz[i]/(intz[i] + delz[i]), 2) + c(-sez[i], sez[i]))
  text(seq(0, 1, by = 0.05)[i] + 0.025, 0.05, colSums(rbind(delz, intz))[i], cex =0.7)
}
#standard error is sqrt(pq/n)
axis(1)
axis(2)
title(xlab = 'Striga habitat suitability')
title(ylab = 'Allele frequency')
dev.off()
```
![alt text](https://github.com/aayudh454/Lasky-Morris-Lab-Sorghum-project/blob/main/Mapping%20SNPs.png)

<div id='id-section6'/>

## Chapter 6: Finding nearby SNPs

Creating a vector of correlation coefficients for gwas resequencing snp vs everything else nearby snps in gbs set. We are looking to see if the resequencing snp is ‘Tagged’ by a gbs snp, so that the gbs snp could be used to infer accessions alleles at the resequencing gwas hit. Moreover, we want to know if your gwas snp is in strong LD with any gbs snp.

```
setwd("~/Library/CloudStorage/OneDrive-UniversityofVermont/PENN STATE/GBS data processing /Bigdata_cgr5_699")
library(SNPRelate)
library(gdsfmt)

genofile <- snpgdsOpen('SNPs_imp.recode.gds')
sample.id <- read.gdsn(index.gdsn(genofile, 'sample.id'))
snp.id <- read.gdsn(index.gdsn(genofile, 'snp.id'))
snp.position <- read.gdsn(index.gdsn(genofile, 'snp.position'))
chromosome.id <- read.gdsn(index.gdsn(genofile, 'snp.chromosome'))
snp.allele <- read.gdsn(index.gdsn(genofile, 'snp.allele'))
snp.rs.id <- read.gdsn(index.gdsn(genofile, 'snp.rs.id'))



data_sampleID <- as.data.frame(sample.id)
data_sampleID_1 <- data.frame(do.call("rbind", strsplit(as.character(data_sampleID$sample.id), ".", fixed = TRUE)))
data_sampleID_1 <- as.data.frame(data_sampleID_1$X1)
dim(data_sampleID_1)
names(data_sampleID_1)[1]<-paste("Accession")

gbs_sampleID <-cbind (data_sampleID,data_sampleID_1)

setwd("~/Library/CloudStorage/OneDrive-UniversityofVermont/PENN STATE/eGWAS_revised list/sample accessions ")
list.files()
Reseq_accessions <- read.csv("1. ReseqGWAS_traits_revissed_predsALL.csv", header = T)

Reseq_GBS_accessions<- merge(Reseq_accessions,gbs_sampleID, by="Accession")


mySNPs <- snpgdsSelectSNP(genofile, Reseq_GBS_accessions$sample.id, maf = 0.05, missing.rate = 0)
length(mySNPs)

mySNPmatrix <- snpgdsGetGeno(genofile, sample.id = Reseq_GBS_accessions$sample.id, snp.id = mySNPs, with.id = TRUE)
dim(mySNPmatrix$genotype)

rs <- data.frame(chr = chromosome.id[mySNPs], positions = snp.position[mySNPs], allele = snp.allele[mySNPs])
SNPs_bim <- data.frame(paste0('rs_', rs$chr, '_', rs$positions, ';', rs$allele), t(mySNPmatrix$genotype)) 
dim(SNPs_bim)
length(SNPs_bim)
SNPs_bim[1,1]


data_2 <- SNPs_bim[,2:240]

colnames(data_2) <- mySNPmatrix$sample.id

data_3 <- cbind(SNPs_bim[,1],data_2)
names(data_3)[1]<-paste("SNPs")
dim(data_3)
data_3[1:5,1:5]

data_4 <- data.frame(do.call("rbind", strsplit(as.character(data_3$SNPs), ";", fixed = TRUE)))
data_5 <- cbind (data_4,data_3)
data_5[1:5,1:5]
names(data_5)[2]<-paste("allele")
data_6 <- data.frame(do.call("rbind", strsplit(as.character(data_5$X1), "_", fixed = TRUE)))
data_7 <- cbind (data_6,data_5)
data_7[1:10,1:10]
names(data_7)[1]<-paste("rs")
names(data_7)[2]<-paste("chr")
names(data_7)[3]<-paste("ps")

data_7 = within(data_7, rm(rs,X1))
data_7[1:10,1:10]
dim(data_7)
data_chr4<- subset(data_7, chr=="4") ##DO == "4"
dim(data_chr4)
data_chr4[1:10,1:10]

setwd("~/Library/CloudStorage/OneDrive-UniversityofVermont/PENN STATE/eGWAS_revised list/compare SNP in gbs")
write.table(data_chr4, "gbs_chr4_SNPs.csv", sep=",")

library(tidyverse)
#library(ggtext)
library(normentR)
library(dplyr)
library(fuzzyjoin)


new_vector <- as.numeric(as.character(data_chr4$ps))
#adding these 2 colums to the dataframe will give you the +/- 5 kb window
start_5kb <- as.data.frame(new_vector-5000)
stop_5kb <- as.data.frame(new_vector + 5000)

data8 <- cbind(start_5kb,stop_5kb)
names(data8)[1]<-paste("start_5kb")
names(data8)[2]<-paste("stop_5kb")

chr4_final <- cbind (data8,data_chr4)
chr4_final[1:10,1:10]


setwd("~/Library/CloudStorage/OneDrive-UniversityofVermont/PENN STATE/eGWAS_revised list/Annotation")
list.files()
top100SNPs <- read.csv("Reseq_preds_all_top100SNPS.csv", header=T)
head(top100SNPs)

rs_04_49813813 <- top100SNPs[1,]
names(rs_04_49813813)[3]<-paste("position")

library(dplyr)
library(fuzzyjoin)
nerarbySNP <- fuzzy_left_join(rs_04_49813813, chr4_final, by=c("position"="start_5kb", "position"="stop_5kb", "chr"="chr"), 
                              match_fun=list(`>=`, `<=`, `==`)) %>% select(SNPs,position,allele,ps)
head(nerarbySNP)
dim(nerarbySNP)

setwd("~/Library/CloudStorage/OneDrive-UniversityofVermont/PENN STATE/eGWAS_revised list/compare SNP in gbs")
write.table(nerarbySNP, "rs_04_49813813_nerarbySNP_accession_match.csv", sep=",")
```

| SNPs in GBS dataset  | allele |
| -------------------- | ------------- |
| rs_4_49813787;C/T  | C/T  |
| rs_4_49813804;T/G  | T/G  |

<div id='id-section7'/>

## Chapter 7: GWAS using vcftools



**Installing vcftools**

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
