# GWAS and variant callings

## Table of contents    
* [Page 1: 2022-25-06](#id-section1). Chapter 1: Getting resequencing data 

* [Page 2: 2022-30-06](#id-section2). Chapter 2: Running Beagle 

* [Page 3 2022-30-06](#id-section3). Chapter 3: SnpEff

* [Page 4 2022-08-08](#id-section4). Chapter 4: Environmental GWAS 

* [Page 5 2022-17-08](#id-section5). Chapter 5: Mapping SNPs 

* [Page 6 2022-22-08](#id-section6). Chapeter 6: Finding nearby SNPs 

* [Page 7 2022-02-08](#id-section7). Chapeter 7: Indels association 

* [Page 8 2022-12-09](#id-section8). Chapeter 8: GWAS using vcftools

* [Page 9 2022-12-10](#id-section9). Chapter 9: Tajima's D 

* [Page 10 2022-24-10](#id-section10). Chapter 10: Fst estimation of Sorghum Reseq 

* [Page 11 2022-21-11](#id-section11). Chapter11: LD decay

* [Page 12 2022-19-12](#id-section12). Chapter12: Putative promoter analysis

* [Page 13 2023-19-01](#id-section13). Chapter13: ABA colocalization

* [Page 14 2023-09-02](#id-section14). Chapter14: Sampling (closely related but high HS difference)

* [Page 15 2023-23-02](#id-section15). Chapter15: Redundancy analysis (RDA) for drought

* [Page 16 2023-23-03](#id-section16). Chapter16: Enrichment analysis

* [Page 17 2023-28-03](#id-section17). Chapter17: Genic Vs Intergenic enrichment 

* [Page 18 2023-5-05](#id-section18). Chapter18: Python based eGWAS annotation pipeline
  
* [Page 19 2023-23-07](#id-section19). Chapter19: Analysing individual candidate

* [Page 20 2024-08-04](#id-section20). Chapter20: snpeff for Bromus tectorum
------
<div id='id-section1'/>

## Chapter 1: Getting resequencing data (by Aayudh)

#### Getting to roar collab
ssh azd6024@submit.hpc.psu.edu

my folder: cd /storage/group/jrl35/default/aayudh/custom_genome_snpeff

#Getting resequencing data
Login info: ssh -l azd6024 submit.aci.ics.psu.edu

### Data details
Here is a vcf for variants (SNP, MNP, indel independently) for all the resequencing libraries of Sorghum. This includes
Sorghum_1757g_AllChr_1757g.tar.gz
Sorghum_all_1757g_metadata.txt

1. Code to wget
```
wget --user xyz --password xyz https://hagsc.org/restricted_access/sorg_wga/Sorghum_1757g_AllChr_1757g.tar.gz
```
2. Unzip data

```
tar -zxvf Sorghum_1757g_AllChr_1757g.tar.gz
```

**DATA STORAGE IN LASKY LAB**

```
/gpfs/group/jrl35/default/aayudh
```
use recent R version in ICDS
```
module use /gpfs/group/RISE/sw7/modules
module load r/4.1.2
```

Delete a job
```
qdel 38423184 
```

#### Subset a cvf.gz with bcftools

```
bcftools view -r chr1 HG01148_ERR022469.GRCh38.variants.vcf.gz -Oz -o HG01148_ERR022469.GRCh38.chr1.variants.vcf.gz

bcftools query -f '%CHROM,%POS,%REF,%ALT\n' HG01148_ERR022469.GRCh38.chr1.variants.vcf.gz > chr1_variants.csv
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

## Chapter 3: SnpEff 

**Introduction**
SnpEff is an open source tool that annotates variants and predicts their effects on genes by using an interval forest approach. This program takes pre-determined variants listed in a data file that contains the nucleotide change and its position and predicts if the variants are deleterious.

**Rationale**
Running SnpEff provides predicted effect sizes for the observed SNPs in the variant dataset

###switchgrass

1.  unzip the vcf.gz
2. Make database for your ref genome file

```
[azd6024@submit03 switchgrass]$ cp BT_assembly.fa BT_assembly
[azd6024@submit03 switchgrass]$ java -jar snpEff/snpEff.jar build -gff3 -v BT_assembly
```
3. Now snpeff-

```
#!/bin/bash

#PBS -l nodes=1:ppn=8
#PBS -l walltime=12:00:00
#PBS -l pmem=24gb
#PBS -M azd6024@psu.edu
#PBS -A open
#PBS -j oe
#PBS -m abe

# Define variables
WORKING_DIR="/storage/group/jrl35/default/aayudh/switchgrass"
GENOME_FA="BT_assembly.fa"
GENOME_VERSION="BT_assembly"  # Custom name for your genome version
SNPEFF_DIR="$WORKING_DIR/snpEff" # Path to SnpEff installation
SNPEFF_CONFIG="$SNPEFF_DIR/snpEff.config"
VCF_FILE="307BRTE.imputed_filtered.vcf"      # Specified VCF file


# Step 4: Run SnpEff for Variant Annotation
echo "Running SnpEff for variant annotation..."
java -Xmx24g -jar $SNPEFF_DIR/snpEff.jar BT_assembly $VCF_FILE > ${VCF_FILE}.ann.vcf

echo "SnpEff annotation completed."

``` 

```
java -jar snpEff.jar build -gtf22 -v BT_assembly
```




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

NOTE: You need to match the genofile$sample.id format with the metadata column.

```
setwd("~/Library/CloudStorage/OneDrive-UniversityofVermont/PENN STATE/eGWAS_revised list/sample accessions ")
list.files()
Fanna_list <- read.csv("Sorghum_metadata_Coordinates_set.csv")
head(Fanna_list)

Reseq_matadata <- read.csv("Sorghum_all_1757g_metadata.csv")
head(Reseq_matadata)
final_data_1 <- gsub("_","",Reseq_matadata$sample)
head(final_data_1)
Reseq_matadata_1 <- cbind (Reseq_matadata,final_data_1)
head(Reseq_matadata_1)
names(Reseq_matadata_1)[names(Reseq_matadata_1) == 'final_data_1'] <- 'Accession'
head(Reseq_matadata_1)

Reseq_accessions<- merge(Reseq_matadata_1,Fanna_list, by="Accession")
dim(Reseq_accessions)
head(Reseq_accessions)
library(data.table)
Reseq_accessions_1 <- unique(setDT(Reseq_accessions)[order(Accession, -Accession)], by = "Accession")
dim(Reseq_accessions_1)
head(Reseq_accessions_1)

library(tidyr)
MASTER_data_final <- Reseq_accessions_1  %>% drop_na()
dim(MASTER_data_final)
head(MASTER_data_final)

library(raster)
library(sp)
library(rgdal)
library(tidyverse)

setwd("~/OneDrive - University of Vermont/PENN STATE/RAstor data")
list.files()

preds.sorg <- raster(paste0("~/OneDrive - University of Vermont/PENN STATE/RAstor data/preds.sorg.tif"))
preds.sorg

HS_score <- raster::extract(preds.sorg, MASTER_data_final[,c("Lon","Lat")])
ReseqGWAS_traits <- cbind(MASTER_data_final,HS_score)
head(ReseqGWAS_traits)
dim(ReseqGWAS_traits)

setwd("~/Library/CloudStorage/OneDrive-UniversityofVermont/PENN STATE/eGWAS_revised list/sample accessions ")
#write.table(ReseqGWAS_traits, "ReseqGWAS_traits_with_NA.csv", sep=",")


library(tidyr)
ReseqGWAS_traits_1 <- ReseqGWAS_traits %>% drop_na()
dim(ReseqGWAS_traits_1)
head(ReseqGWAS_traits_1)

ReseqGWAS_traits_2 <- unique(setDT(ReseqGWAS_traits_1)[order(Accession, -Accession)], by = "Accession")
dim(ReseqGWAS_traits_2)
head(ReseqGWAS_traits_2)

setwd("~/Library/CloudStorage/OneDrive-UniversityofVermont/PENN STATE/eGWAS_revised list/sample accessions ")
write.table(ReseqGWAS_traits_2, "1_ReseqGWAS_traits_revissed_preds_sorg.csv", sep=",")
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
mySNPs <- snpgdsSelectSNP(genofile, traits$LIB, maf = 0.05, missing.rate = 0.5)
length(mySNPs)

mySNPmatrix <- snpgdsGetGeno(genofile, sample.id = traits$LIB, snp.id = mySNPs, with.id = TRUE)
dim(mySNPmatrix$genotype)

#these lines create the bed file for gemma
rs <- data.frame(chr = chromosome.id[mySNPs], positions = snp.position[mySNPs], allele = snp.allele[mySNPs])
SNPs_bim <- data.frame(paste0('rs_', rs$chr, '_', rs$positions, '_', rs$allele), 'A', 'T', t(mySNPmatrix$genotype)) 
write.table(SNPs_bim, 'mySNPbed.bed', row.names = FALSE, quote = FALSE, col.names = FALSE, sep = ',') 

#kinship matrix for gemma obtained with 'identity by state'
ibs <- snpgdsIBS(genofile, sample.id=traits$LIB, snp.id=mySNPs, maf = 0.05, missing.rate = 0.5)
ibs.matrix <- ibs$ibs
dim(ibs.matrix)

setwd("~/work/preds_all_gwas")
write.table(ibs.matrix, 'ibs_matrix.txt', row.names = FALSE, col.names = FALSE, quote = FALSE, sep = ', ')

#prepare the trait data for gemma, accessions should be in the same order as in the genetic matrices
order <- as.data.frame(ibs$sample.id)
dim(order)
colnames(order) <- c("sample.id")

order.idx <- match(order$sample.id, traits$LIB)
order.idx

ordered <- traits[order.idx,]

setwd("~/work/preds_all_gwas")
write.table(ordered$HS_score, 'HS_score_BimBam.txt', row.names = FALSE, quote = FALSE, col.names = FALSE)
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
#PBS -A open
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
#setwd("~/Library/CloudStorage/OneDrive-UniversityofVermont/PENN STATE/eGWAS_revised list/clean_gemma")
gemma_output <- read.table("Reseq_preds_all_HS_score_out.assoc.txt", header = TRUE, sep = "", dec = ".")
head(gemma_output)
tail(gemma_output)
dim(gemma_output)

rs <- gemma_output[,2]
rs_split <- data.frame(do.call("rbind", strsplit(as.character(gemma_output$rs), "_", fixed = TRUE)))
head(rs_split)

gemma_output[, "chr"] <- rs_split$X2
gemma_output[, "ps"] <- rs_split$X3
allele <- rs_split$X4

gemma_output_clean <- cbind(gemma_output,allele)
head(gemma_output_clean)

write.table(gemma_output_clean, 'Reseq_preds_all_HS_score_out.assoc.clean.txt', row.names = FALSE, quote = FALSE, col.names = TRUE, sep = '\t')

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
setwd("~/Library/CloudStorage/OneDrive-UniversityofVermont/PENN STATE/eGWAS_revised list/manhattan plot")
list.files()

library(tidyverse)
library(ggtext)
library(normentR)
library(ggplot2)

Chr10corr <- read.csv("Reseq_preds_all_chr10_corrected.csv", header=T)
head(Chr10corr)

sig_data <- Chr10corr %>% 
  subset(p_wald < 0.05)
notsig_data <- Chr10corr %>% 
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

##Fdr line
Chr10corr <- read.csv("Reseq_gwas_chr10_corrected.csv", header=T)
head(Chr10corr)
Chr10corr_sub <- subset(Chr10corr, p_FDR < 0.05)
Chr10corr_sub_1 <- Chr10corr_sub[order(Chr10corr_sub$p_FDR),]

fdr_line <- Chr10corr_sub_1[166,4]


setwd("~/Library/CloudStorage/OneDrive-UniversityofVermont/PENN STATE/eGWAS_revised list/manhattan plot")
tiff("Reseq_preds_all_HS_score_ggplot.tiff", width = 11, height = 7, units = 'in', res = 300)
ggplot(gwas_data, aes(x = bp_cum, y = -log10(p_wald), 
                      color = as_factor(chr), size = -log10(p_wald))) +
  geom_hline(yintercept = -log10(fdr_line), color = "red", linetype = "dashed") +
  geom_point(alpha = 0.75) +
  scale_x_continuous(label = axis_set$chr, breaks = axis_set$center) +
  scale_y_continuous(expand = c(0,0), limits = c(0, 12)) +
  scale_color_manual(values = rep(c("#276FBF", "#183059"), unique(length(axis_set$chr)))) +
  scale_size_continuous(range = c(0.5,3)) +
  labs(x = NULL, 
       y = "-log<sub>10</sub>(p)") + 
  theme_minimal() +
  theme( 
    legend.position = "none",text=element_text(size=18, colour = "black", family="Times"),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    axis.title.y = element_markdown(),
    axis.text.x = element_text(size = 18, vjust = 0.5))
dev.off()

png(file = 'Reseq_preds_all_HS_score_ggplot.png', width=1400, height=960, res=300)
ggplot(gwas_data, aes(x = bp_cum, y = -log10(p_wald), 
                      color = as_factor(chr), size = -log10(p_wald))) +
  geom_hline(yintercept = -log10(fdr_line), color = "red", linetype = "dashed") +
  geom_point(alpha = 0.75) +
  scale_x_continuous(label = axis_set$chr, breaks = axis_set$center) +
  scale_y_continuous(expand = c(0,0), limits = c(0, 12)) +
  scale_color_manual(values = rep(c("#276FBF", "#183059"), unique(length(axis_set$chr)))) +
  scale_size_continuous(range = c(0.5,3)) +
  labs(x = NULL, 
       y = "-log<sub>10</sub>(p)") + 
  theme_minimal() +
  theme( 
    legend.position = "none",text=element_text(size=18, colour = "black", family="Times"),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    axis.title.y = element_markdown(),
    axis.text.x = element_text(size = 18, vjust = 0.5))
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
                            match_fun=list(`>=`, `<=`, `==`)) %>% select(Gene_name,ps,rs,af,p_wald,start,stop,direction)


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
-----
<div id='id-section5'/>

## Chapter 5: Mapping SNPs 

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
setwd("~/Library/CloudStorage/OneDrive-UniversityofVermont/PENN STATE/cycles_predicted/candidates")
list.files()
##expressed protein_SNP_rs_09_53552965
SNP_53552965 <- read.csv("SNP_rs_09_53552965_cyclin.csv")
head(SNP_53552965)

setwd("~/Library/CloudStorage/OneDrive-UniversityofVermont/PENN STATE/cycles_predicted")
metadata <- read.csv("1. Cycles_new.csv")
head(metadata)

SNP_53552965_full<- merge(SNP_53552965,metadata, by="LIB")
head(SNP_53552965_full)

SNP_53552965_full$ref_alt <- ifelse(SNP_53552965_full$rs_53552965=="2",'REF','ALT')
head(SNP_53552965_full)

##MAP
library(tidyverse)
library(sf)
library(mapview)
library(ggplot2)

world <- map_data("world")
mid<-mean(SNP_53552965_full$water_veg_avg_log)
setwd("~/Library/CloudStorage/OneDrive-UniversityofVermont/PENN STATE/cycles_predicted/candidates")
tiff("cyclin_legend.tiff", width = 6, height = 6, units = 'in', res = 300)
ggplot() + 
  geom_map(data = world, map = world,aes(long, lat, map_id = region),
           color = "black", fill = "white", size = 0.1) +
  geom_point(data = SNP_53552965_full,aes(Lon, Lat, color = water_veg_avg_log,shape = ref_alt),alpha = 1) +
  theme_classic() +
  coord_sf(xlim = c(-20, 50), ylim = c(-35, 35), expand = FALSE)+
  scale_color_gradient2(midpoint=mid, low="blue", mid="white", high="red", space ="Lab" )+
  scale_shape_manual(values=c(19, 2))+
  labs(x = NULL, y = NULL) +
  theme(legend.position="right",text=element_text(size=16, colour = "black", family="Times New Roman"),
        axis.line = element_line(size=0.5, colour = "black"),
        axis.text.x=element_text(colour="black", size = 16),
        axis.text.y=element_text(colour="black", size = 16)) 
dev.off()
```
**For RDA the REF and ALT**

```
grep -w '5569009' Chr05.txt > rs_5_5569009.csv
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

-----
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

**Now doing the correlations**

```
setwd("~/Library/CloudStorage/OneDrive-UniversityofVermont/PENN STATE/eGWAS_revised list/compare SNP in gbs")
list.files()

gbs_chr4_SNPs <- read.csv("gbs_chr4_SNPs.csv", header = T)
head(gbs_chr4_SNPs)

nearby_SNP <- read.csv("1.rs_04_49813813_nerarbySNP_accession_match.csv", header = T)
head(nearby_SNP)

matched_SNPs <- merge (nearby_SNP,gbs_chr4_SNPs, by = "SNPs")
matched_SNPs  = within(matched_SNPs , rm(allele.x,chr,ps,allele.y))

gbs_SNPs <- t(matched_SNPs)
gbs_SNPs <- as.data.frame(gbs_SNPs)

colnames(gbs_SNPs) <- matched_SNPs$SNPs

gbs_SNPs_1 <- gbs_SNPs [2:240,]
head(gbs_SNPs_1)

#write.table(gbs_SNPs_1, "gbs_SNPs_1.csv", sep=",")

data <- read.csv("gbs_SNPs_1.csv")

gbs_SNPs_2 <- data.frame(do.call("rbind", strsplit(as.character(data$Accession), ".", fixed = TRUE)))
gbs_SNPs_3 <- as.data.frame(gbs_SNPs_2$X1)
dim(gbs_SNPs_3)
names(gbs_SNPs_3)[1]<-paste("Accession")

final_data <- cbind(gbs_SNPs_3,gbs_SNPs_1)
head(final_data)

#write.table(final_data, "gbs_SNPs_2.csv", sep=",")
setwd("~/Library/CloudStorage/OneDrive-UniversityofVermont/PENN STATE/eGWAS_revised list/compare SNP in gbs")
#delete first col in excel
df3 <- read.csv ("gbs_SNPs_2.csv")
head(df3)

df1 <- read.csv("1. ReseqGWAS_traits_revissed_predsALL.csv")
df2 <- read.csv("SNP_rs_04_49813813.csv")
df <- merge(df1,df2, by ="LIB")
df  = within(df , rm(LIB,sample,population,Lat,Lon,HS_score))
head(df)


merged_data <- merge(df,df3,by = "Accession")
head(merged_data)
dim(merged_data)

matrix_format <- merged_data[,2:15]

data_cor <- cor(matrix_format[ , colnames(matrix_format) != "x1"],  # Calculate correlations
                matrix_format$x1)
data_cor 
dim(data_cor)

results <- as.data.frame(data_cor[2:14,1])
names(results)[1]<-paste("rs_04_49813813_C.G")

#write.table(results, "1.correlation_rs_04_49813813_C.G.csv", sep=",")

tiff("nearbySNPs.tiff", width = 10, height = 6, units = 'in', res = 300)
heatmap(t(`row.names<-`(as.matrix(merged_data[-1]), merged_data$Accession)))
dev.off()
```

| SNPs in GBS dataset  | Reseq_rs_4_49813813_C.G (cor) |
| -------------------- | ------------- |
| rs_4_49813787;C/T  | 0.600259327  |
| rs_4_49813804;T/G  | 0.752821624  |

*cor=Pearson correlation coefficient*

-----
<div id='id-section7'/>

## Chapter 7: Indel asscoiation 

unzip indel.gz first

```
#!/bin/bash

#PBS -l nodes=1:ppn=8
#PBS -l walltime=24:00:00
#PBS -l pmem=48gb
#PBS -M azd6024@psu.edu
#PBS -A open
#PBS -j oe

WORKINGDIR=/gpfs/group/jrl35/default/aayudh/results_beagle

cd $WORKINGDIR

gunzip Sorghum_1757g_AllChr.polymorphic.indel.noRepeats.5pctMasked.imputed.vcf.gz
```

**Create gds file**: method = "copy.num.of.ref": to extract and store dosage (0, 1, 2) of the reference allele for all variant sites, including bi-allelic SNPs, multi-allelic SNPs, **indels** and structural variants.

```
#!/usr/bin/env Rscript

#PBS -l nodes=1:ppn=8
#PBS -l walltime=12:00:00
#PBS -l pmem=24gb
#PBS -M azd6024@psu.edu
#PBS -A open
#PBS -j oe

setwd("/gpfs/group/jrl35/default/aayudh/results_beagle")

library(gdsfmt)
library(SNPRelate)

vcf.fn <- "/gpfs/group/jrl35/default/aayudh/results_beagle/Sorghum_1757g_AllChr.polymorphic.indel.noRepeats.5pctMasked.imputed.vcf"

snpgdsVCF2GDS(vcf.fn, "Sorghum_Reseq_INDELS.gds", method="copy.num.of.ref")
```
The total number of samples: 1757 ; The total number of SNPs: 4292419 ; SNP genotypes are stored in SNP-major mode (Sample X SNP).

### Association of INDELS using GEMMA

```
setwd("~/Library/CloudStorage/OneDrive-UniversityofVermont/PENN STATE/eGWAS_revised list/indels")
list.files()

library(gdsfmt)
library(SNPRelate)

genofile <- snpgdsOpen('Sorghum_Reseq_INDELS.gds')
sample.id <- read.gdsn(index.gdsn(genofile, 'sample.id'))
chromosome.id <- read.gdsn(index.gdsn(genofile, 'snp.chromosome'))
snp.position <- read.gdsn(index.gdsn(genofile, 'snp.position'))
snp.allele <- read.gdsn(index.gdsn(genofile, 'snp.allele'))

snpgdsSummary('Sorghum_Reseq_INDELS.gds')

traits <- data.frame(read.csv('1. ReseqGWAS_traits_revissed_predsALL.csv', header=TRUE))
head(traits)

mySNPs <- snpgdsSelectSNP(genofile, traits$LIB, maf = 0.05, missing.rate = 0.5)
length(mySNPs)

mySNPmatrix <- snpgdsGetGeno(genofile, sample.id = traits$LIB, snp.id = mySNPs, with.id = TRUE)
dim(mySNPmatrix$genotype)

#these lines create the bed file for gemma
rs <- data.frame(chr = chromosome.id[mySNPs], positions = snp.position[mySNPs],allele = snp.allele[mySNPs])
SNPs_bim <- data.frame(paste0('rs_', rs$chr, '_', rs$positions, ';', rs$allele), 'A', 'T', t(mySNPmatrix$genotype)) 
dim(SNPs_bim)
SNPs_bim[1:5.1:5]

setwd("~/work/sorgh.preds_gwas")
write.table(SNPs_bim, 'mySNPbed.bed', row.names = FALSE, quote = FALSE, col.names = FALSE, sep = ',') 
```
Now take the ibs.matrix of preds_all SNP and HS_score_BimBam.txt to a same folder along with .bed file. You can copy **gemma-0.98.5-linux-static** there as well. Then run the script- 

```
#!/bin/bash

#PBS -l nodes=1:ppn=8
#PBS -l walltime=12:00:00
#PBS -l pmem=24gb
#PBS -M azd6024@psu.edu
#PBS -A open
#PBS -j oe

WORKINGDIR=/storage/home/azd6024/work/preds_all_gwas/indels
cd $WORKINGDIR

./gemma-0.98.5-linux-static -g mySNPbed.bed -k ibs_matrix.txt -lmm 4 -miss 0.1 -p HS_score_BimBam.txt -o Reseq_preds_all_indels_out
```
Download the output and clean it-
```
###clean gemma output
setwd("~/Library/CloudStorage/OneDrive-UniversityofVermont/PENN STATE/eGWAS_revised list/indels")
list.files()
gemma_output <- read.table("Reseq_preds_all_indels_out.assoc.txt", header = TRUE, sep = "", dec = ".")
head(gemma_output)
tail(gemma_output)
dim(gemma_output)

rs <- gemma_output[,2]
rs_split <- data.frame(do.call("rbind", strsplit(as.character(gemma_output$rs), "_", fixed = TRUE)))
rs_allele <- data.frame(do.call("rbind", strsplit(as.character(rs_split$X3), ";", fixed = TRUE)))
rs_noallele <- data.frame(do.call("rbind", strsplit(as.character(gemma_output$rs), ";", fixed = TRUE)))
head(rs_split)
head(rs_allele)
head(rs_noallele)

gemma_output[, "chr"] <- rs_split$X2
gemma_output[, "ps"] <- rs_allele$X1
gemma_output[, "rs"] <- rs_noallele$X1
allele <- rs_allele$X2

gemma_output_clean <- cbind(gemma_output,allele)
head(gemma_output_clean)
dim(gemma_output_clean)

require(dplyr)
##Move column
gemma_output_final <- gemma_output_clean %>% relocate(allele, .before = n_miss)
head(gemma_output_final)
dim(gemma_output_final)

setwd("~/Library/CloudStorage/OneDrive-UniversityofVermont/PENN STATE/eGWAS_revised list/indels")
#write.table(gemma_output_final, 'Reseq_preds_all_indels_out.assoc_clean.txt', row.names = FALSE, quote = FALSE, col.names = TRUE, sep = '\t')

gemma_output_final <- read.table("Reseq_preds_all_indels_out.assoc_clean.txt", header = TRUE, sep = "", dec = ".")

top_SNPs <- gemma_output_final[with(gemma_output_final,order(p_wald)),]
head(top_SNPs)

Reseq_topSNPS <- top_SNPs[1:100,]
head(Reseq_topSNPS)

write.csv(Reseq_topSNPS, file="Reseq_INDELS_preds_all_top100SNPS.csv", quote = T, eol = "\n", na= "NA")
```

### INDELS Annotations 

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

#setwd("~/Library/CloudStorage/OneDrive-UniversityofVermont/PENN STATE/eGWAS_revised list/indels")
setwd("~/work/preds_all_gwas/indels")
top100SNPs <- read.csv("Reseq_INDELS_preds_all_top100SNPS.csv", header=T)
head(top100SNPs)

library(dplyr)
library(fuzzyjoin)
gwas_annot <- fuzzy_left_join(top100SNPs, Sbicolor_annot, by=c("ps"="start_5kb", "ps"="stop_5kb", "chr"="chr"),
                              match_fun=list(`>=`, `<=`, `==`)) %>% select(Gene_name,ps,rs,af,p_wald,start,stop,direction)

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

write.table(gwas_annot_4, "Annotation_reseq_INDELS_top100_5kb_preds_all.csv", sep=",")
```
Now do in R

```
#-------R part----------------------
setwd("~/Library/CloudStorage/OneDrive-UniversityofVermont/PENN STATE/eGWAS_revised list/indels")
gwas_annot <- read.csv("Annotation_reseq_INDELS_top100_5kb_preds_all.csv", header = T)
head(gwas_annot)

library(tidyr)
gwas_annot <- gwas_annot %>% drop_na()
dim(gwas_annot)
head(gwas_annot)

setwd("~/Library/CloudStorage/OneDrive-UniversityofVermont/PENN STATE/eGWAS_revised list/Annotation")
Sorghum_mdata <- read.csv("Sbicolor_454_v3.1.1.annotation_info.csv", header = T)
head(Sorghum_mdata)

annotation_Loc <- merge(Sorghum_mdata,gwas_annot, by = "pacId")
dim(annotation_Loc)
head(annotation_Loc)
annotation_Loc_1 <- annotation_Loc[order(annotation_Loc$p_wald),]

setwd("~/Library/CloudStorage/OneDrive-UniversityofVermont/PENN STATE/eGWAS_revised list/indels")
write.table(annotation_Loc_1, "1. Annotation_INDELS_top100_preds_all.csv", sep=",")
```
### INDELS-SNP correlation

```
setwd("~/Library/CloudStorage/OneDrive-UniversityofVermont/PENN STATE/eGWAS_revised list/indels")
list.files()

library(gdsfmt)
library(SNPRelate)

genofile <- snpgdsOpen('Sorghum_Reseq_INDELS.gds')
sample.id <- read.gdsn(index.gdsn(genofile, 'sample.id'))
chromosome.id <- read.gdsn(index.gdsn(genofile, 'snp.chromosome'))
snp.position <- read.gdsn(index.gdsn(genofile, 'snp.position'))
snp.allele <- read.gdsn(index.gdsn(genofile, 'snp.allele'))

snpgdsSummary('Sorghum_Reseq_INDELS.gds')

traits <- data.frame(read.csv('1. ReseqGWAS_traits_revissed_predsALL.csv', header=TRUE))
head(traits)

mySNPs <- snpgdsSelectSNP(genofile, traits$LIB, maf = 0.05, missing.rate = 0.5)
length(mySNPs)

mySNPmatrix <- snpgdsGetGeno(genofile, sample.id = traits$LIB, snp.id = mySNPs, with.id = TRUE)
dim(mySNPmatrix$genotype)

#these lines create the bed file for gemma
rs <- data.frame(chr = chromosome.id[mySNPs], positions = snp.position[mySNPs],allele = snp.allele[mySNPs])
SNPs_bim <- data.frame(paste0('rs_', rs$chr, '_', rs$positions, ';', rs$allele), t(mySNPmatrix$genotype)) 
dim(SNPs_bim)
SNPs_bim[1:5.1:5]

data_2 <- SNPs_bim[,2:340]
colnames(data_2) <- mySNPmatrix$sample.id
data_3 <- cbind(SNPs_bim[,1],data_2)
head(data_3)
dim(data_3)
data_3[1,1]
names(data_3)[1]<-paste("Indels")

data_4 <- data.frame(do.call("rbind", strsplit(as.character(data_3$Indels), ";", fixed = TRUE)))

data_5 <- cbind (data_4,data_3)
names(data_5)[1]<-paste("position")
names(data_5)[2]<-paste("allele")
data_5 = within(data_5, rm(Indels))


data_49813789 <- subset(data_5, position=="rs_04_49813789")
data_new <- data_49813789[ , colSums(is.na(data_49813789)) < nrow(data_49813789)] 
data_indel_49813789 <- t(data_new)

Indels_49813789 <- as.data.frame(data_indel_49813789)
setwd("~/Library/CloudStorage/OneDrive-UniversityofVermont/PENN STATE/eGWAS_revised list/indel_snp_corr")
write.table(Indels_49813789, "INDEL_rs_04_49813789.csv", sep=",")

#-----------------------------------------------------------------SNPvsIndel------------------
setwd("~/Library/CloudStorage/OneDrive-UniversityofVermont/PENN STATE/eGWAS_revised list/indels")
list.files()
SNP <- read.csv("SNP_rs_04_49813813.csv")
INDEL <- read.csv("INDEL_rs_04_49813789.csv")

SNP_INDEL_corr<- merge(SNP,INDEL, by="LIB")
head(SNP_INDEL_corr)

cor.test(SNP_INDEL_corr$SNP_04_49813813_C.G,SNP_INDEL_corr$INDEL_04_49813789_TA.T, method="pearson")

setwd("~/Library/CloudStorage/OneDrive-UniversityofVermont/PENN STATE/eGWAS_revised list/indels")
tiff("SNP_INDEL_corr.tiff",width = 5, height = 5, units = 'in', res = 300)
par(family="Times")
par(pty="s")
plot(SNP_INDEL_corr$SNP_04_49813813_C.G,SNP_INDEL_corr$INDEL_04_49813789_TA.T, 
     pch = 19, col="blue",xaxt="n",yaxt="n",ylab = "INDEL:Chr4-49813789 (TA/T)",xlab = "SNP:Chr4-49813813 (C/G)",font = 1, las=1)
abline(lm(SNP_INDEL_corr$INDEL_04_49813789_TA.T~SNP_INDEL_corr$SNP_04_49813813_C.G), 
       col = "black", lty = 2, lwd=2)
axis(side = 1, at=seq(min(SNP_INDEL_corr$SNP_04_49813813_C.G), max(SNP_INDEL_corr$SNP_04_49813813_C.G), by = ((-min(SNP_INDEL_corr$SNP_04_49813813_C.G)+max(SNP_INDEL_corr$SNP_04_49813813_C.G))/2)),las=1,cex.axis = 1.5)
axis(side = 2, at=seq(min(SNP_INDEL_corr$INDEL_04_49813789_TA.T), max(SNP_INDEL_corr$INDEL_04_49813789_TA.T), by = ((-min(SNP_INDEL_corr$INDEL_04_49813789_TA.T)+max(SNP_INDEL_corr$INDEL_04_49813789_TA.T))/2)),las=1,cex.axis = 1.5)
dev.off()

library("ggpubr")
ggscatter(SNP_INDEL_corr, x = "SNP_04_49813813_C.G", y = "INDEL_04_49813789_TA.T", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "pearson",
          xlab = "SNP_Chr4_49813813_C/G", ylab = "INDEL_Chr4_49813789_TA/T")
```

-----
<div id='id-section8'/>

## Chapter 8: GWAS using vcftools


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

see what's in it-
```
zcat reseq_SNPeff_variant_call_SNPs.vcf.gz | head -n 11
```

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

-----
<div id='id-section9'/>

## Chapter 9: Tajima's D

### Tajima's D for reseq Sorghum dataset

Intro: https://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-14-289

Tajima's D is computed as the difference between two measures of genetic diversity: the mean number of pairwise differences and the number of segregating sites, each scaled so that they are expected to be the same in a neutrally evolving population of constant size.

A *negative Tajima's D* signifies an excess of low frequency polymorphisms relative to expectation, indicating population size expansion (e.g., after a bottleneck or a selective sweep)

A *positive Tajima's D* signifies low levels of both low and high frequency polymorphisms, indicating a decrease in population size and/or balancing selection.

```
#!/bin/bash

#PBS -l nodes=1:ppn=8
#PBS -l walltime=12:00:00
#PBS -l pmem=24gb
#PBS -M azd6024@psu.edu
#PBS -A open
#PBS -j oe

WORKINGDIR=/gpfs/group/jrl35/default/aayudh/gwas_reseq
cd $WORKINGDIR

vcftools --gzvcf Sorghum_1757g_AllChr.polymorphic.snp.noRepeats.5pctMasked.imputed.combined.vcf.gz --out tajimasd --TajimaD 1000
```

R part

```
setwd("~/Library/CloudStorage/OneDrive-UniversityofVermont/PENN STATE/Fst_tajimasD_reseq")
list.files()
library(ggplot2)
library(tidyverse)
library(ggtext)
library(normentR)
library(tidyr)

tajimasD <- read.csv("tajimasd_1kb_window.csv")
head(tajimasD)

data <- tajimasD[1:683650,]
data <- data %>% drop_na()
chr = substring(data$CHROM, 4)
data[, "CHROM"] <- chr
tail(data)

mean(data$TajimaD)
median(data$TajimaD)

# Basic box plot
setwd("~/Library/CloudStorage/OneDrive-UniversityofVermont/PENN STATE/Fst_tajimasD_reseq")
tiff("Reseq_tajimasD.tiff", width = 5, height = 5, units = 'in', res = 300)
ggplot(data, aes(x=CHROM, y=TajimaD,color=CHROM)) + geom_boxplot() +
  scale_color_manual(values=c("brown","black","brown","black","brown","black","brown","black","brown","black"))+ 
  labs(title=NULL,x="Chromosome", y = "Tajima's D")+
  scale_y_continuous(limits = c(-2.6, 6.5), breaks = seq(-2.6, 6.5,3))+
  theme_classic() +
  theme(legend.position="none",text=element_text(size=16, colour = "black", family="Times New Roman"),
        axis.line = element_line(size=0.5, colour = "black"),
        axis.text.x=element_text(colour="black", size = 16),
        axis.text.y=element_text(colour="black", size = 16))
dev.off()


setwd("~/Library/CloudStorage/OneDrive-UniversityofVermont/PENN STATE/Fst_tajimasD_reseq")
tiff("tajimasD_histogram.tiff", width = 5, height = 5, units = 'in', res = 300)
hist(data$TajimaD, prob = TRUE, main = NULL, xlab="Tajima's D",col = "lightblue")
x <- seq(min(data$TajimaD), max(data$TajimaD), length = 40)
f <- dnorm(x, mean = mean(data$TajimaD), sd = sd(data$TajimaD))
lines(x, f, col = "red", lwd = 2)
dev.off()
```

-----
<div id='id-section10'/>

## Chapter 10: Fst estimation of Sorghum Reseq 

### Choosing Fst population

1. We have to just take 384 accessions that has geographic info and then identify country

```
setwd("~/Library/CloudStorage/OneDrive-UniversityofVermont/PENN STATE/Fst_tajimasD_reseq")
list.files()
data <- read.csv("Reseq_accessions_Aayudh.csv")
head(data)

library(maps)
data$country <- map.where("world", data$Lon, data$Lat)
head(data)
```

2. Now divide africa into different pop grids

```
setwd("~/Library/CloudStorage/OneDrive-UniversityofVermont/PENN STATE/Fst_tajimasD_reseq")
list.files()
data <- read.csv("Reseq_accessions_country.csv")
head(data)
##Remove a country
#cc <- data[- grep("India", data$country),]
##Remove multiple country
onlyAfrica <- data[!grepl("India|USA|Sri Lanka|Pakistan|China|Australia|Turkey|Syria", data$country),]

world <- map_data("world")
head(world)

World_africa <- world[(world$long > -20) & (world$long < 55), ]
World_africa_1 <- World_africa[(World_africa$lat > -35) & (World_africa$lat < 35), ]

setwd("~/Library/CloudStorage/OneDrive-UniversityofVermont/PENN STATE/Fst_tajimasD_reseq/fst populations")
tiff("Africa_grids.tiff", width = 6, height = 6, units = 'in', res = 300)
ggplot() + geom_map(data = world, map = world,aes(long, lat, map_id = region),
                    color = "black", fill = "white", size = 0.5) +
  geom_point(data = onlyAfrica,size=0.5,color = "red",aes(Lon, Lat),alpha = 1) +
  coord_sf(xlim = c(-20, 50), ylim = c(-35, 35), expand = FALSE) +  
  theme(panel.grid = element_line(color = "#8ccde3",size = 0.75,linetype = 2),
        panel.ontop = TRUE, panel.background = element_rect(color = NA, fill = NA))
dev.off()

#further modify grid

world <- map_data("world")
head(world)

World_africa <- world[(world$long > -20) & (world$long < 55), ]
World_africa_1 <- World_africa[(World_africa$lat > -35) & (World_africa$lat < 35), ]

setwd("~/Library/CloudStorage/OneDrive-UniversityofVermont/PENN STATE/Fst_tajimasD_reseq/fst populations")
tiff("Africa_grids_2degree.tiff", width = 6, height = 6, units = 'in', res = 300)
ggplot() + geom_map(data = World_africa_1, map = World_africa_1,aes(long, lat, map_id = region),
                    color = "black", fill = "white", size = 0.5) +
  geom_point(data = onlyAfrica,size=0.5,color = "red",aes(Lon, Lat),alpha = 1) +
  theme(panel.grid = element_line(color = "#8ccde3",size = 0.15,linetype = 1),
        panel.ontop = TRUE, panel.background = element_rect(color = NA, fill = NA))+
  scale_y_continuous(breaks = seq(-35, 35, by = 2))+
  scale_x_continuous(breaks = seq(-20, 55, by = 2))
dev.off()

###populations
setwd("~/Library/CloudStorage/OneDrive-UniversityofVermont/PENN STATE/Fst_tajimasD_reseq/fst populations")
pop1 <- onlyAfrica[(onlyAfrica$Lat > 12) & (onlyAfrica$Lat < 20) & 
                     (onlyAfrica$Lon > -20) & (onlyAfrica$Lon < -10), ]
write.csv(pop1,"pop1.csv")

pop2 <- onlyAfrica[(onlyAfrica$Lat > 12) & (onlyAfrica$Lat < 20) & 
                     (onlyAfrica$Lon > -10) & (onlyAfrica$Lon < 0), ]
write.csv(pop2,"pop2.csv")

pop3 <- onlyAfrica[(onlyAfrica$Lat > 12) & (onlyAfrica$Lat < 20) & 
                     (onlyAfrica$Lon > -0) & (onlyAfrica$Lon < 10), ]
write.csv(pop3,"pop3.csv")

pop4 <- onlyAfrica[(onlyAfrica$Lat > 12) & (onlyAfrica$Lat < 20) & 
                     (onlyAfrica$Lon > 10) & (onlyAfrica$Lon < 20), ]
write.csv(pop4,"pop4.csv")

pop5 <- onlyAfrica[(onlyAfrica$Lat > 12) & (onlyAfrica$Lat < 20) & 
                     (onlyAfrica$Lon > 20) & (onlyAfrica$Lon < 30), ]
write.csv(pop5,"pop5.csv")

pop6 <- onlyAfrica[(onlyAfrica$Lat > 12) & (onlyAfrica$Lat < 20) & 
                     (onlyAfrica$Lon > 30) & (onlyAfrica$Lon < 40), ]
write.csv(pop6,"pop6.csv")

pop7 <- onlyAfrica[(onlyAfrica$Lat > 12) & (onlyAfrica$Lat < 20) & 
                     (onlyAfrica$Lon > 40) & (onlyAfrica$Lon < 50), ]
write.csv(pop7,"pop7.csv")

pop8 <- onlyAfrica[(onlyAfrica$Lat > 2) & (onlyAfrica$Lat < 12) & 
                     (onlyAfrica$Lon > -10) & (onlyAfrica$Lon < 0), ]
write.csv(pop8,"pop8.csv")

pop9 <- onlyAfrica[(onlyAfrica$Lat > 2) & (onlyAfrica$Lat < 12) & 
                     (onlyAfrica$Lon > -0) & (onlyAfrica$Lon < 10), ]
write.csv(pop9,"pop9.csv")

pop10 <- onlyAfrica[(onlyAfrica$Lat > 2) & (onlyAfrica$Lat < 12) & 
                      (onlyAfrica$Lon > 10) & (onlyAfrica$Lon < 20), ]
write.csv(pop10,"pop10.csv")

pop11 <- onlyAfrica[(onlyAfrica$Lat > 2) & (onlyAfrica$Lat < 12) & 
                      (onlyAfrica$Lon > 20) & (onlyAfrica$Lon < 30), ]
write.csv(pop11,"pop11.csv")

pop12 <- onlyAfrica[(onlyAfrica$Lat > 2) & (onlyAfrica$Lat < 12) & 
                      (onlyAfrica$Lon > 30) & (onlyAfrica$Lon < 40), ]
write.csv(pop12,"pop12.csv")

pop13 <- onlyAfrica[(onlyAfrica$Lat > 2) & (onlyAfrica$Lat < 12) & 
                      (onlyAfrica$Lon > 40) & (onlyAfrica$Lon < 50), ]
write.csv(pop13,"pop13.csv")

pop14 <- onlyAfrica[(onlyAfrica$Lat > -10) & (onlyAfrica$Lat < 2) & 
                      (onlyAfrica$Lon > 20) & (onlyAfrica$Lon < 30), ]
write.csv(pop14,"pop14.csv")

pop15 <- onlyAfrica[(onlyAfrica$Lat > -10) & (onlyAfrica$Lat < 2) & 
                      (onlyAfrica$Lon > 30) & (onlyAfrica$Lon < 40), ]
write.csv(pop15,"pop15.csv")

pop16 <- onlyAfrica[(onlyAfrica$Lat > -20) & (onlyAfrica$Lat < -10) & 
                      (onlyAfrica$Lon > 20) & (onlyAfrica$Lon < 30), ]

pop17 <- onlyAfrica[(onlyAfrica$Lat > -20) & (onlyAfrica$Lat < -10) & 
                      (onlyAfrica$Lon > 30) & (onlyAfrica$Lon < 40), ]

pop18 <- onlyAfrica[(onlyAfrica$Lat > -30) & (onlyAfrica$Lat < -20) & 
                      (onlyAfrica$Lon > 10) & (onlyAfrica$Lon < 20), ]

pop19 <- onlyAfrica[(onlyAfrica$Lat > -30) & (onlyAfrica$Lat < -20) & 
                      (onlyAfrica$Lon > 20) & (onlyAfrica$Lon < 30), ]

pop20 <- onlyAfrica[(onlyAfrica$Lat > -30) & (onlyAfrica$Lat < -20) & 
                      (onlyAfrica$Lon > 30) & (onlyAfrica$Lon < 40), ]

pop21 <- onlyAfrica[(onlyAfrica$Lat > -40) & (onlyAfrica$Lat < -30) & 
                      (onlyAfrica$Lon > 20) & (onlyAfrica$Lon < 30), ]

### combined pop

pop16_comb <- rbind(pop16, pop17)
write.csv(pop16_comb,"pop16.csv")

pop17_comb <- rbind(pop18, pop21)
write.csv(pop17_comb,"pop17.csv")

pop18_comb <- rbind(pop19, pop20)
write.csv(pop18_comb,"pop18.csv")
```

3. Same for Asia as well

```
setwd("~/Library/CloudStorage/OneDrive-UniversityofVermont/PENN STATE/Fst_tajimasD_reseq")
list.files()
library(maps)
library(tidyverse)
library(sf)
library(mapview)
library(ggplot2)

data <- read.csv("Reseq_accessions_country.csv")

Asia <- data[(data$Lat > 0) & (data$Lat < 40), ]
SSA <- Asia[(Asia$Lon > 60) & (Asia$Lon < 100), ]

world <- map_data("world")
head(world)

setwd("~/Library/CloudStorage/OneDrive-UniversityofVermont/PENN STATE/Fst_tajimasD_reseq")
tiff("SouthAsia_grids.tiff", width = 6, height = 6, units = 'in', res = 300)
ggplot() + geom_map(data = world, map = world,aes(long, lat, map_id = region),
                    color = "black", fill = "white", size = 0.5) +
  geom_point(data = SSA,size=0.5,color = "red",aes(Lon, Lat),alpha = 1) +
  coord_sf(xlim = c(60, 100), ylim = c(0, 40), expand = FALSE)+  
  theme(panel.grid = element_line(color = "#8ccde3",size = 0.75,linetype = 2),
        panel.ontop = TRUE, panel.background = element_rect(color = NA, fill = NA))
dev.off()

###populations
setwd("~/Library/CloudStorage/OneDrive-UniversityofVermont/PENN STATE/Fst_tajimasD_reseq/fst populations")
###Divide based on east and west

SSA_west <- SSA[(SSA$Lat > 5) & (SSA$Lat < 35) & 
                  (SSA$Lon > 70) & (SSA$Lon < 76.5), ]
write.csv(SSA_west,"pop19.csv")

SSA_east <- SSA[(SSA$Lat > 5) & (SSA$Lat < 35) & 
                  (SSA$Lon > 76.5) & (SSA$Lon < 85), ]
write.csv(SSA_east,"pop20.csv")

```
### Estimate Fst with VCF tools

Upload the text files contaning accessions of each pop without the column name, then run vcftools. 

```
#!/bin/bash

#PBS -l nodes=1:ppn=8
#PBS -l walltime=12:00:00
#PBS -l pmem=24gb
#PBS -M azd6024@psu.edu
#PBS -A open
#PBS -j oe

WORKINGDIR=/storage/home/azd6024/work/fst_reseq/fst
cd $WORKINGDIR

vcftools --gzvcf Sorghum_reseq_SNP.vcf.gz --weir-fst-pop pop1.txt --weir-fst-pop pop2.txt --weir-fst-pop pop3.txt --weir-fst-pop pop4.txt --weir-fst-pop pop5.txt --weir-fst-pop pop6.txt --weir-fst-pop pop7.txt --weir-fst-pop pop8.txt --weir-fst-pop pop9.txt --weir-fst-pop pop10.txt --weir-fst-pop pop11.txt --weir-fst-pop pop12.txt --weir-fst-pop pop13.txt --weir-fst-pop pop14.txt --weir-fst-pop pop15.txt --weir-fst-pop pop16.txt --weir-fst-pop pop17.txt --weir-fst-pop pop18.txt --weir-fst-pop pop19.txt --weir-fst-pop pop20.txt --fst-window-size 1000 --out reseq20pop1kb_Fst
```

### Find Fst of your gene of interest 
```
setwd("~/Library/CloudStorage/OneDrive-UniversityofVermont/PENN STATE/Fst_tajimasD_reseq")
list.files()
library(ggplot2)
library(tidyverse)
library(ggtext)
library(normentR)
library(tidyr)

Fst <- read.csv("FST_20pop_reseq_1kbWindow.csv")
head(Fst)

data <- Fst[1:652922,]
data <- data %>% drop_na()
chr = substring(data$CHROM, 4)
data[, "CHROM"] <- chr
##Weir's $F_ST$ can give values <0 so we reset to 0
data$WEIGHTED_FST[which(data$WEIGHTED_FST<0)]=0
head(data)


chr4_Fst <- subset(data, CHROM=="04")
head(chr4_Fst)
Chr4_49813813_50kb_Fst <- subset(chr4_Fst, BIN_START >49763813 & BIN_START<49863813)

Chr4_49813813_50kb_Fst$Fst_round <- round(Chr4_49813813_50kb_Fst$WEIGHTED_FST, digits = 2)
head(Chr4_49813813_50kb_Fst)

Chr4_49813813_50kb_Fst$BIN_START=Chr4_49813813_50kb_Fst$BIN_START/1E6
Chr4_49813813_50kb_Fst$BIN_END=Chr4_49813813_50kb_Fst$BIN_END/1E6

setwd("~/Library/CloudStorage/OneDrive-UniversityofVermont/PENN STATE/Fst_tajimasD_reseq")
tiff("Chr4_49813813_50kb_Fst.tiff", width = 5, height = 5, units = 'in', res = 300)
par(family="Times")
par(pty="s")
with(Chr4_49813813_50kb_Fst, plot(BIN_START, Fst_round, pch=20,lwd=4,xaxt="n",col="black", yaxt = "n", main=NA, ylab= expression(F[ST]), xlab="Mb Chromosome"))
axis(side = 1, at=seq(min(Chr4_49813813_50kb_Fst$BIN_START), max(Chr4_49813813_50kb_Fst$BIN_START), by = (max(Chr4_49813813_50kb_Fst$BIN_START)-min(Chr4_49813813_50kb_Fst$BIN_START))/3),las=1,cex.axis = 1)
axis(side = 2, at=seq(min(Chr4_49813813_50kb_Fst$Fst_round), max(Chr4_49813813_50kb_Fst$Fst_round), by = 0.11),las=1,cex.axis = 1)
abline(v = 49.813813, col = "red", lwd = 2, lty = 2)
dev.off()
```

### Find locus with high Fst

```
setwd("~/Library/CloudStorage/OneDrive-UniversityofVermont/PENN STATE/Fst_tajimasD_reseq")
list.files()
library(ggplot2)
library(tidyverse)
library(ggtext)
library(normentR)
library(tidyr)

###Top 10 high fst locus
Fst <- read.csv("FST_20pop_reseq_1kbWindow.csv")
head(Fst)

data <- Fst[1:652922,]
data <- data %>% drop_na()
chr = substring(data$CHROM, 4)
data[, "CHROM"] <- chr
##Weir's $F_ST$ can give values <0 so we reset to 0
data$WEIGHTED_FST[which(data$WEIGHTED_FST<0)]=0
head(data)
top10_highFst <- data[with(data,order(-WEIGHTED_FST)),]
top10_highFst <- top10_highFst[1:10,1:5]
top10_highFst_round <- round(top10_highFst$WEIGHTED_FST, digits = 3)
top10_highFst[, "WEIGHTED_FST"] <- top10_highFst_round
head(top10_highFst)

#write.table(top10_highFst, "top10_highFst.csv", sep=",")

top10_highFst <- read.csv("top10_highFst.csv")
dim(top10_highFst)
head(top10_highFst)

#export table
library("gridExtra")
tiff("top10Fst_locus.tiff", width = 8, height = 5, units = 'in', res = 300)
par(family="Times")
grid.table(top10_highFst)
dev.off()

setwd("~/Library/CloudStorage/OneDrive-UniversityofVermont/PENN STATE/eGWAS_revised list/Annotation")
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

#adding these 2 colums to the dataframe will give you the +/- 2.5 kb window
Sbicolor_annot$start_5kb <- Sbicolor_annot$start - 2500
Sbicolor_annot$stop_5kb <- Sbicolor_annot$stop + 2500
head(Sbicolor_annot)

setwd("~/Library/CloudStorage/OneDrive-UniversityofVermont/PENN STATE/Fst_tajimasD_reseq")
list.files()

library(dplyr)
library(fuzzyjoin)
gwas_annot <- fuzzy_left_join(top10_highFst, Sbicolor_annot, by=c("BIN_START"="start_5kb", "BIN_START"="stop_5kb", "chr"="chr"),
                              match_fun=list(`>=`, `<=`, `==`)) %>% select(Gene_name,Chromosome,BIN_START,BIN_END,N_VARIANTS,WEIGHTED_FST,start,stop,direction)


gwas_annot_split <- data.frame(do.call("rbind", strsplit(as.character(gwas_annot$Gene_name), ";", fixed = TRUE)))
pacId <- sub('pacid=', '', gwas_annot_split$X3)
head(gwas_annot_split)

gwas_annot_1 <- cbind (pacId,gwas_annot)
gwas_annot_1 = within(gwas_annot_1, rm(Gene_name))
head(gwas_annot_1)
gwas_annot_2 <- gwas_annot_1[- grep("ancestorIdentifier", gwas_annot_1$pacId),]
head(gwas_annot_2)

library(data.table)
gwas_annot_3 <- unique(setDT(gwas_annot_2)[order(BIN_START, -BIN_START)], by = "BIN_START")
dim(gwas_annot_3)
head(gwas_annot_3)

library(tidyr)
gwas_annot <- gwas_annot_3 %>% drop_na()
dim(gwas_annot)
head(gwas_annot)

setwd("~/Library/CloudStorage/OneDrive-UniversityofVermont/PENN STATE/eGWAS_revised list/Annotation")
Sorghum_mdata <- read.csv("Sbicolor_454_v3.1.1.annotation_info.csv", header = T)
head(Sorghum_mdata)

setwd("~/Library/CloudStorage/OneDrive-UniversityofVermont/PENN STATE/Fst_tajimasD_reseq")
annotation_Loc <- merge(Sorghum_mdata,gwas_annot, by = "pacId")
dim(annotation_Loc)
head(annotation_Loc)

write.table(annotation_Loc, "1. Annotation_top10_highFst_final.csv", sep=",")
```

-----
<div id='id-section11'/>

## Chapter 11: LD Decay 

### Creating genotype file and mapping file

```
#!/usr/bin/env Rscript

#PBS -l nodes=1:ppn=8
#PBS -l walltime=12:00:00
#PBS -l pmem=24gb
#PBS -M azd6024@psu.edu
#PBS -A open
#PBS -j oe
#PBS -m abe

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
mySNPs <- snpgdsSelectSNP(genofile, traits$LIB, maf = 0.05, missing.rate = 0.5)
length(mySNPs)

mySNPmatrix <- snpgdsGetGeno(genofile, sample.id = traits$LIB, snp.id = mySNPs, with.id = TRUE)
dim(mySNPmatrix$genotype)

#these lines create the bed file for gemma
rs <- data.frame(chr = chromosome.id[mySNPs], positions = snp.position[mySNPs])
SNPs_bim <- data.frame(paste0('rs_', rs$chr, '_', rs$positions), t(mySNPmatrix$genotype))

SNPs_bim[1:4,1:4]
dim(SNPs_bim)
length(SNPs_bim)

df <- SNPs_bim[,2:340]
df [1:4,1:4]

#sommer requires data as -1,0,1
data_2 <- df - 1
colnames(data_2) <- mySNPmatrix$sample.id
data_2 [1:4,1:4]

Accessions <- as.data.frame(SNPs_bim[,1])
names(Accessions)[1]<-paste("Accessions")

data_3 <- cbind(Accessions,data_2)
data_3[1:4,1:4]

# transpose all but the first column (name)
n <- data_3$Accessions
df.aree <- as.data.frame(t(data_3[,-1]))

colnames(df.aree) <- n
df.aree[1:4,1:4]

geno_reseq <- tibble::rownames_to_column(df.aree, "Accessions")

geno_reseq[1:4,1:4]
dim(geno_reseq)

##createing genofile
geno_reseq_HBP=df.aree[geno_reseq$Accessions,]
geno_reseq_HBP[1:4,1:4]
dim(geno_reseq_HBP)

##Now SNP info extraction
map_reseq <- as.data.frame(SNPs_bim$paste0..rs_...rs.chr..._...rs.positions.)

rs_split <- data.frame(do.call("rbind", strsplit(as.character(map_reseq$`SNPs_bim$paste0..rs_...rs.chr..._...rs.positions.`), "_", fixed = TRUE)))
head(rs_split)

rs_split[, "X1"] <- map_reseq$`SNPs_bim$paste0..rs_...rs.chr..._...rs.positions.`

names(rs_split)[1]<-paste("Locus")
names(rs_split)[2]<-paste("LG")
names(rs_split)[3]<-paste("Position")

head(rs_split)
rs_split$Position <- as.numeric(as.character(rs_split$Position)) 
rs_split[, "Position"]=rs_split$Position/1E6

rs_split$LG <- as.numeric(as.character(rs_split$LG))
head(rs_split)

map_reseq <- rs_split
head(map_reseq)
dim(map_reseq)
summary(map_reseq)

mapCP <- map_reseq
dim(mapCP)
#### just chromosome 4 (change according to your snp hit)
mapCP <- mapCP[mapCP$LG == '4',]
head(mapCP)
summary(mapCP)

##taking 1mb in both direction from the SNP
mapCP_5mb <- subset(mapCP, Position >48.81 & Position<50.81)
head(mapCP_5mb)
tail(mapCP_5mb)
dim(mapCP_5mb)

setwd("/storage/home/azd6024/work/linkage_D")
write.csv(mapCP_5mb,"mapCP_5mb.csv")

##genofile
CPgeno <- geno_reseq_HBP
CPgeno [1:4,1:4]
dim(CPgeno)

#find what column label is what number in the genofile. You have to make sure the genofile also
#in same length as the map file
head(mapCP_5mb)
#take top of the head SNP name
match("rs_04_48810669",names(geno_reseq_HBP))
tail(mapCP_5mb)
#take bottom of the tail SNP name
match("rs_04_50809956",names(geno_reseq_HBP))

CPgeno_5mb <- CPgeno[,match("rs_04_48810669",names(geno_reseq_HBP)):match("rs_04_50809956",names(geno_reseq_HBP))]
dim(CPgeno_5mb)

write.csv(CPgeno_5mb,"CPgeno_5mb.csv")
```

### Running Sommer

```
#!/usr/bin/env Rscript

#PBS -l nodes=1:ppn=8
#PBS -l walltime=12:00:00
#PBS -l pmem=24gb
#PBS -M azd6024@psu.edu
#PBS -A open
#PBS -j oe
#PBS -m abe

##LD_decay
# https://rdrr.io/cran/sommer/man/LD.decay.html  

library(sommer)
library(dplyr)
library(stringr)
library(ggplot2)
library(cowplot)
library(data.table)

CPgeno_5mb <- read.csv("CPgeno_5mb.csv")
CPgeno_5mb <- CPgeno_5mb [,2:33457]

mapCP_5mb <- read.csv("mapCP_5mb.csv")
head(mapCP_5mb)
dim(mapCP_5mb)

#running sommer
res <- LD.decay(CPgeno_5mb, mapCP_5mb,silent=FALSE,unlinked=FALSE)
names(res)
```

$resp a list with 3 elements; "by.LG", "all.LG", "LDM". The first element (by.LG) is a list with as
many elements as chromosomes where each contains a matrix with 3 columns, the distance,
the r2 value, and the p-value associated to the chi-square test for disequilibrium. The second
element (all.LG) has a big matrix with distance, r2 values and p-values, for each point from
all chromosomes in a single data.frame. The third element (LDM) is the matrix of linkage
disequilibrium between pairs of markers.
If unlinked is selected the program should return the gamma percentile interchromosomal LD
(r2) for each chromosome and average.

### Creating LDM matrix and plotting

```
LDM <- res$LDM
df <- LDM[[4]]
LDM_matrix <- as.data.frame(df)
dim(LDM_matrix)

write.csv(LDM_matrix, "LDM_matrix_reseq_chr4_2mb.csv")

LDM_matrix <- read.csv("LDM_matrix_reseq_chr4_2mb.csv")
LDM_matrix[1:4,1:4]

#ur SNP of interest
match("rs_04_49813813",names(LDM_matrix))
SNP_centered <- LDM_matrix[,match("rs_04_49813813",names(LDM_matrix))]
SNP_main <- as.data.frame(SNP_centered)
SNP_rows <- LDM_matrix[,1]
SNP_full <- cbind(SNP_rows,SNP_main)
head(SNP_full)

write.csv(SNP_full,"SNP_full.csv")

library(ggplot2)
setwd("~/Library/CloudStorage/OneDrive-UniversityofVermont/PENN STATE/LD")
list.files()
data <- read.csv("1.PBS1_SNP_full.csv")
head(data)

rs_split <- data.frame(do.call("rbind", strsplit(as.character(data$SNP_rows), "_", fixed = TRUE)))
head(rs_split)

rs_split[, "X1"] <- data$SNP_rows

names(rs_split)[1]<-paste("SNP")
names(rs_split)[2]<-paste("chr")
names(rs_split)[3]<-paste("Position")

head(rs_split)
rs_split$Position <- as.numeric(as.character(rs_split$Position)) 
rs_split[, "Position"]=rs_split$Position/1E6
head(rs_split)

SNP_LDM <- cbind(rs_split,data$SNP_centered)
names(SNP_LDM)[4]<-paste("r2")
head(SNP_LDM)

SNP_LDM <- 

library(ggplot2)
setwd("~/Library/CloudStorage/OneDrive-UniversityofVermont/PENN STATE/LD")
tiff("LDM_Reseq_chr4_PBS1_2mb_ggpplot.tiff", width = 5, height = 5, units = 'in', res = 300)
ggplot()+
  geom_point(aes(x=SNP_LDM$Position,y=SNP_LDM$r2),size=0.4,colour="cadetblue")+
  labs(x="Position (Megabases)",y=expression(LD~(r^{2})))+
  geom_smooth(aes(x=SNP_LDM$Position,y=SNP_LDM$r2),se=FALSE, method="loess", span=0.75)+
  #scale_y_continuous(breaks=seq(0,1,.2))+
  scale_x_continuous(breaks=seq(48.8,50.8,0.5))+
  theme_classic() +
  geom_vline(xintercept = 49.813813, linetype="dashed", color = "red", size=1)+
  theme(legend.position="none",text=element_text(size=16, colour = "black", family="Times New Roman"),
        axis.line = element_line(size=0.5, colour = "black"),
        axis.text.x=element_text(colour="black", size = 16),
        axis.text.y=element_text(colour="black", size = 16))
dev.off()
SNP_LDM_r2 <- SNP_LDM[which(SNP_LDM$r2 > .15),]

setwd("~/Library/CloudStorage/OneDrive-UniversityofVermont/PENN STATE/LD")
tiff("LDM_r2_.15_Reseq_chr4_PBS1_2mb_ggpplot.tiff", width = 5, height = 5, units = 'in', res = 300)
ggplot()+
  geom_point(aes(x=SNP_LDM_r2$Position,y=SNP_LDM_r2$r2),size=0.4,colour="cadetblue")+
  labs(x="Position (Megabases)",y=expression(LD~(r^{2})))+
  geom_smooth(aes(x=SNP_LDM_r2$Position,y=SNP_LDM_r2$r2),se=FALSE, method="loess", span=0.75)+
  scale_y_continuous(breaks=seq(0,1,.25))+
  scale_x_continuous(breaks=seq(48.8,50.8,0.5))+
  theme_classic() +
  geom_vline(xintercept = 49.813813, linetype="dashed", color = "red", size=1)+
  theme(legend.position="none",text=element_text(size=16, colour = "black", family="Times New Roman"),
        axis.line = element_line(size=0.5, colour = "black"),
        axis.text.x=element_text(colour="black", size = 16),
        axis.text.y=element_text(colour="black", size = 16))
dev.off()
```
### Plotting with gene length 800kb region

```


list.files()
sbicolor <- read.csv("Sbicolor_454_v3.1.1.gene.csv")
head(sbicolor)
sbicolor_chr4 <- subset(sbicolor, chr == "Chr04")
head(sbicolor_chr4)

annot_300kb <- subset(sbicolor_chr4, sbicolor_chr4$start > 49110181 & sbicolor_chr4$start < 49863813)
head(annot_300kb)
genes <- subset(annot_300kb, gene == "gene")
genes$size <- genes$stop - genes$start
genes$position <- genes$start + (genes$size/2)
genes[, "position"]=genes$position/1E3
genes[, "size"]=genes$size/1E3
head(genes)

write.csv(genes, "geneinfo_49010-49838.csv")

genes <- read.csv("geneinfo_49010-49838.csv")

library(ggplot2)
tiff("genes_pbs1near.tiff", width = 13, height = 3.5, units = 'in', res = 300)
ggplot(genes, aes(x=position,y=size))+
  geom_bar(stat="identity", width=4,color="coral4")+
  scale_x_continuous(breaks=seq(49110,49838,182))+
  scale_y_continuous(breaks=seq(0,23,7.6))+
  labs(x=NULL,y="Gene length (kb)")+
  theme_classic() +
  geom_vline(xintercept = 49813.813, linetype="dashed", color = "red", size=0.5)+
  theme(legend.position="none",text=element_text(size=16, colour = "black", family="Times New Roman"),
        axis.line = element_line(size=0.5, colour = "black"),
        axis.text.x=element_blank(),
        axis.ticks.x = element_blank(),
        axis.line.x = element_blank(),
        axis.text.y=element_text(colour="black", size = 16))
dev.off()

###Now more300kb far
#004G156732=49513181-49516616
#004G156800=#49519315-49526466
##49511622	49512961	+	ID=Sobic.004G156666.v3.2;Name=Sobic.004G156666
#49541926
haplotype_block <- subset(SNP_LDM_r2,  SNP_LDM_r2$Position > 49110 & SNP_LDM_r2$Position < 49838.813) 
head(haplotype_block)

library(ggplot2)
setwd("~/Library/CloudStorage/OneDrive-UniversityofVermont/PENN STATE/LD")
tiff("700_kbfromSNP.tiff", width = 13, height = 3.5, units = 'in', res = 300)
ggplot()+
  geom_point(aes(x=haplotype_block$Position,y=haplotype_block$r2),size=1,colour="cadetblue")+
  labs(x="Position (kb)",y=expression(LD~(r^{2})))+
  geom_smooth(aes(x=haplotype_block$Position,y=haplotype_block$r2),se=FALSE, method="loess", span=0.75)+
  #scale_y_continuous(breaks=seq(0,1,.2))+
  scale_x_continuous(breaks=seq(49110,49838,182))+
  theme_classic() +
  geom_vline(xintercept = 49813.813, linetype="dashed", color = "red", size=0.5)+
  theme(legend.position="none",text=element_text(size=16, colour = "black", family="Times New Roman"),
        axis.line = element_line(size=0.5, colour = "black"),
        axis.text.x=element_text(colour="black", size = 16),
        axis.text.y=element_text(colour="black", size = 16),
        panel.border = element_rect(colour = "black", fill=NA, size=1))
dev.off()
```

-----
<div id='id-section12'/>

## Chapter 12: Putative promoter analysis 

### Make a data file first

```
setwd("~/Library/CloudStorage/OneDrive-UniversityofVermont/PENN STATE/eGWAS_revised list/manhattan plot")
list.files()
Chr10corr <- read.csv("Reseq_preds_all_chr10_corrected.csv", header=T)
head(Chr10corr)

data_p0.0005<- subset(Chr10corr, p_wald < 0.0005)
write.csv(data_p0.0005,"Entire_genome_Reseq_preds_all_gemma_p0.0005.csv")

setwd("~/Library/CloudStorage/OneDrive-UniversityofVermont/PENN STATE/Putative promoter regions")
sbicolor <- read.csv("Sbicolor_454_v3.1.1.gene.csv")
head(sbicolor)
dim(sbicolor)
##sbicolor_chr10<- subset(sbicolor, chr=="Chr10"); just need to know end of chr10 as we need to remove supers
sbicolor <- sbicolor [1:425313,]

##+stand
genes_plus <- subset(sbicolor, direction=="+")
genes_minus <- subset(sbicolor, direction=="-")

##CDS
genes_plus_CDS <- subset(genes_plus, gene=="CDS")
genes_minus_CDS <- subset(genes_minus, gene=="CDS")

##3kb region
genes_plus_CDS$start_3kb <- genes_plus_CDS$start - 3000
genes_plus_CDS$stop_3kb <- genes_plus_CDS$stop + 3000

chr_new = substring(genes_plus_CDS$chr, 5)
#replace column
genes_plus_CDS[, "chr"] <- chr_new
head(genes_plus_CDS)

#write.csv(genes_plus_CDS,"Entire_genome_plus_strand_CDS.csv")

###Now the gene_minus starnd 
##3kb region
genes_minus_CDS$start_3kb <- genes_minus_CDS$start - 3000
genes_minus_CDS$stop_3kb <- genes_minus_CDS$stop + 3000

chr_new = substring(genes_minus_CDS$chr, 5)
#replace column
genes_minus_CDS[, "chr"] <- chr_new
head(genes_minus_CDS)

##write.csv(genes_minus_CDS,"Entire_genome_minus_strand_CDS.csv")
```


##this part should be in server-

```
#!/usr/bin/env Rscript

#PBS -l nodes=1:ppn=8
#PBS -l walltime=12:00:00
#PBS -l pmem=24gb
#PBS -M azd6024@psu.edu
#PBS -A open
#PBS -j oe
#PBS -m abe

setwd("~/work/promoter_analysis")

genes_minus_CDS <- read.csv("Entire_genome_minus_strand_CDS.csv")

entire_genome <- read.csv("Entire_genome_Reseq_preds_all_gemma_p0.0005.csv")

library(dplyr)
library(fuzzyjoin)
gwas_annot <- fuzzy_left_join(entire_genome, genes_minus_CDS, by=c("ps"="start_3kb", "ps"="stop_3kb", "chr"="chr"),
                              match_fun=list(`>=`, `<=`, `==`)) %>% select(Gene_name,ps,rs,af,p_wald,start,stop,direction)

genes_minus_CDS_pwald <- na.omit(gwas_annot)
genes_minus_CDS_pwald$snpsite <- genes_minus_CDS_pwald$stop - genes_minus_CDS_pwald$ps
head(genes_minus_CDS_pwald)

write.csv(genes_minus_CDS_pwald, "1.Entire_genome_MINUS_promoter.csv")

###PLUS STRAND----------------------------
genes_plus_CDS <- read.csv("Entire_genome_plus_strand_CDS.csv")

entire_genome <- read.csv("Entire_genome_Reseq_preds_all_gemma_p0.0005.csv")

library(dplyr)
library(fuzzyjoin)
gwas_annot <- fuzzy_left_join(entire_genome, genes_plus_CDS, by=c("ps"="start_3kb", "ps"="stop_3kb", "chr"="chr"),
                              match_fun=list(`>=`, `<=`, `==`)) %>% select(Gene_name,ps,rs,af,p_wald,start,stop,direction)

genes_plus_CDS_pwald <- na.omit(gwas_annot)
genes_plus_CDS_pwald$snpsite <- genes_plus_CDS_pwald$ps - genes_plus_CDS_pwald$start
head(genes_plus_CDS_pwald)

write.csv(genes_plus_CDS_pwald, "1.Entire_genome_PLUS_promoter.csv")

```
### Calculate quantile and make the plot

```
setwd("~/Library/CloudStorage/OneDrive-UniversityofVermont/PENN STATE/cycles_predicted/promoter")
list.files()
library(ggplot2)
library(tidyverse)
library(ggtext)
library(normentR)

list.files()
gene_plus <- read.csv("1.water_stress_veg_PLUS_promoter.csv")
library(data.table)
gene_plus_1 <- unique(setDT(gene_plus)[order(rs, -rs)], by = "rs")
summary(gene_plus_1)
gene_plus_1$pval <- -log10(gene_plus_1$p_wald)
gene_plus_1 <- gene_plus_1[order(gene_plus_1$snpsite),]
gene_plus_1$gene_ength <- gene_plus_1$stop - gene_plus_1$start
head(gene_plus_1)
dim(gene_plus_1)

#write.csv(gene_plus_1,"1.Entire_genome_PLUS_promoter_no_overlap_water_veg.csv", row.names = FALSE)

setwd("~/Library/CloudStorage/OneDrive-UniversityofVermont/PENN STATE/cycles_predicted/promoter")
gene_plus <-read.csv("1.Entire_genome_PLUS_promoter_no_overlap_water_veg.csv")
head(gene_plus)
###Figure out window size
##i +3 because it's taking first 4 rows and doing this for all dataset

dddf <- NULL
## column number 10 as we want snpsite length to figure what bp would be good
for (i in seq(1,nrow(gene_plus),by =3)){
  window <- gene_plus[i:(i +3),10]
  dddf <- rbind(dddf,window)}

data <- as.data.frame(dddf)
head(data)
data$size <- data$V1-data$V4
tail(data)
dim(data)
data <- data[1:387,]
mean(data$size)

###WINDOW1
##Loop to calculate quantile 
#use column number where logp is, it's 11 here.
qnt_w1 <- NULL

for (i in seq(1,nrow(gene_plus),by =3)){
  qres <- quantile(gene_plus[i:(i +3),11], probs = 0.9)
  qnt_w1 <- rbind(qnt_w1,qres)}

##Loop to center snpsite 
#use column number where snpsite is, it's 10 here.
ps_w1 <- NULL
##Column 11 is logp
for (i in seq(1,nrow(gene_plus),by =3)){
  snpcenter <- (-(gene_plus[i,10] - gene_plus[i+3,10])/2) + gene_plus[i,10]
  ps_w1 <- rbind(ps_w1,snpcenter)}

w1 <- as.data.frame(cbind(qnt_w1,ps_w1[1:387,]))
w1[, "V2"] <- (w1$V2/1000)
head(w1)

###w2
##Loop to calculate quantile 
qnt_w2 <- NULL

for (i in seq(2,nrow(gene_plus),by =3)){
  qres <- quantile(gene_plus[i:(i +3),11], probs = 0.9)
  qnt_w2 <- rbind(qnt_w2,qres)}

##Loop to center snpsite 
ps_w2 <- NULL

for (i in seq(2,nrow(gene_plus),by =3)){
  snpcenter <- (-(gene_plus[i,10] - gene_plus[i+3,10])/2) + gene_plus[i,10]
  ps_w2 <- rbind(ps_w2,snpcenter)}

w2 <- as.data.frame(cbind(qnt_w2,ps_w2[1:386,]))
w2[, "V2"] <- (w2$V2/1000)
head(w2)

###w3
##Loop to calculate quantile 
qnt_w3 <- NULL

for (i in seq(3,nrow(gene_plus),by =3)){
  qres <- quantile(gene_plus[i:(i +3),11], probs = 0.9)
  qnt_w3 <- rbind(qnt_w3,qres)}

##Loop to center snpsite 
ps_w3 <- NULL

for (i in seq(3,nrow(gene_plus),by =3)){
  snpcenter <- (-(gene_plus[i,10] - gene_plus[i+3,10])/2) + gene_plus[i,10]
  ps_w3 <- rbind(ps_w3,snpcenter)}

w3 <- as.data.frame(cbind(qnt_w3,ps_w3[1:386,]))
w3[, "V2"] <- (w3$V2/1000)
head(w3)



plus_qnt <- rbind(w1,w2,w3)
plus_qnt <- plus_qnt[order(plus_qnt$V2),]
head(plus_qnt)
dim(plus_qnt)

plus <- ggplot(plus_qnt, aes(x=V2, y=`90%`)) + 
  geom_point() +
  geom_segment(data=plus_qnt,aes(x=plus_qnt$V2, xend=plus_qnt$V2,y=3,yend=plus_qnt$`90%`),
               color = '#B4AF46',alpha=0.8,size=0.5) + 
  labs(x=NULL,y=NULL)+
  scale_x_continuous(breaks=seq(-3,5,1))+
  scale_y_continuous(breaks=seq(3,11,2))+
  annotate(geom="text", x=-0.35, y=11, label="TSS",color="#B4464B")+
  theme_classic() +
  geom_vline(xintercept = 0, linetype="dashed", color = "#B4464B", size=1)+
  theme(legend.position="none",text=element_text(size=16, colour = "black", family="Times New Roman"),
        axis.line = element_line(size=0.5, colour = "black"),
        axis.text.x=element_text(colour="black", size = 16),
        axis.text.y=element_text(colour="black", size = 16))
plus


###minus strand plot----------------------------------------------------------------------------------------
setwd("~/Library/CloudStorage/OneDrive-UniversityofVermont/PENN STATE/cycles_predicted/promoter")
list.files()
gene_minus <- read.csv("1.water_stress_veg_MINUS_promoter.csv")
gene_minus_1 <- unique(setDT(gene_minus)[order(rs, -rs)], by = "rs")
summary(gene_minus_1)
gene_minus_1$pval <- -log10(gene_minus_1$p_wald)
gene_minus_1 <- gene_minus_1[order(gene_minus_1$snpsite),]
gene_minus_1$gene_ength <- gene_minus_1$stop - gene_minus_1$start
head(gene_minus_1)
dim(gene_minus_1)

write.csv(gene_minus_1,"1.Entire_genome_minus_promoter_no_overlap.csv", row.names = FALSE)

gene_minus <-read.csv("1.Entire_genome_minus_promoter_no_overlap.csv")
head(gene_minus)

###Figure out window size
dddf <- NULL

for (i in seq(1,nrow(gene_minus),by =3)){
  window <- gene_minus[i:(i +3),10]
  dddf <- rbind(dddf,window)}

data <- as.data.frame(dddf)
head(data)
data$size <- data$V1-data$V4
tail(data)
dim(data)
data <- data[1:425,]
mean(data$size)

###WINDOW1
##Loop to calculate quantile 
qnt_w1 <- NULL

for (i in seq(1,nrow(gene_minus),by =3)){
  qres <- quantile(gene_minus[i:(i +3),11], probs = 0.9)
  qnt_w1 <- rbind(qnt_w1,qres)}

##Loop to center snpsite 
ps_w1 <- NULL

for (i in seq(1,nrow(gene_minus),by =3)){
  snpcenter <- (-(gene_minus[i,10] - gene_minus[i+3,10])/2) + gene_minus[i,10]
  ps_w1 <- rbind(ps_w1,snpcenter)}

w1 <- as.data.frame(cbind(qnt_w1,ps_w1[1:425,]))
w1[, "V2"] <- (w1$V2/1000)
head(w1)

###w2
##Loop to calculate quantile 
qnt_w2 <- NULL

for (i in seq(2,nrow(gene_minus),by =3)){
  qres <- quantile(gene_minus[i:(i +3),11], probs = 0.9)
  qnt_w2 <- rbind(qnt_w2,qres)}

##Loop to center snpsite 
ps_w2 <- NULL

for (i in seq(2,nrow(gene_minus),by =3)){
  snpcenter <- (-(gene_minus[i,10] - gene_minus[i+3,10])/2) + gene_minus[i,10]
  ps_w2 <- rbind(ps_w2,snpcenter)}

w2 <- as.data.frame(cbind(qnt_w2,ps_w2[1:425,]))
w2[, "V2"] <- (w2$V2/1000)
head(w2)

###w3
##Loop to calculate quantile 
qnt_w3 <- NULL

for (i in seq(3,nrow(gene_minus),by =3)){
  qres <- quantile(gene_minus[i:(i +3),11], probs = 0.9)
  qnt_w3 <- rbind(qnt_w3,qres)}

##Loop to center snpsite 
ps_w3 <- NULL

for (i in seq(3,nrow(gene_minus),by =3)){
  snpcenter <- (-(gene_minus[i,10] - gene_minus[i+3,10])/2) + gene_minus[i,10]
  ps_w3 <- rbind(ps_w3,snpcenter)}

w3 <- as.data.frame(cbind(qnt_w3,ps_w3[1:424,]))
w3[, "V2"] <- (w3$V2/1000)
head(w3)


minus_qnt <- rbind(w1,w2,w3)
minus_qnt <- minus_qnt[order(minus_qnt$V2),]
head(minus_qnt)
dim(minus_qnt)


minus <- ggplot(minus_qnt, aes(x=V2, y=`90%`)) + 
  geom_point() +
  geom_segment(data=minus_qnt,aes(x=minus_qnt$V2, xend=minus_qnt$V2,y=3,yend=minus_qnt$`90%`),
               color = '#B4AF46',alpha=0.8,size=0.5) + 
  labs(x=NULL,y=NULL)+
  scale_x_reverse(breaks=seq(-4,4,1))+
  annotate(geom="text", x=-0.35, y=11, label="TTS",color="#B4464B")+
  theme_classic() +
  geom_vline(xintercept = 0, linetype="dashed", color = "#B4464B", size=1)+
  theme(legend.position="none",text=element_text(size=16, colour = "black", family="Times New Roman"),
        axis.line.y = element_blank(),
        axis.text.x=element_text(colour="black", size = 16),
        axis.text.y=element_blank(),
        axis.ticks.y = element_blank())
minus

library(tidyverse)
library(gridExtra)
library(grid)
library(gridtext)

setwd("~/Library/CloudStorage/OneDrive-UniversityofVermont/PENN STATE/cycles_predicted/promoter")
tiff("Water_veg_promoter_plot.tiff", width = 8, height = 6, units = 'in', res = 300)
library(gridExtra)
grid.arrange(plus, minus, ncol = 2)
dev.off()


```

-----
<div id='id-section13'/>

## Chapter 13: ABA analysis 

### Finding colocalized SNPs

```
setwd("~/Library/CloudStorage/OneDrive-UniversityofVermont/PENN STATE/Awal_aayudh")
list.files()
ABA <- read.csv("Annotation (ABA).csv")
head(ABA)

sbicolor <- read.csv("Sbicolor_454_v3.1.1.gene.csv")
head(sbicolor)

rs_split <- data.frame(do.call("rbind", strsplit(as.character(sbicolor$Gene_name), ";", fixed = TRUE)))
head(rs_split)

rs_split1 <- data.frame(do.call("rbind", strsplit(as.character(rs_split$X2), "=", fixed = TRUE)))
head(rs_split1)

rs_split2 <- data.frame(do.call("rbind", strsplit(as.character(rs_split1$X2), ".", fixed = TRUE)))
head(rs_split2)

sbicolor$sobicID = paste(rs_split2$X1, rs_split2$X2, sep=".")
head(sbicolor)

sbicolor = within(sbicolor, rm(phytozomev12,Gene_name))
head(sbicolor)

setwd("~/Library/CloudStorage/OneDrive-UniversityofVermont/PENN STATE/Awal_aayudh")
#write.csv(sbicolor,"Sbicolor_454_v3.1.1.gene_sobicIDs.csv")

Sorghum_sobics <- read.csv("Sbicolor_454_v3.1.1.gene_sobicIDs.csv")
head(Sorghum_sobics)

Sorghum_sobics_onlyGene<- subset(Sorghum_sobics, gene=="gene")

ABA_start_stop <- merge(ABA, Sorghum_sobics_onlyGene, by ="sobicID")
head(ABA_start_stop)

ABA_start_stop$start_10kb <- ABA_start_stop$start - 10000
ABA_start_stop$stop_10kb <- ABA_start_stop$stop + 10000

chr_new = substring(ABA_start_stop$chr, 4)
#replace column
ABA_start_stop[, "chr"] <- as.numeric(chr_new)
head(ABA_start_stop)

summary(ABA_start_stop)

###Water stress Vegetative GWAS
setwd("~/Library/CloudStorage/OneDrive-UniversityofVermont/PENN STATE/cycles_predicted/manhattan plot")
list.files()
Water_stress_veg <- read.csv("Reseq_water_veg_avg_fdr.csv", header=T)
head(Water_stress_veg)
Water_stress_veg <- Water_stress_veg[order(Water_stress_veg$p_wald),]
head(Water_stress_veg)

Water_stress_veg_sig <- subset(Water_stress_veg , p_wald < 0.001)
head(Water_stress_veg_sig)

##Sig SNPs gemma output
setwd("~/Library/CloudStorage/OneDrive-UniversityofVermont/PENN STATE/Awal_aayudh")
#write.csv(Water_stress_veg_sig ,"Water_stress_veg_sig_p_001.csv")

Water_stress_veg_sig <- read.csv("Water_stress_veg_sig_p_001.csv")
head(Water_stress_veg_sig)
summary(Water_stress_veg_sig)

library(dplyr)
library(fuzzyjoin)
gwas_annot <- fuzzy_left_join(Water_stress_veg_sig, ABA_start_stop, by=c("ps"="start_10kb", "ps"="stop_10kb", "chr"="chr"),
                              match_fun=list(`>=`, `<=`, `==`)) %>% select(sobicID,ps,rs,af,p_wald,start,stop,direction,Pathway..species.)

head(gwas_annot)

gwas_annot <- na.omit(gwas_annot)
head(gwas_annot)
dim(gwas_annot)
gwas_annot$logp <- -log10(gwas_annot$p_wald)

setwd("~/Library/CloudStorage/OneDrive-UniversityofVermont/PENN STATE/Awal_aayudh")
write.csv(gwas_annot,"ABA.gene_sobicIDs_SNP_cyclesoutput.csv")

gwas_annot <- read.csv("ABA.gene_sobicIDs_SNP_cyclesoutput.csv")
head(gwas_annot)

setwd("~/Library/CloudStorage/OneDrive-UniversityofVermont/PENN STATE/eGWAS_revised list/Annotation")
Sorghum_mdata <- read.csv("Sbicolor_454_v3.1.1.annotation_info.csv", header = T)
head(Sorghum_mdata)

names(Sorghum_mdata)[2]<-paste("sobicID")

annotation_Loc <- merge(Sorghum_mdata,gwas_annot, by = "sobicID")
dim(annotation_Loc)
head(annotation_Loc)

library(data.table)
annotation_Loc_1 <- unique(setDT(annotation_Loc)[order(rs, -rs)], by = "rs")
dim(annotation_Loc_1)
head(annotation_Loc_1)
annotation_Loc_2 <- annotation_Loc_1[order(annotation_Loc_1$p_wald),]
setwd("~/Library/CloudStorage/OneDrive-UniversityofVermont/PENN STATE/Awal_aayudh")
#write.csv(annotation_Loc_2,"1.Cycles_ABA_final.csv")

gwas_annot <- read.csv("1.Cycles_ABA_final.csv")
head(gwas_annot)

##arrange chromosome wise

rs_split <- data.frame(do.call("rbind", strsplit(as.character(gwas_annot$rs), "_", fixed = TRUE)))
head(rs_split)

gwas_annot$chr <- rs_split$X2

gwas_annot_2 <- gwas_annot[order(gwas_annot$chr),]
head(gwas_annot_2)

###Other plotting style
library(tidyverse)
library(ggtext)
library(normentR)
library(ggplot2)

sig_data <- gwas_annot_2 %>% 
  subset(p_wald < 0.05)
notsig_data <- gwas_annot_2 %>% 
  subset(p_wald >= 0.0001) %>%
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

fdr_line <- 0.001


setwd("~/Library/CloudStorage/OneDrive-UniversityofVermont/PENN STATE/Awal_aayudh")
tiff("ABA_cycles_manhattan.tiff", width = 5, height = 5, units = 'in', res = 300)
ggplot(gwas_data, aes(x = bp_cum, y = -log10(p_wald), color = gwas_data$Pathway..species.
                      , size = -log10(p_wald))) +
  geom_hline(yintercept = -log10(fdr_line), color = "red", linetype = "dashed") +
  geom_point() +
  scale_y_continuous(breaks=seq(3.0,7.6,1.5))+
  scale_x_continuous(label = axis_set$chr, breaks = axis_set$center) +
  scale_color_manual(values = rep(c('slategrey','#00008B', '#8B8B00', "goldenrod"), unique(length(axis_set$chr)))) +
  scale_size_continuous(range = c(0.5,3)) +
  labs(x = NULL, y = NULL) + 
  theme_minimal() +
  theme( legend.position = "none",text=element_text(size=16, colour = "black",family="Times New Roman"),
         axis.line = element_line(size=0.5, colour = "black"),
         panel.grid.major.x = element_blank(),
         panel.grid.minor.x = element_blank(),
         panel.grid.major.y = element_blank(),
         panel.grid.minor.y = element_blank(),
         axis.text.x = element_text(size = 12, vjust = 0.5, family="Times New Roman"),
         panel.border = element_rect(color = "black",fill = NA,size = 1))+
  geom_text_repel(aes(bp_cum, -log10(p_wald),label=gene.symbol),max.overlaps = Inf)
dev.off()

```

### Zoomed manhattan of CKA2 and HAB1
```
setwd("~/Library/CloudStorage/OneDrive-UniversityofVermont/PENN STATE/cycles_predicted/manhattan plot")
list.files()
Water_stress_veg <- read.csv("Reseq_water_veg_avg_fdr.csv", header=T)
head(Water_stress_veg)
Water_stress_veg <- Water_stress_veg[order(Water_stress_veg$p_wald),]
head(Water_stress_veg)

##CKA2

data_chr1<- subset(Water_stress_veg, chr=="1")
head(data_chr1)


plus_100kb <- 6246059 + 100000 
minus_100kb <- 6246059 - 100000 

Chr1_6246059_100kb <- subset(data_chr1, ps >minus_100kb & ps<plus_100kb)

setwd("~/Library/CloudStorage/OneDrive-UniversityofVermont/PENN STATE/Awal_aayudh")
#write.csv(Chr1_6246059_100kb, "CKA2_100kb.csv")

data<- read.csv("CKA2_100kb.csv", header =T)
head(data)


data$logp <- as.numeric(-log10(data$p_wald))
data$ps_kb <- data$ps/1000
summary(data)

head(data)
data <- na.omit(data)
  
#
setwd("~/Library/CloudStorage/OneDrive-UniversityofVermont/PENN STATE/Awal_aayudh")
tiff("CKA2_zoomed_100kb.tiff", width = 5, height = 5, units = 'in', res = 300)
par(family="Times")
par(pty="s")
with(data, plot(ps_kb, logp, pch=20,lwd=4,xaxt="n",yaxt="n"))
# Add colored points: red if padj<0.05, orange of log2FC>1, green if both)
with(subset(data, logp > 6), points(ps_kb, logp, pch=20, col="red"))
with(subset(data, logp < 6), points(ps_kb, logp, pch=20, col="black"))
abline(h = 6, col = "red", lwd = 2, lty = 3) 
abline(v = 6233.999, col = "blue", lwd = 2, lty = 3) 
abline(v = 6240.632, col = "blue", lwd = 2, lty = 3) 
axis(side = 1, at=seq(min(data$ps_kb), max(data$ps_kb), by = 
                        (max(data$ps_kb)-min(data$ps_kb))/3),las=1,cex.axis = 1)
axis(side = 2, at=seq(0, 7.6, by = 2.5),las=1,cex.axis = 1)
dev.off()


##HAB1
data_chr9<- subset(Water_stress_veg, chr=="9")
head(data_chr9)


plus_100kb <- 55867122 + 100000 
minus_100kb <- 55867122 - 100000 

chr9_55867122_100kb <- subset(data_chr9, ps >minus_100kb & ps<plus_100kb)

setwd("~/Library/CloudStorage/OneDrive-UniversityofVermont/PENN STATE/Awal_aayudh")
#write.csv(chr9_55867122_100kb, "HAB1_100kb.csv")

data1<- read.csv("HAB1_100kb.csv", header =T)
head(data1)

data1$logp <- as.numeric(-log10(data1$p_wald))
data1$ps_kb <- data1$ps/1000
summary(data1)

head(data1)
data1 <- na.omit(data1)

setwd("~/Library/CloudStorage/OneDrive-UniversityofVermont/PENN STATE/Awal_aayudh")
tiff("HAB1_zoomed_100kb.tiff", width = 5, height = 5, units = 'in', res = 300)
par(family="Times")
par(pty="s")
with(data1, plot(ps_kb, logp, pch=20,lwd=4,xaxt="n",yaxt="n"))
# Add colored points: red if padj<0.05, orange of log2FC>1, green if both)
with(subset(data1, logp > 6), points(ps_kb, logp, pch=20, col="red"))
with(subset(data1, logp < 6), points(ps_kb, logp, pch=20, col="black"))
abline(h = 6, col = "red", lwd = 2, lty = 3) 
abline(v = 55861.006, col = "blue", lwd = 2, lty = 3) 
abline(v = 55863.087, col = "blue", lwd = 2, lty = 3) 
axis(side = 1, at=seq(min(data1$ps_kb), max(data1$ps_kb), by = (max(data1$ps_kb)-min(data1$ps_kb))/3),las=1,cex.axis = 1)
axis(side = 2, at=seq(0, 7.6, by = 2.5),las=1,cex.axis = 1)
dev.off()
```
### Promoter analysis

```
setwd("~/Library/CloudStorage/OneDrive-UniversityofVermont/PENN STATE/Awal_aayudh")
list.files()

ABA <- read.csv("1.Cycles_ABA_final.csv")
head(ABA)
dim(ABA)
rs_split <- data.frame(do.call("rbind", strsplit(as.character(ABA$rs), "_", fixed = TRUE)))
head(rs_split)
ABA$chr <- as.numeric(rs_split$X2)
head(ABA)
summary(ABA)



ABA_minus <- subset(ABA, direction == "-")
head(ABA_minus)
ABA_minus_1 <- data.frame(ABA_minus[,1:4], ABA_minus$sobicID, ABA_minus$chr)
head(ABA_minus_1)
names(ABA_minus_1)[6]<-paste("chr")
names(ABA_minus_1)[5]<-paste("sobicID")
head(ABA_minus_1)

genes_minus_CDS <- read.csv("Entire_genome_minus_strand_CDS.csv")
head(genes_minus_CDS)

rs_split <- data.frame(do.call("rbind", strsplit(as.character(genes_minus_CDS$Gene_name), ".", fixed = TRUE)))
rs_split1 <- data.frame(do.call("rbind", strsplit(as.character(rs_split$X1), "=", fixed = TRUE)))

genes_minus_CDS$sobicID <- paste(rs_split1$X2, rs_split$X2, sep=".")
head(genes_minus_CDS)

merge_minus <- merge(genes_minus_CDS,ABA_minus_1, by="sobicID")
head(merge_minus)
library(data.table)
ABA_minus_final <- unique(setDT(merge_minus)[order(rs, -rs)], by = "rs")
ABA_minus_final$snpsite <- ABA_minus_final$stop - ABA_minus_final$ps
head(ABA_minus_final)


ABA_minus_final$logp <- -log10(ABA_minus_final$p_wald)
ABA_minus_final[, "snpsite"] <- (ABA_minus_final$snpsite/1000)

setwd("~/Library/CloudStorage/OneDrive-UniversityofVermont/PENN STATE/Awal_aayudh")
write.csv(ABA_minus_final,"ABA_minus_promoter.csv")


####plus
ABA_plus <- subset(ABA, direction == "+")
head(ABA_plus)
ABA_plus_1 <- data.frame(ABA_plus[,1:4], ABA_plus$sobicID,ABA_plus$chr)
names(ABA_plus_1)[6]<-paste("chr")
names(ABA_plus_1)[5]<-paste("sobicID")
head(ABA_plus_1)


genes_plus_CDS <- read.csv("Entire_genome_plus_strand_CDS.csv")
head(genes_plus_CDS)

rs_split <- data.frame(do.call("rbind", strsplit(as.character(genes_plus_CDS$Gene_name), ".", fixed = TRUE)))
rs_split1 <- data.frame(do.call("rbind", strsplit(as.character(rs_split$X1), "=", fixed = TRUE)))

genes_plus_CDS$sobicID <- paste(rs_split1$X2, rs_split$X2, sep=".")
head(genes_plus_CDS)

merge_plus <- merge(genes_plus_CDS,ABA_plus_1, by="sobicID")
head(merge_plus)
library(data.table)
ABA_plus_final <- unique(setDT(merge_plus)[order(rs, -rs)], by = "rs")
ABA_plus_final$snpsite <- ABA_plus_final$stop - ABA_plus_final$ps
head(ABA_plus_final)

ABA_plus_final$logp <- -log10(ABA_plus_final$p_wald)
ABA_plus_final[, "snpsite"] <- (ABA_plus_final$snpsite/1000)

write.csv(ABA_plus_final,"ABA_plus_promoter.csv")

library(tidyverse)
library(ggtext)
library(normentR)
library(ggplot2)

plus <- ggplot(ABA_plus_final, aes(x=snpsite, y=logp)) + 
  geom_point() +
  geom_segment(data=ABA_plus_final,aes(x=snpsite, xend=snpsite,y=0,yend=logp),
               color = '#B4AF46',alpha=0.8,size=0.5) + 
  labs(x=NULL,y=NULL)+
  scale_x_continuous(breaks=seq(-15,19,5))+
  scale_y_continuous(breaks=seq(3,11,2))+
  annotate(geom="text", x=-2.05, y=11, label="TSS",color="#B4464B")+
  theme_classic() +
  geom_vline(xintercept = 0, linetype="dashed", color = "#B4464B", size=1)+
  theme(legend.position="none",text=element_text(size=16, colour = "black", family="Times New Roman"),
        axis.line = element_line(size=0.5, colour = "black"),
        axis.text.x=element_text(colour="black", size = 16),
        axis.text.y=element_text(colour="black", size = 16))
plus




minus <- ggplot(ABA_minus_final,aes(x=snpsite, y=logp)) + 
  geom_point() +
  geom_segment(data=ABA_minus_final,aes(x=snpsite, xend=snpsite,y=0,yend=logp),
            color = '#B4AF46',alpha=0.8,size=0.5) +
  labs(x=NULL,y=NULL)+
  scale_x_reverse(breaks=seq(-14,15,5))+
  annotate(geom="text", x=-1.55, y=11, label="TTS",color="#B4464B")+
  theme_classic() +
  geom_vline(xintercept = 0, linetype="dashed", color = "#B4464B", size=1)+
  theme(legend.position="none",text=element_text(size=16, colour = "black", family="Times New Roman"),
        axis.line.y = element_blank(),
        axis.text.x=element_text(colour="black", size = 16),
        axis.text.y=element_blank(),
        axis.ticks.y = element_blank())
minus

library(tidyverse)
library(gridExtra)
library(grid)
library(gridtext)

setwd("~/Library/CloudStorage/OneDrive-UniversityofVermont/PENN STATE/Awal_aayudh")
tiff("ABA_promoter_plot.tiff", width = 6, height = 6, units = 'in', res = 300)
library(gridExtra)
grid.arrange(plus, minus, ncol = 2)
dev.off()

```
-----
<div id='id-section14'/>

## Chapter 14: Sampling (closely related but high HS difference) 

### Getting distance matrix and NJ tree

```
setwd("~/Library/CloudStorage/OneDrive-UniversityofVermont/PENN STATE/GBS data processing /Hu_GBS_2019")
library(SNPRelate)
library(gdsfmt)
genofile <- snpgdsOpen('SNPs_imp.recode.gds')
sample.id <- read.gdsn(index.gdsn(genofile, 'sample.id'))
length(sample.id)

GBS_sample <- as.data.frame(sample.id)
  
rs_split <- data.frame(do.call("rbind", strsplit(as.character(GBS_sample$sample.id), ".", fixed = TRUE)))
head(rs_split)

GBS_sample$Accessions <- rs_split$X1
head(GBS_sample)

setwd("~/Library/CloudStorage/OneDrive-UniversityofVermont/PENN STATE/joel sampling")
list.files()

landraces <- read.csv("Accessions_Georef_Final.csv")  
head(landraces)
names(landraces)[2]<-paste("Lat")
names(landraces)[3]<-paste("Lon")

library(maps)
landraces$Country <- map.where("world", landraces$Lon, landraces$Lat)
head(landraces)

rs_split <- data.frame(do.call("rbind", strsplit(as.character(landraces$Accession_Info), "_", fixed = TRUE)))
head(rs_split)

landraces$Accessions <- paste0(rs_split$X1,rs_split$X2)

sampleinfo <- merge(landraces,GBS_sample,merge = "sample.id")
library(data.table)
sampleinfo <- unique(setDT(sampleinfo)[order(sample.id, -sample.id)], by = "sample.id")
head(sampleinfo)
sampleinfo <- sampleinfo[3:377,] ##removed india and sri lanka

sampleinfo <- rbind(sampleinfo[1:17,],sampleinfo[19:375,]) ##removed pakistan

setwd("~/Library/CloudStorage/OneDrive-UniversityofVermont/PENN STATE/joel sampling")
write.csv(sampleinfo,"1.GBS_sampleID_joel.csv", row.names = FALSE)

#LD based SNP pruning, instead of snpgdsSelectSNP you can use snpgdsLDpruning to get unliked SNPs
#You can look into specific LD thresholds for Sorghum
set.seed(1000)
snpset = snpgdsLDpruning(genofile,maf=0.05,missing.rate=0.5,ld.threshold=0.4, slide.max.bp=50000)
snp.id = unlist(snpset)

setwd("~/Library/CloudStorage/OneDrive-UniversityofVermont/PENN STATE/joel sampling")
#here's a dataframe with some traits and sample.id matching genofile
sampleinfo <- data.frame(read.csv('1.GBS_sampleID_joel.csv', header=TRUE))
head(sampleinfo)
sampleinfo$sample.id

#you can use either an ibs or dissimilarity matrix for the nj tree
ibs  =  snpgdsIBS(genofile,sample.id=sampleinfo$sample.id,snp.id=snp.id,maf=0.05,missing.rate=0.5,num.thread=4,verbose=TRUE)
dis  =  snpgdsDiss(genofile,sample.id=sampleinfo$sample.id,snp.id=snp.id,maf=0.05,missing.rate=0.5,num.thread=4,verbose=TRUE)

#install.packages("ape", dependencies = T)
library(ape)
#devtools::install_github("KlausVigo/phangorn")
library(phangorn)
library(phyclust)
library(phyloch)
ibs.matrix <- ibs$ibs
dim(ibs.matrix)
#sometimes the order of the ibs matrix ends up not matching your trait dataframe, so you need to reorder
order <- as.data.frame(ibs$sample.id)
dim(order)
colnames(order) <- c("sample.id")

order.idx <- match(order$sample.id,sampleinfo$sample.id)
order.idx
length(order.idx)

ordered_ibs <- sampleinfo[order.idx,]

row.names(ibs.matrix) <- ordered_ibs$Accessions
colnames(ibs.matrix) <- ordered_ibs$Accessions

#this was me comparing  
#the distance matrix from ibs
mm <- 1-ibs.matrix

matrix_data <- as.data.frame(mm)

library('ggplot2')
library('ape')
library('phangorn')
library('dplyr')
library('ggtree')
library('phylobase')
setwd("~/Library/CloudStorage/OneDrive-UniversityofVermont/PENN STATE/joel sampling")
write.table(matrix_data, "1.GBS_distance_matrix.csv", sep=",")

pdf(file = 'ibs_nj_plain.pdf', width = 40, height = 40)
my_nj <- ape::nj(mm)
plot(my_nj)
dev.off()

pdf(file = 'ibs_nj_phangorn.pdf', width = 10, height = 40)
my_upgma <- phangorn::upgma(mm)
par(family="Times")
plot(my_upgma, cex=0.7)
dev.off()
```
### Getting genetic distance of the closely related pairs
```
setwd("~/Library/CloudStorage/OneDrive-UniversityofVermont/PENN STATE/joel sampling")
list.files()

data <- read.csv("1.GBS_distance_matrix.csv")
head(data)

row.names(data) <- data[,1]
data <- data[,2:378]


#------looping--------------------

dddf_1 <- NULL

for(i in 1:ncol(data)) { 
  dddf <- NULL
  dddf$accession_1 <- print(colnames(data[i]))
  dddf$accession_2 = rownames(data)[which(data[,colnames(data[i])] == min(data[,i][which(data[,i]>0)]))]
  dddf$genedist = min(data[,i][which(data[,i]>0)])
  
  dddf_1 <- rbind(dddf_1, dddf)
}
setwd("~/Library/CloudStorage/OneDrive-UniversityofVermont/PENN STATE/joel sampling/map")
write.csv(dddf_1,"pairs_closest_genetic_distance.csv", row.names = FALSE)
```
### Get HS score difference of the closely related pairs

```
setwd("~/Library/CloudStorage/OneDrive-UniversityofVermont/PENN STATE/joel sampling/map")
list.files()
pairs <- read.csv("pairs_closest_genetic_distance.csv")
head(pairs)

setwd("~/Library/CloudStorage/OneDrive-UniversityofVermont/PENN STATE/joel sampling")
sampleinfo <- data.frame(read.csv('1.GBS_sampleID_joel.csv', header=TRUE))
head(sampleinfo)

names(sampleinfo)[1]<-paste("accession_1")

ac1 <- pairs[,1]
ac1 <- as.data.frame(ac1)
names(ac1)[1]<-paste("accession_1")

order.idx <- match(ac1$accession_1, sampleinfo$accession_1)
order.idx
ordered <- sampleinfo[order.idx,]
head(ordered)

names(ordered )[3]<-paste("Lat")
names(ordered )[4]<-paste("Lon")
head(ordered)
library(raster)
library(sp)
library(rgdal)
library(tidyverse)

setwd("~/OneDrive - University of Vermont/PENN STATE/RAstor data")
list.files()

preds.all <- raster(paste0("~/OneDrive - University of Vermont/PENN STATE/RAstor data/preds.all.tif"))
pairs$accession1_HS <- raster::extract(preds.all, ordered[,c("Lon","Lat")])
head(pairs)

###HS score for accession 2
setwd("~/Library/CloudStorage/OneDrive-UniversityofVermont/PENN STATE/joel sampling")
sampleinfo <- data.frame(read.csv('1.GBS_sampleID_joel.csv', header=TRUE))
head(sampleinfo)

names(sampleinfo)[1]<-paste("accession_2")

ac2 <- pairs[,2]
ac2 <- as.data.frame(ac2)
names(ac2)[1]<-paste("accession_2")

order.idx_2 <- match(ac2$accession_2, sampleinfo$accession_2)
order.idx_2
ordered_2 <- sampleinfo[order.idx_2,]
head(ordered_2)

names(ordered_2)[3]<-paste("Lat")
names(ordered_2)[4]<-paste("Lon")
head(ordered_2)
library(raster)
library(sp)
library(rgdal)
library(tidyverse)

setwd("~/OneDrive - University of Vermont/PENN STATE/RAstor data")
list.files()

preds.all <- raster(paste0("~/OneDrive - University of Vermont/PENN STATE/RAstor data/preds.all.tif"))
pairs$accession2_HS <- raster::extract(preds.all, ordered_2[,c("Lon","Lat")])

head(pairs)
pairs <- na.omit(pairs)
pairs$HS_diff <- pairs$accession1_HS - pairs$accession2_HS
head(pairs)
pairs[,"HS_diff"] <- abs(pairs$HS_diff)
pairs_final <- pairs[with(pairs,order(-HS_diff)),]
head(pairs_final)
pairs_w_HS <- subset(pairs_final, HS_diff > 0.3)
head(pairs_w_HS)

setwd("~/Library/CloudStorage/OneDrive-UniversityofVermont/PENN STATE/joel sampling")
write.csv(pairs_w_HS, "1.Pairs_with_HS score0.3.csv", row.names = FALSE)

```

### histogram to see the distribution of the data

```
setwd("~/Library/CloudStorage/OneDrive-UniversityofVermont/PENN STATE/joel sampling")
data <- read.csv("1.Pairs_with_HS score0.3.csv")

ac1 <- data.frame(data$accession_1,data$accession1_HS)
names(ac1 )[1]<-paste("Accessions")
names(ac1 )[2]<-paste("HS_score")
ac2 <- data.frame(data$accession_2,data$accession2_HS)
names(ac2 )[1]<-paste("Accessions")
names(ac2 )[2]<-paste("HS_score")

ac1_2 <- rbind(ac1,ac2)

library(data.table)
ac1_2 <- unique(setDT(ac1_2)[order(Accessions, -Accessions)], by = "Accessions")


#check histogram_HS_score
setwd("~/Library/CloudStorage/OneDrive-UniversityofVermont/PENN STATE/joel sampling")
tiff("Joel_sampling_77pairs.tiff", width = 5, height = 5, units = 'in', res = 300)
hist(ac1_2$HS_score, prob = TRUE, main = NA, xlab="HS Score",col = "lightblue")
x <- seq(min(ac1_2$HS_score), max(ac1_2$HS_score), length = 40)
f <- dnorm(x, mean = mean(ac1_2$HS_score), sd = sd(ac1_2$HS_score))
lines(x, f, col = "red", lwd = 2)
dev.off()
```

### Creating Map and Haversine_distance

```
setwd("~/Library/CloudStorage/OneDrive-UniversityofVermont/PENN STATE/joel sampling")
data <- read.csv("1.Pairs_with_HS score0.3.csv")
head(data)

ac1 <- data.frame(data$accession_1,data$accession1_HS)
names(ac1 )[1]<-paste("Accessions")
names(ac1 )[2]<-paste("HS_score")
ac2 <- data.frame(data$accession_2,data$accession2_HS)
names(ac2 )[1]<-paste("Accessions")
names(ac2 )[2]<-paste("HS_score")

sampleinfo <- data.frame(read.csv('1.GBS_sampleID_joel.csv', header=TRUE))
head(sampleinfo)

ac1_cord <- merge(ac1,sampleinfo, by ="Accessions")
head(ac1_cord)
ac1_cord_order <- match(ac1$Accessions, ac1_cord$Accessions)
ac1_cord_1 <- ac1_cord [ac1_cord_order,]
head(ac1_cord_1)
names(ac1_cord_1)[4]<-paste("lat1")
names(ac1_cord_1)[5]<-paste("lon1")

ac2_cord <- merge(ac2,sampleinfo, by ="Accessions")
head(ac2_cord)
ac2_cord_order <- match(ac2$Accessions, ac2_cord$Accessions)
ac2_cord_1 <- ac2_cord [ac2_cord_order,]
head(ac2_cord_1)
names(ac2_cord_1)[4]<-paste("lat2")
names(ac2_cord_1)[5]<-paste("lon2")

data1 <- cbind(ac1_cord_1,ac2_cord_1)
setwd("~/Library/CloudStorage/OneDrive-UniversityofVermont/PENN STATE/joel sampling/map")
write.csv(data1,"twocords_csv",row.names = FALSE)

all <- rbind(ac1_cord,ac2_cord)
head(all)

library(data.table)
all <- unique(setDT(all)[order(Accessions, -Accessions)], by = "Accessions")

setwd("~/Library/CloudStorage/OneDrive-UniversityofVermont/PENN STATE/joel sampling/map")
library(maps)
library(ggplot2)
library(wesanderson)

mid<-mean(all$HS_score)

base_world <- map_data("world")
setwd("~/Library/CloudStorage/OneDrive-UniversityofVermont/PENN STATE/joel sampling/map")
tiff("77pairs.tiff", width = 6, height = 6, units = 'in', res = 300)
ggplot() + 
  geom_map(data = base_world, map = base_world,aes(long, lat, map_id = region),
           color = "black", fill = "white", size = 0.5) +
  geom_point(data = all,size=2,aes(longitude, latitude,color = HS_score),alpha = 1) + #set the color outside of `aes`
  theme(text = element_text(size=20), legend.position="none") +
  coord_sf(xlim = c(-20, 50), ylim = c(-35, 35), expand = FALSE)+
  scale_color_gradient2(midpoint=mid, low="blue", mid="white", high="red", space ="Lab" )+
  theme_classic() +
  theme(legend.position="right",text=element_text(size=16, colour = "black", family="Times New Roman"),
        axis.line = element_line(size=0.5, colour = "black"),
        axis.text.x=element_text(colour="black", size = 16),
        axis.text.y=element_text(colour="black", size = 16)) 
dev.off()

###Getting distance
setwd("~/Library/CloudStorage/OneDrive-UniversityofVermont/PENN STATE/joel sampling/map")
list.files()

data1 <- read.csv("twocords_csv")
head(data1)

library(geosphere)

dddf_1 <- NULL

for(i in 1:nrow(data1)) { 
  dddf <- NULL
  dddf$dist <- distHaversine(c(data1[i,5], data1[i,4]), c(data1[i,11], data1[i,10]), r=6378137) 
  
  dddf_1 <- rbind(dddf_1, dddf)
}

df <- as.numeric(dddf_1,row.names = FALSE)
head(df)
fd <- as.data.frame(df)

final <- cbind(data,fd)
names(final)[6]<-paste("Haversine_distance")
head(final)

summary(final)

PI514390
write.csv(final, "Haversine_distance.csv", row.names = FALSE)

```
### Correlatios netween HS diff vs genetic distance of closely related pairs
```
setwd("~/Library/CloudStorage/OneDrive-UniversityofVermont/PENN STATE/joel sampling")
list.files()

data <- read.csv("1.Pairs_with_HS score0.3.csv")
head(data)

library(ggplot2)
setwd("~/Library/CloudStorage/OneDrive-UniversityofVermont/PENN STATE/joel sampling/correlations")

tiff("Closely related_cor.tiff", width = 5, height = 5, units = 'in', res = 300)
ggplot()+
  geom_point(aes(x=data$genedist,y=data$HS_diff),size=2,colour="black")+
  labs(x="Genetic Distance",y="HS Score diff.")+
  geom_smooth(aes(x=data$genedist,y=data$HS_diff),se=FALSE, method="lm", span=0.75)+
  theme_classic() +
  scale_x_continuous(breaks=seq(0,0.21,0.07))+
  scale_y_continuous(breaks=seq(0.30,0.77,0.15))+
  annotate(geom="text", x=0.20, y=0.77, label="Cor = 0.154",color="#B4464B")+
  theme(legend.position="none",text=element_text(size=16, colour = "black", family="Times New Roman"),
        axis.line = element_line(size=0.5, colour = "black"),
        axis.text.x=element_text(colour="black", size = 16),
        axis.text.y=element_text(colour="black", size = 16),
        panel.border = element_rect(colour = "black", fill=NA, size=1))
dev.off()
```


### Creating random pairs and HS diff vs genetic distance

```
setwd("~/Library/CloudStorage/OneDrive-UniversityofVermont/PENN STATE/joel sampling")
sampleinfo <- data.frame(read.csv('1.GBS_sampleID_joel.csv', header=TRUE))
head(sampleinfo)

data <- read.csv("1.GBS_distance_matrix.csv")
head(data)

data1 <- as.data.frame(data$PI218112)

data_acc <- as.data.frame(sampleinfo[,1])

dddf_1 <- NULL
for(i in 87) { 
  dddf <- NULL
  dddf$accession_1 <- (data_acc[sample(nrow(data_acc), i), ])
  dddf$accession_2 <- (data_acc[sample(nrow(data_acc), i), ])
  
  dddf_1 <- rbind(dddf_1, dddf)
}

random_pairs <- as.data.frame(dddf)
head(random_pairs)

random_pairs$genedist <- (data1[sample(nrow(data1), 87), ])

head(random_pairs)

names(sampleinfo)[1]<-paste("accession_1")

df1 <- merge(random_pairs, sampleinfo, by ="accession_1")
names(df1)[5]<-paste("Lat")
names(df1)[6]<-paste("Lon")
head(df1)
library(raster)
library(sp)
library(rgdal)
library(tidyverse)

setwd("~/OneDrive - University of Vermont/PENN STATE/RAstor data")
list.files()

preds.all <- raster(paste0("~/OneDrive - University of Vermont/PENN STATE/RAstor data/preds.all.tif"))
random_pairs$accession1_HS <- raster::extract(preds.all, df1[,c("Lon","Lat")])

###HS score for accession 2
setwd("~/Library/CloudStorage/OneDrive-UniversityofVermont/PENN STATE/joel sampling")
sampleinfo <- data.frame(read.csv('1.GBS_sampleID_joel.csv', header=TRUE))
head(sampleinfo)

names(sampleinfo)[1]<-paste("accession_2")

df2 <- merge(random_pairs, sampleinfo, by ="accession_2")
names(df2)[6]<-paste("Lat")
names(df2)[7]<-paste("Lon")
head(df2)
library(raster)
library(sp)
library(rgdal)
library(tidyverse)

setwd("~/OneDrive - University of Vermont/PENN STATE/RAstor data")
list.files()
preds.all <- raster(paste0("~/OneDrive - University of Vermont/PENN STATE/RAstor data/preds.all.tif"))
random_pairs$accession2_HS <- raster::extract(preds.all, df2[,c("Lon","Lat")])

head(random_pairs)
pairs <- na.omit(random_pairs)
pairs$HS_diff <- pairs$accession1_HS - pairs$accession2_HS
head(pairs)
pairs_final <- pairs[with(pairs,order(HS_diff)),]
head(pairs_final)
pairs_final[,"HS_diff"] <- abs(pairs_final$HS_diff)

pairs_final <- pairs_final[1:77,]
head(pairs_final)

setwd("~/Library/CloudStorage/OneDrive-UniversityofVermont/PENN STATE/joel sampling/random pair")
write.csv(pairs_final,"random_pairs_final.csv",row.names = FALSE)
```


### Correlation between all possible pairs from the matrix and its HS diff vs genetic distance

```
setwd("~/Library/CloudStorage/OneDrive-UniversityofVermont/PENN STATE/joel sampling")
###All pairs
data <- read.csv("1.GBS_distance_matrix.csv")
head(data)

row.names(data) <- data[,1]
data <- data[,2:378]

df <- as.matrix(sapply(data, as.numeric))
row.names(df) <- row.names(data)

library(tidyr)
all_pairs <- crossing(var1 = row.names(data), var2 = row.names(data))

fd1 <- as.data.frame(all_pairs)
head(fd1)


library(phylin)
library(geometry)
##https://rdrr.io/cran/phylin/man/extract.var.html 
fd1$genedist <- extract.val(df, fd1) 
head(fd1)

library(dplyr)
All_pairs_values <- filter(fd1, pair.data > 0)
names(All_pairs_values)[1]<-paste("accession_1")
names(All_pairs_values)[2]<-paste("accession_2")
head(All_pairs_values)

###get HS scores
setwd("~/Library/CloudStorage/OneDrive-UniversityofVermont/PENN STATE/joel sampling")
sampleinfo <- data.frame(read.csv('1.GBS_sampleID_joel.csv', header=TRUE))
head(sampleinfo)

names(sampleinfo)[1]<-paste("accession_1")

ac1 <- All_pairs_values[,1]
ac1 <- as.data.frame(ac1)
names(ac1)[1]<-paste("accession_1")

order.idx <- match(ac1$accession_1, sampleinfo$accession_1)
order.idx
ordered <- sampleinfo[order.idx,]
head(ordered)

names(ordered )[3]<-paste("Lat")
names(ordered )[4]<-paste("Lon")

library(raster)
library(sp)
library(rgdal)
library(tidyverse)

setwd("~/OneDrive - University of Vermont/PENN STATE/RAstor data")
list.files()

preds.all <- raster(paste0("~/OneDrive - University of Vermont/PENN STATE/RAstor data/preds.all.tif"))
All_pairs_values$accession1_HS <- raster::extract(preds.all, ordered[,c("Lon","Lat")])
head(All_pairs_values)



###HS score for accession 2
setwd("~/Library/CloudStorage/OneDrive-UniversityofVermont/PENN STATE/joel sampling")
sampleinfo <- data.frame(read.csv('1.GBS_sampleID_joel.csv', header=TRUE))
head(sampleinfo)

names(sampleinfo)[1]<-paste("accession_2")

ac2 <- All_pairs_values[,2]
ac2 <- as.data.frame(ac2)
names(ac2)[1]<-paste("accession_2")

order.idx_2 <- match(ac2$accession_2, sampleinfo$accession_2)
order.idx_2
ordered_2 <- sampleinfo[order.idx_2,]
head(ordered_2)

names(ordered_2)[3]<-paste("Lat")
names(ordered_2)[4]<-paste("Lon")
head(ordered_2)
library(raster)
library(sp)
library(rgdal)
library(tidyverse)

setwd("~/OneDrive - University of Vermont/PENN STATE/RAstor data")
list.files()
preds.all <- raster(paste0("~/OneDrive - University of Vermont/PENN STATE/RAstor data/preds.all.tif"))
All_pairs_values$accession2_HS <- raster::extract(preds.all, ordered_2[,c("Lon","Lat")])

head(All_pairs_values)
pairs <- na.omit(All_pairs_values)
pairs$HS_diff <- pairs$accession1_HS - pairs$accession2_HS
pairs[,"HS_diff"] <- abs(pairs$HS_diff)
pairs <- filter(pairs, HS_diff > 0)
head(pairs)
pairs_final <- pairs[with(pairs,order(-HS_diff)),]
head(pairs_final)

setwd("~/Library/CloudStorage/OneDrive-UniversityofVermont/PENN STATE/joel sampling/all_pairs")
tiff("all_pairs.tiff", width = 5, height = 5, units = 'in', res = 300)
ggplot()+
  geom_point(aes(x=pairs_final$genedist,y=pairs_final$HS_diff),size=1,colour="black")+
  labs(x="Genetic Distance",y="HS Score diff.")+
  geom_smooth(aes(x=pairs_final$genedist,y=pairs_final$HS_diff),se=FALSE, method="lm", span=0.75)+
  theme_classic() +
  scale_x_continuous(breaks=seq(0.05,0.35,0.10))+
  scale_y_continuous(breaks=seq(0,0.90,0.30))+
  annotate(geom="text", x=0.34, y=0.89, label="Cor = 0.049",color="#B4464B")+
  theme(legend.position="none",text=element_text(size=16, colour = "black", family="Times New Roman"),
        axis.line = element_line(size=0.5, colour = "black"),
        axis.text.x=element_text(colour="black", size = 16),
        axis.text.y=element_text(colour="black", size = 16),
        panel.border = element_rect(colour = "black", fill=NA, size=1))
dev.off()

cor.test(pairs_final$genedist,pairs_final$HS_diff, method="pearson")

```
``
-----
<div id='id-section15'/>

## Chapter 15: Redundancy analysis (RDA) for drought

### Getting climate variables

```
setwd("~/Library/CloudStorage/OneDrive-UniversityofVermont/PENN STATE/acessions and cords/env variables")
cycles <- read.csv("1. Cycles_new.csv")
head(cycles)
##Add bioclim variables

data <- data.frame(cycles[,1:4],cycles[7:8])

setwd("~/Library/CloudStorage/OneDrive-UniversityofVermont/PENN STATE/bioclim variables/wc2.1_30s_bio")
list.files()
library(raster)
library(sp)
library(rgdal)
library(tidyverse)

Annual_Mean_Temperature <- raster(paste0("~/Library/CloudStorage/OneDrive-UniversityofVermont/PENN STATE/bioclim variables/wc2.1_30s_bio/wc2.1_30s_bio_1.tif"))
data$Annual_Mean_Temperature <- raster::extract(Annual_Mean_Temperature, data[,c("Lon","Lat")])

Max_Temperature_of_Warmest_Month <- raster(paste0("~/Library/CloudStorage/OneDrive-UniversityofVermont/PENN STATE/bioclim variables/wc2.1_30s_bio/wc2.1_30s_bio_5.tif"))
data$Max_Temperature_of_Warmest_Month <- raster::extract(Max_Temperature_of_Warmest_Month, data[,c("Lon","Lat")])

Mean_Temp_of_Wettest_Quarter <- raster(paste0("~/Library/CloudStorage/OneDrive-UniversityofVermont/PENN STATE/bioclim variables/wc2.1_30s_bio/wc2.1_30s_bio_8.tif"))
data$Mean_Temp_of_Wettest_Quarter <- raster::extract(Mean_Temp_of_Wettest_Quarter , data[,c("Lon","Lat")])

Mean_Temp_Driest_Quarter <- raster(paste0("~/Library/CloudStorage/OneDrive-UniversityofVermont/PENN STATE/bioclim variables/wc2.1_30s_bio/wc2.1_30s_bio_9.tif"))
data$Mean_Temp_Driest_Quarter <- raster::extract(Mean_Temp_Driest_Quarter, data[,c("Lon","Lat")])

MeanTemp_Warmest_Quarter <- raster(paste0("~/Library/CloudStorage/OneDrive-UniversityofVermont/PENN STATE/bioclim variables/wc2.1_30s_bio/wc2.1_30s_bio_10.tif"))
data$MeanTemp_Warmest_Quarter <- raster::extract(MeanTemp_Warmest_Quarter, data[,c("Lon","Lat")])

Annual_Precipitation <- raster(paste0("~/Library/CloudStorage/OneDrive-UniversityofVermont/PENN STATE/bioclim variables/wc2.1_30s_bio/wc2.1_30s_bio_12.tif"))
data$Annual_Precipitation <- raster::extract(Annual_Precipitation, data[,c("Lon","Lat")])

Precipitation_Wettest_Month <- raster(paste0("~/Library/CloudStorage/OneDrive-UniversityofVermont/PENN STATE/bioclim variables/wc2.1_30s_bio/wc2.1_30s_bio_13.tif"))
data$Precipitation_Wettest_Month <- raster::extract(Precipitation_Wettest_Month, data[,c("Lon","Lat")])

Precipitation_Driest_Month <- raster(paste0("~/Library/CloudStorage/OneDrive-UniversityofVermont/PENN STATE/bioclim variables/wc2.1_30s_bio/wc2.1_30s_bio_14.tif"))
data$Precipitation_Driest_Month <- raster::extract(Precipitation_Driest_Month, data[,c("Lon","Lat")])
p <- (data$Precipitation_Driest_Month + 1)
data$Precipitation_Driest_Month_log <- log(p)

Precipitation_Seasonality <- raster(paste0("~/Library/CloudStorage/OneDrive-UniversityofVermont/PENN STATE/bioclim variables/wc2.1_30s_bio/wc2.1_30s_bio_15.tif"))
data$Precipitation_Seasonality <- raster::extract(Precipitation_Seasonality, data[,c("Lon","Lat")])

Precipitation_Wettest_Quarter <- raster(paste0("~/Library/CloudStorage/OneDrive-UniversityofVermont/PENN STATE/bioclim variables/wc2.1_30s_bio/wc2.1_30s_bio_16.tif"))
data$Precipitation_Wettest_Quarter <- raster::extract(Precipitation_Wettest_Quarter, data[,c("Lon","Lat")])

Precipitation_Driest_Quarter <- raster(paste0("~/Library/CloudStorage/OneDrive-UniversityofVermont/PENN STATE/bioclim variables/wc2.1_30s_bio/wc2.1_30s_bio_17.tif"))
data$Precipitation_Driest_Quarter <- raster::extract(Precipitation_Driest_Quarter, data[,c("Lon","Lat")])
q <- (data$Precipitation_Driest_Quarter +1)
data$Precipitation_Driest_Quarter_log <- log(q)

Precipitation_Warmest_Quarter <- raster(paste0("~/Library/CloudStorage/OneDrive-UniversityofVermont/PENN STATE/bioclim variables/wc2.1_30s_bio/wc2.1_30s_bio_18.tif"))
data$Precipitation_Warmest_Quarter <- raster::extract(Precipitation_Warmest_Quarter, data[,c("Lon","Lat")])

### soil moisture capacity derived from soil texture  - http://webmap.ornl.gov/wcsdown/wcsdown.jsp?dg_id=545_1
setwd("~/Library/CloudStorage/OneDrive-UniversityofVermont/PENN STATE/bioclim variables/soil data")
soil_m <-raster('~/Library/CloudStorage/OneDrive-UniversityofVermont/PENN STATE/bioclim variables/soil data/wrtext.asc')

data$soil_m<- raster::extract(soil_m, data[,c("Lon","Lat")])
data$soil_m_log <- log10(data$soil_m)

AI <- raster(paste0("~/OneDrive - University of Vermont/PENN STATE/Aridity_index_gwas/ai_v3_yr.tif"))
Global_AI <- raster::extract(AI, data[,c("Lon","Lat")])
data$Aridity_Index <- ((Global_AI)*0.0001)

head(data)

setwd("~/Library/CloudStorage/OneDrive-UniversityofVermont/PENN STATE/bioclim variables")


library(maps)
data$country <- map.where("world", data$Lon, data$Lat)
head(data)
data <- na.omit (data)
data_1 <- data %>% relocate(country, .before = water_veg_avg_log)
setwd("~/Library/CloudStorage/OneDrive-UniversityofVermont/PENN STATE/bioclim variables")
write.csv(data_1, "1.Drought_variables.csv", row.names = FALSE)
data <-  read.csv("Drought_bioclim_variables_30sec.csv")
head(data)

setwd("~/Library/CloudStorage/OneDrive-UniversityofVermont/PENN STATE/bioclim variables")
png("clim_var_histo.png", width = 12, height = 10, units = 'in', res = 300)
par(mfrow=c(4,4))

hist(data$water_veg_avg_log, prob = TRUE, main = NA, xlab="log (Water stress veg.)",col = "ivory4")
x <- seq(min(data$water_veg_avg_log), max(data$water_veg_avg_log), length = 40)
f <- dnorm(x, mean = mean(data$water_veg_avg_log), sd = sd(data$water_veg_avg_log))
lines(x, f, col = "red", lwd = 2)

hist(data$water_repro_adj_avg_log, prob = TRUE, main = NA, xlab="log (Water stress veg.)",col = "ivory4")
x <- seq(min(data$water_repro_adj_avg_log), max(data$water_repro_adj_avg_log), length = 40)
f <- dnorm(x, mean = mean(data$water_repro_adj_avg_log), sd = sd(data$water_repro_adj_avg_log))
lines(x, f, col = "red", lwd = 2)

hist(data$Aridity_Index , prob = TRUE, main = NA, xlab="Aridity_Index",col = "ivory4")
x <- seq(min(data$Aridity_Index ), max(data$Aridity_Index ), length = 40)
f <- dnorm(x, mean = mean(data$Aridity_Index ), sd = sd(data$Aridity_Index ))
lines(x, f, col = "red", lwd = 2)

hist(data$Annual_Mean_Temperature, prob = TRUE, main = NA, xlab="Annual_Mean_Temperature",col = "lightblue")
x <- seq(min(data$Annual_Mean_Temperature), max(data$Annual_Mean_Temperature), length = 40)
f <- dnorm(x, mean = mean(data$Annual_Mean_Temperature), sd = sd(data$Annual_Mean_Temperature))
lines(x, f, col = "red", lwd = 2)

hist(data$Max_Temperature_of_Warmest_Month, prob = TRUE, main = NA, xlab="Max_Temperature_of_Warmest_Month",col = "lightblue")
x <- seq(min(data$Max_Temperature_of_Warmest_Month), max(data$Max_Temperature_of_Warmest_Month), length = 40)
f <- dnorm(x, mean = mean(data$Max_Temperature_of_Warmest_Month), sd = sd(data$Max_Temperature_of_Warmest_Month))
lines(x, f, col = "red", lwd = 2)

hist(data$Mean_Temp_of_Wettest_Quarter, prob = TRUE, main = NA, xlab="Mean_Temp_of_Wettest_Quarter",col = "lightblue")
x <- seq(min(data$Mean_Temp_of_Wettest_Quarter), max(data$Mean_Temp_of_Wettest_Quarter), length = 40)
f <- dnorm(x, mean = mean(data$Mean_Temp_of_Wettest_Quarter), sd = sd(data$Mean_Temp_of_Wettest_Quarter))
lines(x, f, col = "red", lwd = 2)

hist(data$Mean_Temp_Driest_Quarter, prob = TRUE, main = NA, xlab="Mean_Temp_Driest_Quarter",col = "lightblue")
x <- seq(min(data$Mean_Temp_Driest_Quarter), max(data$Mean_Temp_Driest_Quarter), length = 40)
f <- dnorm(x, mean = mean(data$Mean_Temp_Driest_Quarter), sd = sd(data$Mean_Temp_Driest_Quarter))
lines(x, f, col = "red", lwd = 2)

hist(data$MeanTemp_Warmest_Quarter, prob = TRUE, main = NA, xlab="MeanTemp_Warmest_Quarter",col = "lightblue")
x <- seq(min(data$MeanTemp_Warmest_Quarter), max(data$MeanTemp_Warmest_Quarter), length = 40)
f <- dnorm(x, mean = mean(data$MeanTemp_Warmest_Quarter), sd = sd(data$MeanTemp_Warmest_Quarter))
lines(x, f, col = "red", lwd = 2)

hist(data$Annual_Precipitation, prob = TRUE, main = NA, xlab="Annual_Precipitation",col = "lightblue")
x <- seq(min(data$Annual_Precipitation), max(data$Annual_Precipitation), length = 40)
f <- dnorm(x, mean = mean(data$Annual_Precipitation), sd = sd(data$Annual_Precipitation))
lines(x, f, col = "red", lwd = 2)

hist(data$Precipitation_Wettest_Month, prob = TRUE, main = NA, xlab="Precipitation_Wettest_Month",col = "lightblue")
x <- seq(min(data$Precipitation_Wettest_Month), max(data$Precipitation_Wettest_Month), length = 40)
f <- dnorm(x, mean = mean(data$Precipitation_Wettest_Month), sd = sd(data$Precipitation_Wettest_Month))
lines(x, f, col = "red", lwd = 2)

hist(data$Precipitation_Driest_Month_log, prob = TRUE, main = NA, xlab="log (Precipitation Driest Month)",col = "lightblue")
x <- seq(min(data$Precipitation_Driest_Month_log), max(data$Precipitation_Driest_Month_log), length = 40)
f <- dnorm(x, mean = mean(data$Precipitation_Driest_Month_log), sd = sd(data$Precipitation_Driest_Month_log))
lines(x, f, col = "red", lwd = 2)

hist(data$Precipitation_Seasonality, prob = TRUE, main = NA, xlab="Precipitation_Seasonality",col = "lightblue")
x <- seq(min(data$Precipitation_Seasonality), max(data$Precipitation_Seasonality), length = 40)
f <- dnorm(x, mean = mean(data$Precipitation_Seasonality), sd = sd(data$Precipitation_Seasonality))
lines(x, f, col = "red", lwd = 2)

hist(data$Precipitation_Warmest_Quarter, prob = TRUE, main = NA, xlab="Precipitation_Warmest_Quarter",col = "lightblue")
x <- seq(min(data$Precipitation_Warmest_Quarter), max(data$Precipitation_Warmest_Quarter), length = 40)
f <- dnorm(x, mean = mean(data$Precipitation_Warmest_Quarter), sd = sd(data$Precipitation_Warmest_Quarter))
lines(x, f, col = "red", lwd = 2)

hist(data$Precipitation_Wettest_Quarter, prob = TRUE, main = NA, xlab="Precipitation_Wettest_Quarter",col = "lightblue")
x <- seq(min(data$Precipitation_Wettest_Quarter), max(data$Precipitation_Wettest_Quarter), length = 40)
f <- dnorm(x, mean = mean(data$Precipitation_Wettest_Quarter), sd = sd(data$Precipitation_Wettest_Quarter))
lines(x, f, col = "red", lwd = 2)

hist(data$Precipitation_Driest_Quarter_log, prob = TRUE, main = NA, xlab="log (Precipitation_Driest_Quarter)",col = "lightblue")
x <- seq(min(data$Precipitation_Driest_Quarter_log), max(data$Precipitation_Driest_Quarter_log), length = 40)
f <- dnorm(x, mean = mean(data$Precipitation_Driest_Quarter_log), sd = sd(data$Precipitation_Driest_Quarter_log))
lines(x, f, col = "red", lwd = 2)

hist(data$soil_m_log, prob = TRUE, main = NA, xlab="log (Soil moisture)",col = "khaki1")
x <- seq(min(data$soil_m_log), max(data$soil_m_log), length = 40)
f <- dnorm(x, mean = mean(data$soil_m_log), sd = sd(data$soil_m_log))
lines(x, f, col = "red", lwd = 2)

dev.off()
```
### Getting LD prunned snpfile

```
#!/usr/bin/env Rscript

#PBS -l nodes=1:ppn=4
#PBS -l walltime=12:00:00
#PBS -l pmem=48gb
#PBS -M azd6024@psu.edu
#PBS -A open
#PBS -j oe
#PBS -m abe

setwd("/gpfs/group/jrl35/default/aayudh/gwas_reseq")

library(gdsfmt)
library(SNPRelate)

genofile <- snpgdsOpen('Sorghum_1757g_2ndtry.gds')
sample.id <- read.gdsn(index.gdsn(genofile, 'sample.id'))
chromosome.id <- read.gdsn(index.gdsn(genofile, 'snp.chromosome'))
snp.position <- read.gdsn(index.gdsn(genofile, 'snp.position'))

set.seed(1000)
snpset = snpgdsLDpruning(genofile,maf=0.05,missing.rate=0.5,ld.threshold=0.4, slide.max.bp=50000)
snp.id = unlist(snpset)

setwd("~/work/rda_drought")
traits <- data.frame(read.csv('1.Drought_variables.csv', header=TRUE))
head(traits)

#mySNPs <- snpgdsSelectSNP(genofile, traits$LIB, maf = 0.05, missing.rate = 0.5)
#length(mySNPs)

mySNPmatrix <- snpgdsGetGeno(genofile, sample.id = traits$LIB, snp.id = snp.id, with.id = TRUE)
dim(mySNPmatrix$genotype)

#these lines create the bed file for gemma
rs <- data.frame(chr = chromosome.id[snp.id], positions = snp.position[snp.id])
SNPs_bim <- data.frame(paste0('rs_', rs$chr, '_', rs$positions), t(mySNPmatrix$genotype))
dim(SNPs_bim)
length(SNPs_bim)
SNPs_bim[1:10,1:10]

data_2 <- SNPs_bim[,2:ncol(SNPs_bim)]
colnames(data_2) <- mySNPmatrix$sample.id
data_3 <- cbind(SNPs_bim[,1],data_2)
dim(data_3)
data_3[1:10,1:10]

library(data.table)
data_snp <- transpose(data_3)
colnames(data_snp) <- rownames(data_3)
rownames(data_snp) <- colnames(data_3)
colnames(data_snp) <- data_snp[1,]
dim(data_snp)
df <- data_snp[2:nrow(data_snp),1:ncol(data_snp)]
dim(df)
df[1:10,1:10]

df <- df[,1:ncol(df)]

write.csv(df,"snpLDprunned.csv")
```
### Check correlation matrix and create drought.rda

#### RDA resources: https://popgen.nescent.org/2018-03-27_RDA_GEA.html 

```
#!/usr/bin/env Rscript

#PBS -l nodes=1:ppn=4
#PBS -l walltime=12:00:00
#PBS -l pmem=192gb
#PBS -M azd6024@psu.edu
#PBS -A open
#PBS -j oe
#PBS -m abe

setwd("~/work/rda_drought")
library(psych)    # Used to investigate correlations among predictors
library(permute)
library(lattice)
library(vegan)

env <- read.csv("1.Drought_variables.csv")
pred <- env

pred <- pred[,6:21]
colnames(pred) <- c("WSvg","WSrepro","AMT","MTWM","MTWQ","MTDQ","MTWaQ","AP","PWM","PDM","PS","PWQ","PDQ","PWaQ","SM","AI")

#based on the cor matrix we are removing MTWaQ, PWQ, PDQ, AP
pred = within(pred, rm(MTWaQ, PWQ, PDQ, AP))
tiff("cor.matrix.tiff", width = 13, height = 7, units = 'in', res = 300)
pairs.panels(pred, scale=T)
dev.off()

setwd("/gpfs/group/jrl35/default/aayudh/rda_drought")
df <- read.csv("snpLDprunned.csv", row.names = 1)
gen.imp <- apply(df, 2, function(x) replace(x, is.na(x), as.numeric(names(which.max(table(x))))))

drought.rda <- rda(gen.imp ~ ., data=pred, scale=T)

#setwd("/gpfs/group/jrl35/default/aayudh/rda_drought")
setwd("~/work/rda_drought")
saveRDS(drought.rda, file="drought_rda.RData")

ddf <-  drought.rda$CCA$v
RDA_values<- as.data.frame(ddf)
write.csv(RDA_values,"RDA_values.csv")
```
### Get SNP candidates for drought.rda

```
#!/usr/bin/env Rscript

#PBS -l nodes=1:ppn=4
#PBS -l walltime=12:00:00
#PBS -l pmem=192gb
#PBS -M azd6024@psu.edu
#PBS -A open
#PBS -j oe
#PBS -m abe

setwd("~/work/rda_drought")
library(psych)    # Used to investigate correlations among predictors
library(permute)
library(lattice)
library(vegan)    # Used to run RDA

#setwd("~/Library/CloudStorage/OneDrive-UniversityofVermont/PENN STATE/Drought_multivar_envGWAS/RDA/test")
env <- read.csv("1.Drought_variables.csv")
pred <- env
pred <- pred[,6:21]
colnames(pred) <- c("WSvg","WSrepro","AMT","MTWM","MTWQ","MTDQ","MTWaQ","AP","PWM","PDM","PS","PWQ","PDQ","PWaQ","SM","AI")

#based on the cor matrix we are removing MTWaQ, PWQ, PDQ, AP
pred = within(pred, rm(MTWaQ, PWQ, PDQ, AP))

setwd("/gpfs/group/jrl35/default/aayudh/rda_drought")
df <- read.csv("snpLDprunned.csv", row.names = 1)
gen.imp <- apply(df, 2, function(x) replace(x, is.na(x), as.numeric(names(which.max(table(x))))))

drought.rda <- rda(gen.imp ~ ., data=pred, scale=T)

load.rda <- scores(drought.rda, choices=c(1:3), display="species")

outliers <- function(x,z){
  lims <- mean(x) + c(-1, 1) * z * sd(x)     # find loadings +/-z sd from mean loading     
  x[x < lims[1] | x > lims[2]]               # locus names in these tails
}

cand1 <- outliers(load.rda[,1],3) # 38
cand2 <- outliers(load.rda[,2],3) # 69
cand3 <- outliers(load.rda[,3],3) # 34

ncand <- length(cand1) + length(cand2) + length(cand3)
ncand

cand1 <- cbind.data.frame(rep(1,times=length(cand1)), names(cand1), unname(cand1))
cand2 <- cbind.data.frame(rep(2,times=length(cand2)), names(cand2), unname(cand2))
cand3 <- cbind.data.frame(rep(3,times=length(cand3)), names(cand3), unname(cand3))

colnames(cand1) <- colnames(cand2) <- colnames(cand3) <- c("axis","snp","loading")

cand <- rbind(cand1, cand2, cand3)
cand$snp <- as.character(cand$snp)

foo <- matrix(nrow=(ncand), ncol=12)  # 12 columns for 16 predictors
colnames(foo) <- c("WSvg","WSrepro","AMT","MTWM","MTWQ","MTDQ","PWM","PDM","PS","PWaQ","SM","AI")

for (i in 1:length(cand$snp)) {
  nam <- cand[i,2]
  snp.gen <- gen.imp[,nam]
  foo[i,] <- apply(pred,2,function(x) cor(x,snp.gen))
  }

cand <- cbind.data.frame(cand,foo)
head(cand)

#Investigate the candidates= if you get any then do this analysis
length(cand$snp[duplicated(cand$snp)])

foo <- cbind(cand$axis, duplicated(cand$snp))
table(foo[foo[,1]==1,2]) # no duplicates on axis 1
table(foo[foo[,1]==2,2]) #  no duplicates on axis 2
table(foo[foo[,1]==3,2]) # 6 duplicates on axis 3
cand <- cand[!duplicated(cand$snp),]

setwd("~/work/rda_drought")
write.csv(cand, "1.cand.csv")
```

### Plotting 
```
setwd("~/Library/CloudStorage/OneDrive-ThePennsylvaniaStateUniversity/RDA_drought")
list.files()
library(psych)    # Used to investigate correlations among predictors
library(permute) 
library(lattice)
library(vegan)  

drought.rda <- readRDS("drought_rda.RData")
env <- read.csv("1.RDA_Drought_variables.csv")
pred <- env
pred <- pred[,7:18]
colnames(pred) <- c("WSvg","WSrepro","AMT","MTWM","MTWQ","MTDQ","PWM","PDM","PS","PWaQ","SM","AI")
#based on the cor matrix we are removing MTWaQ, PWQ, PDQ, AP
#pred = within(pred, rm(MTWaQ, PWQ, PDQ, AP))


setwd("~/Library/CloudStorage/OneDrive-ThePennsylvaniaStateUniversity/RDA_drought")
tiff("cor.matrix.tiff", width = 13, height = 7, units = 'in', res = 300)
pairs.panels(pred, scale=T)
dev.off()

ddf <-  drought.rda$CCA$v
RDA_values<- as.data.frame(ddf)

write.csv(RDA_values,"RDA_values.csv")

RDA_values <- read.csv("RDA_values.csv")

png("screeplot.png", width = 4, height = 7, units = 'in', res = 300)
screeplot(drought.rda)
dev.off()

png("RDA1_2.png", width = 13, height = 7, units = 'in', res = 300)
par(mfrow=c(1,2))
plot(drought.rda, scaling=3)          # default is axes 1 and 2
plot(drought.rda, choices = c(1, 3), scaling=3)  # axes 1 and 3
dev.off()

levels(env$Continent) <- c("America","Asia","East Africa","Southern Africa","West Africa")
Continent <- env$Continent
bg <- c("#ff7f00","#1f78b4","#ffff33","#a6cee3","#33a02c") 

# axes 1 & 2
png("rda1_2_continent.png", width = 13, height = 8, units = 'in', res = 300)
par(mfrow=c(1,2))
plot(drought.rda, type="n", scaling=3)
points(drought.rda, display="species", pch=20, cex=0.7, col="gray32", scaling=3)           # the SNPs
points(drought.rda, display="sites", pch=21, cex=1.3, col="gray32", scaling=3, bg=bg[Continent]) # the wolves
text(drought.rda, scaling=3, display="bp", col="#0868ac", cex=1)                           # the predictors
legend("bottomright", legend=levels(Continent), bty="n", col="gray32", pch=21, cex=1, pt.bg=bg)

# axes 1 & 3
plot(drought.rda, type="n", scaling=3, choices=c(1,3))
points(drought.rda, display="species", pch=20, cex=0.7, col="gray32", scaling=3, choices=c(1,3))
points(drought.rda, display="sites", pch=21, cex=1.3, col="gray32", scaling=3, bg=bg[Continent], choices=c(1,3))
text(drought.rda, scaling=3, display="bp", col="#0868ac", cex=1, choices=c(1,3))
#legend("bottomright", legend=levels(Continent), bty="n", col="gray32", pch=21, cex=1, pt.bg=bg)
dev.off()

load.rda <- scores(drought.rda, choices=c(1:3), display="species")

png("Loadings on RDA.png", width = 13, height = 7, units = 'in', res = 300)
par(mfrow=c(1,4))
screeplot(drought.rda)
hist(load.rda[,1], main="Loadings on RDA1")
hist(load.rda[,2], main="Loadings on RDA2")
hist(load.rda[,3], main="Loadings on RDA3") 
dev.off()

cand <- read.csv("1.cand.csv")
#Investigate the candidates= if you get any then do this analysis
length(cand$snp[duplicated(cand$snp)])  

foo <- cbind(cand$axis, duplicated(cand$snp)) 
table(foo[foo[,1]==1,2]) # no duplicates on axis 1
table(foo[foo[,1]==2,2]) #  no duplicates on axis 2
table(foo[foo[,1]==3,2]) # 6 duplicates on axis 3
cand <- cand[!duplicated(cand$snp),] 

##--------------------
cand <- read.csv("1.cand.csv")
head(cand)

#Next, we’ll see which of the predictors each candidate SNP is most strongly correlated with:
for (i in 1:length(cand$snp)) {
  bar <- cand[i,]
  cand[i,16] <- names(which.max(abs(bar[4:15]))) # gives the variable
  cand[i,17] <- max(abs(bar[4:15]))              # gives the correlation
}

colnames(cand)[16] <- "predictor"
colnames(cand)[17] <- "correlation"

gg <- table(cand$predictor) 
write.table(gg , file = "olstab.csv", sep = ",", quote = FALSE, row.names = F)

gg <- read.csv("olstab.csv")
head(gg)
gg1 <- gg[order(-gg$Freq),]
gg1$percent <- (gg1$Freq/sum(gg1$Freq))*100
gg1$percent <- round(gg1$percent, digits = 2)


library("gridExtra")
tiff("predictor.tiff", width = 3, height = 6, units = 'in', res = 300)
par(family="Times")
grid.table(gg1)
dev.off()

#Plot the SNPs

sel <- cand$snp
env <- cand$predictor
env[env=="WSvg"] <- '#1f78b4'
env[env=="WSrepro"] <- '#a6cee3'
env[env=="AMT"] <- '#6a3d9a'
env[env=="MTWM"] <- '#e31a1c'
env[env=="MTWQ"] <- '#33a02c'
env[env=="MTDQ"] <- '#ffff33'
env[env=="PWM"] <- '#ff7f00'
env[env=="PDM"] <- 'cadetblue'
env[env=="PS"] <- '#fb9a99'
env[env=="PWaQ"] <- 'blue'
env[env=="SM"] <- 'green'
env[env=="AI"] <- 'skyblue'

col.pred <- rownames(drought.rda$CCA$v) # pull the SNP names
col.pred = as.factor(col.pred)

for (i in 1:length(sel)) {           # color code candidate SNPs
  foo <- match(sel[i],col.pred)
  col.pred[foo] <- env[i]
}

col.pred[grep("chr",col.pred)] <- '#f1eef6' # non-candidate SNPs
empty <- col.pred
empty[grep("#f1eef6",empty)] <- rgb(0,1,0, alpha=0) # transparent
empty.outline <- ifelse(empty=="#00FF0000","#00FF0000","gray32")
bg <- c('#1f78b4','#a6cee3','#6a3d9a','#e31a1c','#33a02c','#ffff33',
        '#ff7f00','cadetblue','#fb9a99','blue','green','skyblue')


png("RDA_SNPbased.png", width = 13, height = 7, units = 'in', res = 300)
par(mfrow=c(1,2))
# axes 1 & 2 
plot(drought.rda, type="n", scaling=3, xlim=c(-1,1), ylim=c(-1,1))
points(drought.rda, display="species", pch=21, cex=1, col="gray32", bg=col.pred, scaling=3)
points(drought.rda, display="species", pch=21, cex=1, col=empty.outline, bg=empty, scaling=3)
text(drought.rda, scaling=3, display="bp", col="#0868ac", cex=1)
legend("bottomright", legend=c("WSvg","WSrepro","AMT","MTWM","MTWQ","MTDQ","PWM","PDM","PS","PWaQ","SM","AI"), bty="n", col="gray32", pch=21, cex=1, pt.bg=bg)

# axes 1 & 3
plot(drought.rda, type="n", scaling=3, xlim=c(-1,1), ylim=c(-1,1), choices=c(1,3))
points(drought.rda, display="species", pch=21, cex=1, col="gray32", bg=col.pred, scaling=3, choices=c(1,3))
points(drought.rda, display="species", pch=21, cex=1, col=empty.outline, bg=empty, scaling=3, choices=c(1,3))
text(drought.rda, scaling=3, display="bp", col="#0868ac", cex=1, choices=c(1,3))
#legend("bottomright", legend=c("WSvg","WSrepro","AMT","MTWM","MTWQ","MTDQ","PWM","PDM","PS","PWaQ","SM","AI"), bty="n", col="gray32", pch=21, cex=1, pt.bg=bg)
dev.off()
```

-----
<div id='id-section16'/>

## Chapter 16: Enrichment Analysis

### Extracting variants from snpeff SNP version

```
#!/bin/bash

#PBS -l nodes=1:ppn=8
#PBS -l walltime=12:00:00
#PBS -l pmem=24gb
#PBS -M azd6024@psu.edu
#PBS -A open
#PBS -j oe
#PBS -m abe

WORKINGDIR=/gpfs/group/jrl35/default/aayudh/snpeff/snp
cd $WORKINGDIR

grep missense_variant Sorghum_1757g_AllChr.polymorphic.snp.noRepeats.5pctMasked.ann.vcf > non_synonymous_missense_variant.txt

grep initiator_codon_variant Sorghum_1757g_AllChr.polymorphic.snp.noRepeats.5pctMasked.ann.vcf > non_synonymous_initiator_codon_variant.txt

grep stop_retained_variant Sorghum_1757g_AllChr.polymorphic.snp.noRepeats.5pctMasked.ann.vcf > non_synonymous_stop_retained_variant.txt

grep synonymous_variant Sorghum_1757g_AllChr.polymorphic.snp.noRepeats.5pctMasked.ann.vcf > synonymous_variant.txt

grep start_retained_variant Sorghum_1757g_AllChr.polymorphic.snp.noRepeats.5pctMasked.ann.vcf > synonymous_start_retained.txt

grep stop_retained_variant Sorghum_1757g_AllChr.polymorphic.snp.noRepeats.5pctMasked.ann.vcf > synonymous_stop_retained_variant.txt
```
### Making syno and non-syno datset

```
#!/usr/bin/env Rscript

#PBS -l nodes=1:ppn=8
#PBS -l walltime=12:00:00
#PBS -l pmem=24gb
#PBS -M azd6024@psu.edu
#PBS -A open
#PBS -j oe
#PBS -m abe

setwd("/gpfs/group/jrl35/default/aayudh/snpeff/snp/compare_synoVSnon_syn")
data <- read.delim("synonymous_variant.txt", header=FALSE)
dim(data)
data1 <- data[,1:8]
data1 = within(data1, rm(V3,V6,V7))
head(data1)
names(data1)[1]<-paste("chr")
names(data1)[2]<-paste("ps")
names(data1)[3]<-paste("REF")
names(data1)[4]<-paste("ALT")
names(data1)[5]<-paste("type")
head(data1)

rs_split <- data.frame(do.call("rbind", strsplit(as.character(data1$type), ";", fixed = TRUE)))
head(rs_split)

data1[, "type"] <- rs_split$X11
head(data1)
rs_split1 <- data.frame(do.call("rbind", strsplit(as.character(data1$type), "|", fixed = TRUE)))
head(rs_split1)
data1[, "type"] <- rs_split1$X2

head(data1)
write.csv(data1, "synonymous_variant.txt_clean.csv")
```

```
#!/usr/bin/env Rscript

#PBS -l nodes=1:ppn=8
#PBS -l walltime=12:00:00
#PBS -l pmem=24gb
#PBS -M azd6024@psu.edu
#PBS -A open
#PBS -j oe
#PBS -m abe

setwd("/gpfs/group/jrl35/default/aayudh/snpeff/snp/compare_synoVSnon_syn")
data <- read.delim("non_synonymous_missense_variant.txt", header=FALSE)
dim(data)
data1 <- data[,1:8]
data1 = within(data1, rm(V3,V6,V7))
head(data1)
names(data1)[1]<-paste("chr")
names(data1)[2]<-paste("ps")
names(data1)[3]<-paste("REF")
names(data1)[4]<-paste("ALT")
names(data1)[5]<-paste("type")
head(data1)

rs_split <- data.frame(do.call("rbind", strsplit(as.character(data1$type), ";", fixed = TRUE)))
head(rs_split)

data1[, "type"] <- rs_split$X11
head(data1)
rs_split1 <- data.frame(do.call("rbind", strsplit(as.character(data1$type), "|", fixed = TRUE)))
head(rs_split1)
data1[, "type"] <- rs_split1$X2

head(data1)
write.csv(data1, "non_synonymous_missense_variant_clean.csv")
```

### Merging snpeff calls with envGWAS hits

```
setwd("~/Library/CloudStorage/OneDrive-UniversityofVermont/PENN STATE/eGWAS_revised list/manhattan plot")
list.files()
Chr10corr <- read.csv("Reseq_preds_all_chr10_corrected.csv", header=T)
head(Chr10corr)

setwd("~/Library/CloudStorage/OneDrive-UniversityofVermont/PENN STATE/SNPEFF/syn and non-syn")
list.files()

syn <- read.csv("synonymous_variant.txt_clean.csv")
head(syn)
syn[, "chr"] = substring(syn$chr, 4)
rs <- data.frame(paste0('rs_', syn$chr, '_', syn$ps)) 
syn[, "X"] = rs
names(syn)[1]<-paste("rs")
head(syn)
syn <- syn[,1:6]
syn_merged_HS <- merge(syn,Chr10corr, by = "rs")
head(syn_merged_HS)

setwd("~/Library/CloudStorage/OneDrive-UniversityofVermont/PENN STATE/SNPEFF/syn and non-syn")
#write.csv(syn_merged_HS,"syn_merged_HS.csv")

non_syn <- read.csv("non_synonymous_missense_variant_clean.csv")
head(non_syn)

head(non_syn)
non_syn[, "chr"] = substring(non_syn$chr, 4)
rs <- data.frame(paste0('rs_', non_syn$chr, '_', non_syn$ps)) 
non_syn[, "X"] = rs
names(non_syn)[1]<-paste("rs")
head(non_syn)
non_syn_merged_HS <- merge(non_syn,Chr10corr, by = "rs")
head(non_syn_merged_HS)
#write.csv(non_syn_merged_HS,"non_syn_merged_HS.csv")

setwd("~/Library/CloudStorage/OneDrive-UniversityofVermont/PENN STATE/SNPEFF/syn and non-syn")
syn <- read.csv("1.syn_merged_HS.csv")
syn$log_pwald = -log10(syn$p_wald)
head(syn)

syn1 <- subset(syn, syn$type == "synonymous_variant")

non_syn <- read.csv("1. non_syn_merged_HS.csv")
non_syn$log_pwald = -log10(non_syn$p_wald)
head(non_syn)

non_syn1 <- subset(non_syn, non_syn$type == "missense_variant")

data <- rbind(syn1,non_syn1)
head(data)

## https://bookdown.org/ybrandvain/Applied-Biostats/perm1.html 

SoilSciGuylabs <- c("Non-synonymous", "Synonymous")

setwd("~/Library/CloudStorage/OneDrive-UniversityofVermont/PENN STATE/SNPEFF/syn and non-syn")
tiff("synVSnon-syn.tiff", width = 5, height = 5, units = 'in', res = 300)
library(ggplot2)
library(ggforce)
ggplot(data, aes(x = type, y =  log_pwald, color = type, shape = type))+
  geom_sina(show.legend = FALSE)+ # if you don't have ggforce you can use geom_jitter
  stat_summary(color = "black",show.legend = FALSE)+
  labs(x=NULL,y="-log10 (pvalue)")+
scale_x_discrete(labels= SoilSciGuylabs)+
scale_color_manual(values=c("#999999", "#E69F00"))+
  theme_classic() +
  theme(legend.position="none",text=element_text(size=16, colour = "black", family="Times New Roman"),
        axis.line = element_line(size=0.5, colour = "black"),
        axis.text.x=element_text(colour="black", size = 16),
        axis.text.y=element_text(colour="black", size = 16))
dev.off()
```

### Chromosome rearrangement

```
setwd("~/Library/CloudStorage/OneDrive-UniversityofVermont/PENN STATE/SNPEFF/syn and non-syn")
list.files()

syn <- read.csv("1.syn_merged_HS.csv")
head(syn)

syn1 <- subset(syn, syn$type == "synonymous_variant")
dim(syn1)
tail(syn1)

non_syn <- read.csv("1. non_syn_merged_HS.csv")
non_syn1 <- subset(non_syn, non_syn$type == "missense_variant")

syn_nonsyn <- rbind(syn1,non_syn1)
dim(syn_nonsyn)

##arrange chromosome wise
chr1 <- subset(syn_nonsyn, syn_nonsyn$chr == "1")
chr1 <- unique(setDT(chr1)[order(ps, -ps)], by = "ps")
head(chr1)

chr2 <- subset(syn_nonsyn, syn_nonsyn$chr == "2")
chr2 <- unique(setDT(chr2)[order(ps, -ps)], by = "ps")
head(chr2)

chr3 <- subset(syn_nonsyn, syn_nonsyn$chr == "3")
chr3 <- unique(setDT(chr3)[order(ps, -ps)], by = "ps")
head(chr3)

chr4 <- subset(syn_nonsyn, syn_nonsyn$chr == "4")
chr4 <- unique(setDT(chr4)[order(ps, -ps)], by = "ps")
head(chr4)

chr5 <- subset(syn_nonsyn, syn_nonsyn$chr == "5")
chr5 <- unique(setDT(chr5)[order(ps, -ps)], by = "ps")
head(chr5)

chr6 <- subset(syn_nonsyn, syn_nonsyn$chr == "6")
chr6 <- unique(setDT(chr6)[order(ps, -ps)], by = "ps")
head(chr6)

chr7 <- subset(syn_nonsyn, syn_nonsyn$chr == "7")
chr7 <- unique(setDT(chr7)[order(ps, -ps)], by = "ps")
head(chr7)

chr8 <- subset(syn_nonsyn, syn_nonsyn$chr == "8")
chr8 <- unique(setDT(chr8)[order(ps, -ps)], by = "ps")
head(chr8)

chr9 <- subset(syn_nonsyn, syn_nonsyn$chr == "9")
chr9 <- unique(setDT(chr9)[order(ps, -ps)], by = "ps")
head(chr9)

chr10 <- subset(syn_nonsyn, syn_nonsyn$chr == "10")
chr10 <- unique(setDT(chr10)[order(ps, -ps)], by = "ps")
head(chr10)

chr1_10 <- rbind(chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10)
dim(chr1_10)
tail(chr1_10)
head(chr1_10)

setwd("~/Library/CloudStorage/OneDrive-UniversityofVermont/PENN STATE/SNPEFF/syn and non-syn")
write.csv(chr1_10,"1.Syn_Nonsyn_chr_rearranged.csv",row.names = FALSE)
```

### Enrichment analysis with MAF corrections

```
setwd("~/Library/CloudStorage/OneDrive-ThePennsylvaniaStateUniversity/snp enrichment")
list.files()

data1 <- read.csv("1.Syn_Nonsyn_chr_rearranged.csv")
data <- data1[,1:9]
head(data)

data$logp <- -log10(data$p_wald)

Actual_fold_enrichment_synTOnonsyn <- quantile(subset(data, type == "synonymous_variant")$logp, probs = 0.95)/quantile(subset(data, type == "missense_variant")$logp, probs = 0.95)
Actual_fold_enrichment_synTOnonsyn

dddf <- NULL

for(i in 1:100) { 
  random_position <- data[sample(nrow(data), 1), ]
  f <- rbind(data[which(data$rs == random_position$rs, arr.ind=TRUE):nrow(data),],data[1:(which(data$rs == random_position$rs, arr.ind=TRUE)-1),])
  data[,"type"] <- f$type 
  fold_enrichment <- quantile(subset(data, type == "synonymous_variant")$logp, probs = 0.95)/quantile(subset(data, type == "missense_variant")$logp, probs = 0.95)
  dddf <- rbind(dddf,fold_enrichment)}

synBYnonsyn  <- as.data.frame(dddf)
#write.csv(synBYnonsyn, "synBYnonsyn.csv", row.names = FALSE)

synBYnonsyn$ratio <- paste("synBYnonsyn")

#fold_enrich <- rbind(nonsynBYsyn,synBYnonsyn)

#write.csv(fold_enrich, "1.fold_enrichment.csv", row.names = FALSE)

###Now miniAF correction

setwd("~/Library/CloudStorage/OneDrive-ThePennsylvaniaStateUniversity/snp enrichment")
list.files()

data1 <- read.csv("1.Syn_Nonsyn_chr_rearranged.csv")
data <- data1[,1:9]
head(data)

##Minor allele frequency correction, whatever aboive 0.5 af just deduct 0.5.
data_below_.5 <- subset(data, af < .5)
data_above_.5 <- subset(data, af > .5)

data_above_.5$af <- data_above_.5$af -0.5 

data_cortd <- rbind(data_below_.5,data_above_.5)

data_cortd <- data_cortd[with(data_cortd,order(chr)),]

##minaf 0-0.1
minaf_0.1 <- subset(data_cortd, af < .1) 
minaf_0.1$logp <- -log10(minaf_0.1$p_wald)

Actual_fold_enrichment_minaf_0.1 <- quantile(subset(minaf_0.1, type == "synonymous_variant")$logp, probs = 0.95)/quantile(subset(minaf_0.1, type == "missense_variant")$logp, probs = 0.95)
Actual_fold_enrichment_minaf_0.1

dddf <- NULL

for(i in 1:100) { 
  random_position <- minaf_0.1[sample(nrow(minaf_0.1), 1), ]
  f <- rbind(minaf_0.1[which(minaf_0.1$rs == random_position$rs, arr.ind=TRUE):nrow(minaf_0.1),],minaf_0.1[1:(which(minaf_0.1$rs == random_position$rs, arr.ind=TRUE)-1),])
  minaf_0.1[,"type"] <- f$type 
  fold_enrichment <- quantile(subset(minaf_0.1, type == "synonymous_variant")$logp, probs = 0.95)/quantile(subset(minaf_0.1, type == "missense_variant")$logp, probs = 0.95)
  dddf <- rbind(dddf,fold_enrichment)}

synBYnonsyn_minaf_0.1  <- as.data.frame(dddf)
setwd("~/Library/CloudStorage/OneDrive-ThePennsylvaniaStateUniversity/snp enrichment/minaf_correction")
write.csv(synBYnonsyn_minaf_0.1, "synBYnonsyn_minaf_0.1.csv", row.names = FALSE)

synBYnonsyn_minaf_0.1$ratio <- paste("synBYnonsyn_minaf_0.1")

##minaf 0.1-0.2
minaf_0.2 <- subset(data_cortd, af < 0.2 & af > 0.1) 
minaf_0.2$logp <- -log10(minaf_0.2$p_wald)

Actual_fold_enrichment_minaf_0.2 <- quantile(subset(minaf_0.2, type == "synonymous_variant")$logp, probs = 0.95)/quantile(subset(minaf_0.2, type == "missense_variant")$logp, probs = 0.95)
Actual_fold_enrichment_minaf_0.2

dddf <- NULL

for(i in 1:100) { 
  random_position <- minaf_0.2[sample(nrow(minaf_0.2), 1), ]
  f <- rbind(minaf_0.2[which(minaf_0.2$rs == random_position$rs, arr.ind=TRUE):nrow(minaf_0.2),],minaf_0.2[1:(which(minaf_0.2$rs == random_position$rs, arr.ind=TRUE)-1),])
  minaf_0.2[,"type"] <- f$type 
  fold_enrichment <- quantile(subset(minaf_0.2, type == "synonymous_variant")$logp, probs = 0.95)/quantile(subset(minaf_0.2, type == "missense_variant")$logp, probs = 0.95)
  dddf <- rbind(dddf,fold_enrichment)}

synBYnonsyn_minaf_0.2  <- as.data.frame(dddf)
setwd("~/Library/CloudStorage/OneDrive-ThePennsylvaniaStateUniversity/snp enrichment/minaf_correction")
write.csv(synBYnonsyn_minaf_0.2, "synBYnonsyn_minaf_0.2.csv", row.names = FALSE)

synBYnonsyn_minaf_0.2$ratio <- paste("synBYnonsyn_minaf_0.2")


##minaf 0.2-0.3
minaf_0.3 <- subset(data_cortd, af < 0.3 & af > 0.2) 
minaf_0.3$logp <- -log10(minaf_0.3$p_wald)

Actual_fold_enrichment_minaf_0.3 <- quantile(subset(minaf_0.3, type == "synonymous_variant")$logp, probs = 0.95)/quantile(subset(minaf_0.3, type == "missense_variant")$logp, probs = 0.95)
Actual_fold_enrichment_minaf_0.3

dddf <- NULL

for(i in 1:100) { 
  random_position <- minaf_0.3[sample(nrow(minaf_0.3), 1), ]
  f <- rbind(minaf_0.3[which(minaf_0.3$rs == random_position$rs, arr.ind=TRUE):nrow(minaf_0.3),],minaf_0.3[1:(which(minaf_0.3$rs == random_position$rs, arr.ind=TRUE)-1),])
  minaf_0.3[,"type"] <- f$type 
  fold_enrichment <- quantile(subset(minaf_0.3, type == "synonymous_variant")$logp, probs = 0.95)/quantile(subset(minaf_0.3, type == "missense_variant")$logp, probs = 0.95)
  dddf <- rbind(dddf,fold_enrichment)}

synBYnonsyn_minaf_0.3  <- as.data.frame(dddf)
setwd("~/Library/CloudStorage/OneDrive-ThePennsylvaniaStateUniversity/snp enrichment/minaf_correction")
write.csv(synBYnonsyn_minaf_0.3, "synBYnonsyn_minaf_0.3.csv", row.names = FALSE)

synBYnonsyn_minaf_0.3$ratio <- paste("synBYnonsyn_minaf_0.3")



##minaf 0.3-0.4
minaf_0.4 <- subset(data_cortd, af < 0.4 & af > 0.3) 
minaf_0.4$logp <- -log10(minaf_0.4$p_wald)

Actual_fold_enrichment_minaf_0.4 <- quantile(subset(minaf_0.4, type == "synonymous_variant")$logp, probs = 0.95)/quantile(subset(minaf_0.4, type == "missense_variant")$logp, probs = 0.95)
Actual_fold_enrichment_minaf_0.4

dddf <- NULL

for(i in 1:100) { 
  random_position <- minaf_0.4[sample(nrow(minaf_0.4), 1), ]
  f <- rbind(minaf_0.4[which(minaf_0.4$rs == random_position$rs, arr.ind=TRUE):nrow(minaf_0.4),],minaf_0.4[1:(which(minaf_0.4$rs == random_position$rs, arr.ind=TRUE)-1),])
  minaf_0.4[,"type"] <- f$type 
  fold_enrichment <- quantile(subset(minaf_0.4, type == "synonymous_variant")$logp, probs = 0.95)/quantile(subset(minaf_0.4, type == "missense_variant")$logp, probs = 0.95)
  dddf <- rbind(dddf,fold_enrichment)}

synBYnonsyn_minaf_0.4  <- as.data.frame(dddf)
setwd("~/Library/CloudStorage/OneDrive-ThePennsylvaniaStateUniversity/snp enrichment/minaf_correction")
write.csv(synBYnonsyn_minaf_0.4, "synBYnonsyn_minaf_0.4.csv", row.names = FALSE)

synBYnonsyn_minaf_0.4$ratio <- paste("synBYnonsyn_minaf_0.4")


##minaf 0.4-0.5
minaf_0.5 <- subset(data_cortd, af < 0.5 & af > 0.4) 
minaf_0.5$logp <- -log10(minaf_0.5$p_wald)

Actual_fold_enrichment_minaf_0.5 <- quantile(subset(minaf_0.5, type == "synonymous_variant")$logp, probs = 0.95)/quantile(subset(minaf_0.5, type == "missense_variant")$logp, probs = 0.95)
Actual_fold_enrichment_minaf_0.5

dddf <- NULL

for(i in 1:100) { 
  random_position <- minaf_0.5[sample(nrow(minaf_0.5), 1), ]
  f <- rbind(minaf_0.5[which(minaf_0.5$rs == random_position$rs, arr.ind=TRUE):nrow(minaf_0.5),],minaf_0.5[1:(which(minaf_0.5$rs == random_position$rs, arr.ind=TRUE)-1),])
  minaf_0.5[,"type"] <- f$type 
  fold_enrichment <- quantile(subset(minaf_0.5, type == "synonymous_variant")$logp, probs = 0.95)/quantile(subset(minaf_0.5, type == "missense_variant")$logp, probs = 0.95)
  dddf <- rbind(dddf,fold_enrichment)}

synBYnonsyn_minaf_0.5  <- as.data.frame(dddf)
setwd("~/Library/CloudStorage/OneDrive-ThePennsylvaniaStateUniversity/snp enrichment/minaf_correction")
write.csv(synBYnonsyn_minaf_0.5, "synBYnonsyn_minaf_0.5.csv", row.names = FALSE)

synBYnonsyn_minaf_0.5$ratio <- paste("synBYnonsyn_minaf_0.5")


###bind all

fold_enrich <- rbind(synBYnonsyn,synBYnonsyn_minaf_0.1,synBYnonsyn_minaf_0.2,synBYnonsyn_minaf_0.3,synBYnonsyn_minaf_0.4,synBYnonsyn_minaf_0.5)

setwd("~/Library/CloudStorage/OneDrive-ThePennsylvaniaStateUniversity/snp enrichment/minaf_correction")
write.csv(fold_enrich, "1.fold_enrichment_MAF_correction_with_main.csv", row.names = FALSE)

df <- read.csv("1.fold_enrichment_MAF_correction_with_main.csv") 
head(df)
##plotiing
library(ggplot2)
library(tidyverse)
library(ggtext)
library(normentR)
library(tidyverse)
library(gridExtra)
library(grid)
library(gridtext)
library(palmerpenguins)

setwd("~/Library/CloudStorage/OneDrive-ThePennsylvaniaStateUniversity/snp enrichment/minaf_correction")
tiff("Fold_enrich_HS_100perm_MAFcor.tiff", width = 5, height = 6, units = 'in', res = 300)
ggplot(df, aes(x=X95., y=ratio, color=ratio))+
  geom_point(size=1)+
  labs(x="Fold enrichment",y=NULL)+
  scale_x_continuous(breaks=seq(0.964,1.040,0.025))+
  scale_color_manual(values = c("#B4AF46", "#B4464B","#B47846","#7846B4","#82B446","deeppink"))+
  geom_vline(xintercept = 1, linetype="dashed", color = "slategrey", size=1)+
  annotate("pointrange", x = Actual_fold_enrichment_synTOnonsyn , y = 1, ymin = 1.010, ymax = 1.008,colour = "#4682B4",size = 2)+
  annotate("pointrange", x = Actual_fold_enrichment_minaf_0.1 , y = 2, ymin = 1.010, ymax = 1.008,colour = "#4682B4",size = 2)+
  annotate("pointrange", x = Actual_fold_enrichment_minaf_0.2 , y = 3, ymin = 2.010, ymax = 2.008,colour = "#4682B4",size = 2)+
  annotate("pointrange", x = Actual_fold_enrichment_minaf_0.3 , y = 4, ymin = 1.010, ymax = 1.008,colour = "#4682B4",size = 2)+
  annotate("pointrange", x = Actual_fold_enrichment_minaf_0.4 , y = 5, ymin = 2.010, ymax = 2.008,colour = "#4682B4",size = 2)+
  annotate("pointrange", x = Actual_fold_enrichment_minaf_0.5 , y = 6, ymin = 2.010, ymax = 2.008,colour = "#4682B4",size = 2)+
  theme_classic()+
  theme(legend.position="none",text=element_text(size=16, colour = "black", family="Times New Roman"),
        axis.line = element_line(size=0.5, colour = "black"),
        axis.text.x=element_text(colour="black", size = 16),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=1))
dev.off()
```

-----
<div id='id-section17'/>

## Chapter 17: Genic VS Intergenic enrichment

### Extract all the genic region variant 

What are genic variants? - Synonymous, non-synonumous, initiator codon, UTRs, non-coding intron and splice variants. Rest everything is intergenic.

```
#!/bin/bash

#PBS -l nodes=1:ppn=8
#PBS -l walltime=12:00:00
#PBS -l pmem=24gb
#PBS -M azd6024@psu.edu
#PBS -A open
#PBS -j oe
#PBS -m abe

WORKINGDIR=/gpfs/group/jrl35/default/aayudh/snpeff/snp
cd $WORKINGDIR

grep missense_variant Sorghum_1757g_AllChr.polymorphic.snp.noRepeats.5pctMasked.ann.vcf > non_synonymous_missense_variant.txt

grep initiator_codon_variant Sorghum_1757g_AllChr.polymorphic.snp.noRepeats.5pctMasked.ann.vcf > non_synonymous_initiator_codon_variant.txt

grep stop_retained_variant Sorghum_1757g_AllChr.polymorphic.snp.noRepeats.5pctMasked.ann.vcf > non_synonymous_stop_retained_variant.txt

grep synonymous_variant Sorghum_1757g_AllChr.polymorphic.snp.noRepeats.5pctMasked.ann.vcf > synonymous_variant.txt

grep start_retained_variant Sorghum_1757g_AllChr.polymorphic.snp.noRepeats.5pctMasked.ann.vcf > synonymous_start_retained.txt

grep stop_retained_variant Sorghum_1757g_AllChr.polymorphic.snp.noRepeats.5pctMasked.ann.vcf > synonymous_stop_retained_variant.txt

grep 3_prime_UTR_variant Sorghum_1757g_AllChr.polymorphic.snp.noRepeats.5pctMasked.ann.vcf > 3_prime_UTR_variant.txt

grep 5_prime_UTR_premature_start_codon_gain_variant Sorghum_1757g_AllChr.polymorphic.snp.noRepeats.5pctMasked.ann.vcf > 5_prime_UTR_premature_start_codon_gain_variant.txt

grep 5_prime_UTR_variant Sorghum_1757g_AllChr.polymorphic.snp.noRepeats.5pctMasked.ann.vcf > 5_prime_UTR_variant.txt

grep non_coding_transcript_exon_variant Sorghum_1757g_AllChr.polymorphic.snp.noRepeats.5pctMasked.ann.vcf > non_coding_transcript_exon_variant.txt

grep splice_acceptor_variant Sorghum_1757g_AllChr.polymorphic.snp.noRepeats.5pctMasked.ann.vcf > splice_acceptor_variant.txt

grep splice_donor_variant Sorghum_1757g_AllChr.polymorphic.snp.noRepeats.5pctMasked.ann.vcf > splice_donor_variant.txt

grep splice_region_variant Sorghum_1757g_AllChr.polymorphic.snp.noRepeats.5pctMasked.ann.vcf > splice_region_variant.txt
```

### Clean each variant file

```
#!/usr/bin/env Rscript

#PBS -l nodes=1:ppn=8
#PBS -l walltime=12:00:00
#PBS -l pmem=24gb
#PBS -M azd6024@psu.edu
#PBS -A open
#PBS -j oe
#PBS -m abe

setwd("/gpfs/group/jrl35/default/aayudh/snpeff/snp")
data <- read.delim("3_prime_UTR_variant.txt", header=FALSE)
dim(data)
data1 <- data[,1:8]
data1 = within(data1, rm(V3,V6,V7))
head(data1)
names(data1)[1]<-paste("chr")
names(data1)[2]<-paste("ps")
names(data1)[3]<-paste("REF")
names(data1)[4]<-paste("ALT")
names(data1)[5]<-paste("type")
head(data1)

rs_split <- data.frame(do.call("rbind", strsplit(as.character(data1$type), ";", fixed = TRUE)))
head(rs_split)

data1[, "type"] <- rs_split$X11
head(data1)
rs_split1 <- data.frame(do.call("rbind", strsplit(as.character(data1$type), "|", fixed = TRUE)))
head(rs_split1)
data1[, "type"] <- rs_split1$X2

head(data1)
setwd("/gpfs/group/jrl35/default/aayudh/snpeff/snp/genic")
write.csv(data1, "3_prime_UTR_variant_clean.csv", row.names = FALSE)
```

### Combine all genic and rest name it intergenic

```
setwd("~/Library/CloudStorage/OneDrive-ThePennsylvaniaStateUniversity/eGWAS_revised list/manhattan plot")
list.files()
Chr10corr <- read.csv("Reseq_preds_all_chr10_corrected.csv", header=T)
head(Chr10corr)
reseq_HS <- Chr10corr[order(Chr10corr$p_wald),]
head(reseq_HS)

setwd("~/Library/CloudStorage/OneDrive-ThePennsylvaniaStateUniversity/SNPEFF/syn and non-syn")
list.files()

syn <- read.csv("synonymous_variant.txt_clean.csv")
head(syn)
syn[, "chr"] = substring(syn$chr, 4)
rs <- data.frame(paste0('rs_', syn$chr, '_', syn$ps)) 
syn[, "X"] = rs
names(syn)[6]<-paste("rs")
syn$geneic_intergenic <- paste('genic')
head(syn)

non_syn <- read.csv("non_synonymous_missense_variant_clean.csv")
head(non_syn)
non_syn[, "chr"] = substring(non_syn$chr, 4)
rs <- data.frame(paste0('rs_', non_syn$chr, '_', non_syn$ps))
non_syn[, "X"] = rs
names(non_syn)[6]<-paste("rs")
non_syn$geneic_intergenic <- paste('genic')
head(non_syn)


setwd("~/Library/CloudStorage/OneDrive-ThePennsylvaniaStateUniversity/SNPEFF/syn and non-syn")
initiator_codon <- read.csv("initiator_codon_variant_clean.csv")
head(initiator_codon)

initiator_codon[, "chr"] = substring(initiator_codon$chr, 4)
rs <- data.frame(paste0('rs_', initiator_codon$chr, '_', initiator_codon$ps)) 
initiator_codon[, "X"] = rs
names(initiator_codon)[6]<-paste("rs")
initiator_codon <- subset(initiator_codon, type == "initiator_codon_variant")
initiator_codon$geneic_intergenic <- paste('genic')
head(initiator_codon)

setwd("~/Library/CloudStorage/OneDrive-ThePennsylvaniaStateUniversity/snp enrichment/genic vs non-genic")
list.files()

prime_3_UTR_variant <- read.csv("3_prime_UTR_variant_clean.csv")
head(prime_3_UTR_variant)

prime_3_UTR_variant[, "chr"] = substring(prime_3_UTR_variant$chr, 4)
rs <- data.frame(paste0('rs_', prime_3_UTR_variant$chr, '_', prime_3_UTR_variant$ps)) 
prime_3_UTR_variant[, "X"] = rs
names(prime_3_UTR_variant)[6]<-paste("rs")
prime_3_UTR_variant <- subset(prime_3_UTR_variant, type == "3_prime_UTR_variant")
prime_3_UTR_variant$geneic_intergenic <- paste('genic')
head(prime_3_UTR_variant)

setwd("~/Library/CloudStorage/OneDrive-ThePennsylvaniaStateUniversity/snp enrichment/genic vs non-genic")
list.files()

prime_5_UTR_start <- read.csv("5_prime_UTR_premature_start_codon_gain_variant_clean.csv")
head(prime_5_UTR_start)

prime_5_UTR_start[, "chr"] = substring(prime_5_UTR_start$chr, 4)
rs <- data.frame(paste0('rs_', prime_5_UTR_start$chr, '_', prime_5_UTR_start$ps)) 
prime_5_UTR_start[, "X"] = rs
names(prime_5_UTR_start)[6]<-paste("rs")
prime_5_UTR_start <- subset(prime_5_UTR_start, type == "5_prime_UTR_premature_start_codon_gain_variant")
prime_5_UTR_start$geneic_intergenic <- paste('genic')
head(prime_5_UTR_start)

prime_5_UTR <- read.csv("5_prime_UTR_variant_clean.csv")
head(prime_5_UTR)

prime_5_UTR[, "chr"] = substring(prime_5_UTR$chr, 4)
rs <- data.frame(paste0('rs_', prime_5_UTR$chr, '_', prime_5_UTR$ps)) 
prime_5_UTR[, "X"] = rs
names(prime_5_UTR)[6]<-paste("rs")
prime_5_UTR <- subset(prime_5_UTR, type == "5_prime_UTR_variant")
prime_5_UTR$geneic_intergenic <- paste('genic')
head(prime_5_UTR)

setwd("~/Library/CloudStorage/OneDrive-ThePennsylvaniaStateUniversity/snp enrichment/genic vs non-genic")
list.files()
non_coding <- read.csv("non_coding_transcript_exon_variant_clean.csv")
head(non_coding)

non_coding[, "chr"] = substring(non_coding$chr, 4)
rs <- data.frame(paste0('rs_', non_coding$chr, '_', non_coding$ps)) 
non_coding[, "X"] = rs
names(non_coding)[6]<-paste("rs")
non_coding <- subset(non_coding, type == "non_coding_transcript_exon_variant")
non_coding$geneic_intergenic <- paste('genic')
head(non_coding)

splice_acceptor <- read.csv("splice_acceptor_variant_clean.csv")
head(splice_acceptor)

splice_acceptor[, "chr"] = substring(splice_acceptor$chr, 4)
rs <- data.frame(paste0('rs_', splice_acceptor$chr, '_', splice_acceptor$ps)) 
splice_acceptor[, "X"] = rs
names(splice_acceptor)[6]<-paste("rs")
splice_acceptor <- subset(splice_acceptor, type == "splice_acceptor_variant&intron_variant")
splice_acceptor$geneic_intergenic <- paste('genic')
head(splice_acceptor)


splice_donor <- read.csv("splice_donor_variant_clean.csv")
head(splice_donor)

splice_donor[, "chr"] = substring(splice_donor$chr, 4)
rs <- data.frame(paste0('rs_', splice_donor$chr, '_', splice_donor$ps)) 
splice_donor[, "X"] = rs
names(splice_donor)[6]<-paste("rs")
splice_donor <- subset(splice_donor, type == "splice_donor_variant&intron_variant")
splice_donor$geneic_intergenic <- paste('genic')
head(splice_donor)

splice_region <- read.csv("splice_region_variant_clean.csv")
head(splice_region)

splice_region[, "chr"] = substring(splice_region$chr, 4)
rs <- data.frame(paste0('rs_', splice_region$chr, '_', splice_region$ps)) 
splice_region[, "X"] = rs
names(splice_region)[6]<-paste("rs")
#splice_region <- subset(splice_region, type == "splice_region_variant")
splice_region$geneic_intergenic <- paste('genic')
head(splice_region)

###add all dataset 
genic <- rbind (syn,non_syn,initiator_codon,prime_3_UTR_variant,prime_5_UTR_start,prime_5_UTR,non_coding,
                splice_acceptor,splice_donor,splice_region)
head(genic)
dim(genic)
setwd("~/Library/CloudStorage/OneDrive-ThePennsylvaniaStateUniversity/snp enrichment/genic vs non-genic")
#write.csv(genic,"1.genic.csv", row.names = FALSE)

all <- merge(reseq_HS,genic, by = "rs",all.x = TRUE)
head(all)

all$geneic_intergenic <- ifelse(is.na(all$geneic_intergenic), "intergenic", all$geneic_intergenic)
all = within(all, rm(chr.y,ps.y,REF,ALT))
names(all)[2]<-paste("chr")
names(all)[3]<-paste("ps")
head(all)
dim(all)
setwd("~/Library/CloudStorage/OneDrive-ThePennsylvaniaStateUniversity/snp enrichment/genic vs non-genic")
#write.csv(all, "1.reseq_HS_genicVSnongenic.csv", row.names = FALSE)

data <- read.csv("1.reseq_HS_genicVSnongenic.csv")
head(data)

data$logp <- -log10(data1$p_wald)

#chr rearrangement
dddf <- NULL

for(i in 1:10) { 
  chr <- subset(data, chr == i)[with(subset(data, chr == i),order(ps)),]
  dddf <- rbind(dddf,chr)}

data <- dddf 
head(data)

setwd("~/Library/CloudStorage/OneDrive-ThePennsylvaniaStateUniversity/snp enrichment/genic vs non-genic")
write.csv(data, "1.reseq_HS_genicVSnongenic_chrwise.csv", row.names = FALSE)

```

### SNP genic/intergenic enrichment

```
setwd("~/Library/CloudStorage/OneDrive-ThePennsylvaniaStateUniversity/snp enrichment/genic vs non-genic")
list.files()

data <- read.csv("1.reseq_HS_genicVSnongenic_chrwise.csv")
head(data)

Actual_fold_enrichment_genicTOintergenic <- quantile(subset(data, geneic_intergenic == "genic")$logp, probs = 0.95)/quantile(subset(data, geneic_intergenic == "intergenic")$logp, probs = 0.95)
Actual_fold_enrichment_genicTOintergenic

dddf <- NULL

for(i in 1:10) { 
  random_position <- data[sample(nrow(data), 1), ]
  f <- rbind(data[which(data$rs == random_position$rs, arr.ind=TRUE):nrow(data),],data[1:(which(data$rs == random_position$rs, arr.ind=TRUE)-1),])
  data[,"geneic_intergenic"] <- f$geneic_intergenic 
  fold_enrichment <- quantile(subset(data, geneic_intergenic == "genic")$logp, probs = 0.95)/quantile(subset(data, geneic_intergenic == "intergenic")$logp, probs = 0.95)
  dddf <- rbind(dddf,fold_enrichment)}

genicTOintergenic  <- as.data.frame(dddf)
write.csv(genicTOintergenic, "genicTOintergenic.csv", row.names = FALSE)

genicTOintergenic$ratio <- paste("genicTOintergenic")

#fold_enrich <- rbind(nonsynBYsyn,genicTOintergenic)

#write.csv(fold_enrich, "1.fold_enrichment.csv", row.names = FALSE)

###Now miniAF correction

##Minor allele frequency correction, whatever aboive 0.5 af just deduct 0.5.
data_below_.5 <- subset(data, af < .5)
data_above_.5 <- subset(data, af > .5)

data_above_.5$af <- data_above_.5$af -0.5 

data_cortd <- rbind(data_below_.5,data_above_.5)

data_cortd <- data_cortd[with(data_cortd,order(chr)),]
head(data_cortd)
##minaf 0-0.1
minaf_0.1 <- subset(data_cortd, af < .1) 
minaf_0.1 <- minaf_0.1[with(minaf_0.1,order(chr)),]
head(minaf_0.1)

Actual_fold_enrichment_minaf_0.1 <- quantile(subset(minaf_0.1, geneic_intergenic == "genic")$logp, probs = 0.95)/quantile(subset(minaf_0.1, geneic_intergenic == "intergenic")$logp, probs = 0.95)
Actual_fold_enrichment_minaf_0.1

dddf <- NULL

for(i in 1:100) { 
  random_position <- minaf_0.1[sample(nrow(minaf_0.1), 1), ]
  f <- rbind(minaf_0.1[which(minaf_0.1$rs == random_position$rs, arr.ind=TRUE):nrow(minaf_0.1),],minaf_0.1[1:(which(minaf_0.1$rs == random_position$rs, arr.ind=TRUE)-1),])
  minaf_0.1[,"geneic_intergenic"] <- f$geneic_intergenic 
  fold_enrichment <- quantile(subset(minaf_0.1, geneic_intergenic == "genic")$logp, probs = 0.95)/quantile(subset(minaf_0.1, geneic_intergenic == "intergenic")$logp, probs = 0.95)
  dddf <- rbind(dddf,fold_enrichment)}

genicTOintergenic_minaf_0.1  <- as.data.frame(dddf)
setwd("~/Library/CloudStorage/OneDrive-ThePennsylvaniaStateUniversity/snp enrichment/minaf_correction")
write.csv(genicTOintergenic_minaf_0.1, "genicTOintergenic_minaf_0.1.csv", row.names = FALSE)

genicTOintergenic_minaf_0.1$ratio <- paste("genicTOintergenic_minaf_0.1")

##minaf 0.1-0.2
minaf_0.2 <- subset(data_cortd, af < 0.2 & af > 0.1) 
minaf_0.2 <- minaf_0.2[with(minaf_0.2,order(chr)),]
head(minaf_0.2)

Actual_fold_enrichment_minaf_0.2 <- quantile(subset(minaf_0.2, geneic_intergenic == "genic")$logp, probs = 0.95)/quantile(subset(minaf_0.2, geneic_intergenic == "intergenic")$logp, probs = 0.95)
Actual_fold_enrichment_minaf_0.2

dddf <- NULL

for(i in 1:100) { 
  random_position <- minaf_0.2[sample(nrow(minaf_0.2), 1), ]
  f <- rbind(minaf_0.2[which(minaf_0.2$rs == random_position$rs, arr.ind=TRUE):nrow(minaf_0.2),],minaf_0.2[1:(which(minaf_0.2$rs == random_position$rs, arr.ind=TRUE)-1),])
  minaf_0.2[,"geneic_intergenic"] <- f$geneic_intergenic 
  fold_enrichment <- quantile(subset(minaf_0.2, geneic_intergenic == "genic")$logp, probs = 0.95)/quantile(subset(minaf_0.2, geneic_intergenic == "intergenic")$logp, probs = 0.95)
  dddf <- rbind(dddf,fold_enrichment)}

genicTOintergenic_minaf_0.2  <- as.data.frame(dddf)
setwd("~/Library/CloudStorage/OneDrive-ThePennsylvaniaStateUniversity/snp enrichment/minaf_correction")
write.csv(genicTOintergenic_minaf_0.2, "genicTOintergenic_minaf_0.2.csv", row.names = FALSE)

genicTOintergenic_minaf_0.2$ratio <- paste("genicTOintergenic_minaf_0.2")


##minaf 0.2-0.3
minaf_0.3 <- subset(data_cortd, af < 0.3 & af > 0.2) 
minaf_0.3 <- minaf_0.3[with(minaf_0.3,order(chr)),]
head(minaf_0.3)

Actual_fold_enrichment_minaf_0.3 <- quantile(subset(minaf_0.3, geneic_intergenic == "genic")$logp, probs = 0.95)/quantile(subset(minaf_0.3, geneic_intergenic == "intergenic")$logp, probs = 0.95)
Actual_fold_enrichment_minaf_0.3

dddf <- NULL

for(i in 1:100) { 
  random_position <- minaf_0.3[sample(nrow(minaf_0.3), 1), ]
  f <- rbind(minaf_0.3[which(minaf_0.3$rs == random_position$rs, arr.ind=TRUE):nrow(minaf_0.3),],minaf_0.3[1:(which(minaf_0.3$rs == random_position$rs, arr.ind=TRUE)-1),])
  minaf_0.3[,"geneic_intergenic"] <- f$geneic_intergenic 
  fold_enrichment <- quantile(subset(minaf_0.3, geneic_intergenic == "genic")$logp, probs = 0.95)/quantile(subset(minaf_0.3, geneic_intergenic == "intergenic")$logp, probs = 0.95)
  dddf <- rbind(dddf,fold_enrichment)}

genicTOintergenic_minaf_0.3  <- as.data.frame(dddf)
setwd("~/Library/CloudStorage/OneDrive-ThePennsylvaniaStateUniversity/snp enrichment/minaf_correction")
write.csv(genicTOintergenic_minaf_0.3, "genicTOintergenic_minaf_0.3.csv", row.names = FALSE)

genicTOintergenic_minaf_0.3$ratio <- paste("genicTOintergenic_minaf_0.3")

##minaf 0.3-0.4
minaf_0.4 <- subset(data_cortd, af < 0.4 & af > 0.3) 
minaf_0.4 <- minaf_0.4[with(minaf_0.4,order(chr)),]
head(minaf_0.4)

Actual_fold_enrichment_minaf_0.4 <- quantile(subset(minaf_0.4, geneic_intergenic == "genic")$logp, probs = 0.95)/quantile(subset(minaf_0.4, geneic_intergenic == "intergenic")$logp, probs = 0.95)
Actual_fold_enrichment_minaf_0.4

dddf <- NULL

for(i in 1:100) { 
  random_position <- minaf_0.4[sample(nrow(minaf_0.4), 1), ]
  f <- rbind(minaf_0.4[which(minaf_0.4$rs == random_position$rs, arr.ind=TRUE):nrow(minaf_0.4),],minaf_0.4[1:(which(minaf_0.4$rs == random_position$rs, arr.ind=TRUE)-1),])
  minaf_0.4[,"geneic_intergenic"] <- f$geneic_intergenic 
  fold_enrichment <- quantile(subset(minaf_0.4, geneic_intergenic == "genic")$logp, probs = 0.95)/quantile(subset(minaf_0.4, geneic_intergenic == "intergenic")$logp, probs = 0.95)
  dddf <- rbind(dddf,fold_enrichment)}

genicTOintergenic_minaf_0.4  <- as.data.frame(dddf)
setwd("~/Library/CloudStorage/OneDrive-ThePennsylvaniaStateUniversity/snp enrichment/minaf_correction")
write.csv(genicTOintergenic_minaf_0.4, "genicTOintergenic_minaf_0.4.csv", row.names = FALSE)

genicTOintergenic_minaf_0.4$ratio <- paste("genicTOintergenic_minaf_0.4")


##minaf 0.4-0.5
minaf_0.5 <- subset(data_cortd, af > 0.4) 
minaf_0.5 <- minaf_0.5[with(minaf_0.5,order(chr)),]
head(minaf_0.5)

Actual_fold_enrichment_minaf_0.5 <- quantile(subset(minaf_0.5, geneic_intergenic == "genic")$logp, probs = 0.95)/quantile(subset(minaf_0.5, geneic_intergenic == "intergenic")$logp, probs = 0.95)
Actual_fold_enrichment_minaf_0.5

dddf <- NULL

for(i in 1:100) { 
  random_position <- minaf_0.5[sample(nrow(minaf_0.5), 1), ]
  f <- rbind(minaf_0.5[which(minaf_0.5$rs == random_position$rs, arr.ind=TRUE):nrow(minaf_0.5),],minaf_0.5[1:(which(minaf_0.5$rs == random_position$rs, arr.ind=TRUE)-1),])
  minaf_0.5[,"geneic_intergenic"] <- f$geneic_intergenic 
  fold_enrichment <- quantile(subset(minaf_0.5, geneic_intergenic == "genic")$logp, probs = 0.95)/quantile(subset(minaf_0.5, geneic_intergenic == "intergenic")$logp, probs = 0.95)
  dddf <- rbind(dddf,fold_enrichment)}

genicTOintergenic_minaf_0.5  <- as.data.frame(dddf)
setwd("~/Library/CloudStorage/OneDrive-ThePennsylvaniaStateUniversity/snp enrichment/minaf_correction")
write.csv(genicTOintergenic_minaf_0.5, "genicTOintergenic_minaf_0.5.csv", row.names = FALSE)

genicTOintergenic_minaf_0.5$ratio <- paste("genicTOintergenic_minaf_0.5")


###bind all

fold_enrich <- rbind(genicTOintergenic,genicTOintergenic_minaf_0.1,genicTOintergenic_minaf_0.2,genicTOintergenic_minaf_0.3,genicTOintergenic_minaf_0.4,genicTOintergenic_minaf_0.5)

#setwd("~/Library/CloudStorage/OneDrive-ThePennsylvaniaStateUniversity/snp enrichment/minaf_correction")
write.csv(fold_enrich, "1.fold_enrichment_MAF_correction_with_main.csv", row.names = FALSE)

setwd("~/Library/CloudStorage/OneDrive-ThePennsylvaniaStateUniversity/snp enrichment/genic vs non-genic")
list.files()
df <- read.csv("1.fold_enrichment_MAF_correction_with_main.csv") 
head(df)
##plotiing
library(ggplot2)
library(tidyverse)
library(ggtext)
library(normentR)
library(tidyverse)
library(gridExtra)
library(grid)
library(gridtext)
library(palmerpenguins)

setwd("~/Library/CloudStorage/OneDrive-ThePennsylvaniaStateUniversity/snp enrichment/genic vs non-genic")
tiff("1.genicVSnon-genic_HS_100perm_MAFcor.tiff", width = 5, height = 6, units = 'in', res = 300)
ggplot(df, aes(x=X95., y=ratio, color=ratio))+
  geom_point(size=1)+
  labs(x="Fold enrichment",y=NULL)+
  scale_x_continuous(breaks=seq(0.914,1.081,0.055))+
  scale_color_manual(values = c("#B4AF46", "#B4464B","#B47846","#7846B4","#82B446","deeppink"))+
  geom_vline(xintercept = 1, linetype="dashed", color = "slategrey", size=1)+
  annotate("pointrange", x = Actual_fold_enrichment_genicTOintergenic , y = 1, ymin = 1.010, ymax = 1.008,colour = "#4682B4",size = 2)+
  annotate("pointrange", x = Actual_fold_enrichment_minaf_0.1 , y = 2, ymin = 1.010, ymax = 1.008,colour = "#4682B4",size = 2)+
  annotate("pointrange", x = Actual_fold_enrichment_minaf_0.2 , y = 3, ymin = 2.010, ymax = 2.008,colour = "#4682B4",size = 2)+
  annotate("pointrange", x = Actual_fold_enrichment_minaf_0.3 , y = 4, ymin = 1.010, ymax = 1.008,colour = "#4682B4",size = 2)+
  annotate("pointrange", x = Actual_fold_enrichment_minaf_0.4 , y = 5, ymin = 2.010, ymax = 2.008,colour = "#4682B4",size = 2)+
  annotate("pointrange", x = Actual_fold_enrichment_minaf_0.5 , y = 6, ymin = 2.010, ymax = 2.008,colour = "#4682B4",size = 2)+
  theme_classic()+
  theme(legend.position="none",text=element_text(size=16, colour = "black", family="Times New Roman"),
        axis.line = element_line(size=0.5, colour = "black"),
        axis.text.x=element_text(colour="black", size = 16),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=1))
dev.off()
```

-----
<div id='id-section18'/>

## Chapter 18: Python based annotations

### Create a python script

**Sbicolor_annot_5kb.csv consist of columns - "chr","gene","start","stop","direction","sobicID","start_5kb","stop_5kb"*
*top1K_RDA1Sq_nlyAfrica.csv consist of columns - "SNPS","RDA1","RDA2","RDA1_squared","RDA2_squared","sq_RDA1plus2","chr","ps"**

```
#!/usr/bin/env python
# coding: utf-8

# In[7]:


import pandas as pd
import multiprocess

# read in the first CSV file with two columns
df1 = pd.read_csv('Sbicolor_annot_5kb.csv')


# In[8]:


# read in the second CSV file with one column
df2 = pd.read_csv('top1K_RDA1Sq_nlyAfrica.csv')


# In[9]:


import multiprocessing

# create a function that filters df1 for a single PS value and returns the filtered DataFrame
def filter_df1_for_ps(ps, df1, df2):
    mask = (df1['chr'] == df2.loc[df2['ps'] == ps, 'chr'].iloc[0]) & (df1['start_5kb'] <= ps) & (df1['stop_5kb'] >= ps)
    filtered_df1 = df1[mask]
    if not filtered_df1.empty:
        filtered_df1['SNPS'] = df2.loc[df2['ps'] == ps, 'SNPS'].iloc[0]
        filtered_df1['RDA1_squared'] = df2.loc[df2['ps'] == ps, 'RDA1_squared'].iloc[0]
        filtered_df1['ps'] = ps
    return filtered_df1

# define the number of CPU cores to use for parallelization
num_cores = multiprocessing.cpu_count()

# create a pool of worker processes
pool = multiprocessing.Pool(num_cores)

# apply the filter_df1_for_ps function to each PS value in df2 in parallel
results = [pool.apply_async(filter_df1_for_ps, args=(ps, df1, df2)) for ps in df2['ps']]

# combine the results into a single DataFrame
merged_df = pd.concat([result.get() for result in results if not result.get().empty], ignore_index=True)
merged_df.to_csv('1.RDA1_suared_nlyAfrica_Annotations.csv', index=False)
```

### Now to execute the .py script use this

qsub this .sh file now.
```
#!/bin/bash

#PBS -l nodes=1:ppn=8
#PBS -l walltime=12:00:00
#PBS -l pmem=24gb
#PBS -M azd6024@psu.edu
#PBS -A open
#PBS -j oe
#PBS -m abe

# Load any required modules
source ~/.bashrc
module load anaconda3
conda activate my_env

# Change to the directory where the code and datafiles are located
cd /gpfs/group/jrl35/default/aayudh/africa_drought_RDA

# Run the Python script
python3 annotation.py
```

### Use the Python based output to R
```
##Annotation-------
#do python script 
setwd("~/Library/CloudStorage/OneDrive-ThePennsylvaniaStateUniversity/RDA_drought/RDA_only_Africa/annotations")
list.files()

gwas_annot <- read.csv("1.RDA1_suared_nlyAfrica_Annotations.csv")
head(gwas_annot)

library(data.table)
gwas_annot2 <- unique(setDT(gwas_annot)[order(SNPS, -SNPS)], by = "SNPS")
dim(gwas_annot2)

setwd("~/Library/CloudStorage/OneDrive-ThePennsylvaniaStateUniversity/acessions and cords")
Sorghum_mdata <- read.csv("Sbicolor_454_v3.1.1.annotation_info.csv", header = T)
head(Sorghum_mdata)

names(Sorghum_mdata)[2]<-paste("sobicID")

annotation_Loc <- merge(Sorghum_mdata,gwas_annot2, by = "sobicID")
dim(annotation_Loc)
head(annotation_Loc)

library(data.table)
annotation_Loc_1 <- unique(setDT(annotation_Loc)[order(SNPS, -SNPS)], by = "SNPS")
dim(annotation_Loc_1)
head(annotation_Loc_1)
annotation_Loc_2 <- annotation_Loc_1[order(-annotation_Loc_1$RDA1_squared),]
head(annotation_Loc_2)

library(data.table)
new_order <- c("SNPS", "RDA1_squared", "sobicID", "arabi.defline", "rice.defline", 
               setdiff(names(annotation_Loc_2), c("rs", "p_wald", "sobicID", "arabi.defline", "rice.defline")))
setcolorder(annotation_Loc_2, new_order)

head(annotation_Loc_2)

setwd("~/Library/CloudStorage/OneDrive-ThePennsylvaniaStateUniversity/univariate_drought/annotations")
list.files()
write.csv(annotation_Loc_2,"1.RDA1_squared_nlyAfrica_annotations_final", row.names = FALSE)

### before this step you need to manually edit the annot list

setwd("~/Library/CloudStorage/OneDrive-ThePennsylvaniaStateUniversity/univariate_drought/annotations")
list.files()
top10 <- read.csv("top10_Wsrepro_annotations.csv")
#export table
library("gridExtra")
tiff("top10_WSrepro_annotations.tiff", width = 10, height = 5, units = 'in', res = 300)
par(family="Arial")
grid.table(top10)
dev.off()
```

-----
<div id='id-section19'/>

## Chapter 19: Individual candidate analysis

### Expression profiles and Co-expression network 

Go to https://plantnexus.ohio.edu/e/?plant=sorghum. Then change the sobic format to SORBI_3002G129900 format.  

![alt text](https://github.com/aayudh454/Lasky-Morris-Lab-Sorghum-project/blob/main/expression_plantnexus.png)

For co expression network, start from the global option, play around with mutual rank filter.

![alt text](https://github.com/aayudh454/Lasky-Morris-Lab-Sorghum-project/blob/main/Co-expression_network.png)

-----
<div id='id-section20'/>

## Chapter 20: snpeff for Bromus tectorum

##### Getting to roar collab
ssh azd6024@submit.hpc.psu.edu

my folder: cd /storage/group/jrl35/default/aayudh/custom_genome_snpeff

### Step 1: Database Build

1) Inslide the snpeff directory make this path **/data/BT_assembly**. Put **sequences.fa** which is your refernce sequence (just use 1-7 chromosomes). then **genes.gtf**. Now we had to cobvert gff to gtf--

```
##gff to gtf (install gffread)
/storage/group/jrl35/default/aayudh/cheatgrass/gffread/gffread filtered_genes.gff -T -o filtered_genes.gtf
```

2) from sequence.fa to cds.fa--

```
/storage/group/jrl35/default/aayudh/cheatgrass/gffread/gffread -w /storage/group/jrl35/default/aayudh/cheatgrass/cds_protein/BT_cds.fa -g /storage/group/jrl35/default/aayudh/cheatgrass/cds_protein/filtered_sequences.fa /storage/group/jrl35/default/aayudh/cheatgrass/cds_protein/filtered_genes.gff
```

3) make a config file  

```
vi snpEffect.config
# BT_assembly genome
BT_assembly.genome : BT_assembly
BT_assembly.reference : /storage/group/jrl35/default/aayudh/custom_genome_snpeff/snpEff/data/BT_assembly/sequences.fa
```
4) Now make a script to execute the database build

```
#!/bin/bash

#PBS -l nodes=1:ppn=8
#PBS -l walltime=12:00:00
#PBS -l pmem=24gb
#PBS -M azd6024@psu.edu
#PBS -A open
#PBS -j oe
#PBS -m abe

# Define variables
SNPEFF_DIR="/storage/group/jrl35/default/aayudh/custom_genome_snpeff/snpEff"
CONFIG_FILE="$SNPEFF_DIR/snpEffect.config"

# Build the snpEff database for BT_assembly
java -Xmx4g -jar $SNPEFF_DIR/snpEff.jar build -c $CONFIG_FILE -noCheckProtein -gtf22 -v BT_assembly
```

### Step 2: Annnotation

1) First, in your file change bt to chr exactly the way its in yhour refrence file and database. 

```
#!/bin/bash

#PBS -l nodes=1:ppn=8
#PBS -l walltime=12:00:00
#PBS -l pmem=24gb
#PBS -M azd6024@psu.edu
#PBS -A open
#PBS -j oe
#PBS -m abe

cd /storage/group/jrl35/default/aayudh/custom_genome_snpeff/vcf_gz_modufy

zcat 307BRTE.imputed_filtered.vcf.gz | \
awk 'BEGIN{OFS="\t"}
    /^#/ {print; next}
    {split($1,a,"chr");
     if(length(a)==2 && a[2]+0>0)
        $1="Bt" sprintf("%02d", a[2]);
     print}' | \
gzip > chrBt_modified_307BRTE.imputed_filtered.vcf.gz
```

2) Convert to a vcf file

```
#!/bin/bash

#PBS -l nodes=1:ppn=8
#PBS -l walltime=12:00:00
#PBS -l pmem=24gb
#PBS -M azd6024@psu.edu
#PBS -A open
#PBS -j oe
#PBS -m abe

cd /storage/group/jrl35/default/aayudh/custom_genome_snpeff/vcf_gz_modufy

zcat chrBt_modified_307BRTE.imputed_filtered.vcf.gz > chr_modified_variant_filtered.vcf
```

3) Now run the annnotation script.

```
#!/bin/bash

#PBS -l nodes=1:ppn=8
#PBS -l walltime=12:00:00
#PBS -l pmem=24gb
#PBS -M azd6024@psu.edu
#PBS -A open
#PBS -j oe
#PBS -m abe

# Define variables
SNPEFF_DIR="/storage/group/jrl35/default/aayudh/custom_genome_snpeff/snpEff"
CONFIG_FILE="$SNPEFF_DIR/snpEffect.config"

# Build the snpEff database for BT_assembly
#java -Xmx4g -jar $SNPEFF_DIR/snpEff.jar build -c $CONFIG_FILE -noCheckProtein -gtf22 -v BT_assembly

java -Xmx4g -jar $SNPEFF_DIR/snpEff.jar ann BT_assembly -c $CONFIG_FILE /storage/group/jrl35/default/aayudh/custom_genome_snpeff/vcf_gz_modufy/chr_modified_variant_filtered.vcf > /storage/group/jrl35/default/aayudh/custom_genome_snpeff/1_BT_variants_annotated.vcf
```

### Step 3: Find "HIGH" effect variant

Find the 'HIGH' variant--

```
#!/bin/bash

#PBS -l nodes=1:ppn=8
#PBS -l walltime=12:00:00
#PBS -l pmem=24gb
#PBS -M azd6024@psu.edu
#PBS -A open
#PBS -j oe
#PBS -m abe

grep -v "^#" 1_BT_variants_annotated.vcf | awk '{if ($8 ~ /HIGH/) print $0}' | wc -l

grep -v "^#" 1_BT_variants_annotated.vcf | awk '{if ($8 ~ /HIGH/) print $0}' > 2_HIGH_BT_variants_.vcf
```

take the header 1_BT_variants_annotated.vcf----

```
grep "^#" 1_BT_variants_annotated.vcf > header.txt
cat 2_HIGH_BT_variants.vcf >> header.txt

mv header.txt 3_HIGH_BT_variants.txt

```

modify high variant

```
data <- read.csv("New_out.csv", header=FALSE, check.names=FALSE)

affected_rows <- grep("GT", data$V9, invert = TRUE, value = FALSE)

for (row in affected_rows) {
  # Shift columns V9 to V316 one position to the left
  # Now we are considering V317 as the last column to shift data into
  data[row, 8:316] <- data[row, 9:317]
}
data$V317[affected_rows] <- NA

data1 <- data [,1:316]

df <- read.csv("out.csv", header=FALSE, check.names=FALSE)
df = df[, 1:316]

names(data1) <- as.character(unlist(df[1, ]))

library(tidyr)
library(dplyr)

# Assuming you know or guess a max number of splits, e.g., 10 here
# Replace `10` with the actual number of splits you expect
data1 <- data1 %>%
  separate(INFO, into = paste0("INFO_", 1:10), sep = "\\|", fill = "right")

# Renaming columns using dplyr
data1 <- data1 %>%
  rename(Impact = INFO_3, Annotation = INFO_2, Allele = INFO_1,Gene_Name = INFO_4, Gene_ID = INFO_5,)

data2 <- data1 %>%
  filter(Impact == "HIGH")
```
count high variant

```
# Initialize a vector to store the counts
counts <- integer()

# Loop through columns 19 to 325 (assuming column 325 exists, adjust if needed)
for(i in 19:325) {
  # Count how many times "1|1" appears in each column
  counts[i - 18] <- sum(data2[,i] == "1|1", na.rm = TRUE)
}

# Create a data frame with the results
# Extracting names of the columns 19 to 325 for the names in our table
column_names <- names(data2)[19:325]

# Create the final table
count_table <- data.frame(Genotype_ID = column_names, No_HIGH_variant = counts)

# Print the table
print(count_table)
```
email

```
echo "Please find the attached CSV file." | mail -s "High Variant Count per Genotype CSV" -a /storage/group/jrl35/default/aayudh/bromus_tectorum_snpeff/1.high_variant_count_per_genotype.csv aayudhdas@gmail.com

```
visualization- 

```
setwd("~/Library/CloudStorage/OneDrive-ThePennsylvaniaStateUniversity/snpeff_bromus")
list.files()

data <- read.csv("1.High_variant_rawDATA.csv")
head(data)
# Assuming your data is stored in a data.frame or data.table named 'filtered_data'
number_of_unique_genes <- length(unique(data$gene))

# Print the number of unique genes
print(number_of_unique_genes)


#write.csv(filtered_data, "GEB0017_130_SNPS_Intron_filtered_data.csv", row.names = FALSE)
library(ggplot2)
png("pie_chart.png", width = 10, height = 6, units = 'in', res = 300)
ggplot(data, aes(x = "", fill = Annotation)) +
  geom_bar(width = 1, stat = "count") +
  geom_text(stat = "count", aes(label = stat(count)), position = position_stack(vjust = 0.5)) +
  coord_polar("y") +
  labs(title = paste(nrow(data), " HIGH variant"), y = NULL, x = NULL) +
  theme_classic(base_family = "Arial") +
  theme_void() +
  theme(plot.title = element_text(hjust = 0.5))  # Center the title
dev.off()  

df <- read.csv("1.high_variant_count_per_genotype.csv")
head(df)

png("bar_chart.png", width = 6, height = 6, units = 'in', res = 300)
ggplot(df, aes(x = reorder(factor(Genotype_ID), No_HIGH_variant), y = No_HIGH_variant, fill = No_HIGH_variant)) +
  geom_col() +
  scale_fill_gradient(low = "blue", high = "red") +
  labs(x = "Genotype ID", y = "Number of HIGH Variants", fill = "No. HIGH Variant") +
  coord_flip() + # This inverts the plot
  theme_classic(base_family = "Arial") +
  theme(axis.text.x = element_text(angle = 0)) # Adjust the angle if needed for clarity
dev.off()

df1 <- df[order(-df$No_HIGH_variant),][1:50,]

png("bar_top50_chart.png", width = 6, height = 6, units = 'in', res = 300)
ggplot(df1, aes(x = reorder(factor(Genotype_ID), No_HIGH_variant), y = No_HIGH_variant, fill = No_HIGH_variant)) +
  geom_col() +
  scale_fill_gradient(low = "blue", high = "red") +
  labs(x = "Genotype ID", y = "Number of HIGH Variants", fill = "No. HIGH Variant") +
  coord_flip() + # This inverts the plot
  theme_classic(base_family = "Arial") +
  theme(axis.text.x = element_text(angle = 0))# Adjust the angle if needed for clarity
dev.off()
```
