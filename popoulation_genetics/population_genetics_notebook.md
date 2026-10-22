# Population Genomics Notebook

**Course:** Intro to Ecological Genomics - Fall 2026

**Name:** Walt Regan-Loomis

------------------------------------------------------------------------

## \*\* - Template

-   

**Working Directory**

`/gpfs1/home/w/r/wreganlo/eco_genomics_2026/population_genomics`

**Input Files:**

`none`

**Output Files:**

`/gpfs1/home/w/r/wreganlo/eco_genomics_2026/transcriptomics/population_genomics_notebook.md`

`/gpfs1/home/w/r/wreganlo/eco_genomics_2026/transcriptomics/myscripts/10.22.2026.genetic.diversity`

**Programs and dependencies:**

-   `R version 4.5.1`

-   `R-Studio`

**Scripts:**

`none`

**Code:**

``` r
```

**Table:**

| Col1 | Col2 | Col3 |
|------|------|------|
|      |      |      |
|      |      |      |
|      |      |      |

**Image:**

**Notes / Observations:**

-   

**Next Steps:**

------------------------------------------------------------------------

## 10.20.2026 - VCF Messing Around

-   Setting up and copying files into my data

-   using plink to run analysis

**Working Directory**

`/gpfs1/home/w/r/wreganlo/eco_genomics_2026/population_genomics`

**Input Files:**

`~/projects/eco_genomics_2026/popoulation_genetics/mydata/Lg_LD_0.2eigenval`

`~/projects/eco_genomics_2026/popoulation_genetics/mydata/Lg_LD_0.2eigenvec`

**Output Files:**

`/gpfs1/home/w/r/wreganlo/eco_genomics_2026/population_genetics/population_genomics_notebook.md`

`/gpfs1/home/w/r/wreganlo/eco_genomics_2026/population_genetics/myscripts/10.20.2026_VCF_PCA.R`

`/gpfs1/home/w/r/wreganlo/eco_genomics_2026/population_genetics/myscripts/10.20.2026_VCF_admixture.R`

**Programs and dependencies:**

-   `R version 4.5.1`

-   `R-Studio`

-   `PLINK v1.9.0-b.7.7 64-bit (22 Oct 2024)`

**Scripts:**

`DATA_CIR=/users/w/r/wreganlo/projects/eco_genomics_2026/popoulation_genetics/mydata/`

`VCF=$DATA_DIR/Lottia_gigantea_flt3_sub0.01.vcf.gz`

**Code:**

``` r
##We should start organising our results file a bit, lets make a sub directory called PopulationStructure
mkdir -p ~/path/to/my/directory/PopulationGenomics/myresults/PopulationStructure 

# where the VCF lives (make sure this is the right path!)
DATA_DIR=/path/to/my/directory/PopulationGenomics/mydata
VCF=$DATA_DIR/Lottia_gigantea_pl_flt3.renamed.vcf.gz

# check the variable holds what you expect, and that the file is really there
echo $VCF
ls -lh $VCF

cd $DATA_DIR
pwd ##check you are in the right place

##We should start organising our results file a bit, lets make a sub directory called PopulationStructure
mkdir -p ~/path/to/my/directory/PopulationGenomics/myresults/PopulationStructure 

# where the VCF lives (make sure this is the right path!)
DATA_DIR=/path/to/my/directory/PopulationGenomics/mydata
VCF=$DATA_DIR/Lottia_gigantea_pl_flt3.renamed.vcf.gz

# check the variable holds what you expect, and that the file is really there
echo $VCF
ls -lh $VCF

cd $DATA_DIR
pwd ##check you are in the right place

##We should start organising our results file a bit, lets make a sub directory called PopulationStructure
mkdir -p ~/path/to/my/directory/PopulationGenomics/myresults/PopulationStructure 

# where the VCF lives (make sure this is the right path!)
DATA_DIR=/path/to/my/directory/PopulationGenomics/mydata
VCF=$DATA_DIR/Lottia_gigantea_pl_flt3.renamed.vcf.gz

# check the variable holds what you expect, and that the file is really there
echo $VCF
ls -lh $VCF

cd $DATA_DIR
pwd ##check you are in the right place
```

**Code for Admixture:**

``` bash
# keep a copy of the original .bim file (only run this line once! see below)
cp Lg_LD_0.2.bim Lg_LD_0.2.bim.orig

# read the original, replace column 1 (the chromosome) with 0, and write the new .bim
# OFS="\t" keeps the columns tab-separated
awk 'BEGIN{OFS="\t"} {$1="0"; print}' Lg_LD_0.2.bim.orig > Lg_LD_0.2.bim

# compare before and after
head -n 3 Lg_LD_0.2.bim.orig
head -n 3 Lg_LD_0.2.bim

# keep a copy of the original .bim file (only run this line once! see below)
cp Lg_LD_0.2.bim Lg_LD_0.2.bim.orig

# read the original, replace column 1 (the chromosome) with 0, and write the new .bim
# OFS="\t" keeps the columns tab-separated
awk 'BEGIN{OFS="\t"} {$1="0"; print}' Lg_LD_0.2.bim.orig > Lg_LD_0.2.bim

# compare before and after
head -n 3 Lg_LD_0.2.bim.orig
head -n 3 Lg_LD_0.2.bim

##Running as Bash Job:

#!/bin/bash
#SBATCH --job-name=admixture
#SBATCH --partition=general
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=24:00:00
#SBATCH --output=./logs/%x.%j.out
#SBATCH --error=./logs/%x.%j.err
#SBATCH --mail-type=END
#SBATCH --mail-user=NETID@uvm.edu

cd ~/popgen/class2
module load admixture

INPUT="Lg_LD_0.2.bed"
THREADS=${SLURM_CPUS_PER_TASK}   # use as many threads as CPUs we asked for

rm -f cv_errors.tsv
for K in $(seq 1 10); do
    admixture --cv -j${THREADS} "${INPUT}" ${K} | tee "admixture_K${K}.out"
    CV=$(grep "CV error" "admixture_K${K}.out" | awk '{print $NF}')
    echo -e "${K}\t${CV}" >> cv_errors.tsv
done
```

**Image:**

![](~/projects/eco_genomics_2026/popoulation_genetics/myresults/population_structure/variance_0.2.png)

![](~/projects/eco_genomics_2026/popoulation_genetics/myresults/population_structure/Lg_0.2_pca.png)

![](~/projects/eco_genomics_2026/popoulation_genetics/myresults/population_structure/Lg_0.2_pca_map.png)

![](~/projects/eco_genomics_2026/popoulation_genetics/myresults/population_structure/Lg_0.2_pca_latitude.png)

![](~/projects/eco_genomics_2026/popoulation_genetics/myresults/population_structure/cv_error.png)

![](~/projects/eco_genomics_2026/popoulation_genetics/myresults/population_structure/admixture_K1-5.png)

**Notes / Observations:**

-   

**Next Steps:**

-   More in class learning

------------------------------------------------------------------------

## 10.22.2026 - Genetic Diversity

-   Working more with the data and finding genetic diversity metrics (nucleotide diversity (π), Tajima’s D and heterozygosity measure)

-   Made Pi folder in myresults

**Working Directory**

`/gpfs1/home/w/r/wreganlo/eco_genomics_2026/population_genomics`

**Input Files:**

`~/projects/eco_genomics_2026/popoulation_genetics/mydata/pi/BOD.windowed.pi`

`~/projects/eco_genomics_2026/popoulation_genetics/mydata/pi/PGR.windowed.pi`

**Output Files:**

`/gpfs1/home/w/r/wreganlo/eco_genomics_2026/population_genetics/population_genomics_notebook.md`

**Programs and dependencies:**

-   `R version 4.5.1`

-   `R-Studio`

-   `VCF Tools 0.1.16`

**Scripts:**

`DATA_CIR=/users/w/r/wreganlo/projects/eco_genomics_2026/popoulation_genetics/mydata/`

`VCF=$DATA_DIR/Lottia_gigantea_flt3_sub0.01.vcf.gz`

**Code:**

``` bash
cd ~/path/to/my/files/PopulationGenomics/mydata

module load vcftools
vcftools --version   # note the version in your notebook

# the 1% subset from Class 1 (make sure this is the right path!)
# We will assign our data directory and VCF now so we dont have to move around to find this. We will now call #$VCF when needed
DATA_DIR=~</path/to/my/files>/PopulationGenomics/mydata ###remember change whats in the <>
VCF=$DATA_DIR/Lottia_gigantea_pl_flt3_noMAF_sub0.01.vcf.gz

# check which columns hold the sample ID and the site
head sample_metadata.csv

# skip the header line (tail -n +2), keep columns 1 and 2 (cut),
# and swap the comma between them for a tab (tr)
tail -n +2 sample_metadata.csv | cut -d ',' -f 1,2 | tr ',' '\t' > popmap.txt
head popmap.txt ##check it looks correct

# check which columns hold the sample ID and the site
head sample_metadata.csv

# skip the header line (tail -n +2), keep columns 1 and 2 (cut),
# and swap the comma between them for a tab (tr)
tail -n +2 sample_metadata.csv | cut -d ',' -f 1,2 | tr ',' '\t' > popmap.txt
head popmap.txt ##check it looks correct

mkdir -p pops ##makes subfolder in mydata (where you should be, if unsure use pwd)

# take each site name once (cut column 2, then sort -u), and for each one
# write the IDs of its individuals (column 1) to pops/SITE.txt
cut -f 2 popmap.txt | sort -u | while read -r POP; do
    awk -v pop="$POP" -F '\t' '$2 == pop {print $1}' popmap.txt > pops/${POP}.txt
done

# how many individuals are there per site?
wc -l pops/*.txt
```

**Images:**

![](myresults/pi/pi_by_site.png)

![](myresults/pi/pi_vs_latitude.png)

**Notes:**

![](images/Code Cheat Sheet for Day 4 Genetic Diversity.png)
