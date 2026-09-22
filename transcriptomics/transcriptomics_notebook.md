# Transcriptomics Notebook

**Course:** Intro to Ecological Genomics - Fall 2026

**Name:** Walt Regan-Loomis

------------------------------------------------------------------------

## 9.15.2026 - Setting up lab notebook and learning markdown

-   Setting up transcriptomics notebook

-   Learn how to take notes in markdown

-   Push notes to github

**Working Directory**

`/gpfs1/home/w/r/wreganlo/eco_genomics_2026/transcriptomics`

**Input Files:**

`none`

**Output Files:**

`/gpfs1/home/w/r/wreganlo/eco_genomics_2026/transcriptomics/transcriptomics_notebook.md`

**Programs and dependencies:**

-   `R version 4.5.1`

-   `R-Studio`

**Scripts:**

`none`

**Code:**

``` r
print("Hello World")
```

**Table:**

| Col1 | Col2 | Col3 |
|------|------|------|
|      |      |      |
|      |      |      |
|      |      |      |

**Image:**

![](images/Markdown_Cheat_Sheet){width="360"}

**Notes / Observations:**

-   Very nice graph/ interpretation

**Next Steps:**

Lets dive into project!

------------------------------------------------------------------------

## 9.17.2026 - introduce study system:

**Working Directory**

`/gpfs1/home/w/r/wreganlo/eco_genomics_2026/transcriptomics`

**Input Files:**

`none`

**Output Files:**

`/gpfs1/home/w/r/wreganlo/eco_genomics_2026/transcriptomics/transcriptomics_notebook.md`

**Programs and dependencies:**

-   `R version 4.5.1`

-   `R-Studio`

**Scripts:**

`none`

**Code:**

```         
cd /gpfs1/cl/biol3990 #change directory
ll #look in directory
zcat #un cat something (open zip)
-n #number of lines
head #number of lines from the top you want to see
wc -l #number of lines in file (word count)
```

**Image:**

![](images/9.17.2026.terminal.screenshot.png)

**Notes / Observations:**

-   

**Next Steps:**

------------------------------------------------------------------------

## 9.22.2026 - Day 3 of Transcriptomics

**Working Directory**

`/gpfs1/home/w/r/wreganlo/eco_genomics_2026/transcriptomics`

**Input Files:**

`none`

**Output Files:**

`/gpfs1/home/w/r/wreganlo/eco_genomics_2026/transcriptomics/transcriptomics_notebook.md`

**Programs and dependencies:**

-   `R version 4.5.1`

-   `R-Studio`

**Scripts:**

`none`

**Code:**

```{r}
git status
git add .gitignore
git commit -m "Ignore results directory"
git push

touch ./mydata/test.txt
git status

#General Flow:
1. git pull
2. Do work
3. git status
4. git add .
5. git commit -m "meaningful message"
6. git push

#Diagnosing:
git status
git branch
git remote -v





```

**Notes:**

-   Today we set up our r working environment and copied some data.

-   We started looking at the data and exploring it with the commands mean and made a histogram to look at how the reads turned out (looking a frequency)

-   We then defined our model and started working with DESeq2

-   After that we transformed the data a bit and made a cluster tree to look for outliers

-   Then we made a PCA to visualize global gene expression patterns

    -   The PCAs look specifically at F0, F2, F4, and F11

    -   Arranged all togehter

    -   Saved them

**Photos:**

![](myresults/PCA_allGens.png)

------------------------------------------------------------------------
