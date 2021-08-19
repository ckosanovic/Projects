Methods in Bioinformatics Final Exam (Project \#2)
================

\#Christina Kosanovic (matr.923764) and Sampath Reddy (matr.937686).

\#Exam: July 22, 2019.

In this project, we worked on a human gene expression profiling RNA-Seq
data set composed of 36 samples from 6 tissues. We utilized specific
packages collected under Bioconductor, a dedicated repository/project:
which provides tools for the analysis of high-throughput genomic data.
Thus, the first order of business is to install these packages (edgeR,
limma, DESeq2) into R.

``` r
BiocManager::install(c("edgeR"))
```

    ## Bioconductor version 3.12 (BiocManager 1.30.16), R 4.0.3 (2020-10-10)

    ## Warning: package(s) not installed when version(s) same as current; use `force = TRUE` to
    ##   re-install: 'edgeR'

    ## Installation paths not writeable, unable to update packages
    ##   path: C:/Program Files/R/R-4.0.3/library
    ##   packages:
    ##     boot, class, cluster, codetools, foreign, KernSmooth, lattice, MASS,
    ##     Matrix, mgcv, nlme, nnet, spatial, survival

    ## Old packages: 'gert', 'later', 'RCurl', 'rtracklayer', 'xfun', 'XML'

``` r
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("limma")
```

    ## Bioconductor version 3.12 (BiocManager 1.30.16), R 4.0.3 (2020-10-10)

    ## Warning: package(s) not installed when version(s) same as current; use `force = TRUE` to
    ##   re-install: 'limma'

    ## Installation paths not writeable, unable to update packages
    ##   path: C:/Program Files/R/R-4.0.3/library
    ##   packages:
    ##     boot, class, cluster, codetools, foreign, KernSmooth, lattice, MASS,
    ##     Matrix, mgcv, nlme, nnet, spatial, survival
    ## Old packages: 'gert', 'later', 'RCurl', 'rtracklayer', 'xfun', 'XML'

``` r
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("DESeq2")
```

    ## Bioconductor version 3.12 (BiocManager 1.30.16), R 4.0.3 (2020-10-10)

    ## Warning: package(s) not installed when version(s) same as current; use `force = TRUE` to
    ##   re-install: 'DESeq2'

    ## Installation paths not writeable, unable to update packages
    ##   path: C:/Program Files/R/R-4.0.3/library
    ##   packages:
    ##     boot, class, cluster, codetools, foreign, KernSmooth, lattice, MASS,
    ##     Matrix, mgcv, nlme, nnet, spatial, survival
    ## Old packages: 'gert', 'later', 'RCurl', 'rtracklayer', 'xfun', 'XML'

Our datasets consist of the expression values (read counts) for each
sample: - 18829 (high quality) proteing coding genes and - 9438 (high
quality) non protein coding RNAs (6520 lincRNA, 1780 snRNAs and 1137
miRNAs)

The total dataset is comprised of 3 files (tables):

Counts.csv: gene expression values (read counts) 28266 human genes in
the 36 replicates

Annot.csv: annotation (gene symbol and class) for the 28266 genes

Design.csv: design of the RNAseq experiment (tissue, individual and sex
associated with each biological replicate)

We downloaded the files and imported them into our environment using the
read.table() function, subsequently setting them in our working
directory. We set parameters for row names, header and print separator.
By using the library() function we were able to access the edgeR package
to begin our data analysis.

``` r
data <- read.table(file = "./data.csv", row.names = 1, header = 1, sep = "\t")
conditions <- read.table(file = "./conditions.csv", row.names = 1, header = 1, sep = "\t")
annot <- read.table(file = "./annot.csv", row.names = 1, header = 1, sep = "\t")
library(edgeR)
```

    ## Loading required package: limma

Since our files have not been ordered, this is the first task in data
handling. We ordered the files (beginning with the “data.csv”" file) by
using the order function. Both columns and rows of each data set were
ordered in decreasing fashion. We opted to keep the original file names.

``` r
ordered_column <-order(colnames(data), decreasing = T)
ordered_row<-order(row.names(data), decreasing = T)
data<-data[ordered_row,ordered_column]

ordered_column <-order(colnames(annot), decreasing = T)
ordered_row<-order(row.names(annot), decreasing = T)
annot<-annot[ordered_row,ordered_column]

ordered_column <-order(colnames(conditions), decreasing = T)
ordered_row<-order(row.names(conditions), decreasing = T)
conditions<-conditions[ordered_row,ordered_column]
```

\#COMMON PART: (1) Read the data into a dgeList object

edgeR provides a platform for differential expression analysis of
RNA-seq expression profiles with biological replication. It stores data
in a simple list-based data object called a DGEList which can be
manipulated like any other list in R. Since a table of counts is already
available as a matrix or a data.frame then a DGEList object can be made
(using counts = data) and a grouping factor (using group = conditions)
can be added at the same time. The main elements of the DGEList object
are a matrix (counts), containing the integer counts and a data.frame
(samples) with information about the samples or libraries. The
data.frame samples will also have a corresponding column group to
identify the group membership of each sample.

``` r
data_DGEList<-DGEList(counts=data, group=conditions$tissue)
```

2.  Keep only genes that are likely to be expressed (i.e genes that have
    more than 10 reads in at least 2 replicates).

Since this is a differential expression analysis we want to get rid of
genes with very low counts across all libraries. Thus, prior to further
analysis, we must filter our data sets. Genes not expressed in all
samples for any of the conditions will be filtered out and thus we set
our cuttoff due to the suggestion of more than 10 reads in at least 2
replicates. Since edgeR prefers to use CPM rather than raw counts, we
have to be careful to select for counts as CPM will only select the
highly expressed genes. To accomplish this, we started with a logical
test to assign the counts in our data set greater than 10 to a new
matrix we named “cond”. Then we were able to select for all of the genes
that met this new logical condition, denoted as “1” for TRUE. By using
the apply() function, we were able to aggregate the sums of the TRUE
conditions from each row. Then, by using the cbind() function we were
able to add these results to the new matrix as a new column. Finally, we
could add our second condition (at least 2 replicates) by creating a
logical test for the rows in our new matrix with a value of 2 or greater
against the results of our previous logical test, found in the last
column (column 37). The results of the logical tests were renamed into a
subset of filtered data “exp” to denote expressed genes. After
filtering, the result is 19077 obs.(out of 28266) of the 36 variables.

``` r
cond<-data>10 #filtering data to the suggested more than 10 reads
cond<-cbind(cond, apply(cond,1,sum)) #selecting for genes that meet this condition
exp<-cond[,37]>=2 #selecting for second condition of at least 2 replicates
data_exp<-data[exp,] #aggregating the results of the logical test into a new subset
dim(data_exp) #returns dimensions of the data set
```

    ## [1] 19077    36

``` r
data_DGEList<-DGEList(counts=data_exp, group=conditions$tissue) #modified DGEList using expressed genes as the counts.
```

3.  Perform normalization (both for library size and composition)

The next crucial step in data analysis is performing normalization of
the data in order to make it comparable. The best practices for RNA-seq
data processing requires normalization for both library size and
composition. We can utlize the function, Calculate Normalization Factors
(calcNormFactors) to scale the raw library sizes. calcNormFactors uses
TMM normalisation for the actual differential expression analysis. The
main aim in TMM normalization is to account for variation in library
size between samples of interest, accounting for the fact that some
extremely differentially expressed genes would negatively impact the
normalization procedure. TMM (Trimmed Mean of M-values) is calculated
where a trimmed mean is the average after removing the upper and lower %
of the data (removing most/lowest expressed genes and genes with
highest/lowest log ratios from a gene set). It is an assumption of TMM
that the majority of the genes are not differentially expressed.

``` r
data_DGEList<-calcNormFactors(data_DGEList) #calculate normalization factors
print(data_DGEList)
```

    ## An object of class "DGEList"
    ## $counts
    ##                 GTEX.NPJ7.0011.R8a.SM.2I3G2 GTEX.NPJ7.0011.R5a.SM.33HBK
    ## ENSG00000273493                           2                           6
    ## ENSG00000273492                         126                         190
    ## ENSG00000273487                          72                         398
    ## ENSG00000273481                           4                          16
    ## ENSG00000273476                           2                           6
    ##                 GTEX.NPJ7.0011.R4a.SM.2I3GJ GTEX.NPJ7.0011.R2a.SM.2I3GF
    ## ENSG00000273493                           6                           2
    ## ENSG00000273492                         231                          83
    ## ENSG00000273487                         246                         270
    ## ENSG00000273481                           4                          14
    ## ENSG00000273476                          10                           3
    ##                 GTEX.NPJ7.0011.R1a.SM.3GACT GTEX.NPJ7.0011.R10A.SM.2I3E5
    ## ENSG00000273493                           5                            2
    ## ENSG00000273492                         152                          194
    ## ENSG00000273487                         155                          292
    ## ENSG00000273481                           2                            5
    ## ENSG00000273476                          10                            5
    ##                 GTEX.N7MT.0011.R8a.SM.2I5GU GTEX.N7MT.0011.R5a.SM.2I3G6
    ## ENSG00000273493                           6                           8
    ## ENSG00000273492                         328                         608
    ## ENSG00000273487                         263                         286
    ## ENSG00000273481                          10                          16
    ## ENSG00000273476                          13                          29
    ##                 GTEX.N7MT.0011.R4a.SM.2I3G9 GTEX.N7MT.0011.R2a.SM.2I3GI
    ## ENSG00000273493                           3                           6
    ## ENSG00000273492                         292                          80
    ## ENSG00000273487                         181                         246
    ## ENSG00000273481                           2                           8
    ## ENSG00000273476                           7                          12
    ##                 GTEX.N7MT.0011.R1a.SM.5SI7S GTEX.N7MT.0011.R10A.SM.2I3E1
    ## ENSG00000273493                           2                            2
    ## ENSG00000273492                         315                          526
    ## ENSG00000273487                         217                          350
    ## ENSG00000273481                           9                            9
    ## ENSG00000273476                           3                            8
    ##                 GTEX.15CHQ.0011.R8b.SM.6EU32 GTEX.15CHQ.0011.R5b.SM.6AJAN
    ## ENSG00000273493                            1                            2
    ## ENSG00000273492                          262                          374
    ## ENSG00000273487                          189                          227
    ## ENSG00000273481                            0                            0
    ## ENSG00000273476                            4                            8
    ##                 GTEX.15CHQ.0011.R4a.SM.686ZX GTEX.15CHQ.0011.R2b.SM.6AJAD
    ## ENSG00000273493                            4                            0
    ## ENSG00000273492                          260                          126
    ## ENSG00000273487                          270                          206
    ## ENSG00000273481                           16                           12
    ## ENSG00000273476                            7                            3
    ##                 GTEX.15CHQ.0011.R1a.SM.6EU31 GTEX.15CHQ.0011.R10b.SM.6AJBV
    ## ENSG00000273493                           10                            14
    ## ENSG00000273492                          145                           343
    ## ENSG00000273487                          238                           301
    ## ENSG00000273481                            6                             2
    ## ENSG00000273476                            4                             9
    ##                 GTEX.145LS.0011.R8b.SM.5PNXA GTEX.145LS.0011.R5a.SM.5SI65
    ## ENSG00000273493                            6                            3
    ## ENSG00000273492                          169                          223
    ## ENSG00000273487                          287                          327
    ## ENSG00000273481                           27                            7
    ## ENSG00000273476                            6                           14
    ##                 GTEX.145LS.0011.R4b.SM.5S2UT GTEX.145LS.0011.R2a.SM.5PNZI
    ## ENSG00000273493                            8                            4
    ## ENSG00000273492                          159                          121
    ## ENSG00000273487                          202                          253
    ## ENSG00000273481                            6                            8
    ## ENSG00000273476                           17                            8
    ##                 GTEX.145LS.0011.R1b.SM.5PNUP GTEX.145LS.0011.R10a.SM.5PNUQ
    ## ENSG00000273493                            3                             7
    ## ENSG00000273492                           86                           268
    ## ENSG00000273487                          196                           415
    ## ENSG00000273481                            7                             0
    ## ENSG00000273476                           17                            14
    ##                 GTEX.13X6I.0011.R8a.SM.5PNZF GTEX.13X6I.0011.R5a.SM.5PNWW
    ## ENSG00000273493                            3                            1
    ## ENSG00000273492                          324                          322
    ## ENSG00000273487                          269                          256
    ## ENSG00000273481                           21                           13
    ## ENSG00000273476                            5                           28
    ##                 GTEX.13X6I.0011.R4b.SM.5PNU9 GTEX.13X6I.0011.R2b.SM.5PNWQ
    ## ENSG00000273493                            5                            1
    ## ENSG00000273492                          303                          158
    ## ENSG00000273487                          245                          226
    ## ENSG00000273481                            6                           10
    ## ENSG00000273476                            6                            2
    ##                 GTEX.13X6I.0011.R1b.SM.5PNZC GTEX.13X6I.0011.R10a.SM.5PNWI
    ## ENSG00000273493                            4                             3
    ## ENSG00000273492                          230                           447
    ## ENSG00000273487                          263                           353
    ## ENSG00000273481                           12                             7
    ## ENSG00000273476                           13                            13
    ##                 GTEX.13RTJ.0011.R8b.SM.5O9DL GTEX.13RTJ.0011.R5a.SM.5P9HR
    ## ENSG00000273493                           10                            9
    ## ENSG00000273492                          288                          230
    ## ENSG00000273487                          149                          498
    ## ENSG00000273481                           11                           14
    ## ENSG00000273476                            4                           35
    ##                 GTEX.13RTJ.0011.R4b.SM.5PNX1 GTEX.13RTJ.0011.R2a.SM.5PNW9
    ## ENSG00000273493                            7                            5
    ## ENSG00000273492                          210                          112
    ## ENSG00000273487                          292                          159
    ## ENSG00000273481                           24                            8
    ## ENSG00000273476                            7                            3
    ##                 GTEX.13RTJ.0011.R1a.SM.5O9D9 GTEX.13RTJ.0011.R10b.SM.5O9CW
    ## ENSG00000273493                            8                            12
    ## ENSG00000273492                           95                           201
    ## ENSG00000273487                          460                           424
    ## ENSG00000273481                           30                             4
    ## ENSG00000273476                           16                             6
    ## 19072 more rows ...
    ## 
    ## $samples
    ##                                               group lib.size norm.factors
    ## GTEX.NPJ7.0011.R8a.SM.2I3G2            Hypothalamus 34188351    1.0298659
    ## GTEX.NPJ7.0011.R5a.SM.33HBK Caudate (basal ganglia) 61484598    1.0069764
    ## GTEX.NPJ7.0011.R4a.SM.2I3GJ                Amygdala 38818445    0.9710285
    ## GTEX.NPJ7.0011.R2a.SM.2I3GF        Substantia nigra 46030528    0.8550101
    ## GTEX.NPJ7.0011.R1a.SM.3GACT             Hippocampus 32146447    1.0441424
    ## 31 more rows ...

4.  Perform a MDS (\~PCA) plot of the data

Given that this is an RNA-Seq experiment, the aim is to perform a
differential expression (DE). To perform DE with a package such as
edgeR, we can draw an MDS plot using the plotMDS() function (from the
limma package). This plot is a ‘summary’ of the data attained through a
process of reduction that can transform the large number of variables
into a lesser number that are uncorrelated (i.e. the ‘principal
components’). PCA finds the principal components of data. In this way,
data is measured on its principal components rather than on a normal x-y
axis. The principal components are essentially the underlying structure
in the data. They are the directions where there is the most
variance/where the data is most spread out. In summary, we are able to
visualize the data thanks to linear dimension reduction. From plotMDS,
results a scatterplot where the tissues with similar gene expression
profiles are clustered closer together while those with different gene
expression profiles will be located further apart.

``` r
plotMDS(data_DGEList,col=rep(c("red","green","blue","yellow","black","purple")),table(levels(conditions$tissue)),labels = conditions$tissue, cex = 0.5, main = "Gene Distribution Profiles of Tissues") #selecting cex = 0.5 for the numeric vector of plot symbol expansions.
```

![](Methods_BioInf_final_report_Kosanovic_Reddy_files/figure-gfm/unnamed-chunk-7-1.png)<!-- -->

5.  Select 2 different (and meaningful) biological conditions and
    perform a differential expression analysis using the exactTest()
    function.

Based on our MDS plot we can use what we understand about PCA to select
genes that are differently expressed according to the plot. Therefore,
based on their positions on the plot, we selected the tissues “Amygdala”
and “Caudate (basal ganglia)”. We utilize the exactTest function to
compute genewise exact tests for differences in the means between two
groups of negative-binomially distributed counts. A negative binomial
distribution is required for RNA-seq to perform DE as the standard
binomial distribution does not effectively take into account the extra
degree of biological variance inherent in RNA-seq experiments. We must
first use the estimateDisp() function and plotBCV() function to fit the
model and calculate the dispersion. We have the option to estimate
different types of dispersion, in this case we use tagwise. The tagwise
dispersion method calculates gene-specific dispersions as opposed to a
common “one-size-fits-all” dispersion value for all genes.

From here we can see if the dispersion of the data is consistent with
the positions on the MDS plot. Ultimately, we are calculating a Fisher’s
exact test which is a statistical test used to determine if there are
nonrandom associations between two categorical variables. We generate a
table of values ordered by p-values and FDR (false discover rates).
Within this table of ordered FDR values, those genes that fall below the
FDR cut-off are considered DE (differentially expressed). We input this
data into deTa and set the dispersion parameter as “tagwise”.

``` r
disp_data<-estimateDisp(data_DGEList)
```

    ## Using classic mode.

``` r
disp_plot<-plotBCV(disp_data) #Plot the genewise biological coefficient of variation (BCV)
```

![](Methods_BioInf_final_report_Kosanovic_Reddy_files/figure-gfm/unnamed-chunk-8-1.png)<!-- -->

``` r
exactTest(disp_data)
```

    ## An object of class "DGEExact"
    ## $table
    ##                       logFC    logCPM    PValue
    ## ENSG00000273493 -0.58725138 -2.734636 0.3034764
    ## ENSG00000273492  0.01442936  2.401448 0.9719018
    ## ENSG00000273487  0.07710504  2.550377 0.7287672
    ## ENSG00000273481 -0.10855502 -2.019798 0.8923620
    ## ENSG00000273476  0.70423742 -1.962671 0.1663791
    ## 19072 more rows ...
    ## 
    ## $comparison
    ## [1] "Amygdala"                "Caudate (basal ganglia)"
    ## 
    ## $genes
    ## NULL

``` r
deTa<-exactTest(disp_data, pair=c("Amygdala", "Caudate (basal ganglia)"), dispersion = "tagwise")
Ta=topTags(deTa,n=nrow(deTa),p.value = 0.01)
```

6.  Create a topTags type of edgeR object containing the list of
    differentially expressed genes (DEGs). DEGs should have a FDR &lt;=
    0.01

We implement the topTags function which gives us a table of the Top
Differentially Expressed Genes/Tags. It extracts the most differentially
expressed genes (or sequence tags) from a test object, ranked either by
p-value or by absolute log-fold-change. We set the p-value parameter
equal to 0.01, to select only those lower or equal to the cut-off value.
Using this function, we find that 1609 genes among a total of 28266 were
differentially expressed.

``` r
Diff.Exp.Gen<- topTags(deTa,nrow(deTa), p.value = 0.01) #returns most differentially expressed genes
nrow(Diff.Exp.Gen) #returns the number of rows of the matrix
```

    ## [1] 1609

7.  Assign all the genes into one of the 4 possible classes: DE\_UP
    (FDR&lt;=0.01 and logFC&gt;0), DE\_DOWN (FDR&lt;=0.01 and
    logFC&lt;0), notDE\_UP (FDR&lt;=0.01 and logFC&gt;0), notDE\_DOWN
    (FDR&lt;=0.01 and logFC&lt;0), and then do a boxplot of the logFC of
    the genes belonging to each class.

To accurately display the distribution of the expression levels of our
genes into 4 different categories, we must perform a differential
expression test for all of these “classes”. In this case, we perform the
same test as in the previous step, utilizing the topTags() function,
however we will not indicate a p-value threshold. Now we can access all
of the genes, not just those that pass the statistical test. To keep it
simple, we will names this category “genes”.

In order to classify these genes into their respective class (DE\_UP,
DE\_DOWN, notDE\_UP, notDE\_DOWN), we must use a method to separate
genes based on FDR and logFC values. Using the accessor function ($), we
can retrieve the tabulated list within our “genes” data. We input this
into “genes\_table”.

From this data, we can categorize it using a series of ifelse()
functions. In R, an if-else statement tells the program to run one block
of code if the conditional statement is TRUE, and a different block of
code if it is FALSE. This effectively arranges the four different
classes based on parameters we set, FDR and logFC values. We called the
resulting categorical data “classes”. To display this in an organized
fashion, we organized the results into a data frame, which contained the
class and logFC for each gene. This we named “results”.

``` r
genes<- topTags(deTa,nrow(deTa)) #test for DE genes without the p-value parameter
genes_table<- genes$table #$ accessor to retrieve gene table list
classes<-paste(ifelse(genes_table$FDR<=0.01,"DE","notDE"),ifelse(genes_table$logFC>0,"UP","DOWN"),sep="_") #arranging the classes by FDR and logFC
results <-data.frame(genes_table$logFC,classes,row.names=row.names(genes_table))
table(results$classes) #to return the values for each class
```

    ## 
    ##    DE_DOWN      DE_UP notDE_DOWN   notDE_UP 
    ##        717        892       8713       8755

Box plot of the logFC of the genes belonging to each class (DE\_UP,
DE\_DOWN, notDE\_UP, notDE\_DOWN).

In R, we can use the boxplot() function to produce a box-and-whisker
plot of the given (grouped) values. The boxplot() function takes in any
number of numeric vectors, drawing a boxplot for each vector. It can
also pass in a list (or data frame) with numeric vectors as its
components.

To use this function, we must input a formula containing logFC and
classes of our genes (logFC\~results$classes). We can make use of
various arguments/tests (width, notch, col) to alter parameters/display
of the boxplot. Now we have a boxplot to visualize the distribution.

``` r
boxplot(results$genes_table.logFC~results$classes,outline=F,varwidth=T,notch=T,col=c("red","green","blue", "yellow"),las=1,cex.axis=1, main="Differentially expressed genes") #Construct plotting parameters: Width of the plot window in inches (otherwise defaults to 4.5), logical test for notch (where we can add a notch with notch = TRUE, col to define color), Specify the size of the tick label numbers/text with cex.axis length = 1.
```

![](Methods_BioInf_final_report_Kosanovic_Reddy_files/figure-gfm/unnamed-chunk-11-1.png)<!-- -->
Results from Part 1:

We filtered the counts to get 19077 obs.(out of 28266). From our gene
set we found 1609 genes to be differentially expressed. We then
classified the genes into their respective classes (DE\_UP, DE\_DOWN,
notDE\_UP, notDE\_DOWN), by separating genes based on FDR and logFC
values. Of the genes from the filtered counts:

DE\_UP (differentially expressed and upregulated):892 genes

DE\_DOWN (differentially expressed and downregulated):717 genes

notDE\_UP (not differentially expressed but upregulated):8755 genes

notDE\_DOWN (not differentially expressed but down regulated): 8713
genes

\#PART 2 (Project \#2)

Identifying genes showing sex specific expression in the 2 tissues
already considered for the “common part”.

First, we filter out the two tissues separately which are selected above
in the common part and we cohere the two tissues using rbind function

``` r
filt_Amygdala <-conditions[conditions$tissue== "Amygdala",]
filt_Caudate <-conditions[conditions$tissue== "Caudate (basal ganglia)",]
filt_tissues <-rbind(filt_Amygdala, filt_Caudate)
head(filt_tissues)
```

    ##                                tissue sex         id
    ## GTEX.NPJ7.0011.R4a.SM.2I3GJ  Amygdala   2  GTEX.NPJ7
    ## GTEX.N7MT.0011.R4a.SM.2I3G9  Amygdala   2  GTEX.N7MT
    ## GTEX.15CHQ.0011.R4a.SM.686ZX Amygdala   1 GTEX.15CHQ
    ## GTEX.145LS.0011.R4b.SM.5S2UT Amygdala   2 GTEX.145LS
    ## GTEX.13X6I.0011.R4b.SM.5PNU9 Amygdala   1 GTEX.13X6I
    ## GTEX.13RTJ.0011.R4b.SM.5PNX1 Amygdala   1 GTEX.13RTJ

Then, we keep only columns related to the filtered tissues i.e.,
Amygdala and caudate

``` r
rownames_cond_tissue<-rownames(filt_tissues)
data_tissue<-data[,rownames_cond_tissue]
rownames_cond_tissue
```

    ##  [1] "GTEX.NPJ7.0011.R4a.SM.2I3GJ"  "GTEX.N7MT.0011.R4a.SM.2I3G9" 
    ##  [3] "GTEX.15CHQ.0011.R4a.SM.686ZX" "GTEX.145LS.0011.R4b.SM.5S2UT"
    ##  [5] "GTEX.13X6I.0011.R4b.SM.5PNU9" "GTEX.13RTJ.0011.R4b.SM.5PNX1"
    ##  [7] "GTEX.NPJ7.0011.R5a.SM.33HBK"  "GTEX.N7MT.0011.R5a.SM.2I3G6" 
    ##  [9] "GTEX.15CHQ.0011.R5b.SM.6AJAN" "GTEX.145LS.0011.R5a.SM.5SI65"
    ## [11] "GTEX.13X6I.0011.R5a.SM.5PNWW" "GTEX.13RTJ.0011.R5a.SM.5P9HR"

``` r
dim(data_tissue) #we get only 12 columns out of 36 after filtering only the columns related to the tissue
```

    ## [1] 28266    12

Then, we order both the tissue data frames with respect to sex and we
also sort the row names of conditions dataframe

``` r
cond_tissues_order<-order(filt_tissues$sex)
cond_tissues_order_sex<-filt_tissues[cond_tissues_order,]
order_rownames_cond_tissues_order_sex<-row.names(cond_tissues_order_sex)
data_ord<-data[,order_rownames_cond_tissues_order_sex]
```

After ordering, we perform a PCA analysis to see the seperation between
the replicates of different sexes

``` r
Expr_genes_tissue1_pseud<-data_ord+0.01 
Log_expr_genes_tissue1_pseud<-log(Expr_genes_tissue1_pseud)
PCA1<-prcomp(t(Log_expr_genes_tissue1_pseud))
plot(PCA1$x, main="PCA of RNAseq data in Caudate and Amygdala, sex separation highlighted", col= rep(c("red", "blue"), rep(6,2)), pch=20, cex=1)
text(PCA1$x, labels = cond_tissues_order_sex$tissue, offset=1, cex=0.5) #in the plot we can see that both the tissues are displayed by different sexes, the red & blue color represents different sexes
```

![](Methods_BioInf_final_report_Kosanovic_Reddy_files/figure-gfm/unnamed-chunk-15-1.png)<!-- -->

We filtered the two tissues so that each condition has at least 10 reads
and 2 replicates for each tissue Again, we filtered the data frames
where we get only the information about one tissue & we do it seperately
for both the tissues

``` r
data_ord_Amygdala <-data_tissue[,1:6]
data_ord_Amygdala_logical<-data_ord_Amygdala>10
data_ord_Amygdala<-data_ord_Amygdala[(apply(data_ord_Amygdala_logical, 1,sum)>=2),]
dim(data_ord_Amygdala) # 17433 genes are likely to be expressed in Amygdala tissue
```

    ## [1] 17433     6

``` r
data_ord_caudate <-data_tissue[,7:12]
data_ord_caudate_logical<-data_ord_caudate >10
data_ord_caudate<-data_ord_caudate[(apply(data_ord_caudate_logical, 1,sum)>=2),]
dim(data_ord_caudate) # 17797 genes are expressed in Caudate tissue 
```

    ## [1] 17797     6

Then we perform the Differential Expression (DE) analysis of genes
according to sex specificity in Amygdala

``` r
Diff.Exp.Gen_Amygdala <-filt_tissues[1:6,]
head(Diff.Exp.Gen_Amygdala)
```

    ##                                tissue sex         id
    ## GTEX.NPJ7.0011.R4a.SM.2I3GJ  Amygdala   2  GTEX.NPJ7
    ## GTEX.N7MT.0011.R4a.SM.2I3G9  Amygdala   2  GTEX.N7MT
    ## GTEX.15CHQ.0011.R4a.SM.686ZX Amygdala   1 GTEX.15CHQ
    ## GTEX.145LS.0011.R4b.SM.5S2UT Amygdala   2 GTEX.145LS
    ## GTEX.13X6I.0011.R4b.SM.5PNU9 Amygdala   1 GTEX.13X6I
    ## GTEX.13RTJ.0011.R4b.SM.5PNX1 Amygdala   1 GTEX.13RTJ

``` r
s_Amygdala<-DGEList(counts = data_ord_Amygdala , group = Diff.Exp.Gen_Amygdala$sex)
s_Amygdala$samples
```

    ##                              group lib.size norm.factors
    ## GTEX.NPJ7.0011.R4a.SM.2I3GJ      2 38812155            1
    ## GTEX.N7MT.0011.R4a.SM.2I3G9      2 30307465            1
    ## GTEX.15CHQ.0011.R4a.SM.686ZX     1 44825349            1
    ## GTEX.145LS.0011.R4b.SM.5S2UT     2 37619120            1
    ## GTEX.13X6I.0011.R4b.SM.5PNU9     1 37302237            1
    ## GTEX.13RTJ.0011.R4b.SM.5PNX1     1 55562421            1

``` r
s_Amygdala<-calcNormFactors(s_Amygdala)
s_Amygdala<-estimateDisp(s_Amygdala)
```

    ## Using classic mode.

``` r
deTa_s_Amygdala<-exactTest(s_Amygdala,dispersion = "tagwise") 
Ta_s_Amygdala<-topTags(deTa_s_Amygdala, n=nrow(deTa_s_Amygdala), p.value = 0.05) #considering genes which have an FDR < 0.5
dim(Ta_s_Amygdala) #there are 17 sex specific genes in Amygdala
```

    ## [1] 17  4

Following this, we perform a differential expression analysis of genes
with respect to sex in Caudate tissue

``` r
Diff.Exp.Gen_caudate <-filt_tissues[7:12,]
head(Diff.Exp.Gen_caudate)
```

    ##                                               tissue sex         id
    ## GTEX.NPJ7.0011.R5a.SM.33HBK  Caudate (basal ganglia)   2  GTEX.NPJ7
    ## GTEX.N7MT.0011.R5a.SM.2I3G6  Caudate (basal ganglia)   2  GTEX.N7MT
    ## GTEX.15CHQ.0011.R5b.SM.6AJAN Caudate (basal ganglia)   1 GTEX.15CHQ
    ## GTEX.145LS.0011.R5a.SM.5SI65 Caudate (basal ganglia)   2 GTEX.145LS
    ## GTEX.13X6I.0011.R5a.SM.5PNWW Caudate (basal ganglia)   1 GTEX.13X6I
    ## GTEX.13RTJ.0011.R5a.SM.5P9HR Caudate (basal ganglia)   1 GTEX.13RTJ

``` r
s_caudate<-DGEList(counts = data_ord_caudate , group = Diff.Exp.Gen_caudate$sex) #considering genes which have an FDR < 0.5
s_caudate$samples
```

    ##                              group lib.size norm.factors
    ## GTEX.NPJ7.0011.R5a.SM.33HBK      2 61478013            1
    ## GTEX.N7MT.0011.R5a.SM.2I3G6      2 53774997            1
    ## GTEX.15CHQ.0011.R5b.SM.6AJAN     1 41571691            1
    ## GTEX.145LS.0011.R5a.SM.5SI65     2 46175212            1
    ## GTEX.13X6I.0011.R5a.SM.5PNWW     1 50218333            1
    ## GTEX.13RTJ.0011.R5a.SM.5P9HR     1 63827675            1

``` r
s_caudate<-calcNormFactors(s_caudate)
s_caudate<-estimateDisp(s_caudate)
```

    ## Using classic mode.

``` r
deTa_s_caudate<-exactTest(s_caudate,dispersion = "tagwise")
Ta_s_caudate<-topTags(deTa_s_caudate, n=nrow(deTa_s_caudate), p.value = 0.05)
dim(Ta_s_caudate) # 23 genes are sex specific in Caudate
```

    ## [1] 23  4

``` r
library(VennDiagram) #command to load venn diagram package
```

    ## Warning: package 'VennDiagram' was built under R version 4.0.5

    ## Loading required package: grid

    ## Loading required package: futile.logger

Then plot a Venn diagram which shows the genes that are sex specific in
both the tissues, create a list accordingly.

``` r
DE_1<-list(rownames(Ta_s_Amygdala), rownames(Ta_s_caudate))
venn.diagram(DE_1, category.names = c("Amygdala", "Caudate"), filename="Venn Diagram 1", main = "Venn Diagram: sex specific DEs between Amygdala and Caudate", col=c("red", "green"), imagetype = "png" )
```

    ## [1] 1

15 genes are found to be sex specific between the tissues Amygdala and
Caudate.

Then, we plot another Venn diagram to show the genes which have sex
specific expression in the tissues

``` r
Sex_Spfc_DE<-c(rownames(Ta_s_Amygdala), rownames(Ta_s_caudate))
DE_2 <-list(Sex_Spfc_DE, rownames(Ta))
venn.diagram(DE_2, category.names = c("Sex_Spfc_DE in Amygdala and/or caudate", "DE between Amygdala and caudate"), filename="Venn Diagram_2", main = "Venn Diagram: DE_genes and Sex_Spfc_DE_genes considering Amygdala and caudate", col=c("blue", "pink"), imagetype = "png" )
```

    ## [1] 1

Conclusions:

The final Venn diagram shows us that there are no genes which are
differentially expressed between two tissues and between two sexes. We
can conclude that only 15 genes have sex specificity in Amygdala and
Caudate. Further, there are no genes which show sex specific expression
between the two tissues.
