# TSLsummerschool2017
A tutorial on genetics that integrates visualization with ggplot2

## Introduction
The goal of this practical is to gain an understanding of how to use natural variation to map a trait of interest. Whether mapping a gene or quantitative trait locus (QTL), the approaches used are identical. We will be using a diverse array of software and scripts that aid genetic analysis, including:
1. R
  * qtl
  * ggplot2
2. QTL Cartographer

## R/qtl
This tutorial is composed of two data sets:
* Oregon Wolfe Barley (OWB) doubled-haploid population
* Steptoe x Morex (SxM) doubled-haploid population

Data is contributed by Matthew Moscou (The Sainsbury Laboratory; unpublished) and Rients Niks (Wageningen University & Research; unpublished).

### Oregon Wolfe Barley doubled-haploid population
The Oregon Wolfe Barley doubled haploid population was created by Dr. Bob Wolfe (Oregon State University) using two barley accessions that had contrasting phenotypes for a range of morphological characteristics: OWB Dominant and OWB Recessive. The names are based on the inheritance of morphological traits that are most commonly associated with each parent. In addition to morphological variation, the population also segregates for a range of disease traits including:
* Barley powdery mildew (*Blumeria graminis* f. sp. *hordei*; *Bgh*)
* Oat crown rust (*Puccinia coronata*; *Pc*)
* Barley leaf rust (*Puccinia hordei*; *Phor*)
* False barley leaf rust (*Puccinia hordei-murini*; *Ph-m*)
* Meadow barley leaf rust (*Puccinia hordei-secalini*; *Ph-s*)
* Wheat leaf rust (alternate host: Ranunculaceae) (*Puccinia persistens* subsp. *triticina*; *Pper*)
* Wheat leaf rust (alternate host: Boraginacec) (*Puccinia triticina*; *Ptri*)

These pathogen range from host pathogens to intermediate host pathogens (Bettgenhaeuser *et al.* (2014)), therefore different phenotypes have been collected based on the pathogen used. These phenotypes include:
* Visible infection sites (VIS)
* Infection frequency (IF)
* Relative latency period (RLP)
* Genotype effect (GE)
* Flecks number (FN)


#### Genetic map
The first step is to open RStudio and load R/qtl. To learn more about R/qtl and all it has to offer, go here [http://www.rqtl.org/](http://www.rqtl.org/). The manual is extremely comprehensive, but below we have a range of examples demonstrating its versatility and use.

```R
library(qtl)
```

Next, we need to import the genetic map, marker information, phenotypic data, chromosome identifiers, marker names, and phenotype names. This information is stored in two different files: a CRO file and MAP file. For more information on the format of these files, see the QTL Cartographer manual [http://statgen.ncsu.edu/qtlcart/#manual](http://statgen.ncsu.edu/qtlcart/#manual). These files have been generated previously so that we can move straight into QTL analysis. See the Advanced User section below for information on how to build these files yourself.

```R
pop.data = read.cross(format="qtlcart", file="qtlcart_OWB.cro", mapfile="qtlcart_OWB.map")
```

The power of interval mapping and composite interval mapping is the use of positions between markers for QTL analysis. To do this, we need to estimate the probability of a genotype at a specific position between markers. Here, we specify a step distance of 1.0 cM and use the Kosambi mapping function (standard distance metric in genetics).

```R
pop.data = calc.genoprob(pop.data, step=1, map.function="kosambi")
```

R/qtl has implemented a quick approach to visualize genetic maps, distributions of phenotypic traits, and the distribution of missing data.

```R
plot(pop.data)
```

This includes visualizing alleles within the ordered genetic map. This figure is not extremely useful for assessing marker x trait association, but it can be used to quickly assess if markers are in a correct order.

```R
geno.image(pop.data)
```

One of the most important tools is `plot.rf`, which we perform a recombination fraction plot on a genetic map. This can be used to assess the quality of the constructed genetic map. Generally, markers should only have linkage with markers nearest to them. There are reasons why this would not be the case, such as segregation distortion or multilocus prezygotic and postzygotic barriers.

```R
plot.rf(pop.data, what="rf")
```

#### Phenotypes
Within R, complex data sets are stored in data frames. Knowing how to access these can be extremely useful at understanding how to visualize data using R. To visualize the types of data within the R/qtl data frame, use the following command:

```R
attributes(pop.data)
```

Now, we can access the names of the phenotypes within this data frame with the following command:

```R
pop.data$pheno[1,]
```

R/qtl has a very focused set of plotting functions. In contrast, ggplot2 is a powerful set of visualization tools that take advantage of the data frame structure in R. Next, we load ggplot2 and use it to visualize histograms for different traits.

```R
library(ggplot2)

ggplot(pop.data$pheno, aes(Bgh5874)) + geom_histogram()
```

This same line of code can be used to visualize the distribution of traits. There are several diverse distribution including bimodal (Bgh5874), skewed (Ph_sIF), and normal (Phor1.2.1RLP).

```R
ggplot(pop.data$pheno, aes(PtriVIS)) + geom_histogram()
ggplot(pop.data$pheno, aes(PtriIF)) + geom_histogram()
ggplot(pop.data$pheno, aes(Ph_mVIS)) + geom_histogram()
ggplot(pop.data$pheno, aes(Ph_mIF)) + geom_histogram()
ggplot(pop.data$pheno, aes(Ph_sVIS)) + geom_histogram()
ggplot(pop.data$pheno, aes(Ph_sIF)) + geom_histogram()
ggplot(pop.data$pheno, aes(PperVIS)) + geom_histogram()
ggplot(pop.data$pheno, aes(PperIF)) + geom_histogram()
ggplot(pop.data$pheno, aes(Phor1.2.1RLP)) + geom_histogram()
ggplot(pop.data$pheno, aes(Phor1.2.1GE)) + geom_histogram()
ggplot(pop.data$pheno, aes(PcFN)) + geom_histogram()
ggplot(pop.data$pheno, aes(BghIFS)) + geom_histogram()
```

For an extended resource in using ggplot, please see Dan MacLean's `ggplotbook` tutorial [https://danmaclean.github.io/ggplotbook/](https://danmaclean.github.io/ggplotbook/).

When performing QTL analysis, it is critical to understand the distribution of your data. The type of statistical model may need to be changed, or alternatively, how you interpret your results may change. Generally, you need to be able to answer the following:
* Qualitative or quantitative or a mixture
* Linear or non-linear metric
* Normal or non-normal distribution
* Degree of missing data

#### QTL analysis
To start, we will perform interval mapping using the barley powdery mildew data set that had a bimodal distribution, which is suggestive of a single gene controlling resistance. While the data is not normally distributed, we can still use it to map the *R* gene locus.

```R
trait1.im = scanone(pop.data, pheno.col=1, model="normal", method="em")
plot(trait1.im)
summary(trait1.im)
```

A more appropriate model to select in this case would be non-parametric. We can change this in the model parameter and rerun the QTL analysis.

```R
trait1.im = scanone(pop.data, pheno.col=1, model="np", method="em")
plot(trait1.im)
```

Next, we need to perform bootstrap analysis to identify linkage that is significant relative to what could happen by change. The standard approach for this is to randomize the data set, perform QTL analysis, and identify the strongest linkage. If we do this 1,000 times and take the 95th quantile, we can control at an equivalent p value of 0.05. 

```R
trait1.im.perm = scanone(pop.data, pheno.col=1, model="np", method="em", n.perm=1000)
summary(trait1.im, perms=trait1.im.perm, alpha=0.05, pvalues=TRUE)

plot(trait1.im)
add.threshold(trait1.im, perms=trait1.im.perm, alpha=0.05, col="red")
```

One of the most powerful approaches at QTL analysis is composite interval mapping. This takes into account major effect loci that contribute to a phenotype and incorporates these into the regression analysis for identifying significant QTL. The main parameter to select is the number of covariate loci (n.marcovar). This can generally be set at 5 loci without need for change.

```R
trait1.cim = cim(pop.data, pheno.col=1, n.marcovar=5, method="em", map.function="kosambi")
plot(trait1.cim)
summary(trait1.cim)
```

### Steptoe x Morex doubled-haploid population
A second data set involves the Steptoe x Morex doubled-haploid population of barley. The population was inoculated with *Blumeria graminis* f. sp. *hordei* isolate CC148. Two phenotypes were collect on the population:
* Infection, which is a measure of mildew susceptibility (scale 0 to 4)
* Necrosis, mildew induces necrosis in this population, which was scored on a scale of 0 to 4 

The goal is to adapt the scripts above to identify the QTL that control these phenotypes and assess the relationship between the two phenotypes. To start, you can load the population with the following command:

```R
pop.data = read.cross(format="qtlcart", file="qtlcart_SxM.cro", mapfile="qtlcart_SxM.map")
```

## Advanced users
Within this tutorial, we have provided pre-generated CRO and MAP files. These are generated using [QTL Cartographer](http://statgen.ncsu.edu/qtlcart/), which is one of the most commonly used software for QTL analysis (Basten et al. (1994)). Our group actively uses QTL Cartographer and has written several scripts that facilitate QTL analysis using QTL Cartographer. These scripts include:
* QKcartographer_preprocess.py, converts tabular flat text files into formatted input files for QTL Cartographer.
* QKcartographer_permutations.py, identifies experiment-wise thresholds based on permuted data sets from QTL Cartographer.
* QKcartographer_visualization.py, generates figures in PNG or postscript format for publication purposes.
* QKcartographer_epistasis.py, parses QTL Cartographer output files and permits curation of significant QTLs. Optional command to generate scripts for epistasis analysis with R/qtl.
* QKcartographer_segregation.py, analyzes the genetic map of a population and generates input files for generating a segregation distortion plot.
* QKcartographer_phenotypes.py, analyzes phenotypic data and generates input files for generating histograms and pairwise plots.

To understand how to use these scripts, following the tutorial available at [QKcartographer](https://github.com/matthewmoscou/QKcartographer).

## References
Basten CJ, Weir BS, Zeng Z-B (1994) Zmap–QTL cartographer. In: Smith C, Gavora JS, Benkel B, Chesnais J, Fairfull W, Gib- son JP, Kennedy BW, Burnside EB (eds) Proceedings of the 5th World Congress on Genetics Applied to Livestock Production: Computing Strategies and Software, Guelph, Ontario, Canada

Bettgenhaeuser, J., Gilbert, M., Ayliffe, M., Moscou, M.J. (2014) Nonhost resistance to rust pathogens – a continuation of continua. ***Frontiers in Plant Science*** **5**:664 [doi:10.3389/fpls.2014.00664](http://doi.org/10.3389/fpls.2014.00664)

