
<!-- README.md is generated from README.Rmd. Please edit that file -->

# GO-a-GO

<!-- badges: start -->

<!-- badges: end -->

GO-a-GO annotates functional terms that are overrepresented in a given
set of gene pairs. The enrichment of Gene Ontology terms is calculated
from a permutation test for overrepresentation of gene pairs that are
associated with a common term.

## Installation

You can install the development version of GO-a-GO from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("ajank/GOaGO")
```

## Example

Let us run GO-a-GO on a set of gene pairs associated with human
GM12878-specific chromatin loops. We can inspect a few rows of the
dataset:

``` r
library(GOaGO)
#> 

tail(genePairsGM12878Specific)
#> Key: <loopID>
#>    loopID chrom1   start1     end1 centroid1 distance_to_TSS1 geneID1 chrom2
#>     <int> <char>    <int>    <int>     <int>            <int>   <int> <char>
#> 1:   9384   chrX 55740000 55750000  55745000                0   10325   chrX
#> 2:   9390   chrX 57610000 57620000  57615000                0  158586   chrX
#> 3:   9418   chrX 73830000 73835000  73832500                0   51132   chrX
#> 4:   9423   chrX 77160000 77170000  77162500                0     538   chrX
#> 5:   9434   chrX 80060000 80070000  80065000                0  254065   chrX
#> 6:   9434   chrX 80060000 80070000  80065000                0  254065   chrX
#>      start2     end2 centroid2 distance_to_TSS2 geneID2 loop_distance pairID
#>       <int>    <int>     <int>            <int>   <int>         <int>  <int>
#> 1: 56750000 56760000  56755000             4016  442454       1010000    806
#> 2: 57940000 57950000  57940000             2933    7789        325000    807
#> 3: 74140000 74145000  74142500              286  340533        310000    808
#> 4: 77350000 77360000  77357500                0    5230        195000    809
#> 5: 80450000 80460000  80455000                0    6451        390000    810
#> 6: 80450000 80460000  80455000                0   79366        390000    811
```

The column names ending with `1` and `2` refer to the first and second
loop anchors, respectively. The essential columns are: `geneID1` and
`geneID2` (gene identifiers from a database of choice), `pairID`
(identifier of a gene pair), and `loopID` (identifier of a chromatin
loop). As some loop anchors overlapped Transcription Start Sites of
multiple genes, the dataset contains all gene combinations for these
loops.

When running GO-a-GO, we should specify that gene identifiers are from
the ENTREZ database, and use the Bioconductor `org.Hs.eg.db` package as
the source of Gene Ontology annotations for human genes, taking all
three subontologies (Biological Process, Molecular Function, and
Cellular Component). We will run 10,000 permutations, and use the
default *p*-value cutoff of 0.05 with Benjamini-Hochberg correction:

``` r
library(org.Hs.eg.db)

# set the number of CPU threads to use
library(BiocParallel)
options(MulticoreParam=MulticoreParam(workers=2))

goago <- GOaGO(genePairsGM12878Specific, keyType = "ENTREZID",
               OrgDb = org.Hs.eg.db, ont = "ALL", numPermutations = 10000L)
#> Warning in uniqueGenePairs(genePairs): removing 6 duplicated gene pair(s)
#> Warning in uniqueGenePairs(genePairs): removing 30 gene pair(s) containing the
#> same gene twice
```

We can inspect the overrepresented GO terms:

``` r
head(goago@result)
#>   ONTOLOGY         ID
#> 1       BP GO:0034728
#> 2       BP GO:0006334
#> 3       BP GO:0040029
#> 4       BP GO:0071824
#> 5       BP GO:0065004
#> 6       BP GO:0019886
#>                                                                         Description
#> 1                                                           nucleosome organization
#> 2                                                               nucleosome assembly
#> 3                                          epigenetic regulation of gene expression
#> 4                                                  protein-DNA complex organization
#> 5                                                      protein-DNA complex assembly
#> 6 antigen processing and presentation of exogenous peptide antigen via MHC class II
#>   Count       Ratio      BgRatio pvalue p.adjust qvalue
#> 1     8 0.010322581 2.807742e-04      0        0      0
#> 2     8 0.010322581 2.061935e-04      0        0      0
#> 3     6 0.007741935 7.201290e-04      0        0      0
#> 4     8 0.010322581 5.343226e-04      0        0      0
#> 5     8 0.010322581 4.640000e-04      0        0      0
#> 6     2 0.002580645 1.987097e-05      0        0      0
```

Note that some of the *p*-values have the value of zero, which means
that for all the randomizations fewer gene pairs were associated with a
given term than in the input dataset.

We can plot the overrepresented GO terms as a dotplot. The x-axis shows
the log<sub>2</sub> fold change of the fraction of gene pairs sharing
the given term; `Ratio` is the fraction in the input gene pairs, while
`BgRatio` is the fraction in permuted ones. Color indicates the adjusted
*p*-value:

``` r
DotPlot(goago)
```

<img src="man/figures/README-dotplot-1.png" width="100%" />

Note that by default only the terms associated with at least 5 gene
pairs are shown; you can change this by setting `minTermPairs` to any
other value.

We can also see the sampling distributions of numbers of gene pairs
sharing each GO term, obtained for the randomized gene pairs. From these
distributions, empirical *p*-values were calculated:

``` r
RidgePlot(goago)
#> Picking joint bandwidth of 0.08
```

<img src="man/figures/README-ridgeplot-1.png" width="100%" />
