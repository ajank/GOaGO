
<!-- README.md is generated from README.Rmd. Please edit that file -->

# GO-a-GO

<!-- badges: start -->

<!-- badges: end -->

GO-a-GO annotates Gene Ontology terms that are enriched in a given set
of gene pairs. The enrichment is calculated from a permutation test for
overrepresentation of gene pairs that are associated with a common term.
Such gene pairs are counted for the original set of gene pairs and
compared against randomized sets in which the structure of the pairs is
preserved, but the gene identities (including the associated terms) are
permuted.

## Installation

You can install the development version of GO-a-GO from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("ajank/GOaGO")
```

## Example

Let us run GO-a-GO on a set of gene pairs associated with chromatin
loops identified in human cell line GM12878. We can inspect a few rows
of the dataset:

``` r
library(GOaGO)
#> 

data("genePairsGM12878")
head(genePairsGM12878)
#>    loopID chrom1    start1      end1 centroid1 geneID1      tss1 strand1 chrom2
#>     <int> <char>     <int>     <int>     <int>  <char>     <int>  <char> <char>
#> 1:      3  chr10 101190000 101195000 101195833    2805 101190530       -  chr10
#> 2:      7  chr10 102100000 102110000 102100714    6319 102106772       +  chr10
#> 3:      7  chr10 102100000 102110000 102100714    6319 102106772       +  chr10
#> 4:     10  chr10 102810000 102815000 102803750   81621 102820999       +  chr10
#> 5:     19  chr10 103600000 103605000 103598750   30819 103603677       -  chr10
#> 6:     23  chr10 103985000 103990000 103987500   83401 103986143       +  chr10
#>       start2      end2 centroid2 geneID2      tss2 strand2
#>        <int>     <int>     <int>  <char>     <int>  <char>
#> 1: 101370000 101375000 101373333   81894 101372701       -
#> 2: 102270000 102280000 102277857   25956 102276754       -
#> 3: 102270000 102280000 102277857   25956 102279595       -
#> 4: 102900000 102905000 102903750    3195 102891061       +
#> 5: 103830000 103835000 103831250   79803 103825124       +
#> 6: 104160000 104165000 104162500    5662 104169041       -
```

The column names ending with `1` and `2` refer to the first and second
loop anchors, respectively. The essential columns are: `geneID1` and
`geneID2` (gene identifiers from a database of choice) and `loopID`
(identifier of a chromatin loop). As some loop anchors overlapped
multiple Transcription Start Sites, possibly of many genes, for these
loops the dataset contains all combinations.

When running GO-a-GO, we should specify that gene identifiers are from
the Entrez database, and use the Bioconductor `org.Hs.eg.db` package as
the source of Gene Ontology annotations for human genes, taking all
three subontologies (Biological Process, Molecular Function, and
Cellular Component). We will run the default 10,000 permutations, and
also use the default *p*-value cutoff of 0.05 with Benjamini-Hochberg
correction:

``` r
library(org.Hs.eg.db)

# set the number of CPU threads to use
library(BiocParallel)
options(MulticoreParam = MulticoreParam(workers = 2))

goago <- GOaGO(genePairsGM12878,
    keyType = "ENTREZID", OrgDb = org.Hs.eg.db, ont = "ALL"
)
#> Warning in uniqueGenePairs(genePairs): removing 509 duplicated gene pair(s)
#> Warning in uniqueGenePairs(genePairs): removing 87 gene pair(s) containing the
#> same gene twice
```

We can inspect the enriched GO terms:

``` r
head(as.data.frame(goago))
#>   ONTOLOGY         ID
#> 1       BP GO:0034728
#> 2       BP GO:0006334
#> 3       BP GO:0071824
#> 4       BP GO:0040029
#> 5       BP GO:0002495
#> 6       BP GO:0019886
#>                                                                         Description
#> 1                                                           nucleosome organization
#> 2                                                               nucleosome assembly
#> 3                                                  protein-DNA complex organization
#> 4                                          epigenetic regulation of gene expression
#> 5           antigen processing and presentation of peptide antigen via MHC class II
#> 6 antigen processing and presentation of exogenous peptide antigen via MHC class II
#>   Count   PairRatio      BgRatio FoldEnrichment pvalue p.adjust qvalue
#> 1    19 0.010900746 1.580608e-04       68.96552      0        0      0
#> 2    19 0.010900746 1.053356e-04      103.48584      0        0      0
#> 3    19 0.010900746 3.423982e-04       31.83646      0        0      0
#> 4    15 0.008605852 4.303500e-04       19.99733      0        0      0
#> 5     3 0.001721170 1.491681e-05      115.38462      0        0      0
#> 6     3 0.001721170 1.187608e-05      144.92754      0        0      0
```

Note that some of the *p*-values have the value of zero, which means
that for all the randomizations fewer gene pairs were associated with a
given term than in the input dataset.

We can plot the enriched GO terms as a dotplot. The X-axis shows the
fold enrichment of the fraction of gene pairs sharing the given term,
calculated as a quotient of `PairRatio` (the fraction in the input gene
pairs) and `BgRatio` (the fraction in permuted ones). Point size
indicates the number of input gene pairs sharing the given term, and
point color – the adjusted *p*-value:

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
#> Picking joint bandwidth of 0.0684
```

<img src="man/figures/README-ridgeplot-1.png" width="100%" />
