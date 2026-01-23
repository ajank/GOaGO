##' Gene pairs associated with GM12878-specific chromatin loops
##'
##' The dataset is based on 3,638 chromatin loops specific to cell line GM12878,
##' i.e. the ones which did not overlap any loop annotated in the other studied
##' cell types. Of these chromatin loops, 659 overlapped at least one gene
##' Transcription Start Site (TSS) at both loop anchors. As some loop anchors
##' overlapped TSSes of multiple genes, the dataset contains all gene
##' combinations for these loops, yielding a total of 811 GM12878-specific gene
##' pairs, of which 775 pairs are unique.
##'
##' @format ## `genePairsGM12878Specific`
##' A data frame with 811 rows and 14 columns:
##' \describe{
##'   \item{loopID}{Loop identifier}
##'   \item{chrom1}{Chromosome of loop anchor 1}
##'   \item{start1}{Start coordinate of loop anchor 1}
##'   \item{end1}{End coordinate of loop anchor 1}
##'   \item{centroid1}{Centroid of loop anchor 1}
##'   \item{distance_to_TSS1}{Distance between loop anchor 1 and the closest TSS of the associated gene}
##'   \item{geneID1}{Entrez identifier of the gene associated to loop anchor 1}
##'   \item{chrom2}{Chromosome of loop anchor 2}
##'   \item{start2}{Start coordinate of loop anchor 2}
##'   \item{end2}{End coordinate of loop anchor 2}
##'   \item{centroid2}{Centroid of loop anchor 2}
##'   \item{distance_to_TSS2}{Distance between loop anchor 2 and the closest TSS of the associated gene}
##'   \item{geneID2}{Entrez identifier of the gene associated to loop anchor 2}
##'   \item{loop_distance}{Genomic distance between the anchor centroids.}
##' }
##' @source <https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE63525>
"genePairsGM12878Specific"
