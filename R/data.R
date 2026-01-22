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
##' @name genePairsGM12878Specific
##' @format A data frame with 811 rows and 15 columns:
##' \describe{
##'   \item{loopID}{Loop identifier}
##'   \item{chrom\code{i}}{Chromosome of loop anchor}
##'   \item{start\code{i}}{Start coordinate of loop anchor}
##'   \item{end\code{i}}{End coordinate of loop anchor}
##'   \item{centroid\code{i}}{Centroid of loop anchor}
##'   \item{distance_to_TSS\code{i}}{Distance between the loop anchor and the closest TSS of the associated gene}
##'   \item{geneID\code{i}}{Entrez identifier of the associated gene}
##'   \item{loop_distance}{Distance between anchor centroids}
##'   \item{pairID}{Gene pair identifier}
##' }
##' where \code{i} is 1 or 2 for loop anchor 1 and 2, respectively.
##' @source <https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE63525>
NULL
