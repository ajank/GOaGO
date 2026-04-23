## Coercion to a data frame

`as.data.frame.GOaGO-result` <- function(x, ...) {
    data.frame(x@result, ...)
}

setAs("GOaGO-result", "data.frame", function(from) {
    `as.data.frame.GOaGO-result`(from)
})

##' Coerce a \code{GOaGO-result} object to a data frame
##'
##' @name as.data.frame
##' @aliases as.data.frame,GOaGO-result-method
##' @docType methods
##' @rdname as.data.frame-methods
##'
##' @usage as.data.frame(x, row.names=NULL, optional=FALSE, ...)
##' @param x The object to coerce.
##' @param row.names,optional,... Not used, inherited from
##'   \code{base::as.data.frame()}.
##' @returns A data frame of the enriched Gene Ontology terms, with the
##'   following columns: \code{ONTOLOGY}, \code{ID}, \code{Description} (all of
##'   the GO term), \code{Count} (number of input gene pairs sharing the given
##'   term), \code{PairRatio} (fraction of input gene pairs sharing the given
##'   term), \code{BgRatio} (fraction of permuted gene pairs sharing the given
##'   term), \code{FoldEnrichment} (quotient of the two fractions),
##'   \code{pvalue}, \code{p.adjust}, \code{qvalue}.
##' @examples
##' library(org.Hs.eg.db)
##' data("genePairsGM12878")
##'
##' genePairsSubset <- subset(genePairsGM12878, chrom1 == "chr11")
##' goago <- GOaGO(genePairsSubset, keyType = "ENTREZID", OrgDb = org.Hs.eg.db)
##' as.data.frame(goago)
setMethod(
    "as.data.frame", signature(x = "GOaGO-result"),
    `as.data.frame.GOaGO-result`
)


## Accessors

##' Accessors and show method for \code{GOaGO-result} objects
##'
##' @name GOaGO-accessors
##' @param object of class \code{GOaGO-result}
##' @returns
##' \code{genePairs} returns a data frame with the input gene pairs, with the
##' columns \code{geneID1}, \code{geneID2} and \code{pairID}.
##'
##' \code{keyType} returns the type of gene identifiers, such as "ENTREZID" or
##' "ENSEMBL".
##'
##' \code{organism} returns the scientific name (i.e. genus and species, or
##' genus and species and subspecies) of the organism.
##'
##' \code{show} displays the object, and returns an invisible \code{NULL}.
##' @examples
##' library(org.Hs.eg.db)
##' data("genePairsGM12878")
##'
##' genePairsSubset <- subset(genePairsGM12878, chrom1 == "chr11")
##' goago <- GOaGO(genePairsSubset, keyType = "ENTREZID", OrgDb = org.Hs.eg.db)
##' show(goago)
##'
##' genePairs(goago)
##' keyType(goago)
##' organism(goago)
NULL

##' @rdname GOaGO-accessors
##' @export
genePairs <- function(object) {
    object@genePairs
}

##' @rdname GOaGO-accessors
##' @export
keyType <- function(object) {
    object@keyType
}

##' @name organism
##' @aliases organism,GOaGO-result-method
##' @rdname GOaGO-accessors
##' @usage organism(object)
##' @export
setMethod(
    "organism", signature(object = "GOaGO-result"),
    function(object) {
        object@organism
    }
)
