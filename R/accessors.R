## Coercion methods

`.as.data.frame.GOaGO-result` <- function(x, ...) {
    data.frame(x@result, ...)
}

setAs("GOaGO-result", "data.frame", function(from) {
    `.as.data.frame.GOaGO-result`(from)
})

##' as.data.frame method for \code{GOaGO-result} instance
##'
##' @name as.data.frame
##' @aliases as.data.frame,GOaGO-result-method
##' @docType methods
##' @rdname as.data.frame-methods
##'
##' @title as.data.frame method
##' @usage as.data.frame(x, row.names=NULL, optional=FALSE, ...)
##' @param x A \code{GOaGO-result} instance to coerce.
##' @param row.names,optional,... Not used. They are inherited from
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
##' library("GOaGO")
##' data("genePairsGM12878")
##'
##' goago <- GOaGO(genePairsGM12878, keyType = "ENTREZID", OrgDb = org.Hs.eg.db)
##' as.data.frame(goago)
setMethod(
    "as.data.frame", signature(x = "GOaGO-result"),
    `.as.data.frame.GOaGO-result`
)


##' as.data.table method for \code{GOaGO-result} instance
##'
##' @name as.data.table
##' @aliases as.data.table,GOaGO-result-method
##' @docType methods
##' @rdname as.data.table-methods
##'
##' @title as.data.table method
##' @usage as.data.table(x, keep.rownames=FALSE, ...)
##' @param x A \code{GOaGO-result} instance to coerce.
##' @param keep.rownames,... Not used. They are inherited from
##'   \code{data.table::as.data.table()}.
##' @returns A data table of the enriched Gene Ontology terms, with the
##'   following columns: \code{ONTOLOGY}, \code{ID}, \code{Description} (all of
##'   the GO term), \code{Count} (number of input gene pairs sharing the given
##'   term), \code{PairRatio} (fraction of input gene pairs sharing the given
##'   term), \code{BgRatio} (fraction of permuted gene pairs sharing the given
##'   term), \code{FoldEnrichment} (quotient of the two fractions),
##'   \code{pvalue}, \code{p.adjust}, \code{qvalue}.
##' @examples
##' library(org.Hs.eg.db)
##' library("GOaGO")
##' data("genePairsGM12878")
##'
##' goago <- GOaGO(genePairsGM12878, keyType = "ENTREZID", OrgDb = org.Hs.eg.db)
##' as.data.table(goago)
`.as.data.table.GOaGO-result` <- function(x, ...) {
    data.table(x@result, ...)
}

setAs("GOaGO-result", "data.table", function(from) {
    `.as.data.table.GOaGO-result`(from)
})

setMethod(
    "as.data.table", signature(x = "GOaGO-result"),
    `.as.data.table.GOaGO-result`
)


## Accessors

##' Key type accessor for \code{GOaGO-result} instance
##'
##' @param object of class \code{GOaGO-result}
##' @returns type of gene identifiers, such as "ENTREZID" or "ENSEMBL"
##' @export
##' @examples
##' library(org.Hs.eg.db)
##' library("GOaGO")
##' data("genePairsGM12878")
##'
##' goago <- GOaGO(genePairsGM12878, keyType = "ENTREZID", OrgDb = org.Hs.eg.db)
##' keyType(goago)
keyType <- function(object) {
    object@keyType
}

##' organism method for \code{GOaGO-result} instance
##'
##' @name organism
##' @aliases organism,GOaGO-result-method
##' @docType methods
##' @rdname organism-methods
##'
##' @title organism method
##' @param object A \code{GOaGO-result} instance.
##' @returns scientific name (i.e. genus and species, or genus and species and
##'   subspecies) of the organism
##' @usage organism(object)
##' @examples
##' library(org.Hs.eg.db)
##' library("GOaGO")
##' data("genePairsGM12878")
##'
##' goago <- GOaGO(genePairsGM12878, keyType = "ENTREZID", OrgDb = org.Hs.eg.db)
##' organism(goago)
setMethod(
    "organism", signature(object = "GOaGO-result"),
    function(object) {
        object@organism
    }
)
