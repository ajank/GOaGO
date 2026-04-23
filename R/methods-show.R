## Show method for objects

##' show method for \code{GOaGO-result} instance
##'
##' @name show
##' @aliases show,GOaGO-result-method
##' @docType methods
##' @rdname show-methods
##'
##' @title show method
##' @usage show(object)
##' @param object A \code{GOaGO-result} instance.
##' @returns message
##' @examples
##' library(org.Hs.eg.db)
##' data("genePairsGM12878")
##'
##' genePairsSubset <- subset(genePairsGM12878, chrom1 == "chr11")
##' goago <- GOaGO(genePairsSubset, keyType = "ENTREZID", OrgDb = org.Hs.eg.db)
##' show(goago)
setMethod(
    "show", signature(object = "GOaGO-result"),
    function(object) {
        cat("#\n")
        cat("# Gene Ontology over-representation test in a set of gene pairs\n")

        cat("#\n")
        cat("#...@organism", "\t", object@organism, "\n")
        cat("#...@ontology", "\t", object@ontology, "\n")
        cat("#...@keyType", "\t", object@keyType, "\n")

        cat("#\n")
        cat(sprintf(
            "#...%d gene pairs considered, %d permutations\n",
            nrow(object@genePairs), object@numPermutations
        ))
        cat(sprintf(
            "#...identified GO terms shared by at least %d %s\n",
            object@minCount,
            ifelse(object@minCount == 1, "gene pair", "gene pairs")
        ))
        cat(sprintf(
            "#...p-values adjusted by '%s' with cutoff < %s\n",
            object@pAdjustMethod, object@pvalueCutoff
        ))

        cat("#\n")
        n <- nrow(object@result)
        cat(sprintf(
            "#...%d enriched %s found\n", n,
            ifelse(n == 1, "term", "terms")
        ))
        if (n > 0) str(object@result)
    }
)
