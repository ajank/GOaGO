## Show method for objects

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
            object@minTermPairs,
            ifelse(object@minTermPairs == 1, "gene pair", "gene pairs")
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
