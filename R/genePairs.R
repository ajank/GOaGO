##' Extract the enriched Gene Ontology terms along with gene pairs sharing them
##'
##' @param object of class \code{GOaGO-result}
##' @param OrgDb OrgDb to map gene identifiers to gene symbols
##' @returns A data frame similar to returned by \code{as.data.frame(object)},
##'   with additional columns \code{pairID}, \code{geneID1} and \code{geneID2},
##'   for each enriched Gene Ontology term providing all the gene pairs sharing
##'   this term, one gene pair in each row. If \code{OrgDb} is provided,
##'   \code{geneSymbol1} and \code{geneSymbol2} will also be added.
##' @seealso \code{\link{as.data.frame,GOaGO-result-method}}
##' @export
##' @examples
##' library(org.Hs.eg.db)
##' library("GOaGO")
##' data("genePairsGM12878")
##'
##' goago <- GOaGO(genePairsGM12878, keyType = "ENTREZID", OrgDb = org.Hs.eg.db)
##' genePairs(goago, OrgDb = org.Hs.eg.db)
genePairs <- function(object, OrgDb = NULL) {
    # first merge the enriched Gene Ontology terms with all the input gene
    # pairs that share the given term (by identifier of the GO term)
    gp <- merge(as.data.table(object@result), object@pairTerms,
        by = "ID", all = TRUE, sort = FALSE, allow.cartesian = TRUE
    )
    # then add the information on gene identifiers for each gene pair
    gp <- merge(gp, object@genePairs, by = "pairID", all.x = TRUE, sort = FALSE)

    # add gene symbols if OrgDb is provided
    if (!is.null(OrgDb)) {
        geneIDs <- unique(with(gp, c(geneID1, geneID2)))
        # suppress the message on 1:1 mapping between keys and columns
        suppressMessages(geneSymbols <- AnnotationDbi::mapIds(OrgDb, geneIDs,
            column = "SYMBOL", keytype = object@keyType
        ))

        gp$geneSymbol1 <- geneSymbols[match(gp$geneID1, geneIDs)]
        gp$geneSymbol2 <- geneSymbols[match(gp$geneID2, geneIDs)]
    }

    # set the column order so that the columns from object@result will be first
    # and column order from object@genePairs will also be preserved
    setcolorder(gp, unique(c(names(object@result), names(object@genePairs))))

    return(gp)
}
