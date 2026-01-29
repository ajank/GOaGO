test_that("`uniqueGenePairs` should throw a warning when loops (gene pairs
containing the same gene twice) are found in the `genePairs` argument", {
    genePairs <- data.table(geneID1 = c("g1", "g3"), geneID2 = c("g2", "g3"))
    expect_warning(
        genePairs <- uniqueGenePairs(genePairs),
        "removing 1 gene pair\\(s\\) containing the same gene twice"
    )
    expect_equal(
        genePairs,
        data.table(geneID1 = "g1", geneID2 = "g2", pairID = 1)
    )
})

test_that("`uniqueGenePairs` should throw a warning when duplicate gene pairs
are found in the `genePairs` argument", {
    genePairs <- data.table(geneID1 = c("g1", "g2"), geneID2 = c("g2", "g1"))
    expect_warning(
        genePairs <- uniqueGenePairs(genePairs),
        "removing 1 duplicated gene pair\\(s\\)"
    )
    expect_equal(
        genePairs,
        data.table(geneID1 = "g1", geneID2 = "g2", pairID = 1)
    )
})


#' Construct a simulated set of gene pairs for testing GOaGO.
#'
#' The result will consist of \code{num_term1_pairs} pairs between genes that
#' have \code{term1}, but do not have \code{term2}, and \code{num_term2_pairs}
#' pairs between genes that have \code{term2}, but do not have \code{term1}. The
#' genes will be sampled from a given genome-wide annotation database.
#'
#' @param OrgDb OrgDb from which the genes will be sampled
#' @param term1 GO term 1
#' @param term2 GO term 2
#' @param num_term1_pairs number of pairs between genes having \code{term1}
#' @param num_term2_pairs number of pairs between genes having \code{term2}
#' @param keyType type of gene identifiers to use, such as "ENTREZID" or "ENSEMBL"
#'
#' @returns A data frame with columns \code{geneID1} and \code{geneID2}.
simulate_gene_pairs <- function(
      OrgDb, term1, term2, num_term1_pairs, num_term2_pairs,
      keyType = "ENTREZID"
) {
    results1 <- AnnotationDbi::select(OrgDb, keys = term1, columns = c(keyType), keytype = "GOALL")
    results2 <- AnnotationDbi::select(OrgDb, keys = term2, columns = c(keyType), keytype = "GOALL")

    candidate_genes1 <- unique(results1[[keyType]])
    candidate_genes2 <- unique(results2[[keyType]])

    genes1 <- setdiff(candidate_genes1, candidate_genes2)
    genes2 <- setdiff(candidate_genes2, candidate_genes1)

    stopifnot(num_term1_pairs > 0)
    stopifnot(num_term2_pairs > 0)
    stopifnot(length(genes1) >= 2 * num_term1_pairs)
    stopifnot(length(genes2) >= 2 * num_term2_pairs)

    genePairs <- data.frame(
        geneID1 = c(
            genes1[seq.int(from = 1, by = 2, length.out = num_term1_pairs)],
            genes2[seq.int(from = 1, by = 2, length.out = num_term2_pairs)]
        ),
        geneID2 = c(
            genes1[seq.int(from = 2, by = 2, length.out = num_term1_pairs)],
            genes2[seq.int(from = 2, by = 2, length.out = num_term2_pairs)]
        )
    )

    return(genePairs)
}

##' Parameters to construct simulated gene pairs for testing
term1 <- "GO:0003677" # DNA binding
term2 <- "GO:0046983" # protein dimerization activity
num_term1_pairs <- 10
num_term2_pairs <- 20
keyType <- "ENTREZID"

##' Note that it would be challenging to test the p-values and enrichment ratios
##' obtained from random permutations, so we do not check the values in columns
##' `BgRatio` and `pvalue`. To ensure that the results are not filtered out by
##' p-value or q-value, we set `pvalueCutoff` and `qvalueCutoff` to permissive
##' values.
##'
##' The GO terms chosen for this test are unlikely to ever become obsolete, but
##' the number of genes having this term (Gene Set Size) will surely change.
##' Hence, we set `minGSSize` and `maxGSSize` so that no term is excluded.

test_that("`GOaGO` should return a `GOaGO-result` instance, with statistics
(Count, Ratio) consistent with the simulated data", {
    testthat::skip_if_not_installed("org.Hs.eg.db")
    genePairs <- simulate_gene_pairs(
        org.Hs.eg.db::org.Hs.eg.db, term1, term2,
        num_term1_pairs, num_term2_pairs, keyType
    )
    goago <- GOaGO(genePairs,
        keyType = keyType, OrgDb = org.Hs.eg.db::org.Hs.eg.db,
        pvalueCutoff = Inf, qvalueCutoff = Inf, minGSSize = 0, maxGSSize = Inf
    )
    expect_s4_class(goago, "GOaGO-result")
    expect_equal(subset(goago@result, ID == term1)$Count, num_term1_pairs)
    expect_equal(subset(goago@result, ID == term2)$Count, num_term2_pairs)
    expect_equal(
        subset(goago@result, ID == term1)$Ratio,
        num_term1_pairs / (num_term1_pairs + num_term2_pairs)
    )
    expect_equal(
        subset(goago@result, ID == term2)$Ratio,
        num_term2_pairs / (num_term1_pairs + num_term2_pairs)
    )
})
