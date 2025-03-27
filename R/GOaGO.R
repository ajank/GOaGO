##' An S4 class to represent the results of GO-a-GO enrichment analysis.
##'
##' @slot result A data frame of the enriched Gene Ontology categories, with the
##'   following columns: \code{ONTOLOGY}, \code{ID}, \code{Description} (all of
##'   the GO term), \code{Count} (number of input gene pairs sharing the given
##'   term), \code{Ratio} (fraction of input gene pairs sharing the given term),
##'   \code{BgRatio} (fraction of permuted gene pairs sharing the given term),
##'   \code{pvalue}, \code{p.adjust}, \code{qvalue}
##' @slot pairTerms A data frame linking the enriched Gene Ontology categories
##'   with the input gene pairs, with the columns \code{pairID} and \code{ID}
##'   (of the GO term).
##' @slot permutedResult A data frame with the columns \code{ID} (of the GO
##'   term) and \code{Count}, keeping the numbers of permuted gene pairs sharing
##'   the term as obtained in every random permutation.
##' @export

setClass("GOaGO-result",
    representation=representation(
        result         = "data.frame",
        pairTerms      = "data.frame",
        permutedResult = "data.frame"
        # pvalueCutoff   = "numeric",
        # pAdjustMethod  = "character",
        # qvalueCutoff   = "numeric",
        # minTermPairs
        # number of permutations
        # organism       = "character",
        # ontology       = "character",
        # gene           = "character",
        # keytype        = "character",
        # universe       = "character",
        # gene2Symbol    = "character",
        # geneSets       = "list",
        # termsim        = "matrix",
        # method         = "character",
        # dr             = "list"
        )
    )


##' Gene Ontology enrichment analysis in a set of gene pairs.
##'
##' Given a data frame of gene pairs, this function will return the enriched
##' Gene Ontology categories after FDR control.
##'
##' @param genePairs a data frame with columns \code{geneID1} and \code{geneID2}
##'   containing gene identifiers.
##' @param OrgDb OrgDb
##' @param keyType keytype of input gene
##' @param ont One of "BP", "MF", and "CC" subontologies, or "ALL" for all
##'   three.
##' @param minTermPairs cutoff for number of pairs that share a GO term for this
##'   term to be considered
##' @param numPermutations number of permutations performed in the enrichment
##'   test
##' @param universe a set of background genes. If missing, all the genes from
##'   all the gene pairs will be used as background.
##' @param pvalueCutoff adjusted p-value cutoff on enrichment tests to report
##' @param pAdjustMethod one of "holm", "hochberg", "hommel", "bonferroni",
##'   "BH", "BY", "fdr", "none"
##' @param qvalueCutoff q-value cutoff on enrichment tests to report as
##'   significant. Tests must pass (i) \code{pvalueCutoff} on unadjusted
##'   p-values, (ii) \code{pvalueCutoff} on adjusted p-values, and (iii)
##'   \code{qvalueCutoff} on q-values to be reported.
##' @param minGSSize minimal size of genes annotated for testing
##' @param maxGSSize maximal size of genes annotated for testing
##' @returns A \code{GOaGO} instance.
##' @seealso \code{\link{GOaGO-result-class}}
##' @export

GOaGO <- function(genePairs, OrgDb, keyType = "ENTREZID", ont = "MF",
    minTermPairs = 1, numPermutations = 10000L, universe,
    pvalueCutoff = 0.05, pAdjustMethod = "BH", qvalueCutoff = 0.2, minGSSize = 10, maxGSSize = 500) {
    # ensure that gene identifiers are provided in the input data
    stopifnot("geneID1" %in% colnames(genePairs))
    stopifnot("geneID2" %in% colnames(genePairs))

    # ensure that pair identifiers are unique if provided
    if ("pairID" %in% colnames(genePairs))
        stopifnot(!anyDuplicated(genePairs$pairID))

    # gene universe
    if (missing(universe)) {
        # if not explicitly provided, universe is constituted by all the genes from all the gene pairs
        universe <- unique(c(genePairs$geneID1, genePairs$geneID2))
    }
    else {
        # check if all the genes are within the universe
        stopifnot(genePairs$geneID1 %in% universe)
        stopifnot(genePairs$geneID2 %in% universe)
        # remove duplicates in the universe
        universe <- unique(universe)
    }

    # use enrichGO to fetch the information on gene-to-term associations
    ego <- enrichGO(gene = universe, OrgDb = OrgDb, keyType = keyType, ont = ont, pAdjustMethod = "none", pvalueCutoff = Inf, qvalueCutoff = Inf, minGSSize = minGSSize, maxGSSize = maxGSSize)
    egoResult <- as.data.table(ego@result)

    # GO term universe
    ID_universe <- egoResult$ID
    stopifnot(!anyDuplicated(ID_universe))
    # take only the relevant columns
    go_split <- egoResult[, list(geneID = factor(strsplit(geneID, "/")[[1]], universe)),
        by = list(ID = factor(ID, ID_universe))]
    stopifnot(as.integer(go_split$geneID) == match(go_split$geneID, universe))
    stopifnot(!is.na(go_split$ID))

    # matrix (not a data frame) defining the gene pairs
    genePairsMatrix <- cbind(factor(genePairs$geneID1, universe), factor(genePairs$geneID2, universe))

    # construct sparse gene x term matrix
    geneTermMatrix <- Matrix(FALSE, nrow = length(universe), ncol = length(ID_universe), sparse = TRUE, doDiag = FALSE)
    geneTermMatrix[cbind(as.integer(go_split$geneID), as.integer(go_split$ID))] <- TRUE

    # now do the actual permutation testing
    countGenePairsPerTerm <- function(genePairsMatrix, geneTermMatrix, permuted = FALSE)
    {
        if (permuted)
        {
            permuted_genes <- sample.int(length(universe))
            gene1 <- permuted_genes[genePairsMatrix[, 1]]
            gene2 <- permuted_genes[genePairsMatrix[, 2]]
        }
        else
        {
            gene1 <- genePairsMatrix[, 1]
            gene2 <- genePairsMatrix[, 2]
        }

        bw <- geneTermMatrix[gene1, ] & geneTermMatrix[gene2, ]
        result <- as.integer(colSums(bw))
        return(result)
    }

    # simplified version of the above, without permutations and per-term aggregation
    pairTermsMatrix <- function(genePairsMatrix, geneTermMatrix)
    {
        gene1 <- genePairsMatrix[, 1]
        gene2 <- genePairsMatrix[, 2]
        bw <- geneTermMatrix[gene1, ] & geneTermMatrix[gene2, ]
        return(bw)
    }

    pairCountsPerTerm <- countGenePairsPerTerm(genePairsMatrix, geneTermMatrix)

    # take only the GO terms that are associated to at least minTermPairs gene pairs
    sel <- which(pairCountsPerTerm >= minTermPairs)
    ID_universe_reduced <- ID_universe[sel]
    pairCountsPerTerm_reduced <- pairCountsPerTerm[sel]
    geneTermMatrix_reduced <- geneTermMatrix[, sel]

    # TODO handle properly the case numPermutations == 0
    # associate actual and permuted gene pairs with their GO terms
    countPermutedGenePairsPerTerm <- function(i)
    {
        # if (i %% 1000 == 0) message(i)
        permuted_lapg <- countGenePairsPerTerm(genePairsMatrix, geneTermMatrix_reduced, permuted = T)
        return(permuted_lapg)
    }
    permutedPairCountsPerTerm_reduced <- do.call(rbind,
        bplapply(seq_len(numPermutations), countPermutedGenePairsPerTerm))

    # add empirical p-value
    pairs <- data.table(internalID = seq_along(pairCountsPerTerm_reduced), count = pairCountsPerTerm_reduced)
    permuted_pairs <- data.table(internalID = as.vector(col(permutedPairCountsPerTerm_reduced)), count = as.vector(permutedPairCountsPerTerm_reduced))

    dt <- rbind(data.table(permuted_pairs, control = TRUE), data.table(pairs, control = FALSE))
    dt <- dt[, list(
        BgRatio = mean(count[control == TRUE]) / nrow(genePairsMatrix),
        pvalue = sum(count[control == TRUE] >= count[control == FALSE]) / sum(control == TRUE)
        ), by = "internalID"]
    stopifnot(dt$internalID == seq_along(ID_universe_reduced))

    cols <- intersect(c("ONTOLOGY", "ID", "Description"), colnames(egoResult))
    result <- cbind(egoResult[sel, ..cols], Count = pairCountsPerTerm_reduced,
        Ratio = pairCountsPerTerm_reduced / nrow(genePairsMatrix), BgRatio = dt$BgRatio, pvalue = dt$pvalue)
    result[, p.adjust := p.adjust(pvalue, method = pAdjustMethod)]
    qobj <- tryCatch(qvalue(p = result$pvalue, lambda = 0.05, pi0.method = "bootstrap"),
        error = function(e) NULL)
    if (inherits(qobj, "qvalue")) {
        result[, qvalue := qobj$qvalues]
        sel_signif <- which(result$pvalue < pvalueCutoff & result$p.adjust <= pvalueCutoff & result$qvalue <= qvalueCutoff)
    }
    else {
        result[, qvalue := NA]
        sel_signif <- which(result$pvalue < pvalueCutoff & result$p.adjust <= pvalueCutoff)
    }

    # take only the GO terms that are significantly overrepresented
    result <- result[sel_signif, ]
    ID_universe_reduced_signif <- ID_universe_reduced[sel_signif]
    geneTermMatrix_reduced_signif <- geneTermMatrix_reduced[, sel_signif]
    permutedPairCountsPerTerm_reduced_signif <- permutedPairCountsPerTerm_reduced[, sel_signif]

    # construct sparse pair x terms matrix, convert it to tidy format
    pt <- which(pairTermsMatrix(genePairsMatrix, geneTermMatrix_reduced_signif), arr.ind=TRUE)
    pairTerms <- data.table()
    if ("pairID" %in% colnames(genePairs))
        pairTerms$pairID <- genePairs$pairID[pt[, 1]]
    else
        pairTerms$pairID <- pt[, 1]
    pairTerms$ID <- factor(ID_universe_reduced_signif)[pt[, 2]]

    # convert the results of the permutations to tidy format
    permutedResult <- data.table(
        ID = factor(ID_universe_reduced_signif)[as.vector(col(permutedPairCountsPerTerm_reduced_signif))],
        Count = as.vector(permutedPairCountsPerTerm_reduced_signif))

    # sort the results according to p-value
    result <- result[order(pvalue), ]

    result <- new("GOaGO-result", result = result, pairTerms = pairTerms, permutedResult = permutedResult)
    return(result)
}
