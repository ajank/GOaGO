##' An S4 class to represent the results of GO-a-GO enrichment analysis.
##'
##' @slot result A data frame of the enriched Gene Ontology terms, with the
##'   following columns: \code{ONTOLOGY}, \code{ID}, \code{Description} (all of
##'   the GO term), \code{Count} (number of input gene pairs sharing the given
##'   term), \code{Ratio} (fraction of input gene pairs sharing the given term),
##'   \code{BgRatio} (fraction of permuted gene pairs sharing the given term),
##'   \code{pvalue}, \code{p.adjust}, \code{qvalue}.
##' @slot pvalueCutoff adjusted p-value cutoff on enrichment tests
##' @slot pAdjustMethod p-value adjustment method
##' @slot qvalueCutoff q-value cutoff on enrichment tests
##' @slot minTermPairs cutoff for number of pairs that share a GO term for this
##'   term to be considered
##' @slot numPermutations number of permutations performed in the enrichment
##'   test
##' @slot minGSSize minimal size of genes annotated for testing
##' @slot maxGSSize maximal size of genes annotated for testing
##' @slot organism scientific name of the organism
##' @slot ontology one of "BP", "MF", and "CC" subontologies, or "ALL" for all
##'   three
##' @slot keyType type of gene identifiers, such as "ENTREZID" or "ENSEMBL"
##' @slot genePairs A data frame with the input gene pairs, with the columns
##'   \code{geneID1}, \code{geneID2} and \code{pairID}.
##' @slot pairTerms A data frame linking the enriched Gene Ontology terms with
##'   the input gene pairs, with the columns \code{pairID} and \code{ID} (of the
##'   GO term).
##' @slot permutedResult A data frame with the columns \code{ID} (of the GO
##'   term) and \code{Count}, keeping the numbers of permuted gene pairs sharing
##'   the term as obtained in every random permutation.
##' @slot universe a set of background genes
##' @export
setClass("GOaGO-result",
    slots = c(
        result          = "data.frame",
        pvalueCutoff    = "numeric",
        pAdjustMethod   = "character",
        qvalueCutoff    = "numeric",
        minTermPairs    = "numeric",
        numPermutations = "numeric",
        minGSSize       = "numeric",
        maxGSSize       = "numeric",
        organism        = "character",
        ontology        = "character",
        keyType         = "character",
        genePairs       = "data.frame",
        pairTerms       = "data.frame",
        permutedResult  = "data.frame",
        universe        = "character"
    )
)


##' Extract unique gene pairs from the data frame provided.
##'
##' Given a data frame of gene pairs, this function will return the unique
##' pairs of genes, removing loops (gene pairs containing the same gene twice)
##' and duplicates. Note that gene pair (A, B) is a duplicate of (B, A).
##'
##' @param genePairs a data frame with columns \code{geneID1} and \code{geneID2}
##'   containing gene identifiers; column \code{pairID} will also be used if
##'   provided.
##' @returns A data frame with columns \code{pairID}, \code{geneID1} and
##'   \code{geneID2}. If loops or duplicates were removed, a warning will
##'   alert you. If column \code{pairID} was not provided in \code{genePairs},
##'   an integer vector equal to \code{seq_len(nrow(result))} will be used.
uniqueGenePairs <- function(genePairs) {
    # ensure that gene identifiers are provided in the input data
    stopifnot("geneID1" %in% colnames(genePairs))
    stopifnot("geneID2" %in% colnames(genePairs))

    # ensure that pair identifiers are unique if provided
    pairID_provided <- "pairID" %in% colnames(genePairs)
    if (pairID_provided) {
        stopifnot(!anyDuplicated(genePairs$pairID))
        genePairs <- with(genePairs, data.table(pairID, geneID1, geneID2))
    } else {
        genePairs <- with(genePairs, data.table(geneID1, geneID2))
    }

    n <- nrow(genePairs)
    genePairs[, gid1 := pmin(geneID1, geneID2)]
    genePairs[, gid2 := pmax(geneID1, geneID2)]
    genePairs <- unique(genePairs, by = c("gid1", "gid2"))
    if (nrow(genePairs) != n) {
        warning("removing ", n - nrow(genePairs), " duplicated gene pair(s)")
    }
    genePairs[, gid1 := NULL]
    genePairs[, gid2 := NULL]

    n <- nrow(genePairs)
    genePairs <- genePairs[geneID1 != geneID2, ]
    if (nrow(genePairs) != n) {
        warning(
            "removing ", n - nrow(genePairs),
            " gene pair(s) containing the same gene twice"
        )
    }

    if (!pairID_provided) {
        genePairs[, pairID := .I]
        setcolorder(genePairs, "pairID")
    }

    return(genePairs)
}


##' Gene Ontology enrichment analysis in a set of gene pairs.
##'
##' Given a data frame of gene pairs, this function will return the enriched
##' Gene Ontology terms after FDR control.
##'
##' @param genePairs a data frame with columns \code{geneID1} and \code{geneID2}
##'   containing gene identifiers; column \code{pairID} will also be used if
##'   provided.
##' @param OrgDb OrgDb
##' @param keyType type of gene identifiers, such as "ENTREZID" or "ENSEMBL"
##' @param ont one of "BP", "MF", and "CC" subontologies, or "ALL" for all three
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
##' @returns A \code{GOaGO-result} instance.
##' @seealso \code{\link{GOaGO-result-class}}
##' @export
##' @examples
##' library(org.Hs.eg.db)
##' library("GOaGO")
##' data("genePairsGM12878Specific")
##'
##' goago <- GOaGO(genePairsGM12878Specific,
##'     keyType = "ENTREZID", OrgDb = org.Hs.eg.db
##' )
##'
##' show(goago)
GOaGO <- function(
    genePairs, OrgDb, keyType = "ENTREZID", ont = "MF",
    minTermPairs = 1, numPermutations = 10000, universe, pvalueCutoff = 0.05,
    pAdjustMethod = "BH", qvalueCutoff = 0.2, minGSSize = 10, maxGSSize = 500
) {
    # extract unique gene pairs from the data frame provided
    genePairs <- uniqueGenePairs(genePairs)
    # ensure that gene identifiers are character vectors
    genePairs$geneID1 <- as.character(genePairs$geneID1)
    genePairs$geneID2 <- as.character(genePairs$geneID2)

    # gene universe
    if (missing(universe)) {
        # if not explicitly provided, universe is constituted by all the genes
        # from all the gene pairs
        universe <- c(genePairs$geneID1, genePairs$geneID2)
    } else {
        # ensure that gene identifiers are character vectors
        universe <- as.character(universe)
        # check if all the genes are within the universe
        stopifnot(genePairs$geneID1 %in% universe)
        stopifnot(genePairs$geneID2 %in% universe)
    }
    # remove duplicates from the universe
    universe <- unique(universe)

    # use enrichGO to fetch the information on gene-to-term associations
    ego <- enrichGO(
        gene = universe, OrgDb = OrgDb, keyType = keyType,
        ont = ont, pAdjustMethod = "none", pvalueCutoff = Inf,
        qvalueCutoff = Inf, minGSSize = minGSSize, maxGSSize = maxGSSize
    )
    egoResult <- as.data.table(ego)

    # GO term universe
    ID_universe <- egoResult$ID
    stopifnot(!anyDuplicated(ID_universe))
    # take only the relevant columns
    go_split <- egoResult[,
        list(geneID = factor(strsplit(geneID, "/")[[1]], universe)),
        by = list(ID = factor(ID, ID_universe)),
        showProgress = FALSE
    ]
    stopifnot(as.integer(go_split$geneID) == match(go_split$geneID, universe))
    stopifnot(!is.na(go_split$ID))

    # matrix (not a data frame) defining the gene pairs
    genePairsMatrix <- cbind(
        factor(genePairs$geneID1, universe),
        factor(genePairs$geneID2, universe)
    )

    # construct sparse gene x term matrix
    geneTermMatrix <- Matrix(FALSE,
        nrow = length(universe),
        ncol = length(ID_universe), sparse = TRUE, doDiag = FALSE
    )
    ind <- cbind(as.integer(go_split$geneID), as.integer(go_split$ID))
    geneTermMatrix[ind] <- TRUE

    # a function useful for permutation testing
    countGenePairsPerTerm <- function(genePairsMatrix,
        geneTermMatrix,
        permuted = FALSE
    ) {
        if (permuted) {
            permuted_genes <- sample.int(length(universe))
            gene1 <- permuted_genes[genePairsMatrix[, 1]]
            gene2 <- permuted_genes[genePairsMatrix[, 2]]
        } else {
            gene1 <- genePairsMatrix[, 1]
            gene2 <- genePairsMatrix[, 2]
        }

        bw <- geneTermMatrix[gene1, , drop = FALSE] &
            geneTermMatrix[gene2, , drop = FALSE]
        result <- as.integer(colSums(bw))
        return(result)
    }

    # simplified version of the above, without permutations and per-term
    # aggregation
    pairTermsMatrix <- function(genePairsMatrix, geneTermMatrix) {
        gene1 <- genePairsMatrix[, 1]
        gene2 <- genePairsMatrix[, 2]
        bw <- geneTermMatrix[gene1, , drop = FALSE] &
            geneTermMatrix[gene2, , drop = FALSE]
        return(bw)
    }

    pairCountsPerTerm <- countGenePairsPerTerm(genePairsMatrix, geneTermMatrix)

    # take only the GO terms that are associated to at least minTermPairs gene
    # pairs
    sel <- which(pairCountsPerTerm >= minTermPairs)
    ID_universe_reduced <- ID_universe[sel]
    pairCountsPerTerm_reduced <- pairCountsPerTerm[sel]
    geneTermMatrix_reduced <- geneTermMatrix[, sel, drop = FALSE]

    # now do the permutation test: for each permutation, count permuted gene
    # pairs sharing each GO term
    countPermutedGenePairsPerTerm <- function(i) {
        result <- countGenePairsPerTerm(genePairsMatrix, geneTermMatrix_reduced,
            permuted = TRUE
        )
        return(result)
    }
    # store the resulting counts in term x permutation matrix
    permutedPairCountsPerTerm_reduced <- do.call(
        cbind,
        bplapply(seq_len(numPermutations), countPermutedGenePairsPerTerm)
    )

    # aggregate the results of permutation testing
    cols <- intersect(c("ONTOLOGY", "ID", "Description"), colnames(egoResult))
    result <- cbind(
        egoResult[sel, ..cols],
        Count = pairCountsPerTerm_reduced,
        Ratio = pairCountsPerTerm_reduced / nrow(genePairsMatrix),
        BgRatio = rowMeans(permutedPairCountsPerTerm_reduced) /
            nrow(genePairsMatrix),
        pvalue = rowMeans(permutedPairCountsPerTerm_reduced >=
            pairCountsPerTerm_reduced)
    )

    # adjust p-values and estimate q-values
    result[, p.adjust := p.adjust(pvalue, method = pAdjustMethod)]
    qobj <- tryCatch(qvalue(
        p = result$pvalue, lambda = 0.05,
        pi0.method = "bootstrap"
    ), error = function(e) NULL)
    if (inherits(qobj, "qvalue")) {
        result[, qvalue := qobj$qvalues]
        sel_signif <- which(result$pvalue < pvalueCutoff &
            result$p.adjust <= pvalueCutoff & result$qvalue <= qvalueCutoff)
    } else {
        result[, qvalue := NA]
        sel_signif <- which(result$pvalue < pvalueCutoff &
            result$p.adjust <= pvalueCutoff)
    }

    # take only the GO terms that are significantly overrepresented
    result <- result[sel_signif, ]
    ID_universe_reduced_signif <- ID_universe_reduced[sel_signif]
    geneTermMatrix_reduced_signif <- geneTermMatrix_reduced[, sel_signif,
        drop = FALSE
    ]
    permutedPairCountsPerTerm_reduced_signif <-
        permutedPairCountsPerTerm_reduced[sel_signif, , drop = FALSE]

    # construct sparse pair x term matrix, convert it to tidy format
    pt <- which(pairTermsMatrix(genePairsMatrix, geneTermMatrix_reduced_signif),
        arr.ind = TRUE
    )
    pairTerms <- data.table(
        pairID = genePairs$pairID[pt[, 1]],
        ID = factor(ID_universe_reduced_signif)[pt[, 2]]
    )

    # convert the results of the permutations to tidy format
    ind <- as.vector(row(permutedPairCountsPerTerm_reduced_signif))
    permutedResult <- data.table(
        ID = factor(ID_universe_reduced_signif)[ind],
        Count = as.vector(permutedPairCountsPerTerm_reduced_signif)
    )

    # sort the results according to p-value
    result <- result[order(pvalue), ]

    result <- new("GOaGO-result",
        result = result,
        pvalueCutoff = pvalueCutoff,
        pAdjustMethod = pAdjustMethod,
        qvalueCutoff = qvalueCutoff,
        minTermPairs = minTermPairs,
        numPermutations = numPermutations,
        minGSSize = minGSSize,
        maxGSSize = maxGSSize,
        organism = AnnotationDbi::species(OrgDb),
        ontology = ont,
        keyType = keyType,
        genePairs = genePairs,
        pairTerms = pairTerms,
        permutedResult = permutedResult,
        universe = universe
    )
    return(result)
}
