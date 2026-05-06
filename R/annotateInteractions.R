# find the metadata column containing gene identifiers, following the convention
# of GenomicRanges::makeGRangesFromDataFrame()
.find_geneid_col <- function(metadata_colnames, geneid_column) {
    ind <- which(metadata_colnames %in% geneid_column)
    if (length(ind) == 0) {
        stop("cannot find gene ID column")
    }
    if (length(ind) >= 2) {
        stop("cannot determine gene ID column unambiguously")
    }
    return(ind)
}


##' Convert gene transcripts to Transcription Start Sites
##'
##' Convert a \code{GRanges} object with gene transcripts to a \code{GRanges}
##' object with gene TSSes of length 1 bp, with duplicate rows removed. Assumes
##' that gene identifiers are provided in one of the metadata columns, either as
##' a vector (possibly containing \code{NA} values) or a \code{CharacterList}.
##' The function is idempotent, i.e. can be applied multiple times without
##' changing the result.
##'
##' @param transcripts \code{TxDb} annotation or other \code{GRanges} object
##' @param geneid_column A character vector of recognized names for the metadata
##' column in \code{transcripts} that contains gene identifiers. If none or more
##' than one is found, an error is raised.
##'
##' @returns A \code{GRanges} object with gene identifiers in metadata column
##' \code{geneID} being a vector.
##' @export
##'
##' @examples
##' library(TxDb.Hsapiens.UCSC.hg19.knownGene)
##'
##' # take only the transcripts of coding genes by ensuring that the coding
##' # sequence strand is not NA
##' transcripts <- transcripts(TxDb.Hsapiens.UCSC.hg19.knownGene,
##'     columns = "gene_id", filter = list(cds_strand = c("-", "+"))
##' )
##'
##' convertTranscriptsToTSS(transcripts)
convertTranscriptsToTSS <- function(
    transcripts, geneid_column = c("gene_id", "GENEID", "geneID")
) {
    # find the metadata column with gene IDs
    geneid_col <- .find_geneid_col(names(mcols(transcripts)), geneid_column)
    geneid_values <- mcols(transcripts)[[geneid_col]]

    # convert transcripts to TSSes
    tss <- resize(transcripts, width = 1, fix = "start")
    # remove metadata columns from the resulting GRanges object
    mcols(tss) <- NULL

    # repopulate gene ID column in the metadata, unlisting a CharacterList
    # or removing NAs from an atomic vector as needed
    if (inherits(geneid_values, "List")) {
        ind <- rep(seq_along(tss), vapply(geneid_values, length, integer(1)))
        tss <- tss[ind]
        tss$geneID <- unlist(geneid_values)
    } else {
        stopifnot(is.vector(geneid_values))
        ind <- !is.na(geneid_values)
        tss <- tss[ind]
        tss$geneID <- geneid_values[ind]
    }

    # remove duplicate rows; here the conversion to a data table is to avoid
    # collapsing TSSes of different genes (geneID column is in the metadata)
    tss <- tss[!duplicated(as.data.table(tss)), ]

    return(tss)
}


# use pre-computed TSSes if provided, otherwise use convertTranscriptsToTSS
.determine_tss <- function(transcripts = NULL, tss = NULL) {
    if (is.null(tss)) {
        if (is.null(transcripts)) {
            stop("either `transcripts` or `tss` must be provided")
        }
        tss <- convertTranscriptsToTSS(transcripts)
    } else {
        if (!is.null(transcripts)) {
            stop("`transcripts` and `tss` cannot be both provided")
        }
        stopifnot(width(tss) == 1)
        stopifnot(is.vector(tss$geneID))
    }

    return(tss)
}


# determine keyType from GRanges metadata
.determine_keyType <- function(tss, keyType = NULL) {
    # if a value is provided, use it
    if (!is.null(keyType)) {
        return(keyType)
    }

    keyType_description <- metadata(tss)$genomeInfo$`Type of Gene ID`
    # the above can be NULL if there is no genomeInfo in the metadata etc.
    if (is.null(keyType_description)) {
        keyType_description <- ""
    }

    if (keyType_description == "Entrez Gene ID") {
        return("ENTREZID")
    } else if (keyType_description == "Ensembl gene ID") {
        return("ENSEMBL")
    } else {
        stop(
            "type of gene identifiers cannot be determined from the metadata, ",
            "please provide it as `keyType` argument"
        )
    }
}


##' Associate interaction anchors to the nearest TSSes
##'
##' Interaction anchors are associated to TSSes as follows. Each anchor is
##' associated to all the TSSes it overlaps. If there is no such overlap, then
##' the anchor is associated to all TSSes with the shortest distance to the
##' anchor, if this distance is not larger than \code{maxDistanceToTSS}.
##'
##' Either \code{transcripts} or \code{tss} must be provided, but not both.
##'
##' @param anchors object of class \code{GRanges}
##' @param transcripts \code{TxDb} annotation or other \code{GRanges} object
##' @param tss object of class \code{GRanges} as returned by
##'   \code{\link{convertTranscriptsToTSS}}
##' @param keyType type of gene identifiers, such as "ENTREZID" or "ENSEMBL", if
##'   it cannot be determined from metadata of \code{transcripts} or \code{tss}
##' @param maxDistanceToTSS maximal distance to extend the search for nearest
##'   TSS outside the anchor, or -1 (the default) to skip the extension
##'
##' @returns A data table with columns \code{interactionID} (index of the
##'   anchor), \code{chrom}, \code{start}, \code{end} (coordinates of the
##'   anchor), \code{geneID} (gene identifier from `transcripts` or `tss`),
##'   \code{tss} (TSS position) and \code{strand} (TSS strand).
##' @seealso \code{\link{convertTranscriptsToTSS}}
##' @export
##'
##' @examples
##' library(TxDb.Hsapiens.UCSC.hg19.knownGene)
##'
##' # take only the transcripts of coding genes by ensuring that the coding
##' # sequence strand is not NA
##' transcripts <- transcripts(TxDb.Hsapiens.UCSC.hg19.knownGene,
##'     columns = "gene_id", filter = list(cds_strand = c("-", "+"))
##' )
##'
##' tss <- convertTranscriptsToTSS(transcripts)
##' gr <- GRanges("chr1", IRanges(c(42001, 890001), c(62000, 900000)))
##'
##' # note that anchors are associated to TSSes outside the anchor only if there
##' # are no TSSes overlapping the anchor
##' annotateAnchors(gr, tss, maxDistanceToTSS = 10e3)
##'
##' # this may yield more associations to TSSes outside the anchors
##' annotateAnchors(gr + 10e3, tss)
annotateAnchors <- function(
    anchors, transcripts = NULL, tss = NULL, keyType = NULL,
    maxDistanceToTSS = -1
) {
    # prevent "no visible binding for global variable" NOTEs in R CMD check
    distance_to_tss <- min_distance_to_tss <- NULL

    # for each anchor, find all TSSes not further away than maxDistanceToTSS
    tss <- .determine_tss(transcripts, tss)
    ov <- findOverlaps(anchors, tss, maxgap = maxDistanceToTSS)

    # keep information on anchor coordinates, distance to the nearest TSS,
    # as well as gene ID, position (1-based) and strand of the nearest TSS
    nearest_tss <- data.table(
        interactionID = queryHits(ov),
        chrom = as.vector(seqnames(anchors))[queryHits(ov)],
        start = start(anchors)[queryHits(ov)],
        end = end(anchors)[queryHits(ov)],
        distance_to_tss = ifelse(
            width(pintersect(anchors[queryHits(ov)], tss[subjectHits(ov)])) > 0,
            -1L,
            distance(anchors[queryHits(ov)], tss[subjectHits(ov)])
        ),
        geneID = tss$geneID[subjectHits(ov)],
        tss = start(tss)[subjectHits(ov)],
        strand = as.vector(strand(tss))[subjectHits(ov)]
    )

    # for each anchor, keep only the nearest TSSes (possibly more than one)
    nearest_tss[, min_distance_to_tss := min(distance_to_tss), by = "interactionID"]
    nearest_tss <- nearest_tss[distance_to_tss == min_distance_to_tss, ]
    nearest_tss[, min_distance_to_tss := NULL]
    nearest_tss[, distance_to_tss := NULL]

    # keep seqinfo and keyType from gene TSSes as attributes
    attr(nearest_tss, "seqinfo") <- seqinfo(tss)
    attr(nearest_tss, "keyType") <- .determine_keyType(tss, keyType)

    return(nearest_tss)
}


##' Associate interactions to gene pairs
##'
##' Both interaction anchors are associated to TSSes as described in
##' \code{\link{annotateAnchors}}. Briefly, each anchor is associated to all the
##' TSSes it overlaps, or to all closest TSSes up to \code{maxDistanceToTSS} if
##' there is no such overlap. The annotation of an interaction is a Cartesian
##' product of annotations for both anchors.
##'
##' Either \code{transcripts} or \code{tss} must be provided, but not both.
##'
##' @param interactions object of class \code{Pairs} (of \code{GRanges}) or
##'   \code{GenomicInteractions}
##' @param transcripts \code{TxDb} annotation or other \code{GRanges} object
##' @param tss object of class \code{GRanges} as returned by
##'   \code{\link{convertTranscriptsToTSS}}
##' @param keyType type of gene identifiers, such as "ENTREZID" or "ENSEMBL", if
##'   it cannot be determined from metadata of \code{transcripts} or \code{tss}
##' @param maxDistanceToTSS maximal distance to extend the search for nearest
##'   TSS outside the anchor, or -1 (the default) to skip the extension
##'
##' @returns A data table with columns \code{interactionID} (index of the
##'   interaction), \code{chrom1}, \code{start1}, \code{end1}, \code{chrom2},
##'   \code{start2}, \code{end2} (coordinates of both anchors), \code{geneID1},
##'   \code{geneID2} (gene identifiers from `transcripts` or `tss` for both
##'   anchors), \code{tss1}, \code{tss2} (TSS position for both anchors),
##'   \code{strand1} and \code{strand2} (TSS strand for both anchors).
##' @seealso \code{\link{convertTranscriptsToTSS}},
##'   \code{\link{annotateAnchors}}
##' @export
##'
##' @examples
##' library(TxDb.Hsapiens.UCSC.hg19.knownGene)
##'
##' # take only the transcripts of coding genes by ensuring that the coding
##' # sequence strand is not NA
##' transcripts <- transcripts(TxDb.Hsapiens.UCSC.hg19.knownGene,
##'     columns = "gene_id", filter = list(cds_strand = c("-", "+"))
##' )
##'
##' fpath <- system.file("extdata", "GM12878_loops.bedpe.gz", package = "GOaGO")
##' pairs <- rtracklayer::import(fpath, genome = "hg19")
##' annotateInteractions(pairs, transcripts, maxDistanceToTSS = 10e3)
annotateInteractions <- function(
    interactions, transcripts = NULL, tss = NULL, keyType = NULL,
    maxDistanceToTSS = -1
) {
    # extract the two interacting regions
    if (inherits(interactions, "Pairs")) {
        # core Bioconductor class, used by package `rtracklayer`
        gr1 <- S4Vectors::first(interactions)
        gr2 <- S4Vectors::second(interactions)
    } else if (inherits(interactions, "GenomicInteractions")) {
        # used by packages `GenomicInteractions` and `mariner`
        if (!requireNamespace("GenomicInteractions", quietly = TRUE)) {
            stop("Package \"GenomicInteractions\" must be installed ",
                "to handle GenomicInteractions objects.",
                call. = FALSE
            )
        }
        gr1 <- GenomicInteractions::anchorOne(interactions)
        gr2 <- GenomicInteractions::anchorTwo(interactions)
    } else {
        stop("unknown class of `interactions` object")
    }

    # annotate the anchors by the overlapping or nearest TSSes
    tss <- .determine_tss(transcripts, tss)
    anchor1 <- annotateAnchors(gr1,
        tss = tss, keyType = keyType,
        maxDistanceToTSS = maxDistanceToTSS
    )
    anchor2 <- annotateAnchors(gr2,
        tss = tss, keyType = keyType,
        maxDistanceToTSS = maxDistanceToTSS
    )

    # construct a Cartesian product of annotations for both anchors
    genePairs <- merge(anchor1, anchor2,
        by = "interactionID",
        allow.cartesian = TRUE, suffixes = c("1", "2")
    )
    setkey(genePairs, NULL)

    # keep seqinfo and keyType from gene TSSes as attributes
    attr(genePairs, "seqinfo") <- seqinfo(tss)
    attr(genePairs, "keyType") <- .determine_keyType(tss, keyType)

    return(genePairs)
}
