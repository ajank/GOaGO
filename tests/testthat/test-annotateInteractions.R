library(TxDb.Hsapiens.UCSC.hg19.knownGene)

# take only the transcripts of coding genes by ensuring that the coding
# sequence strand is not NA
transcripts <- transcripts(TxDb.Hsapiens.UCSC.hg19.knownGene,
    columns = "gene_id", filter = list(cds_strand = c("-", "+"))
)

fpath <- system.file("extdata", "GM12878_loops.bedpe.gz", package = "GOaGO")
pairs <- rtracklayer::import(fpath, genome = "hg19")


test_that("`convertTranscriptsToTSS` should return a GRanges object of width 1
with gene identifiers in metadata column `geneID` being a character vector", {
    tss <- convertTranscriptsToTSS(transcripts)
    expect_s4_class(tss, "GRanges")
    expect_setequal(width(tss), 1)
    expect_type(tss$geneID, "character")
})

test_that("`convertTranscriptsToTSS` should be idempotent", {
    tss <- convertTranscriptsToTSS(transcripts)
    expect_equal(tss, convertTranscriptsToTSS(tss))
})


# synthetic transcripts
st <- GRanges("chr1",
    IRanges(c(11, 11, 21, 51), c(20, 40, 40, 70)),
    strand = c("+", "+", "+", "-"),
    gene_id = CharacterList("A", c("A", "B"), c("A", "B"), "C")
)

st2 <- GRanges("chr1",
    IRanges(c(11, 11, 21, 51), c(20, 40, 40, 70)),
    strand = c("+", "+", "+", "-"),
    gene_id = c("A", "A", NA, "C")
)

test_that("`convertTranscriptsToTSS` works as expected on synthetic data", {
    expect_equal(
        convertTranscriptsToTSS(st),
        GRanges("chr1",
            IRanges(c(11, 11, 21, 21, 70)),
            strand = c("+", "+", "+", "+", "-"),
            geneID = c("A", "B", "A", "B", "C")
        )
    )
    expect_equal(
        convertTranscriptsToTSS(st2),
        GRanges("chr1",
            IRanges(c(11, 70)),
            strand = c("+", "-"),
            geneID = c("A", "C")
        )
    )
})


# synthetic anchors: the first one does not overlap a TSS, but is equidistant
# to three; the second one overlaps a TSS; the third one is adjacent to a TSS
sa <- GRanges("chr1", IRanges(c(31, 61, 71), c(60, 90, 100)))

test_that("`annotateAnchors` works as expected on synthetic data", {
    dt <- data.table(
        interactionID = 2, chrom = "chr1", start = 61, end = 90,
        geneID = "C", tss = 70, strand = "-"
    )
    attr(dt, "seqinfo") <- seqinfo(st)
    attr(dt, "keyType") <- "ENTREZID"
    expect_equal(annotateAnchors(sa, st, keyType = "ENTREZID"), dt)

    dt2 <- data.table(
        interactionID = c(1, 1, 1, 2, 3), chrom = "chr1",
        start = c(31, 31, 31, 61, 71), end = c(60, 60, 60, 90, 100),
        geneID = c("A", "B", "C", "C", "C"),
        tss = c(21, 21, 70, 70, 70),
        strand = c("+", "+", "-", "-", "-")
    )
    attr(dt2, "seqinfo") <- seqinfo(st)
    attr(dt2, "keyType") <- "ENTREZID"
    expect_equal(
        annotateAnchors(sa, st, keyType = "ENTREZID", maxDistanceToTSS = 100),
        dt2
    )
})

test_that("`annotateInteractions` works as expected on synthetic data", {
    dt <- data.table(
        interactionID = c(1, 1, 1), chrom1 = "chr1", start1 = c(31, 31, 31),
        end1 = c(60, 60, 60), geneID1 = c("A", "B", "C"), tss1 = c(21, 21, 70),
        strand1 = c("+", "+", "-"), chrom2 = "chr1", start2 = 61, end2 = 90,
        geneID2 = "C", tss2 = 70, strand2 = "-"
    )
    attr(dt, "seqinfo") <- seqinfo(st)
    attr(dt, "keyType") <- "ENTREZID"

    expect_equal(
        annotateInteractions(Pairs(sa[1], sa[2]), st,
            keyType = "ENTREZID",
            maxDistanceToTSS = 100
        ),
        dt
    )
})


test_that("`annotateInteractions` should return the same result for `Pairs` and
`GenomicInteractions` objects", {
    result1 <- annotateInteractions(pairs, transcripts, maxDistanceToTSS = 10e3)

    gi <- GenomicInteractions::GenomicInteractions(
        S4Vectors::first(pairs),
        S4Vectors::second(pairs)
    )
    result2 <- annotateInteractions(gi, transcripts, maxDistanceToTSS = 10e3)
    expect_equal(result1, result2)
})
