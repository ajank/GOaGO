library(data.table)
library(GOaGO)
library(rtracklayer)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)


#
#  Download from the Gene Expression Omnibus (GEO) repository a dataset
#  of 9,448 chromatin loops identified in human cell line GM12878,
#  and save it as a BEDPE file in the `inst/extdata` directory.
#

url <- paste0(
    "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE63nnn/GSE63525/suppl/",
    "GSE63525_GM12878_primary%2Breplicate_HiCCUPS_looplist.txt.gz"
)

dt <- fread(url, sep = "\t", header = TRUE)
dt[, chr1 := paste0("chr", chr1)]
dt[, chr2 := paste0("chr", chr2)]

si <- GenomeInfoDb::Seqinfo(genome = "hg19")
gr1 <- makeGRangesFromDataFrame(dt,
    seqinfo = si, seqnames.field = "chr1", start.field = "x1", end.field = "x2",
    ignore.strand = TRUE, starts.in.df.are.0based = TRUE
)
gr2 <- makeGRangesFromDataFrame(dt,
    seqinfo = si, seqnames.field = "chr2", start.field = "y1", end.field = "y2",
    ignore.strand = TRUE, starts.in.df.are.0based = TRUE
)
pairs <- Pairs(gr1, gr2)

export(pairs, gzfile(file.path("..", "extdata", "GM12878_loops.bedpe.gz"), compression = 9))


#
#  Annotate the loops using gene TSSes from TxDb.Hsapiens.UCSC.hg19.knownGene,
#  and save the resulting data table with gene pairs in the `data` directory.
#

# take only the transcripts of coding genes by ensuring that the coding sequence strand is not NA
transcripts <- transcripts(TxDb.Hsapiens.UCSC.hg19.knownGene,
    columns = "gene_id", filter = list(cds_strand = c("-", "+"))
)

genePairsGM12878 <- annotateInteractions(pairs, transcripts, maxDistanceToTSS = 10e3)

save(genePairsGM12878,
    file = file.path("..", "..", "data", "genePairsGM12878.rda"),
    version = 2, compress = "xz", compression_level = 9
)
