# sort by enrichment ratio, and impose the same ordering on the label levels

sortedResult <- function(goago) {
    dt <- as.data.table(goago@result)
    dt[, label := Description]
    dt <- dt[order(Ratio/BgRatio), ]
    dt[, label := factor(label, unique(label))]
}


##' Dotplot of the overrepresented Gene Ontology terms
##'
##' The x-axis shows the log~2~ fold change of the fraction of gene pairs
##' sharing the given terms. Color indicates the adjusted p-value.
##'
##' @param goago GO-a-GO results of class \code{GOaGO-result}
##' @param minTermPairs take only the GO terms that are associated to at least
##'   the given number of gene pairs
##' @returns A ggplot object that can be further customized using the
##'   \code{ggplot2} package.
##' @export

DotPlot <- function(goago, minTermPairs = 5) {
    dt <- sortedResult(goago)
    dt <- dt[Count >= minTermPairs, ]

    p <- ggplot(dt, aes(x=log2(Ratio/BgRatio), y=label, size=Count, color=p.adjust)) +
        geom_point() +
        ylab(NULL) +
        scale_size(range=c(3, 8))
    return(p)
}


##' Ridgeplot of the sampling distributions for the randomized gene pairs
##'
##' Ridgeplot of the sampling distributions of numbers of gene pairs sharing
##' each Gene Ontology term, obtained for the randomized gene pairs.
##'
##' @param goago GO-a-GO results of class \code{GOaGO-result}
##' @param minTermPairs take only the GO terms that are associated to at least
##'   the given number of gene pairs
##' @returns A ggplot object that can be further customized using the
##'   \code{ggplot2} package.
##' @export

RidgePlot <- function(goago, minTermPairs = 5) {
    dt <- sortedResult(goago)
    dt <- dt[Count >= minTermPairs, ]

    dt2 <- merge(dt[, c("ID", "label")], goago@permutedResult, by = "ID")

    p <- ggplot(dt2, aes(x = Count, y = label)) +
        geom_density_ridges() +
        geom_point(data = dt, aes(color = p.adjust < 0.05)) +
        scale_color_manual(NULL, values = c(`TRUE` = "#e41a1c", `FALSE` = "#666666"),
                           labels = c(`TRUE` = "Number of actual gene pairs\nsharing a Gene Ontology term", `FALSE` = "not significant")) +
        xlab("Number of randomized gene pairs\nsharing a Gene Ontology term") +
        ylab(NULL) +
        guides(color=guide_legend(ncol=1)) +
        theme(legend.position="bottom") +
        NULL
    return(p)
}
