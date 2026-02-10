# Internal function to select only the terms that are associated to at least
# `minTermPairs` gene pairs, and convert the terms (both `ID` and `Description`
# columns) to factors ordered by the column `orderBy`. Finally, select only the
# number of terms (or the terms themselves) specified by `showCategory`.
.sortedResult <- function(
    object,
    minTermPairs = 5,
    showCategory = 10,
    orderBy = "FoldEnrichment",
    decreasing = TRUE
) {
    dt <- as.data.table(object)
    dt <- dt[Count >= minTermPairs, ]

    if (!orderBy %in% colnames(dt)) {
        message("wrong orderBy parameter: ", orderBy)
    }

    # note the reverse order, so that the terms show up as expected
    # (top-to-bottom, reverse to the scale direction) in the final plot
    ind <- order(dt[[orderBy]], decreasing = !decreasing)
    dt[, ID := factor(ID, ID[ind])]
    dt[, Description := factor(Description, Description[ind])]

    if (is.numeric(showCategory)) {
        dt <- head(dt, showCategory)
    } else {
        dt <- dt[ID %in% showCategory, ]
    }

    return(dt)
}

# Internal function to wrap labels at a given wrap length, or use a labeller
# function provided.
.label_format <- function(label_format = 50) {
    if (is.function(label_format)) {
        label_func <- label_format
    } else {
        label_func <- label_wrap_gen(label_format)
    }
    return(label_func)
}


##' Dotplot of the enriched Gene Ontology terms
##'
##' By default, plots the most enriched terms, with fold enrichment on the
##' X-axis, point size indicating the number of gene pairs sharing the given
##' term, and point color -- the adjusted p-value.
##'
##' @param object GO-a-GO results of class \code{GOaGO-result}
##' @param minTermPairs plot only the GO terms that are associated to at least
##'   the given number of gene pairs
##' @param x Variable for X-axis, one of \code{"FoldEnrichment"}, \code{"Ratio"}
##'   and \code{"Count"}.
##' @param color Variable used to color enriched terms, e.g. \code{"pvalue"},
##'   \code{"p.adjust"} or \code{"qvalue"}.
##' @param size Variable used to scale the sizes of points, one of
##'   \code{"FoldEnrichment"}, \code{"Ratio"} and \code{"Count"}.
##' @param showCategory number of terms to display or a vector of terms
##' @param orderBy The order of the Y-axis, one of \code{"FoldEnrichment"},
##'   \code{"Ratio"} and \code{"Count"}.
##' @param decreasing logical. Should the \code{orderBy} order be increasing or
##'   decreasing?
##' @param font.size font size
##' @param label_format a numeric value sets wrap length, alternatively a custom
##'   function to format axis labels. By default wraps names longer than 50
##'   characters.
##' @returns A \code{ggplot} object that can be further customized using the
##'   \code{ggplot2} package.
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
##' DotPlot(goago)
DotPlot <- function(
    object,
    minTermPairs = 5,
    x = "FoldEnrichment",
    color = "p.adjust",
    size = "Count",
    showCategory = 10,
    orderBy = "FoldEnrichment",
    decreasing = TRUE,
    font.size = 12,
    label_format = 50
) {
    dt <- .sortedResult(object,
        minTermPairs = minTermPairs,
        showCategory = showCategory, orderBy = orderBy, decreasing = decreasing
    )

    label_func <- .label_format(label_format)

    p <- ggplot(
        dt, aes(
            x = .data[[x]], y = Description, size = .data[[size]],
            color = .data[[color]]
        )
    ) +
        geom_point() +
        labs(y = NULL) +
        scale_x_continuous(expand = expansion(mult = 0.1)) +
        scale_y_discrete(labels = label_func) +
        theme_dose(font.size)

    return(p)
}


##' Ridgeplot of the sampling distributions for the randomized gene pairs
##'
##' Ridgeplot of the sampling distributions of numbers of gene pairs sharing
##' each enriched Gene Ontology term, obtained for the randomized gene pairs.
##'
##' @param object GO-a-GO results of class \code{GOaGO-result}
##' @param minTermPairs plot only the GO terms that are associated to at least
##'   the given number of gene pairs
##' @param showCategory number of terms to display or a vector of terms
##' @param orderBy The order of the Y-axis, one of \code{"FoldEnrichment"},
##'   \code{"Ratio"} and \code{"Count"}.
##' @param decreasing logical. Should the \code{orderBy} order be increasing or
##'   decreasing?
##' @param font.size font size
##' @param label_format a numeric value sets wrap length, alternatively a custom
##'   function to format axis labels. By default wraps names longer than 50
##'   characters.
##' @returns A \code{ggplot} object that can be further customized using the
##'   \code{ggplot2} package.
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
##' RidgePlot(goago)
RidgePlot <- function(
    object,
    minTermPairs = 5,
    showCategory = 10,
    orderBy = "FoldEnrichment",
    decreasing = TRUE,
    font.size = 12,
    label_format = 50
) {
    dt <- .sortedResult(object,
        minTermPairs = minTermPairs,
        showCategory = showCategory, orderBy = orderBy, decreasing = decreasing
    )

    dt_permuted <- merge(dt[, c("ID", "Description")],
        object@permutedResult,
        by = "ID"
    )

    label_func <- .label_format(label_format)

    col_color <- c(`TRUE` = "#e41a1c", `FALSE` = "black")
    col_fill <- c(`TRUE` = NA, `FALSE` = "gray70")
    lab <- c(`TRUE` = "actual gene pairs", `FALSE` = "randomized gene pairs")

    p <- ggplot(
        dt_permuted,
        aes(x = Count, y = Description, color = FALSE, fill = FALSE)
    ) +
        geom_density_ridges() +
        geom_point(data = dt, aes(x = Count, color = TRUE, fill = TRUE)) +
        scale_color_manual(NULL, values = col_color, labels = lab) +
        scale_fill_manual(NULL, values = col_fill, labels = lab) +
        labs(x = "Number of gene pairs sharing a term", y = NULL) +
        scale_y_discrete(labels = label_func) +
        theme_dose(font.size) +
        theme(legend.position = "bottom")

    return(p)
}
