## Coercion methods

`as.data.frame.GOaGO-result` <- function(x, ...)  {
    data.frame(x@result, ...)
}

setAs("GOaGO-result", "data.frame", function(from)
    `as.data.frame.GOaGO-result`(from))

setMethod("as.data.frame", signature(x = "GOaGO-result"),
    `as.data.frame.GOaGO-result`
)


`as.data.table.GOaGO-result` <- function(x, ...)  {
    data.table(x@result, ...)
}

setAs("GOaGO-result", "data.table", function(from)
    `as.data.table.GOaGO-result`(from))

setMethod("as.data.table", signature(x = "GOaGO-result"),
    `as.data.table.GOaGO-result`
)


## Accessors

setGeneric("keyType", function(object) {
  standardGeneric("keyType")
})

setMethod("keyType", signature(object = "GOaGO-result"),
    function(object) {
        object@keyType
    }
)

# defined in package BiocGenerics
setMethod("organism", signature(object = "GOaGO-result"),
    function(object) {
        object@organism
    }
)
