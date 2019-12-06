## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## ------ parse ------
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
setClass(
  "surfParsing",
  contains = "DataFrame",
  representation = representation(params = "list")
)

setValidity("surfParsing", function(object) {
  TRUE
})

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## ------ drseq ------
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
setClass("drseqData",
         contains = "surfParsing",
         representation = representation(DEXSeqData = "DEXSeqDataSet"))

setValidity("drseqData", function(object) {
  TRUE
})

setClass("drseqResults",
         contains = "DataFrame",
         representation = representation(
           modelFrameBM = "DataFrameList",
           sampleData = "DataFrame",
           dispersionFunction = "List",
           params = "list"))

setValidity("drseqResults", function(object) {
  stopifnot("sample" %in% colnames(object@sampleData) || !ncol(object@sampleData))
  stopifnot(colnames(object$countData) == as.character(object@sampleData$sample))
  stopifnot(all(names(object@modelFrameBM) == names(object@dispersionFunction)))
  TRUE
})

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## ------ associate ------
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
setClass(
  "surfData",
  contains = "DataFrame",
  representation = representation(
    shRNAsample = "DataFrame",
    eCLIPsample = "DataFrame",
    params = "list")
)

setValidity("surfData", function(object) {
  stopifnot("sample" %in% colnames(object@shRNAsample) || !ncol(object@shRNAsample))
  stopifnot("sample" %in% colnames(object@eCLIPsample) || !ncol(object@eCLIPsample))
  TRUE
})

setClass(
  Class = "surfResults",
  contains = "DataFrame",
  representation = representation(
    trainData = "surfData",
    params = "list"
  )
)

setValidity("surfResults", function(object) {
  stopifnot(all(object$size >= object@params$min.size))
  stopifnot(object@params$trim < 0.5)
  TRUE
})

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## ------ daseq ------
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# @importClassesFrom AUCell aucellResults

aucellResults <- AUCell::aucellResults ## temp hack

daseqRanking <- setClass(
  "daseqRanking",
  contains = "aucellResults",
  representation = representation()
)

#' DASeq results
setClass(
  "daseqResults",
  contains = "DataFrame",
  representation = representation(
    targetAUC = "SummarizedExperiment",
    params = "list")
)

setValidity("daseqResults", function(object) {
  stopifnot("condition" %in% colnames(colData(object@targetAUC)))
  stopifnot(object@params$n.resample > 0)
  TRUE
})

