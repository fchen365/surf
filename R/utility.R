## what: various helper functions 
## who: fan chen (fan.chen@wisc.edu) 
## when: 08/12/2018 
## ------------------------------------------------------------

suppressPackageStartupMessages({
  library(tidyverse)
  library(doParallel)
  library(DEXSeq)
})

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## ------ global variable ------
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## SURF pre-defined 
surf.events <- c("SE", "RI", "A3SS", "A5SS", "AFE", "A5U", "IAP", "TAP")
surf.colors <- c('#fb9a99','#e31a1c','#fdbf6f','#ff7f00','#b2df8a','#33a02c','#a6cee3','#1f78b4')
surf.features <- c("up3", "up2", "up1", "bd1", "bd2", "dn1", "dn2", "dn3")
greek.features <- c(bquote(alpha), bquote(beta), bquote(gamma), bquote(delta), bquote(epsilon), bquote(zeta), bquote(eta), bquote(theta))


## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## ------ class ------
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#' DrSeq class
#' 
#' \code{drseqResults} is a  stand-alone object of DrSeq (the analysis module 1 of SURF) results.  
setClass(
  "drseqResults",
  contains = "DataFrame",
  representation = representation( 
    modelFrameBM = "DataFrameList", 
    dispersionFunction = "List"
  )
)

setValidity(
  "drseqResults", function(object) {
    stopifnot(all(names(object@modelFrameBM) == names(object@dispersionFunction)))
    TRUE
  }
)

#' FASeq class
#' 
#' `faseqResults` is a stand-alone object of FASeq (the analysis module 2 of SURF) results.  
setClass(
  Class = "faseqResults",
  contains = "DataFrame"
)

setValidity("faseqResults", function(object) {
  stopifnot("size" %in% colnames(object))
  stopifnot("min.size" %in% names(metadata(object)))
  stopifnot("trim" %in% names(metadata(object)))
  stopifnot(all(object$size >= metadata(object)$min.size))
  stopifnot(metadata(object)$trim < 0.5)
  TRUE
})

#' DASeq class
#' 
#' `daseqResults` is a stand-alone object of DASeq (the discovery module of SURF) results.  
setClass(
  "daseqResults",
  contains = "DataFrame",
  representation = representation(
    AUC = "SummarizedExperiment"
  )
)

setValidity("daseqResults", function(object) {
  # stopifnot(!!nrow(object@rankings))
  stopifnot("condition" %in% names(colData(object@AUC)))
  stopifnot(all(rownames(object) == rownames(object@AUC)))
  TRUE
})

#' SURF class
#' 
#' The \code{surf} class is an all-in-one analytic object used in SURF framework. 
#' A \code{surf} object contains results of DrSeq (analysis module 1), FASeq (analysis module 2), and DASeq (discovery module 1).
#' In addition, it also contain gene (isoform) parts list, parsing from genome annotation. 
#' The gene parts list is essential for ATR event construction. 
#' If the any analysis module or discovery module is performed, the sample metadata will also be recorded. 
#' Each of the components mentioned above can be accessed through a function listed below.
#' 
#' @seealso \link{genePartsList}, \link{drseqResults}, \link{faseqResults}, \link{daseqResults}, \link{sampleData}.
#' @references Chen, F., & Keles, S. (2020). SURF: integrative analysis of a compendium of RNA-seq and CLIP-seq datasets highlights complex governing of alternative transcriptional regulation by RNA-binding proteins. *Genome Biology*, 21(1), 1-23.
setClass(
  "surf",
  contains = "DataFrame",
  representation = representation(
    genePartsList = "DataFrame",
    drseqData = "DEXSeqDataSet", 
    drseqResults = "drseqResults", 
    faseqData = "RangedSummarizedExperiment", 
    faseqResults = "faseqResults", 
    daseqResults = "daseqResults",
    sampleData = "DataFrameList"
  )
)

setValidity("surf", function(object) {
  ## main columns and @genePartsList
  stopifnot(all(object$event_name %in% surf.events))
  stopifnot(all(object$gene_id %in% object@genePartsList$gene_id)) 
  stopifnot(all(object$transcript_id %in% unlist(object@genePartsList$transcript_id))) 
  
  ## @sampleData colnames
  for (nm in lapply(object@sampleData, colnames)) {
    stopifnot(all(c("sample", "condition") %in% nm))
  }
  
  ## RNA-seq @sampleData rownames
  # stopifnot(rownames(object@sampleData$"RNA-seq") == 
  #             colnames(object@drseqData))
  
  ## CLIP-seq @sampleData rownames
  # stopifnot(rownames(object@sampleData$"CLIP-seq") == 
  #             colnames(object@faseqData))
  
  TRUE
})


## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## ------ accessor ------
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## ------ _ genePartsList ------
#' FASeq results
#' @param object a \code{surf} object output by \link{faseq} or \link{faseqFit} or \link{faseqInference}.
#' @return a \code{DataFrame} object.
setGeneric(
  name = "genePartsList",
  def = function(object)
    standardGeneric("genePartsList")
)

#' @rdname surf-class
#' @exportMethod genePartsList
setMethod(
  "genePartsList", "surf",
  function(object) {
    object@genePartsList
  }
)


## ------ _ **seqResults ------

#' DrSeq Results
#' @param object a \code{surf} object output by \link{drseq} or \link{drseqFit} or \link{drseqFilter}.
#' @return a \code{drseqResults} object
setGeneric(
  name = "drseqResults",
  def = function(object)
    standardGeneric("drseqResults")
)

#' @rdname surf-class
#' @exportMethod drseqResults
setMethod(
  "drseqResults", "surf",
  function(object) {
    object@drseqResults[object$event_id, ]
  }
)

#' FASeq Results
#' @param object a \code{surf} object output by \link{faseq} or \link{faseqFit} or \link{faseqInference}.
#' @return a \code{faseqResults} object
setGeneric(
  name = "faseqResults",
  def = function(object)
    standardGeneric("faseqResults")
)

#' @rdname surf-class
#' @exportMethod faseqResults
setMethod(
  "faseqResults", "surf",
  function(object) {
    object@faseqResults
  }
)

#' DASeq Results
#' @param object a \code{surf} object output by \link{daseq}.
#' @return a \code{daseqResults} object.
setGeneric(
  name = "daseqResults",
  def = function(object)
    standardGeneric("daseqResults")
)

#' @rdname surf-class
#' @exportMethod daseqResults
setMethod(
  "daseqResults", "surf",
  function(object) {
    object@daseqResults
  }
)

## ------ _ sampleData ------
#' Sample Data 
#' @name sampleData
#' @param type \code{character}, any subset of "RNA-seq", "CLIP-seq", and "External".
#' @return \code{DataFrame} or \code{DataFrameList} (if \code{length(type)>1}).
#' @export sampleData 
setGeneric(
  name = "sampleData",
  def = function(object, ...)
    standardGeneric("sampleData")
)

#' @rdname SummarizedExperiment-class
#' @exportMethod sampleData
setMethod(
  "sampleData", "SummarizedExperiment",
  function(object) {
    colData(object)
  }
)

#' @rdname daseqResults-class
#' @exportMethod sampleData
setMethod(
  "sampleData", "daseqResults",
  function(object) {
    sampleData(object@AUC)
  }
)

#' @rdname surf-class
#' @param type \code{character}, any subset of "RNA-seq", "CLIP-seq", and "External".
#' @exportMethod sampleData
setMethod(
  "sampleData", "surf", 
  function(object, type = NULL) {
    notfound <- setdiff(type, names(object@sampleData))
    if (length(notfound)) {
      stop("Sample data of ", 
           paste(notfound, collapse = ", "), 
           " not found.")
    }
    if (is.null(type)) type = names(object@sampleData)
    if (length(type) > 1) {
      object@sampleData[type]
    } else if (length(type) == 1) {
      object@sampleData[[type]]
    } 
  }
)

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## ------ methods ------
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## ------ _ subsetByOverlaps ------
setMethod(
  "subsetByOverlaps", 
  signature(x = "surf",
            ranges = "GenomicRanges"),
  function(x, 
           ranges,
           maxgap = 0L,
           minoverlap = 1L,
           type = c("any", "start", "end", "within", "equal"),
           ignore.strand = FALSE) {
    genomicData <- x$genomicData
    overlaps <- findOverlaps(
      query = genomicData,
      subject = ranges,
      maxgap = maxgap,
      minoverlap = minoverlap,
      type = type,
      ignore.strand = ignore.strand
    )
    query[queryHits(overlaps),]
  }
)

## ------ _ findOverlaps ------
setMethod(
  "findOverlaps", 
  signature(query = "surf", 
            subject = "GenomicRanges"),
  function(query,
           subject,
           maxgap = 0L,
           minoverlap = 1L,
           type = c("any", "start", "end", "within", "equal"),
           ignore.strand = FALSE){
    genomicData <- query$genomicData
    overlaps <- findOverlaps(
      query = genomicData,
      subject = subject,
      maxgap = maxgap,
      minoverlap = minoverlap,
      type = type,
      ignore.strand = ignore.strand
    )
    overlaps
  }
)


## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## ------ methods (drseq) ------
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## ------ _ plotDispFunc ------
#' Plot Dispersion Functions
#' 
#' Plot the fitted mean-dispersion functions of DrSeq models. 
#' DrSeq fits up to eight mean-dispersion functions; each corresponds to one ATR event type.
#' This allows DrSeq to better account for the nuance in the over-dispersion presented by different ATR event types.
#' For more details, please refer to SURF paper.
#' 
#' @param object a \code{drseqResults} object.
#' @return a \code{ggplot} object.
#' @details By SURF default, ATR type are colored with the `Paired` palette (see \url{http://colorbrewer2.org}).
setGeneric(
  name = "plotDispFunc",
  def = function(object, ...)
    standardGeneric("plotDispFunc")
)

#' @rdname plotDispFunc
#' @param x.limits \code{numeric(2)}, limits for x-axis.
#' @exportMethod plotDispFunc
setMethod(
  "plotDispFunc",
  "drseqResults",
  function(object, x.limits = c(1e-1, 1e4)) {
    colrs = setNames(surf.colors, surf.events) ## surf colors
    func <- object@dispersionFunction
    g <- ggplot(data.frame(x = 0), aes(x = x))
    for (event in names(func)) {
      g <- g +
        stat_function(fun = func[[event]], color = colrs[event])
    }
    g <- g +
      scale_x_continuous(limits = x.limits, trans = "log10") +
      scale_y_continuous(trans = "log10") +
      scale_color_manual(values = colrs) +
      labs(x = "normalized mean", y = "estimated dispersion",
           color = 'event') +
      theme_bw()
    return(g)
  }
)

#' @rdname plotDispFunc
#' @exportMethod plotDispFunc
setMethod("plotDispFunc", "surf",
          function(object,...) {
            plotDispFunc(drseqResults(object))
          })


## ------ _ Volcano plot ------
#' Volcano plot
#' 
#' Create a volcano plot for the DrSeq results, stratefied by alternative transcriptional regulation (ATR) event types.
#' A volcano plot is a scatter plot of tested units, where log2 fold change is in x-axis, and -log10(p.value) is in y-axis.
#' 
#' @param object a \code{drseqResults} object.
#' @return a \code{ggplot} object.
#' @details By default, ATR type are colored with the `Paired` palette (see \url{http://colorbrewer2.org}).
setGeneric(
  name = "volcano.plot",
  def = function(object, ...)
    standardGeneric("volcano.plot")
)

#' @rdname volcano.plot
#' @param lfc.cutoff \code{numeric(2)}, the range of log2 fold change that is consider NOT significant.
#' @param fdr.cutoff \code{numeric(1)}, significance level of adjusted p-value.
#' @param x.limits \code{numeric(2)}, range of log2 fold change. Any values beyond this range will be projected onto the boundary.
#' @param y.limits \code{numeric(2)}, range of -log10(p.value). Any values beyond this range will be projected onto the boundary.
#' @param remove.portion.grey \code{numeric}, between 0 and 1, the portion of non-significant points to be randomly remove. This is only for speeding up plotting.
#' @param remove.portion.color \code{numeric}, between 0 and 1, the portion of significant points to be randomly remove. This is only for speeding up plotting.
#' @exportMethod volcano.plot
setMethod(
  "volcano.plot",
  "drseqResults",
  function(object, 
           lfc.cutoff = c(-1, 1), 
           fdr.cutoff = 0.01, 
           x.limits = c(-15, 15), 
           y.limits = c(0, 50), 
           remove.portion.grey = 0,
           remove.portion.color = 0,
           colrs = setNames(surf.colors, surf.events)) {
    stopifnot(length(lfc.cutoff) == 2 && lfc.cutoff[1] <= lfc.cutoff[2])
    stopifnot(fdr.cutoff > 0 && fdr.cutoff < 1)
    stopifnot(length(x.limits) == 2 && x.limits[1] <= x.limits[2])
    stopifnot(length(y.limits) == 2 && y.limits[1] <= y.limits[2])
    fdr.cutoff = -log10(fdr.cutoff)
    dat <- data.frame(x = object$logFoldChange, 
                      y = -log10(object$padj), 
                      group = object$event_name) %>%
      dplyr::filter(!is.na(x), !is.na(y)) %>%
      mutate(x = pmax(x, x.limits[1]),
             x = pmin(x, x.limits[2]),
             y = pmin(y, y.limits[2]),
             color = ifelse(x > lfc.cutoff[1] & x < lfc.cutoff[2] | 
                              y < fdr.cutoff, "No Sig.", as.character(group)), 
             color = factor(color, c(surf.events, "No Sig.")),
             size = ifelse(color == "No Sig.", 1, 2))
    ## remove some "No Sig." to reduce plot size
    if (remove.portion.grey < 1 && remove.portion.grey > 0) {
      index <- which(dat$color == "No Sig.")
      remove.index = sample(index, round(length(index) * remove.portion.grey))
      dat <- dat[-remove.index,]
    }
    if (remove.portion.color < 1 && remove.portion.color > 0) {
      index <- which(dat$color != "No Sig.")
      remove.index = sample(index, round(length(index) * remove.portion.color))
      dat <- dat[-remove.index,]
    }
    g <- ggplot(dat, aes(x, y, color = color, size = size)) + 
      geom_vline(xintercept = c(lfc.cutoff[1], lfc.cutoff[2]), 
                 color = "grey40", linetype = 2, alpha = .9) + 
      geom_hline(yintercept = fdr.cutoff, color = "grey40", linetype = 2, alpha = .9) +
      geom_point(alpha = .7) +
      scale_color_manual(values = c(colrs, "No Sig." = "grey60")) + 
      scale_size(range = c(.1, .7)) + 
      labs(x = "log"[2]~"(fold change), knock-down vs. wild-type", 
           y = "-log"[10]~"(adjusted p value)") +
      guides(size = "none") + 
      scale_x_continuous(limits = x.limits) + 
      scale_y_continuous(limits = y.limits) + 
      facet_wrap(~ group, nrow = 2) +
      theme_bw()
    return(g)
  }
)

#' @rdname volcano.plot
#' @exportMethod volcano.plot
setMethod("volcano.plot", "surf",
          function(object,...) {
            volcano.plot(drseqResults(object))
          })

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## ------ methods (faseq) ------
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#' Functional association plot
#' 
#' This function provides a ggplot implementation for functional association (FA) plot.
#' 
#' @param object a \code{surf} object from \link{faseq} or \link{faseqFit} or \link{faseqInfer}.
#' @param plot.event \code{character} vector, event type wanted. In particular, "all" means all event types.
#' @param fdr.cutoff \code{numeric}, theshold for adjusted p-value in functional association test.
#' @return a \code{ggplot} object.
#' @export
fa.plot <- function(object, 
                    plot.event = "all", 
                    trim = metadata(object)$trim,
                    fdr.cutoff = 0.05){
  stopifnot(is(object, "surf"))
  if (any(plot.event == "all")) plot.event = levels(object$event)
  levels <- outer(surf.events, surf.features, paste) %>% t %>% as.vector
  
  ## box plot: feature signals
  trainData <- object[object$included, ]
  groupSize <- table(trainData$event_name, trainData$group) %>% 
    apply(1, paste, collapse = ", ")
  trainData <- trainData[trainData$event_name %in% plot.event,]
  featureSignal <- trainData$featureSignal
  dat1 = data.frame(
    event = rep(trainData$event_name, elementNROWS(featureSignal)),
    feature = unlist(lapply(featureSignal, names), use.names = F),
    group = rep(trainData$group, elementNROWS(featureSignal)),
    signal = unlist(featureSignal, use.names = F)) %>% 
    group_by(event, feature) %>%
    mutate(lower = quantile(signal, trim * 2), 
           upper = quantile(signal, 1 - trim * 2),
           x = factor(feature, surf.features)) %>% 
    dplyr::filter(signal > lower, signal < upper) %>% 
    ungroup %>% 
    mutate(strip = paste0(event, " (", groupSize[event], ")") %>% 
             factor(paste0(names(groupSize), " (", groupSize, ")")))
  g1 <- ggplot(dat1, aes(x, signal, fill = group)) +
    geom_boxplot(color = "grey30", alpha = .9, 
                 outlier.shape = ".", outlier.color = "grey80") + 
    labs(y = "feature signal") + 
    scale_fill_manual(values = c("increase" = "#4d9221", "no change" = "grey50","decrease" = "#c51b7d")) + 
    scale_x_discrete(breaks = surf.features, labels = greek.features) +
    facet_grid(cols = vars(strip), scales = "free_x", space = "free_x") + 
    theme_bw() + 
    theme(axis.title.x = element_blank(), 
          axis.text.x = element_blank(), 
          axis.ticks.x = element_blank())
  
  ## dot-line plot: adjusted p-values
  suppressWarnings({
    dat2 <- dat1 %>% 
      group_by(event, feature) %>% 
      summarise(n = n()) %>%
      left_join(as.data.frame(object), by = c("event", "feature")) %>% 
      ungroup() %>%
      mutate(
        event = factor(event, surf.events),
        logp = - log10(padj), 
        x = factor(feature, surf.features), 
        functional = ifelse(is.na(logp), "not tested", 
                            as.character(functional)) %>%
          factor(c("exclusion", "inclusion", "not tested")), 
        logp = ifelse(is.na(logp), 0, logp)
      ) 
  })
  g2 <- ggplot(dat2, aes(x, logp, color = functional)) +
    geom_hline(yintercept = -log10(fdr.cutoff), color = "grey40", alpha = .9, linetype = 2, show.legend = T) + 
    geom_point(alpha = .9) +
    stat_summary(aes(group = functional), fun.y = sum, geom = "line", alpha = .8) +
    labs(x = "location feature", y = "-log"[10]~"(adjusted p value)", color = "association") +
    scale_color_manual(values = c(exclusion = "#4d9221", inclusion = "#c51b7d", "not tested" = "grey40")) + 
    scale_x_discrete(breaks = surf.features, labels = greek.features) +
    facet_grid(cols = vars(event), scales = "free_x", space = "free_x") + 
    theme_bw() + 
    theme(strip.background = element_blank(), 
          strip.text = element_blank())
  
  g <- ggpubr::ggarrange(g1, g2, ncol = 1, align = "v")
  # ggsave("fa.plot.pdf", g, width = 8, height = 5)
  return(g)
}

# setMethod("fa.plot", signature(object = "surf"), fa.plot)


#' Extract SURF-inferred location features
#'  
#' @param object a \code{surf} object from \link{faseq} or \link{faseqInfer}.
#' @return a \code{GRanges} object of all SURF-inferred location features.
#' @export
inferredFeature = function(object) {
  stopifnot(is(object, "surf"))
  ## check
  stopifnot(all(c("feature", "inferredFeature") %in% colnames(object)))
  
  feature <- c_granges(object$feature, sep = "-")
  n_feature <- elementNROWS(object$feature)
  feature$event_id <- rep(object$event_id, n_feature)
  feature$event_name <- rep(object$event_name, n_feature)
  feature$gene_id <- rep(object$gene_id, n_feature)
  feature$transcript_id <- rep(object$transcript_id, n_feature)
  feature$feature_name <- factor(names(unlist(object$feature, use.names = F)), surf.features)
  feature$functional <- unlist(object$inferredFeature, use.names = F)
  feature <- feature[feature$functional != "none"]
  return(feature)
}

# setMethod("inferredFeature", signature(object = "surf"), inferredFeature)


## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## ------ methods (daseq) ------
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## ------ _ getAUC ------
#' Get AUC
#' 
#' Get AUC measure for each target set (row) in every sample (column).
#' @return a \code{matrix} of AUC, whose rows correspond to target sets and columns correspond to samples.
setGeneric(
  name = "getAUC",
  def = function(object)
    standardGeneric("getAUC")
)

#' @rdname getAUC
#' @exportMethod getAUC
setMethod(
  "getAUC", "SummarizedExperiment",
  function(object) {
    if("AUC" %in% assayNames(object)) {
      SummarizedExperiment::assays(object)[["AUC"]]
    }else{
      stop("Cannot find the 'AUC' assay")
    }
  }
)

#' @rdname getAUC
#' @exportMethod getAUC
setMethod(
  "getAUC", "daseqResults",
  function(object) {
    getAUC(object@AUC)
  }
)

#' @rdname getAUC
#' @exportMethod getAUC
setMethod(
  "getAUC", "surf",
  function(object) {
    getAUC(object@daseqResults)
  }
)


## ------ _ heatmapAUC ------
#' Heatmap of AUC score
#' 
#' This function produces the heatmap of AUC scores.
#' 
#' @param object for \code{surf} object, it should be output from \link{daseq}.
#' @return a \code{ggplot} object.
setGeneric(
  name = "heatmapAUC",
  def = function(object, ...)
    standardGeneric("heatmapAUC")
)


#' @description Row blocking is allowed by providing the \code{group} parameter.
#' @rdname heatmapAUC
#' @param group an optional \code{factor} vector, group of rows.
#' @exportMethod heatmapAUC
setMethod(
  "heatmapAUC", "SummarizedExperiment",
  function(object, group = NULL) {
    
    mat <- getAUC(object)
    sampleData <- sampleData(object)
    
    ## check group input
    if (is.null(group)) {
      group = rep("set", nrow(object))
    } else {
      stopifnot(length(group) != nrow(object)) 
    }
    names(group) <- rownames(mat)
    
    ## cluster rows and columns
    set_cluster = clusterByGroup(mat, group)
    sample_cluster <- clusterByGroup(t(mat), sampleData[colnames(mat), "condition"])
    
    dat <- aggregateAUCbyCondition(object)
    dat$group = group[dat$set]
    dat$condition <- sampleData$condition[dat$sample]
    
    g <- ggplot(dat, aes(sample, set, fill = AUC)) +
      geom_raster() + 
      scale_x_discrete(breaks = sample_cluster) +
      scale_y_discrete(breaks = set_cluster) +
      scale_fill_distiller(palette = "GnBu", direction = 1) + 
      facet_grid(rows = vars(group), cols = vars(condition), 
                 scales = "free", space = "free") + 
      labs(x = "Sample", y = "Set", fill = "AUC") +
      theme_bw() +
      theme(panel.spacing = unit(.2, "lines"), 
            panel.border = element_blank(), 
            panel.grid = element_blank(), 
            axis.text = element_blank(), 
            axis.line = element_blank(),
            axis.ticks = element_blank())
    return(g)
  }
)

#' @rdname heatmapAUC
#' @exportMethod heatmapAUC
setMethod(
  "heatmapAUC", "daseqResults",
  function(object, group = NULL) {
    heatmapAUC(object@targetAUC, group = group)
  }
)


## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## ------ general helper ------ 
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#' Bind list by row into data.frame
#' 
#' Bind by row a \code{list} of \code{data.frame} with the same \code{ncol} into \code{data.frame}. 
#' In particular, if the input is a list of vector-like objects (e.g. numeric, atomic, double, etc), each unit of list will be coerced into row vector.
#' In addition, the function allows to save list names if needed.
#' 
#' @param x list, all elements are vectors of the same length or array of the same column size
#' @param save.names logical or character, if not FALSE, save list names into an new column
list_rbind = function(x, save.names = F) {
  x = as.list(x)
  x = x[!sapply(x, is.null)]
  if (!length(x)) return(data.frame())
  if (!is.null(dim(x[[1]]))) 
    x = lapply(x, as.data.frame, stringsAsFactors = F)
  res = suppressWarnings(dplyr::bind_rows(x)) 
  if (is.null(dim(x[[1]])) || length(dim(x[[1]])) == 1) {
    res = t(res)
    colnames(res) = names(x[[1]])
    res = as.data.frame(res)
    ## rownames are automatically inherited
  } 
  
  if (save.names != FALSE) {
    if (is.null(names(x))) 
      stop("The input list is unnamed. Either set save.names to FALSE or set names for x.")
    list.names <- rep(names(x), sapply(x, nrow)) 
    res = cbind(list.names, res, stringsAsFactors = F)
    names(res)[1] <- ifelse(save.names == TRUE, "list.names", save.names) 
  }
  
  res
}


#' return maximum in absolute value
abs.max = function(x, ...) {
  x[which.max(abs(x), ...)]
}


#' Hierarchical clustering of rows by groups
#' @param x matrix
#' @param group group label of rows, disregard names, will be coerce to factor
#' @return character vector of length \code{nrow(x)}, row names. 
#' If `x` doesn't have row.names, this function names it sequantially (from 1).
clusterByGroup = function(x, group) {
  ## check
  if (nrow(x) < 2 || is.null(dim(x))) {
    warning('x is empty or has just one row.')
    return(rownames(x))
  }
  if (any(is.na(x))) {
    warning('Replace NA\'s with 0\'s.')
    x[is.na(x)] = 0 
  }
  
  if (nrow(x) != length(group)) 
    stop("Uncomfortable dimensions.")
  group = as.factor(group)
  if (is.null(rownames(x))) 
    rownames(x) <- seq_len(nrow(x))
  
  res = c()
  for (gp in levels(group)) {
    x1 = x[group == gp, ]
    hc = hclust(dist(x1))
    res = append(res, hc$labels[hc$order])
  }
  res
}


## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## ---- genomic processor ----
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#' Concatenate a list of GRanges
#' 
#' @param grl a list of [GRanges].
#' @param use.names logical, whether (`TRUE`) to inherit list names.
#' @param sep character, seperator between list names and GRanges names.
#' @param save.names logical or character, if not FALSE, save list names as an attribute (i.e., [mcols()]).s
#' @return a [GRanges] object
c_granges <- function(grl, 
                      use.names = T, 
                      sep = ".", 
                      save.names = F) {
  grl <- as(grl, "CompressedGRangesList")
  gr <- unlist(grl, use.names = F)
  list.names <- rep(names(grl), elementNROWS(grl)) 
  if (use.names) {
    names(gr) <- paste(list.names, names(gr), sep = sep) 
  }
  if (save.names != FALSE) {
    attr.names <- ifelse(save.names == TRUE, 
                         "list.names", 
                         save.names) 
    mcols(gr)[[attr.names]] <- list.names
  }
  
  return(gr)
}

#' Add identifiers of genes and transcripts
#' Fill up missing "gene" and "transcript" in genome annotation
#' If the genome annotation lacks "gene" and "transcript" entries, add to it.
#' @param anno a \code{GRanges} object of genome annotation, from \link{import}.
#' @return a \code{GRanges} object with added `gene_id` and `transcript_id` columns.
addGeneTx <- function(anno) {
  l <- levels(anno$type)
  if (! "transcript" %in% l) {
    tx = range(GRangesList(split(anno, anno$transcript_id)))
    tx = unlist(tx)
    tx$type = "transcript" 
    tx$transcript_id = names(tx)
    tx$gene_id <- sapply(split(anno$gene_id, anno$transcript_id), head, n = 1)
    tx <- unname(tx)
    anno = c(tx,anno)
  }
  if (! "gene" %in% l) {
    gene = range(GRangesList(split(anno, anno$gene_id)))
    gene = unlist(gene)
    gene$type = "gene" 
    gene$gene_id = names(gene) 
    gene <- unname(gene)
    anno = c(gene, anno)
  }
  return(sort(anno))
}

#' Import annotation files
#' 
#' This is a wrapper of \link{rtracklayer::import} with specialty in genome annotation.
#' Specifically, it checks whether "gene" or "transcript" exist in the annotation and fill them in if possible.
#' 
#' @param ... whatever input to \link{rtracklayer::import}.
#' @return a \code{GRanges} object.
import <- function(...) {
  anno <- rtracklayer::import(...) 
  
  ## check annotation attributes 
  if (is.null(anno$type)) 
    stop("Cannot find 'type' attribute in annotation.")
  if (!("exon" %in% anno$type)) 
    stop("Genome annotation does not contain any row of type 'exon'.")
  if (!all(c("gene", "transcript") %in% anno$type)) 
    anno = addGeneTx(anno)
  
  return(anno)
}

#' Find genes that contain multiple transcripts
#' 
#' @param anno data.frame or GRanges, annotation
#' @param min integer, minimum number of tx's
#' @param max integer, maximum number of tx's
#' @return gene_id's
getMultiTxGene = function(anno, min = 2, max = Inf) {
  anno_tx = anno[anno$type == 'transcript']
  cnt_tx = table(anno_tx$gene_id)
  names(cnt_tx)[cnt_tx >= min & cnt_tx <= max]
}

## FUN: get raw read count matrix from a list of BAM files
## INPUT:
##    peaks - standard GRanges object 
##    bamfiles - the full paths to the list of BAM files
##    colname - the colnames of output objects
##    n.cores - no. of cores 
##    is.PE - is paired end data or not
## OUTPUT: list()
##    $countMat - n x p count matrix. n <- length(peaks), p <- length(bamfiles)
##    $fragment_counts - a list of length equal to no. of bamfiles
##    $designInfo - the sequencing depth (no. of paired reads if is.PE = TRUE) for each BAM file
##    $peaks - the same as input
## DEPENDS:
##    chromVAR, GenomicAlignments
# getCountMatrix <- function(peaks, bamfiles, 
#                            colname = seq_along(bamfiles), 
#                            n.cores = 1, is.PE = TRUE){
#   if(length(bamfiles) != length(colname)){
#     stop("BAM File list must have the same length as column names list!")
#   }
#   
#   fragment_counts <-
#     mclapply(
#       as.list(bamfiles),
#       getCounts,
#       peaks,
#       paired = is.PE,
#       by_rg = FALSE,
#       format = "bam",
#       mc.cores = n.cores,
#       mc.preschedule = TRUE
#     )
#   countMat <- matrix(NA, length(peaks), length(bamfiles))
#   colnames(countMat) <- colname
#   seqDepth <- NULL
#   for(k in seq_along(bamfiles)){
#     countMat[, k] <- counts(fragment_counts[[k]])[,1]
#     seqDepth[k] <- fragment_counts[[k]]@colData[1,1]
#   }
#   designInfo <- data.frame(exps = colname, depth = seqDepth)
#   
#   res <- list(countMat = countMat, 
#               fragment_counts = fragment_counts, 
#               designInfo = designInfo,
#               peaks = peaks)
#   return(res)
# }

#' Get URL's for UCSC Genome Browser 
#' @param gr a \code{GRanges} object
#' @return a vector of URL's
getURL = function(gr, hub = NULL) {
  if (is.null(hub) || !length(hub)) {
    hub = "http://genome.ucsc.edu/cgi-bin/hgTracks?hubUrl=ftp://ftp.cs.wisc.edu/pub/users/kelesgroup/fchen/hub.txt"
  } else {
    hub = paste0("http://genome.ucsc.edu/cgi-bin/hgTracks?hubUrl=ftp://ftp.cs.wisc.edu/pub/users/kelesgroup/fchen/",hub,".txt")
  }
  if (!is.null(gr$feature)) {
    gr = unlist(range(gr$feature))
  }
  paste0(hub, "&position=", seqnames(gr), ':', start(gr), '-', end(gr))
}


#' Customized \code{featureCounts}
#' 
#' This is a wrapper of \code{Rsubread::featureCounts} allowing verbose suppression.
#' This redirects the output to `/dev/null`, so it assumes a UNIX-like
#' system.
#'
#' @param ... parameters for \code{Rsubread::featureCounts}.
#' @return See \code{Rsubread::featureCounts} documentation.
featureCounts <- function(..., verbose = F) {
  if (verbose) {
    Rsubread::featureCounts(...)
  } else {
    withr::with_output_sink("/dev/null", Rsubread::featureCounts(...))
  }
}

#' Import transcriptome quantification results
#' 
#' Customized tximport from \code{DESeq2} and `RSEM` transcriptome quantification results.
#' 
#' @param files character, files to transcriptome quantification results
#' @return a list of length four, named as `abundance`, `counts`, `length`, and `countsFromAbundance`
tximport2 = function (files, 
                      type = c("none", "kallisto", "salmon", "sailfish", "rsem"), 
                      txIn = TRUE, txOut = FALSE, 
                      countsFromAbundance = c("no", "scaledTPM", "lengthScaledTPM"), 
                      tx2gene = NULL, reader = read.delim, 
                      geneIdCol, txIdCol, abundanceCol, countsCol, lengthCol, importer,
                      collatedFiles, ignoreTxVersion = FALSE, 
                      quiet = F) 
{
  type <- match.arg(type, c("none", "kallisto", "salmon", "sailfish", 
                            "rsem"))
  countsFromAbundance <- match.arg(countsFromAbundance, c("no", 
                                                          "scaledTPM", "lengthScaledTPM"))
  stopifnot(all(file.exists(files)))
  if (!txIn & txOut) 
    stop("txOut only an option when transcript-level data is read in (txIn=TRUE)")
  if (type == "kallisto") {
    geneIdCol = "gene_id"
    txIdCol <- "target_id"
    abundanceCol <- "tpm"
    countsCol <- "est_counts"
    lengthCol <- "eff_length"
    importer <- reader
  }
  if (type %in% c("salmon", "sailfish")) {
    geneIdCol = "gene_id"
    txIdCol <- "Name"
    abundanceCol <- "TPM"
    countsCol <- "NumReads"
    lengthCol <- "EffectiveLength"
    importer <- function(x) reader(x, comment = "#")
  }
  # if (type == "rsem") {
  #   txIn <- FALSE
  #   geneIdCol <- "gene_id"
  #   abundanceCol <- "FPKM"
  #   countsCol <- "expected_count"
  #   lengthCol <- "effective_length"
  #   importer <- reader
  # }
  if (type == "rsem") {
    # txIn <- FALSE
    geneIdCol <- "gene_id"
    txIdCol <- "transcript_id"
    abundanceCol <- "TPM"
    countsCol <- "expected_count"
    lengthCol <- "effective_length"
    percentCol <- "IsoPct"
    importer <- reader
  }
  if (type == "cufflinks") {
    stop("reading from collated files not yet implemented")
  }
  if (txIn) {
    if (!quiet) message("reading in files")
    for (i in seq_along(files)) {
      if (!quiet) message(i, " ", appendLF = FALSE)
      raw <- as.data.frame(importer(files[i]))
      if ((i == 1) & (type %in% c("salmon", "sailfish")) & 
          !("EffectiveLength" %in% names(raw))) {
        lengthCol <- "Length"
        importer <- function(x) {
          tmp <- reader(x, comment = "#", header = FALSE)
          names(tmp) <- c("Name", "Length", "TPM", "NumReads")
          tmp
        }
        raw <- try(as.data.frame(importer(files[i])), silent = TRUE)
        if (inherits(raw, "try-error")) {
          importer <- function(x) {
            reader(x, comment = "#", 
                   col_names = c("Name", "Length", "TPM", "NumReads"))
          }
          raw <- try(as.data.frame(importer(files[i])))
          if (inherits(raw, "try-error")) 
            stop("tried but couldn't use reader() without error\n  user will need to define the importer() as well")
        }
      }
      if (is.null(tx2gene) & !txOut) {
        if (!geneIdCol %in% names(raw)) {
          if (!quiet) message()
          stop("\n\n  tximport failed at summarizing to the gene-level.\n  Please see 'Solutions' in the Details section of the man page: ?tximport\n\n")
        }
        stopifnot(all(c(lengthCol, abundanceCol) %in% 
                        names(raw)))
        if (i == 1) {
          geneId <- raw[[geneIdCol]]
        }
        else {
          stopifnot(all(geneId == raw[[geneIdCol]]))
        }
      }
      else {
        stopifnot(all(c(lengthCol, abundanceCol) %in% 
                        names(raw)))
        if (i == 1) {
          txId <- raw[[txIdCol]]
        }
        else {
          stopifnot(all(txId == raw[[txIdCol]]))
        }
      }
      if (i == 1) {
        mat <- matrix(nrow = nrow(raw), ncol = length(files))
        rownames(mat) <- raw[[txIdCol]]
        colnames(mat) <- names(files)
        abundanceMatTx <- mat
        countsMatTx <- mat
        lengthMatTx <- mat
        percentMatTx <- mat
      }
      abundanceMatTx[, i] <- raw[[abundanceCol]]
      countsMatTx[, i] <- raw[[countsCol]]
      lengthMatTx[, i] <- raw[[lengthCol]]
      percentMatTx[, i] <- raw[[percentCol]]
    }
    if (!quiet) message("")
    txi <- list(abundance = abundanceMatTx, counts = countsMatTx, 
                length = lengthMatTx, percent = percentMatTx,
                countsFromAbundance = "no")
    if (txOut) {
      return(txi)
    }
    txi[["countsFromAbundance"]] <- NULL
    txiGene <- summarizeToGene(txi, tx2gene, ignoreTxVersion, 
                               countsFromAbundance)
    return(txiGene)
  }
  else {
    if (!quiet) message("reading in files")
    for (i in seq_along(files)) {
      if (!quiet) message(i, " ", appendLF = FALSE)
      raw <- as.data.frame(importer(files[i]))
      stopifnot(all(c(geneIdCol, abundanceCol, lengthCol) %in% 
                      names(raw)))
      if (i == 1) {
        mat <- matrix(nrow = nrow(raw), ncol = length(files))
        rownames(mat) <- raw[[geneIdCol]]
        colnames(mat) <- names(files)
        abundanceMat <- mat
        countsMat <- mat
        lengthMat <- mat
      }
      abundanceMat[, i] <- raw[[abundanceCol]]
      countsMat[, i] <- raw[[countsCol]]
      lengthMat[, i] <- raw[[lengthCol]]
    }
  }
  if (!quiet) message("")
  return(list(abundance = abundanceMat, counts = countsMat, 
              length = lengthMat, countsFromAbundance = "no"))
}


# ## ------ depreciated ------ 
# ## FUN: transform idr result (overlapped peaks) table into GRanges 
# ## INPUT: idr w/ the following columns:
# # [1] "chr1"       "start1"     "stop1"      "sig.value1" "chr2"      
# # [6] "start2"     "stop2"      "sig.value2" "idr.local"  "IDR"
# ## OUTPUT: GRanges object, w/ following mcols()
# ##   idr.local - local idr score
# ##   IDR - global idr score, the smaller the more reproducible
# ## DEPEND: GenomicRanges
# idr2GRanges = function(idr) {
#   ## Check idr 
#   start = apply(cbind(idr$start1, idr$start2), 1, max)
#   stop = apply(cbind(idr$stop1, idr$stop2), 1, min)
#   
#   wrong = (idr$chr1 != idr$chr2) | 
#     (idr$start1 > idr$stop1) | (idr$start2 > idr$stop2) | 
#     (stop - start > 1e+3) | (stop - start < 0)
#   if (sum(wrong) > 0)
#     warning(paste("There are", sum(wrong), "incorrect overlapping peaks being discarded"))
#   idr = idr[!wrong,] 
#   start <- start[!wrong] 
#   stop <- stop[!wrong] 
#   
#   GRanges(idr$chr1, IRanges(start, stop), 
#           idr.local = idr$idr.local, IDR = idr$IDR) 
# }
# 
# ## FUN: add strand info idr (overlapped) peaks, using original peak files
# ## INPUT: 
# ##   peak - idr peaks
# ##   filename - original peak files
# ## OUTPUT: 
# ##   character, of same length as peak, "+", "-", "*" (if ambiguous)
# extCols_nPk <- c(singnalValue = "numeric",
#                  pValue = "numeric",
#                  qValue = "numeric",
#                  peak = "integer")
# getStrand = function(peak, filename_bed) {
#   peak_raw <- GRanges()
#   for (fn in filename_bed) 
#     peak_raw = c(peak_raw,
#                  import(fn, format = "BED", 
#                         extraCols = extCols_nPk))
#   nIdr_plus = countOverlaps(peak, peak_raw[strand(peak_raw) == '+'])
#   nIdr_minus = countOverlaps(peak, peak_raw[strand(peak_raw) == '-'])
#   ifelse(nIdr_plus > nIdr_minus, '+', 
#          ifelse(nIdr_plus < nIdr_minus, '-', '*'))
# }
# 
# 
# ## huber's location estimation
# huber = function(y, ...) {
#   MASS::huber(y = y, ...)['mu']
# }
# 
# ## brute-force choose, this is just a wraper of combn() function in utils 
# bruteCombn = function(set, ...) {
#   if (length(set) == 1) 
#     return(list(set))
#   return(combn(x = set, ...))
# }
