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
#' @references Chen F and Keles S. "SURF: Integrative analysis of a compendium of RNA-seq and CLIP-seq datasets highlights complex governing of alternative transcriptional regulation by RNA-binding proteins."
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

#' Gene Parts List
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
    event <- object[c("event_id", "event_name", 
                      "gene_id", "transcript_id", 
                      "genomicData", "feature")]
    drr <- object@drseqResults
    res <- new(
      "drseqResults", cbind(event, drr), 
      modelFrameBM = DataFrameList(modelFrameBM), 
      dispersionFunction = List(dispersionFunction
      ))
    metadata(res) <- metadata(object)
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
