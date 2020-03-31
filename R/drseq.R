## SURF Analysis Module 1 

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## ------ constructor ------
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#' Prepare DrSeq annotation files
#' 
#' To use the `useMetaFeatures` functionality of [Rsubread::featureCounts], we need the GTF input. 
#' This function is a helper function that produce the *two* annotation files ("event" and "exon") needed.
#' 
#' @param anno_event a `surf` object
#' @param anno.prefix character, prefix of exon/event annotation files for saving. These files are needed by [Rsubread].
#' @param anno.format character, e.g. "gtf", as accepted by [rtracklayer::export].
#' @param remove.overlap.exon logical, remove overlaping exons across genes (default to `FALSE`).
#' @param cores integer, number of available workers, sent to `nthreads` of [Rsubread::featureCounts].
#' @param verbose logical, whether (`TRUE`) to echo progess.
#' @return NULL, the function is a procedure and only output/export results to file system, except for messages and warnings.
#' @export
prepDrseqAnno = function(anno_event, 
                         anno.prefix = "annotation.drseq", 
                         anno.format = "gff2", 
                         remove.overlap.exon = F, 
                         cores = max(1, detectCores()-2), 
                         verbose = T) {
  ## reconstruct genome annotation from @genePartsList
  ## with $type and $gene_id attributes
  isoPL = anno_event@genePartsList
  if (verbose) cat("Recognizing", nrow(isoPL), "gene parts list...\n")
  registerDoParallel(cores)
  anno <- foreach (segment = isoPL$segment, label = isoPL$label) %dopar% {
    gene <- range(segment[!!label])
    gene$type = "gene"
    ## range() and reduce() should be equivalent here 
    mergedSegment <- c_granges(range(split(segment, label)), 
                               use.names = F, 
                               save.names = "exon_label") 
    exon <- mergedSegment[mergedSegment$exon_label != "0"] 
    exon$type = "exon"
    c(gene, exon)
  }
  stopImplicitCluster()
  anno <- c_granges(setNames(anno, isoPL$gene_id), 
                    use.names = F, save.names = "gene_id")
  
  ## ---- "exon" flattened file
  output.exon <- paste0(anno.prefix, '.exon.', anno.format)
  if (verbose) cat("Outputing exon flattened file to", output.exon, "...\n")
  
  ## merge exons by genes
  exon = anno[anno$type == 'exon']
  exon = GRangesList(split(exon, exon$gene_id))
  exon = unlist(reduce(exon))
  ## remove overlaping exons across genes (optional)
  if (remove.overlap.exon) { 
    exon = exon[countOverlaps(exon, exon) < 2]
  }
  
  ## exonic_part
  exon$source = factor('drseq')
  exon$type = factor('exonic_part')
  exon$gene_id = names(exon)
  num = unlist(lapply(rle(exon$gene_id)$lengths, seq_len))
  exon$exonic_part_number = stringr::str_pad(num, 3, pad = "0")
  exon$exonic_part_id = paste(exon$gene_id, exon$exonic_part_number, sep = ":")
  
  ## aggregate_gene
  gene = anno[anno$type == 'gene' & anno$gene_id %in% exon$gene_id]
  names(gene) = gene$gene_id
  mcols(gene) = NULL
  gene$source = factor('drseq')
  gene$type = factor('aggregate_gene')
  gene$gene_id = names(gene)
  
  ## output exon annotation
  gr1 = c(gene, exon)
  gr1 = unlist(GRangesList(split(unname(gr1), gr1$gene_id)), use.names = F)
  rtracklayer::export(gr1, output.exon, anno.format)
  
  ## ---- "event" flattened file
  output.event <- paste(anno.prefix, 'event', anno.format, sep = ".")
  if (verbose) cat("Outputing event flattened file to", output.event, "...\n")
  
  ## remove events outside exons (if some exons are removed)
  if (remove.overlap.exon) {
    hit <- findOverlaps(anno_event, exon, type = "within") ## AFE1/ALE1 may range multiple exons
    hit <- hit[anno_event$gene_id[from(hit)] == exon$gene_id[to(hit)]]
    anno_event <- anno_event[unique(from(hit))] 
  }
  
  ## exonic_part
  event <- c_granges(anno_event$genomicData, sep = "") 
  exonic_part_count <- elementNROWS(anno_event$genomicData)
  event$source = factor('drseq')
  event$type = factor('exonic_part')
  event$event_id = names(event)
  event$gene_id = rep(anno_event$gene_id, exonic_part_count)
  event$transcript_id = rep(anno_event$gene_id, exonic_part_count)
  num = rep(unlist(lapply(rle(anno_event$gene_id)$lengths, seq_len)), exonic_part_count) 
  event$event_number = stringr::str_pad(num, 3, pad = "0")
  # event$exonic_part_id = paste(event$gene_id, event$exonic_part_number, sep = ":")
  
  # ## aggregate_gene
  # gene = unlist(range(anno_event$genomicData))
  # gene$source = factor('drseq')
  # gene$type = factor('aggregate_gene')
  # gene$gene_id = names(gene)
  
  ## output 
  gr2 = c(gene, event)
  gr2 = unlist(GRangesList(split(unname(gr2), gr2$gene_id)), use.names = F)
  rtracklayer::export(gr2, output.event, anno.format)
  
}


#' DrSeq Data 
#' 
#' Create the `drseqData` slot needed by DrSeq analysis, which is a `DEXSeqDataSet`.
#' It is possible to output event read counts for later usage.
#' 
#' @param event a `surf` object from [parseEvent].
#' @param sampleData data.frame, describes the RNA-seq samples, contains at least two columns -- `bam` and `condition`, the row.names represent sample names.
#' @inheritParams DEXSeq::DEXSeqDataSet ## the `design` parameter.
#' @param remove.overlap.exon logical, remove overlaping exons across genes (default to FALSE).
#' @param anno.prefix file names for outputting annotation files.
#' @param anno.format character, e.g. "gtf", as accepted by [rtracklayer::export].
#' @param minMQS;isPairedEnd as defined in [Rsubread::featureCounts]. Note that the default is customized for SURF (see details for more information). 
#' @param cores integer, number of available workers, sent to `nthreads` of `featureCounts`.
#' @param verbose logical, whether (`TRUE`) to echo progess.
#' @param ... additional parameters for [Rsubread::featureCounts]. 
#' @details If you used Illumina HiSeq 2000, set `strandSpecific = 2` (reversed strand).
#' @return a `surf` object, with `drseqData` slot updated.
#' @export
drseqData = function(event, sampleData, 
                     design = ~ sample + exon + condition:exon,
                     remove.overlap.exon = F, 
                     anno.prefix = "drseq.annotation", 
                     anno.format = "gff2", 
                     minMQS = 10, 
                     isPairedEnd = T, 
                     cores = max(1, detectCores()-2), 
                     verbose = F, 
                     ...) {
  # check input sampleData
  sampleData <- as.data.frame(sampleData)
  if (is.null(sampleData$bam)) 
    stop("sampleData must contain \"bam\" column.")
  bam.files <- sampleData$bam
  if (is.null(sampleData$condition)) 
    stop("sampleData must contain \"condition\" column.")
  sampleData$condition <- as.factor(sampleData$condition)
  if (length(levels(sampleData$condition))!=2) 
    stop("The condition must have two levels, for IP and input respectively.")
  if (is.null(rownames(sampleData))) {
    warning("sampleData has no row.names input, created automatically.")
    rownames(sampleData) <- paste0(
      as.character(sampleData$condition), 
      unlist(lapply(rle(as.character(sampleData$condition))$lengths, seq_len))
    )
  }
  
  exon.file <- paste(anno.prefix, 'exon', anno.format, sep = ".")
  event.file <- paste(anno.prefix, 'event', anno.format, sep = ".")
  if (!file.exists(exon.file) || !file.exists(event.file)) {
    cat("Cannot find event/exon annotation files at\n   ", 
        exon.file, "\n   ", event.file, 
        "\nMaking them freshly...\n")
    ## prepare annotation for exons and events 
    prepDrseqAnno(
      anno_event = event, 
      anno.prefix = anno.prefix, anno.format = anno.format,
      remove.overlap.exon = remove.overlap.exon, 
      cores = cores,
      verbose = verbose
    )
  } 
  
  ## count reads on exons
  count.exon = featureCounts(
    files = bam.files, 
    annot.ext = exon.file, 
    isGTFAnnotationFile = T,
    GTF.featureType = "exonic_part", 
    GTF.attrType = "gene_id", 
    allowMultiOverlap = T, 
    minMQS = minMQS,
    isPairedEnd = isPairedEnd, 
    nthreads = cores, 
    verbose = verbose, ...
  )
  colnames(count.exon$counts) = 
    colnames(count.exon$stat)[-1] = 
    rownames(sampleData)
  # cat("Total count:", sum(count.exon$counts), "\n")
  
  ## count reads on events
  count.event = featureCounts(
    files = bam.files, 
    annot.ext = event.file, 
    isGTFAnnotationFile = T,
    GTF.featureType = "exonic_part", 
    GTF.attrType = "event_id", 
    allowMultiOverlap = T, 
    minMQS = minMQS,
    isPairedEnd = isPairedEnd, 
    nthreads = cores, 
    verbose = verbose, ...
  )
  colnames(count.event$counts) = 
    colnames(count.event$stat)[-1] = 
    rownames(sampleData)
  
  ## DEXSeqDataSet
  index <- match(event$event_id, rownames(count.event$counts))
  if (any(is.na(index))) {
    stop("Count results can not match the event annotation.")
  } ## just to be sure
  drd <- DEXSeqDataSet(
    countData = count.event$counts[index,], 
    sampleData = sampleData, 
    design = design, 
    featureID = event$event_id, 
    groupID = event$gene_id
  )
  
  ## replace with gene reads count (still a DEXSeqDataSet)
  cnt.other = count.exon$counts[rowData(drd)$groupID,] - counts(drd)[,seq_along(bam.files)]
  counts(drd)[,length(bam.files)+seq_along(bam.files)] = cnt.other
  rowMin <- rowMin(counts(drd))
  if (any(rowMin < 0)) {
    warning(sum(rowMin < 0), " event(s) prompt to invalid REU coefficients.")
    drd = drd[rowMin >= 0,]
  }
  
  ## update our "surf" object
  event@drseqData = drd
  event@sampleData$"RNA-seq" = DataFrame(
    sample = rownames(sampleData), 
    sampleData[setdiff(names(sampleData), "sample")]
  )
  metadata(event) = c(
    metadata(event), 
    remove.overlap.exon = remove.overlap.exon
  )
  
  return(event)
}


#' DrSeq Fit
#' 
#' Fit DrSeq models using the DEXSeq as the estimation engine.
#' @param drd a [surf] object from `drseqData`.
#' @param cores integer, number of computing workers.
#' @param verbose logical, whether (`TRUE`) to echo progess.
#' @return a [surf] object with (1) `drseqResults` and `sampleData` slot updated and (2) three added columns:
#' \item{eventBaseMean}{base read coverage of the event from RNA-seq data.}
#' \item{padj}{adjusted p-value for differential REU.}
#' \item{logFoldChange}{estimated log2 fold change of REU: log2(condition of the 1st sample / another condition)}
#' @details the drseqResults slot contains extensive DrSeq results. To access the object, use [deseqResults()] function.
#' @export
drseqFit <- function(drd, 
                     fullModel = design(object),
                     reducedModel = ~ sample + exon,
                     fitExpToVar = "condition",
                     cores = max(1, detectCores()-2), 
                     verbose = F) {
  dxd <- drd@drseqData ## "DEXSeqDataSet"
  if (is.null(dxd)) stop("Cannot find drseqData.")
  
  BPPARAM = MulticoreParam(workers = cores)
  dxdl = split(dxd, drd[mcols(dxd)$featureID, "event_name"])
  
  if (verbose) cat("Fitting DrSeq for", paste(names(dxdl), collapse = ", "), "...\n")
  t <- system.time({
    dxrl <- lapply(dxdl, DEXSeq, 
                   fullModel = fullModel, 
                   reducedModel = reducedModel, 
                   fitExpToVar = fitExpToVar,
                   BPPARAM = BPPARAM, 
                   quiet = !verbose)
  })
  if (verbose) cat("Running time:", t[3], "sec. \n")
  
  ## added columns
  dxr = do.call("rbind", dxrl)
  index <- match(drd$event_id, dxr$featureID)
  if (any(is.na(index))) {
    stop("Event ID's do not match between DrSeq and annotation.")
  }
  dxr_sub <- dxr[index, c(3:10, 12)]
  dxr_sub$padj = p.adjust(dxr_sub$pvalue, method = "fdr") ## re-control FDR across all event types
  names(dxr_sub)[1] <- "eventBaseMean"
  names(dxr_sub)[8] <- "logFoldChange" ## log2(condition of the 1st sample / another condition)
  mcols(dxr_sub)$type <- "drseq"
  mcols(dxr_sub)["eventBaseMean", "description"] <- "mean of the counts across samples in each event"
  mcols(dxr_sub)["dispersion", "description"] <- "event dispersion estimate"
  mcols(dxr_sub)["countData", "type"] <- "RNA-seq"
  
  ## representation
  sampleData = dxrl[[1]]@sampleData
  modelFrameBM = lapply(dxrl, slot, "modelFrameBM")
  dispersionFunction = lapply(dxrl, slot, "dispersionFunction")
  metadata = do.call("c", lapply(dxrl, slot, "metadata"))
  
  ## take a subset of drd's column 
  
  drr <- new("drseqResults", cbind(drd, dxr_sub), 
             modelFrameBM = DataFrameList(modelFrameBM), 
             dispersionFunction = List(dispersionFunction))
  drr@metadata <- c(drd@metadata, metadata)
  
  res <- new(
    "surf", 
    cbind(drd, drr[c("eventBaseMean", "padj", "logFoldChange")]), 
    genePartsList = event@genePartsList, 
    drseqData = event@drseqData, 
    drseqResults = drr, ## newly added
    faseqData = event@faseqData, 
    faseqResults = event@faseqResults, 
    daseqResults = event@daseqResults,
    sampleData = event@sampleData
  )
  res@sampleData$"RNA-seq" <- sampleData
  metadata(res) <- metadata(event)
  return(res)
}

#' Screen SURF training samples
#' 
#' Filter DrSeq results for the analysis module 2 of SURF (DASeq).
#' 
#' @param event a [surf] object output by [drseqFit].
#' @param drseq.fdr FDR (BH procedure) adjusted p-value cut-off.
#' @param read.length numeric, RNA-seq read length. Default is 100 bp (e.g., Illumina TruSeq). This is used to adjust event base count, which is then used to select the representative events if replicated.
#' @param min.adjMean numeric, adjusted event base mean threshold.
#' @param filter.overlap.event logical, whether (default to `TRUE`) to select one representitive event from overlapping ones and remove the others.
#' @param verbose logical, whether (default to `FALSE`) to print out basic summary statistics.
#' 
#' @return a `surf` object, with three columns added: 
#' \item{adjMean}{adjusted base mean of the event from RNA-seq data.}
#' \item{group}{group labels of ATR events, `increase` for increased REU upon RBP knock-down, and `decrease` for decreased, and `no change` for no-changed.}
#' \item{included}{logical, indicating whether the event is included into SURF ananlysis module 2.}
#' @export
drseqFilter = function(event, 
                       drseq.fdr = .05, 
                       read.length = 100,
                       min.adjMean = .05, 
                       filter.overlap.event = T,
                       verbose = F) {
  ## check input
  if (!all(c("eventBaseMean", "padj", "logFoldChange") %in% colnames(event))) {
    stop("You should perform DrSeq first.") 
  }
  if (any(c("adjMean", "group", "included") %in% colnames(event))) {
    event$adjMean = NULL
    event$group = NULL
    event$included = NULL
  }
  
  ## event read coverage
  eventLen = sapply(width(event$genomicData), sum)
  adjMean = event$eventBaseMean / (eventLen + read.length - 1) 
  ## group by logFoldChange x padj
  group = rep("no change", nrow(event))
  group[event$padj < drseq.fdr & event$logFoldChange < 0] = "decrease"
  group[event$padj < drseq.fdr & event$logFoldChange > 0] = "increase"
  group = ordered(group, c("decrease", "no change", "increase"))
  
  ## included event id
  id_evtDU = event$event_id[!is.na(event$padj) & 
                              event$padj < drseq.fdr &
                              adjMean > min.adjMean &
                              !is.na(event$logFoldChange) & 
                              abs(event$logFoldChange) > 0]
  id_evtEU = event$event_id[!is.na(event$padj) & 
                              event$padj > .4 & 
                              adjMean > min.adjMean & 
                              !is.na(event$logFoldChange) & 
                              abs(event$logFoldChange) > 0]
  if (verbose) 
    cat("Identified", length(id_evtDU), "DR events and", 
        length(id_evtEU), "ER events.\n")
  id_evt = union(id_evtDU, id_evtEU)
  
  ## pooling duplicated events by max{adjMean}
  ## (i) same ATR event  
  ## (ii) overlapped event body
  if (filter.overlap.event) {
    event1 = event[event$event_id %in% id_evt,]
    adjMean1 <- adjMean[event$event_id %in% id_evt]
    event1 = event1[order(adjMean1, decreasing = T),] ## sort by normalized read coverage
    hits = findOverlaps(event1$genomicData, event1$genomicData)
    hits = hits[from(hits) > to(hits)] ## from back to front, will remove the smaller in adjMean
    hits = hits[event1$event_name[from(hits)] == event1$event_name[to(hits)]] ## same event names 
    if (verbose) cat(length(unique(from(hits))), "overlapping events.\n")
    id_evt = setdiff(id_evt, event1$event_id[from(hits)])
  }
  
  if (!length(id_evt)) stop("No usable SURF instance!")
  included = event$event_id %in% id_evt
  
  ## new columns
  ds = DataFrame(adjMean, group, included)
  mcols(ds)$type = c("RNA-seq", "faseq", "faseq")
  mcols(ds)$description <- c(
    "adjusted base mean of each event",
    "direction of differential regulation (upon knock-down)", 
    "indicator of usable event")
  
  ## wrap into faseqData object
  res <- new(
    "surf", cbind(event, ds), 
    genePartsList = event@genePartsList, 
    drseqData = event@drseqData, 
    drseqResults = event@drseqResults,
    faseqData = event@faseqData, 
    faseqResults = event@faseqResults, 
    daseqResults = event@daseqResults,
    sampleData = event@sampleData
  )
  metadata(res) = c(
    metadata(event), 
    drseq.fdr = drseq.fdr, 
    RNAseq.read.length = read.length, 
    min.adjMean = min.adjMean, 
    filter.overlap.event = filter.overlap.event
  )
  
  if (verbose) {
    cat(length(id_evt), "events included for SURF analysis.\n", 
        "The distribution of AS/ATI/APA events identified:\n")
    print(table(res[res$included,c("event_name", "group")]))
  }
  
  return(res)
}

#' DrSeq
#' 
#' Perform the differential REU (DrSeq) test in a single command. 
#' This function is a wrapper that calls the necessary functions in order for DrSeq.
#' @inheritParams drseqData
#' @param ... parameters for [Rsubread::featureCounts].
#' @inheritParams drseqFit
#' @inheritParams drseqFilter
#' @return a `surf` object
#' @references Chen F and Keles S. "SURF: Integrative analysis of a compendium of RNA-seq and CLIP-seq datasets highlights complex governing of alternative transcriptional regulation by RNA-binding proteins."
#' @export
drseq <- function(event, 
                  
                  ## data
                  sampleData, 
                  design = ~ sample + exon + condition:exon,
                  remove.overlap.exon = F, 
                  anno.prefix = "drseq.annotation", 
                  anno.format = "gff2", 
                  minMQS = 10, 
                  isPairedEnd = T, 
                  ..., 
                  
                  ## fit
                  fullModel = design,
                  reducedModel = ~ sample + exon,
                  fitExpToVar = "condition",
                  
                  ## filter
                  drseq.fdr = .05, 
                  read.length = 100,
                  min.adjMean = .05, 
                  filter.overlap.event = T, 
                  
                  cores = max(1, detectCores()-2), 
                  verbose = F) {
  if (verbose) 
    cat("Preparing DrSeq count dataset...\n")
  event <- drseqData(
    event, sampleData,
    design = design, 
    remove.overlap.exon = remove.overlap.exon, 
    anno.prefix = anno.prefix, 
    anno.format = anno.format, 
    minMQS = minMQS, 
    isPairedEnd = isPairedEnd, 
    cores = cores, 
    verbose = verbose, ...
  )
  
  if (verbose) 
    cat("Fitting DrSeq models...\n")
  event <- drseqFit(
    event, 
    fullModel = fullModel,
    reducedModel = reducedModel,
    fitExpToVar = fitExpToVar,
    cores = cores, 
    verbose = verbose
  )
  
  if (verbose) 
    cat("Annotating/filtering DrSeq results...\n")
  event <- drseqFilter(
    event,
    drseq.fdr = drseq.fdr,
    read.length = read.length,
    min.adjMean = min.adjMean,
    filter.overlap.event = filter.overlap.event,
    verbose = verbose
  )
  
  return(event)
}

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## ------ methods ------
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
#' @param fdr.cutoff \code{numeric}, significance level of adjusted p-value.
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



