
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## ------ drseq ------
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## SURF Analysis Module 1

#' Prepare DrSeq annotation files
#'
#' To use the \code{useMetaFeatures} functionality of \link{Rsubread::featureCounts}, we need the GTF input.
#' This function helps to produce the *two* annotation files ("event" and "exon") needed by \link{drseqCount}.
#'
#' @param anno_event a \code{surf} object
#' @param anno.prefix \code{character}, prefix of exon/event annotation files for saving. These files are needed by \code{Rsubread::featureCounts}.
#' @param anno.format \code{character}, e.g. "gtf", as accepted by \code{rtracklayer::export}.
#' @param remove.overlap.exon \code{logical}, remove overlaping exons across genes (default to \code{FALSE}).
#' @param cores \code{integer}, number of available workers, sent to \code{nthreads} of \code{Rsubread::featureCounts}.
#' @param verbose \code{logical}, whether (\code{TRUE}) to echo progess.
#' @return \code{NULL}, the function is a procedure and only output/export results to file system, except for messages and warnings.
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
    mergedSegment <- c_granges(range(S4Vectors::split(segment, label)),
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
  exon = GRangesList(S4Vectors::split(exon, exon$gene_id))
  exon = unlist(GenomicRanges::reduce(exon))
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
  gr1 = unlist(GRangesList(S4Vectors::split(unname(gr1), gr1$gene_id)), use.names = F)
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
  gr2 = unlist(GRangesList(S4Vectors::split(unname(gr2), gr2$gene_id)), use.names = F)
  rtracklayer::export(gr2, output.event, anno.format)
  
}


#' Construct DrSeq Data
#'
#' This function creates the \code{drseqData} slot needed by DrSeq analysis, which is a \link{DEXSeqDataSet}.
#' This function requires two DrSeq annotation files created by \link{prepDrseqAnno}. Two files share the same prefix and end with ".exon.<file type>" and ".event.<file.type>" respective.
#' If annotation files are missing, this function can create them freshly, which might take some time.
#' This function counts the RNA-seq reads on ATR events and genes using \code{featureCounts}.
#'
#' @param event a \code{surf} object from \link{parseEvent}.
#' @param sampleData \code{data.frame}, describes the RNA-seq samples, contains at least two columns -- `bam` and `condition`, the row.names represent sample names.
#' @inheritParams DEXSeq::DEXSeqDataSet ## the \code{design} parameter.
#' @param remove.overlap.exon \code{logical}, remove overlaping exons across genes (default to \code{FALSE}).
#' @param anno.prefix \code{character}, file names for outputting annotation files.
#' @param anno.format \code{character}, e.g. "gtf", as accepted by \link{rtracklayer::export}.
#' @param minMQS;isPairedEnd as defined in \link{Rsubread::featureCounts}. Note that the default is customized for SURF (see details for more information).
#' @param cores \code{integer}, number of available workers, sent to \code{nthreads} of \code{featureCounts}.
#' @param verbose \code{logical}, whether (\code{TRUE}) to echo progess.
#' @param ... additional parameters for \code{Rsubread::featureCounts}.
#' @details If you used Illumina HiSeq 2000, set \code{strandSpecific = 2} (reversed strand).
#' @return a \code{surf} object, with \code{drseqData} slot updated.
#' @export
drseqCount = function(event, sampleData,
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
  cat(fill = TRUE)
  
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
  metadata(event) = metadata(event)
  metadata(event)$remove.overlap.exon = remove.overlap.exon
  
  return(event)
}


#' DrSeq Fit
#'
#' Fit DrSeq models using the DEXSeq as the estimation engine.
#' @param drd a \code{surf} object from \link{drseqCount}.
#' @param cores \code{integer}, number of computing workers.
#' @param verbose \code{logical}, whether (\code{TRUE}) to echo progess.
#' @return a \code{surf} object with (1) \code{drseqResults} and \code{sampleData} slot updated and (2) three added columns:
#' \item{eventBaseMean}{base read coverage of the event from RNA-seq data.}
#' \item{padj}{adjusted p-value for differential REU.}
#' \item{logFoldChange}{estimated log2 fold change of REU: log2(condition of the 1st sample / another condition)}
#' @details the \code{drseqResults} slot contains extensive DrSeq results. To access the object, use \link{deseqResults} function.
#' @export
drseqFit <- function(drd,
                     fullModel = design(drd@drseqData),
                     reducedModel = ~ sample + exon,
                     fitExpToVar = "condition",
                     cores = max(1, detectCores()-2),
                     verbose = F) {
  dxd <- drd@drseqData ## "DEXSeqDataSet"
  if (is.null(dxd)) stop("Cannot find drseqData.")
  
  BPPARAM = MulticoreParam(workers = cores)
  dxdl = S4Vectors::split(dxd, drd[mcols(dxd)$featureID, "event_name"])
  
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
  dxr_sub <- as(dxr[index, c(3:10, 12)], "DataFrame")
  dxr_sub$padj <- p.adjust(dxr_sub$pvalue, method = "fdr") ## re-control FDR across all event types
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
  
  drr <- new("drseqResults", 
             cbind(as(drd, "DataFrame"), as(dxr_sub, "DataFrame")),
             modelFrameBM = DataFrameList(modelFrameBM),
             dispersionFunction = List(dispersionFunction))
  drr@metadata <- c(drd@metadata, metadata)
  
  res <- new(
    "surf",
    cbind(as(drd, "DataFrame"), 
          as(drr[c("eventBaseMean", "padj", "logFoldChange")], "DataFrame")),
    genePartsList = drd@genePartsList,
    drseqData = drd@drseqData,
    drseqResults = drr, ## newly added
    faseqData = drd@faseqData,
    faseqResults = drd@faseqResults,
    daseqResults = drd@daseqResults,
    sampleData = drd@sampleData
  )
  res@sampleData$"RNA-seq" <- sampleData
  metadata(res) <- metadata(drd)
  return(res)
}

#' Screen SURF training samples
#'
#' Filter DrSeq results for the analysis module 2 of SURF (DASeq).
#'
#' @param event a \code{surf} object output by \link{drseqFit}.
#' @param drseq.fdr \code{numeric}, FDR (BH procedure) adjusted p-value cut-off.
#' @param read.length \code{numeric}, RNA-seq read length. Default is 100 bp (e.g., Illumina TruSeq). This is used to adjust event base count, which is then used to select the representative events if replicated.
#' @param min.adjMean \code{numeric}, adjusted event base mean threshold.
#' @param filter.overlap.event \code{logical}, whether (default to \code{TRUE}) to select one representitive event from overlapping ones and remove the others.
#' @param verbose \code{logical}, whether (default to \code{FALSE}) to print out basic summary statistics.
#'
#' @return a \code{surf} object, with three columns added:
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
    "surf", 
    cbind(as(event, "DataFrame"), ds),
    genePartsList = event@genePartsList,
    drseqData = event@drseqData,
    drseqResults = event@drseqResults,
    faseqData = event@faseqData,
    faseqResults = event@faseqResults,
    daseqResults = event@daseqResults,
    sampleData = event@sampleData
  )
  metadata(res) = metadata(event)
  metadata(res)$drseq.fdr = drseq.fdr
  metadata(res)$RNAseq.read.length = read.length
  metadata(res)$min.adjMean = min.adjMean
  metadata(res)$filter.overlap.event = filter.overlap.event
  
  if (verbose) {
    cat(length(id_evt), "events included for SURF analysis.\n",
        "The distribution of AS/ATI/APA events identified:\n")
    print(table(data.frame(res[res$included, c("event_name", "group")])))
  }
  
  return(res)
}

#' DrSeq
#'
#' Perform the differential REU (DrSeq) test in a single command.
#' This function is a wrapper that calls the necessary functions in order for DrSeq.
#' @inheritParams drseqCount
#' @param ... parameters for \link{Rsubread::featureCounts}.
#' @inheritParams drseqFit
#' @inheritParams drseqFilter
#' @return a \code{surf} object
#' @references Chen, F., & Keles, S. (2020). SURF: integrative analysis of a compendium of RNA-seq and CLIP-seq datasets highlights complex governing of alternative transcriptional regulation by RNA-binding proteins. *Genome Biology*, 21(1), 1-23.
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
  timer <- system.time({
    event <- drseqCount(
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
  })
  if (verbose) cat("Run time (RNA-seq read counting):", timer[3], "sec.\n")
  
  # if (verbose) cat("Fitting DrSeq models...\n")
  timer <- system.time({
    event <- drseqFit(
      event,
      fullModel = fullModel,
      reducedModel = reducedModel,
      fitExpToVar = fitExpToVar,
      cores = cores,
      verbose = verbose
    )
  })
  if (verbose) cat("Run time (model fitting):", timer[3], "sec.\n")
  
  if (verbose) cat("Annotating/referencing DrSeq results...\n")
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
## ------ faseq ------
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## SURF Analysis Module 2

#' Construct FASeq Data Set
#'
#' This function uantifies feature signals for location features using CLIP-seq data.
#' You align CLIP-seq reads to the genome and provide FASeq with the resulting bam files.
#' We will take care of the rest.
#'
#' @param event a \code{surf} object.
#' @param sampleData \code{data.frame}, must contain two columns -- "bam" and "condition" (for "IP" and "input", where "IP" should come first), whose \code{row.names} represent the sample names. "bam" is the file name of CLIP-seq bam. "condition" will be coerced to factor, whose first level will be treated as IP, and the second level as input.
#' @param signal.type \code{character}, indicate the type of feature signal wanted, support "TPM" for Transcripts Per Kilobase Million, "FPKM" for Fragments Per Kilobase Million (for paired-end reads) and Reads Per Kilobase Million (for single-end reads), and "raw.count" for raw.count read counts
#' @param FUN.aggregate \code{function}, used for aggreating signals within \code{conditon}, default to mean.
#' @param cores \code{integer}, number of available workers, sent to \code{nthreads} of \link{featureCounts}
#' @param verbose \code{logical}, whether (default to \code{TRUE}) to echo progess
#' @param minMQS;minOverlap;isPairedEnd;... parameters for \link{featureCounts}. \code{minMQS} is default to 10, and \code{minOverlap} is default to 12 (25% of the typical read length of CLIP-seq (~50bp)), and \code{isPairedEnd} is default to \code{TRUE}.
#' @return a \code{surf} object, with (i) one column \code{featureSignal} added, (2) \code{faseqData} slot updated, and (3) \code{sampleData} slot updated.
#' @details If your sequencing platform is Illumina HiSeq 2000, set \code{strandSpecific = 2}.
#' @keywords feature signal, CLIP-seq
#' @references \url{https://www.rna-seqblog.com/rpkm-fpkm-and-tpm-clearly-explained/}
#' @export
faseqCount <- function(event, sampleData,
                       signal.type = "FPKM",
                       FUN.aggregate = "mean",
                       minMQS = 10,
                       minOverlap = 12,
                       isPairedEnd = T,
                       cores = max(1, detectCores()-2),
                       verbose = F, ...) {
  ## check @sampleData
  if (length(signal.type) != 1 || !signal.type %in% c("TPM", "FPKM", "raw.count")) {
    stop("Please specify signal type.")
  }
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
      unlist(lapply(rle(as.character(sampleData$condition))$lengths, seq_len)))
  }
  
  
  if (any("featureSignal" %in% colnames(event))) {
    event$featureSignal = NULL
  }
  
  ## featureCounts
  count.feature = featureCounts(
    bam.files,
    annot.ext = ensembldb::toSAF(event$feature),
    useMetaFeatures = F,
    allowMultiOverlap = T,
    minOverlap = minOverlap,
    minMQS = minMQS,
    isPairedEnd = isPairedEnd,
    nthreads = cores,
    verbose = verbose, ...
  )
  colnames(count.feature$counts) =
    colnames(count.feature$stat)[-1] =
    rownames(sampleData)
  
  ## prepare @sampleData
  CLIPseqSampleData = DataFrame(
    sample = rownames(sampleData),
    sampleData[setdiff(names(sampleData), "sample")],
    depth = colSums(count.feature$stat[rownames(sampleData)])
  )
  mcols(CLIPseqSampleData)$type = "input"
  mcols(CLIPseqSampleData)$description = ""
  
  ## construct "faseqData", CLIP-seq read count
  if (verbose) cat("Generating signals from individual samples...\n")
  featureCount <- count.feature$counts
  rownames(featureCount) <- NULL
  rd <- c_granges(event$feature, sep = "-")
  n_feature <- elementNROWS(event$feature)
  rd$event_id <- rep(event$event_id, n_feature)
  rd$event_name <- rep(event$event_name, n_feature)
  rd$gene_id <- rep(event$gene_id, n_feature)
  rd$transcript_id <- rep(event$transcript_id, n_feature)
  rd$feature_name <- factor(unlist(sapply(event$feature, names)), surf.features)
  faseqData = SummarizedExperiment(
    assays = List(count = featureCount),
    rowData = rd,
    colData = CLIPseqSampleData,
  )
  
  ## transform count into signal
  if (signal.type == "raw.count") {
    signal <- count.feature$counts
    f <- 1
  } else {
    rgs <- unlist(event$feature)
    if (signal.type == "FPKM") {
      f <- colSums(count.feature$stat[-1]) * 1e-6 # "per million" scaling factor
      rpk <- t(t(count.feature$counts / f)) # reads per million (RPM)
      signal <- rpkm <- rpk / width(rgs) * 1e3
    }
    if (signal.type == "TPM") {
      rpk <- count.feature$counts / width(rgs) * 1e3 # reads per kilobase (RPK).
      f <- colSums(rpk) * 1e-6 # "per million" scaling factor
      signal <- tpm <- t(t(rpk) / f)
    }
  }
  rownames(signal) <- names(unlist(event$feature))
  
  ## aggregate by condition -- mean, then take the difference: IP - input
  if (verbose) cat("Contrasting IP from SMInput...\n")
  # colnames(signal) = sampleData$condition
  # signal <- signal %>%
  #   melt(varnames = c("feature", "condition")) %>%
  #   group_by(feature, condition) %>%
  #   summarise(aggregate = FUN.aggregate(value)) %>%
  #   group_by(feature) %>%
  #   summarise(contrast = - diff(aggregate))
  signal.IP = signal[,sampleData$condition == sampleData$condition[1], drop = F]
  if (FUN.aggregate == "mean") {
    signal.IP <- rowMeans(signal.IP)
  } else {
    signal.IP <- apply(signal.IP, 1, FUN.aggregate)
  }
  signal.input = signal[,sampleData$condition == sampleData$condition[nrow(sampleData)], drop = F]
  if (FUN.aggregate == "mean") {
    signal.input <- rowMeans(signal.input)
  } else {
    signal.input <- apply(signal.input, 1, FUN.aggregate)
  }
  signal.contrast <- signal.IP - signal.input
  names(signal.contrast) <- names(unlist(event$feature, use.names = F))
  featureSignal = relist(signal.contrast, event$feature)
  
  ## annotate new column's attribute
  ds = DataFrame(featureSignal)
  mcols(ds)$type <- "CLIP-seq"
  mcols(ds)$description <- "normalized CLIP-seq feature signals (contrasted)"
  
  ## add scaling factor to sampleData
  CLIPseqSampleData <- cbind(CLIPseqSampleData, sizeFactor = f)
  mcols(CLIPseqSampleData)["sizeFactor", "type"] = "intermediate"
  mcols(CLIPseqSampleData)["sizeFactor", "description"] = "\"per million\" scaling factor"
  
  res <- new(
    "surf", 
    cbind(as(event, "DataFrame"), ds),
    genePartsList = event@genePartsList,
    drseqData = event@drseqData,
    drseqResults = event@drseqResults,
    faseqData = faseqData, ## newly added
    faseqResults = event@faseqResults,
    daseqResults = event@daseqResults,
    sampleData = event@sampleData
  )
  res@sampleData$"CLIP-seq" <- CLIPseqSampleData
  metadata(res) = metadata(event)
  metadata(res)$signal.type = signal.type
  metadata(res)$FUN.aggregate = FUN.aggregate
  return(res)
}

#' Perform the functional association test
#'
#' This is a learning one unit of SURF.
#' It trains a GLM model for the functional association of one RBP with one ATR event.
#'
#' @param data \code{data.frame}, contains training data for one RBP and one event type
#' @param min.size \code{integer}, the minimum size of "reliable" training set, default to `60`.
#' @param trim \code{numeric}, the percentile used to trim the training data. This is useful in producing a more robust estimation of functional association.
fat = function(data, min.size = 60, trim = 0.025) {
  feature = intersect(colnames(data), c("up3", "up2", "up1", "bd1", "bd2", "dn1", "dn2", "dn3"))
  res = data.frame()
  # inc vs ctrl
  coef.inc = t(sapply(feature, function(f) {
    sub = data[data$group != "decrease", c("group", f)]
    sub = sub[sub[[f]] > quantile(sub[[f]], trim) &
                sub[[f]] < quantile(sub[[f]], 1 - trim),]
    fit.glm = arm::bayesglm(paste("group ~", f), binomial(link = "logit"), sub)
    coef = coef(summary(fit.glm))[2,]
    coef[4] = pnorm(coef[3], lower.tail = F)
    names(coef)[4] = "p.value"
    coef
  }))
  res = rbind(res, cbind(
    feature = feature,
    size = sum(data$group == "increase"),
    as.data.frame(coef.inc),
    functional = "exclusion"))
  
  # dec vs ctrl
  coef.dec = t(sapply(feature, function(f) {
    sub = data[data$group != "increase", c("group", f)]
    sub = sub[sub[[f]] > quantile(sub[[f]], trim) &
                sub[[f]] < quantile(sub[[f]], 1 - trim),]
    fit.glm = arm::bayesglm(paste("group ~", f), binomial(link = "logit"), sub)
    coef = coef(summary(fit.glm))[2,]
    coef[4] = pnorm(coef[3], lower.tail = T)
    names(coef)[4] = "p.value"
    coef
  }))
  res = rbind(res, cbind(
    feature = feature,
    size = sum(data$group == "decrease"),
    as.data.frame(coef.dec),
    functional = "inclusion"))
  
  res = dplyr::filter(res, size >= min.size)
  return(res)
}

#' Functional Association using Sequencing data
#'
#' This function performs functiona association test (FAT).
#' The null hypothesis of FAT is that there is no association between feature signals and differential ATR.
#'
#' @param event a \code{surf} obejcet output by \link{faseqCount}.
#' @inheritParams fat
#' @param verbose \code{logical}, whether to print out progress information.
#' @return a \code{surf} object with \code{faseqResults} slot updated.
#' @export
faseqFit <- function(event,
                     min.size = 60,
                     trim = 0.025,
                     verbose = F) {
  ## check
  stopifnot(all(c("group", "included", "featureSignal") %in% colnames(event)))
  
  ## format data for SURF, group by event_name
  dat = event[event$included, c("group","featureSignal")]
  event_name = event[event$included, "event_name"]
  dat = S4Vectors::split(dat, event_name)
  dat <- dat[!sapply(dat, is.null) & !!sapply(dat, nrow)]
  if (verbose) cat("Testing location features for", length(dat), "events:",
                   paste(names(dat), collapse = " "), "\n")
  
  ## after slipt, ncol() becomes the same within each event type,
  ## thus can coerce into data.frame
  dat <- lapply(dat, function(x) {
    x$featureSignal <- list_rbind(x$featureSignal)
    data.frame(group = x$group, x$featureSignal)
  })
  
  testing <- lapply(dat, fat,
                    min.size = min.size,
                    trim = trim)
  
  res <- list_rbind(testing, save.names = "event")
  res$event <- factor(res$event, surf.events)
  res$padj = p.adjust(res$p.value, method = "fdr")
  res <- DataFrame(res)
  mcols(res)$type = "faseq"
  mcols(res)$description = c(
    "event type/category",
    "positional feature",
    "number of events (sample size in FA test)",
    "estimated feature main effect",
    "estiamted feature standard error",
    "standardized Z value (Gaussian)",
    "p value",
    "inferred regulating function",
    "adjusted p value (BH)"
  )
  if (nrow(res)) {
    rownames(res) <- paste0(res$event, "-", res$feature, ":", res$functional)
  }
  
  metadata(res) <- list(min.size = min.size,
                        trim = trim)
  faseqResults <- new("faseqResults", res)
  
  event@faseqResults <- faseqResults
  metadata(event)$min.size = min.size
  metadata(event)$trim = trim
  
  return(event)
}

#' Functional association inference
#'
#' Inference the functionality of individual locaiton features,
#' where the RBP is likely to interact and regulate the corresponding ATR event.
#' A location feature is inferred as function associated if
#' (i) it is included in FAT, and
#' (ii) the corresponding FAT is significant (padj < cut.off), and
#' (iii) it has strong binding signal (featureSignal > cut.off).
#'
#' @param event a \code{surf} object from \link{faseq} or \link{faseqFit}.
#' @param fdr.cutoff \code{numeric}, significance cutoff for the adjusted p-values of FAT.
#' @param signal.cutoff \code{numeric}, threshold cut-off for the eCLIP signals, default to 20. Set this to 0 if dont wnat to filter those location with low eCLIP signals of the RBP.
#' @return a \code{surf} object, with one added \code{inferredFeature} column (inclusion/exclusion/none).
#' @export
faseqInfer = function(event,
                      fdr.cutoff = 0.05,
                      signal.cutoff = 20) {
  if (any("inferredFeature" %in% colnames(event))) {
    event$inferredFeature = NULL
  }
  
  far <- event@faseqResults
  signal <- unlist(event$featureSignal, use.names = F)
  event_name <- rep(event$event_name, elementNROWS(event$feature))
  group <- rep(event$group, elementNROWS(event$feature))
  included <- rep(event$included, elementNROWS(event$feature))
  testFeature <- paste0(event_name, "-", names(signal))
  inferred <- rep("none", length(signal))
  
  ## infer inclusion features
  subfar <- far[far$padj < fdr.cutoff &
                  far$functional == "inclusion",]
  sigFeature <- paste0(subfar$event, "-", subfar$feature)
  inferred[testFeature %in% sigFeature &
             signal > signal.cutoff &
             group == "decrease" &
             included] <- "inclusion"
  
  ## infer exclusion features
  subfar <- far[far$padj < fdr.cutoff &
                  far$functional == "exclusion",]
  sigFeature <- paste0(subfar$event, "-", subfar$feature)
  inferred[testFeature %in% sigFeature &
             signal > signal.cutoff &
             group == "increase" &
             included] <- "exclusion"
  
  ## recode into factor & relist
  inferred <- factor(inferred,
                     c("inclusion", "exclusion", "none"))
  inferred <- FactorList(relist(inferred, event$featureSignal))
  
  ds <- DataFrame(inferredFeature = inferred)
  mcols(ds)$type <- "faseq"
  mcols(ds)$description <- "inferred functionality of location features"
  
  res <- new(
    "surf", 
    cbind(as(event, "DataFrame"), ds),
    genePartsList = event@genePartsList,
    drseqData = event@drseqData,
    drseqResults = event@drseqResults,
    faseqData = event@faseqData,
    faseqResults = event@faseqResults,
    daseqResults = event@daseqResults,
    sampleData = event@sampleData
  )
  metadata(res) = metadata(event)
  metadata(res)$faseq.fdr = fdr.cutoff
  metadata(res)$signal.cutoff = signal.cutoff
  return(res)
}

#' DASeq
#'
#' Perform the functional association analysis (DASeq) in a single command.
#' This function is a wrapper that calls the necessary functions in order for DASeq.
#'
#' @inheritParams faseqCount
#' @param ... parameters for \link{Rsubread::featureCounts}.
#' @inheritParams faseqFit
#' @inheritParams faseqInfer
#' @return a \code{surf} object DASeq results updated.
#' @references Chen, F., & Keles, S. (2020). SURF: integrative analysis of a compendium of RNA-seq and CLIP-seq datasets highlights complex governing of alternative transcriptional regulation by RNA-binding proteins. *Genome Biology*, 21(1), 1-23.
#' @export
faseq <- function(event,
                  ## data
                  sampleData,
                  signal.type = "FPKM",
                  FUN.aggregate = "mean",
                  minMQS = 10,
                  minOverlap = 12,
                  isPairedEnd = T,
                  cores = max(1, detectCores()-2),
                  ...,
                  
                  ## fit
                  min.size = 100,
                  trim = 0.025,
                  
                  ## inference
                  fdr.cutoff = 0.05,
                  signal.cutoff = 20,
                  
                  verbose = F) {
  
  if (verbose)
    cat("Counting CLIP-seq reads (FASeq data)...\n")
  event <- faseqCount(
    event, sampleData,
    signal.type = signal.type,
    FUN.aggregate = FUN.aggregate,
    minMQS = minMQS,
    minOverlap = minOverlap,
    isPairedEnd = isPairedEnd,
    cores = cores,
    verbose = verbose, ...
  )
  
  if (verbose)
    cat("Performing functional association test (FAT)...\n")
  event <- faseqFit(
    event,
    min.size = min.size,
    trim = trim,
    verbose = verbose
  )
  
  if (verbose)
    cat("Inferencing functional association...\n")
  event <- faseqInfer(
    event,
    fdr.cutoff = fdr.cutoff,
    signal.cutoff = signal.cutoff
  )
  
  return(event)
}


## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## ------ daseq ------
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## SURF Discovery Module 1

#' Get target set based on \code{inferredFeature}
#' @return list of character, sets of target identifiers.
getTargetSet <- function(event,
                         id_column = "transcript_id",
                         verbose = F) {
  targeted <- sapply(event$inferredFeature != "none", any)
  event_targeted <- event[targeted,]
  id <- split(event_targeted[[id_column]],
              paste0(event_targeted$event_name))
  lapply(id, unique)
}

#' Get control sets
#' @return list of character, find a control set for each one of targetSets.
getControlSet <- function(event, targetSets,
                          id_column = "transcript_id",
                          verbose = F) {
  isoPL <- event@genePartsList
  lapply(targetSets, function(targetSet) {
    targeted <- selectMethod("%in%", "Vector")(isoPL[[id_column]], targetSet)
    any_targeted <- sapply(targeted, any)
    unlist(isoPL[[id_column]][any_targeted], use.names = F)
  })
}

#' Build gene/transcript rankings for each sample
#'
#' Builds the "rankings" for each sample: expression-based ranking for all the genes/transcripts in each sample
#' The genes/transcripts with same expression value are shuffled. Therefore, genes/transcripts with expression '0' are randomly sorted at the end of the ranking.
#' These "rankings" can be seen as a new representation of the original dataset. Once they are calculated, they can be saved for future analyses.
#'
#' @param exprMat Expression matrix (genes/trnascripts as rows, samples as columns)
#' The expression matrix can also be provided as one of the Bioconductor classes:
#' \itemize{
#' \item \link{RangedSummarizedExperiment} and derived classes:
#' The matrix will be obtained through assay(exprMatrix),
#' -which will extract the first assay (usually the counts)-
#' or the assay name given in 'assayName'
#' \item \link[Matrix]{dgCMatrix-class}:
#' Sparse matrix
#' \item \code{ExpressionSet}:
#' The matrix will be obtained through exprs(exprMatrix)
#' }
#' @param cores \code{integer}, number of computing workers.
#' @inheritParams AUCell::AUCell_buildRankings
#' @return a \code{SummarizedExperiment} object.
#' @export
getRankings <- function(exprMat,
                        plotStats = F,
                        cores = max(1, detectCores() - 2),
                        verbose = F) {
  
  
  ## rank exprMat
  rankings <- AUCell::AUCell_buildRankings(
    as.matrix(exprMat),
    plotStats = plotStats,
    nCores = cores,
    verbose = verbose
  )
  names(dimnames(assays(rankings, withDimnames=FALSE)$ranking)) =
    c("genomic feature", "sample")
  SummarizedExperiment(assays(rankings))
}

#' Calculate AUC
#' @param set a `list` of sets (or signatures) to test. The sets should be provided as `GeneSet`, `GeneSetCollection` or `character` list.
#' @param cores integer, number of computing cores to use.
#' @inheritParams AUCell::AUCell_calcAUC
#' @return a \code{SummarizedExperiment} object
calculateAUC <- function(set, rankings,
                         cores = 1,
                         verbose = F, ...) {
  if (!is.list(set) || is.null(names(set)))
    stop("set must be a named list.")
  if (!length(set))
    return(SummarizedExperiment())
  AUC <- AUCell::AUCell_calcAUC(
    set,
    AUCell::aucellResults(rankings),
    nCores = cores,
    verbose = verbose, ...
  )
  names(dimnames(assays(AUC, withDimnames=FALSE)$AUC)) = c("set", "sample")
  SummarizedExperiment(assays(AUC))
}

#' Aggregate AUC by sample condition
#'
#' @param object a \code{SummarizedExperiment} output from \link{calculateAUC}
#' @param FUN.aggregate function, used for aggreating AUC within `conditon`, default to \code{mean}.
#' @return a \code{data.frame} with 4 columns: \code{set}, \code{condition.1}, \code{condition.2}, and \code{diff}.
aggregateAUCbyCondition <- function(object, sampleData,
                                    FUN.aggregate = "median") {
  conditions <- levels(sampleData$condition)
  getAUC(object) %>%
    reshape2::melt(value.name = "AUC") %>%
    dplyr::mutate(set = as.character(set),
                  sample = as.character(sample)) %>%
    left_join(data.frame(sampleData), by = "sample") %>%
    group_by(set, condition) %>%
    summarise(AUC = match.fun(FUN.aggregate)(AUC)) %>%
    ungroup() %>%
    pivot_wider(id_cols = set,
                names_from = condition,
                values_from = "AUC") %>%
    # dplyr::select(set, conditions[1], conditions[2]) %>%
    dplyr::mutate(diff = .[[2]] - .[[3]])
  
  # AUC <- getAUC(object)
  # set <- rownames(AUC)
  # condition <- sampleData(object)[set, "condition"]
  # aggrAUC <- apply(getAUC(auc), 1, function(x) {
  #   a <- aggregate(x, list(condition = condition), "mean")
  #   setNames(a$x, a$condition)
  # })
  # res <- cbind(set = set, t(aggrAUC))
  # res$diff = res[[2]] - res[[3]] ## is this necessary?
  #
}

#' Differential activity (via AUC)
#'
#' Detect differential activity using the AUC measure and RNA-seq quantification.
#' This unit is a helper of \code{daseq}.
#'
#' @param targetSet \code{character} vector, set of targeted units.
#' @param controlSet \code{character} vector, control set for contrast. If `NULL` (default), the full set of elements in \code{rankings} will be used.
#' @param rankings a \code{SummarizedExperiment} object, with row correponding to transcript or gene and column corresponding to samples. Each element is a expression measure (e.g., TPM).
#' @param n.sample \code{integer}, number of times of controlSet sampling, default to 1000.
#' @param cores \code{integer}, number of computing cores.
#' @param verbose \code{logical}, whether to print out progress report.
#' @inheritParams calculateAUC
#' @return a \code{data.frame} of DASeq results.
diffAUC <- function(targetSet,
                    controlSet = NULL,
                    rankings,
                    sampleData,
                    n.sample = 1000,
                    cores = 1,
                    verbose = F, ...) {
  conditions = levels(sampleData$condition)
  
  ## calculate targetSet AUC
  auc <- calculateAUC(list("NONAME" = targetSet),
                      rankings,
                      cores = 1,
                      verbose = F, ...)
  auc.obs <- aggregateAUCbyCondition(auc, sampleData)
  
  ## sample control sets
  if (is.null(controlSet)) controlSet = rownames(rankings)
  controlSetSamples = lapply(seq_len(n.sample), sample,
                             x = controlSet,
                             size = length(targetSet))
  names(controlSetSamples) = seq_len(n.sample)
  
  ## calculate control AUC dist'n
  auc <- calculateAUC(controlSetSamples,
                      rankings,
                      cores = cores,
                      verbose = F, ...)
  auc.null <- aggregateAUCbyCondition(auc, sampleData)
  
  
  res <- data.frame(
    base = weighted.mean(auc.obs[2:3], table(sampleData$condition)),
    auc.obs[2:3],
    # p.cond1 = 2 * min(mean(auc.obs[[2]] > auc.null[[2]]),
    #                   mean(auc.obs[[2]] < auc.null[[2]])),
    # p.cond2 = 2 * min(mean(auc.obs[[3]] > auc.null[[3]]),
    #                   mean(auc.obs[[3]] < auc.null[[3]])),
    background = mean(auc.null$diff),
    stat = (auc.obs$diff - mean(auc.null$diff)) / sd(auc.null$diff),
    p.value = mean(abs(auc.null$diff) > abs(auc.obs$diff))
  )
  # names(res)[4:5] = paste0("PD.", conditions)
  return(res)
}

#' Differential activity analysis (DASeq)
#'
#' Detect differential activity using the AUC measure and RNA-seq quantification.
#' For more details about methodology, see Details.
#' This function can be used as part of SURF or as well a stand-along analysis.
#' For the former, input a \code{surf} object to \code{event}.
#' For the latter, input \code{targetSets} and optionally \code{controlSets}.
#'
#' @param rankings a \code{SummarizedExperiment} object from \link{getRankings}.
#' @param sampleData \code{data.frame}, external samples, which contain `condition` column, whose row names must match the column names of \code{rankings}.
#' @param event a \code{surf} object from \link{faseqInference} or \link{faseq}.
#' @param target.type \code{character(1)}, either "transcript" or "gene".
#' @inheritParams diffAUC
## @example inst/examples/example_AUCell_buildRankings.R
#' @return a \code{daseqResults} object if \code{targetSets} was given, a \code{surf} object if \code{event} was give.
#' @references Chen, F., & Keles, S. (2020). SURF: integrative analysis of a compendium of RNA-seq and CLIP-seq datasets highlights complex governing of alternative transcriptional regulation by RNA-binding proteins. *Genome Biology*, 21(1), 1-23.
#' @export
daseq <- function(event = NULL, 
                  rankings, sampleData,
                  target.type = "transcript",
                  targetSets = NULL,
                  controlSets = NULL,
                  n.sample = 1000,
                  cores = max(1, detectCores() - 2),
                  verbose = F, ...) {
  ## check rankings and sampleData
  if (ncol(rankings) != nrow(sampleData) ||
      any(colnames(rankings) != rownames(sampleData))) {
    stop("colnames(rankings) and rownames(sampleData) must match.")
  }
  if (is.null(sampleData$condition)){
    stop("sampleData must contain a \"condition\" column indicating two groups to compare; please add.")
  }
  sampleData$condition <- as.factor(sampleData$condition)
  if (nlevels(sampleData$condition) != 2) {
    stop("The condition must have two levels.")
  }
  
  ## standardize sampleData
  externalSampleData = DataFrame(
    sample = rownames(sampleData),
    sampleData[setdiff(names(sampleData), "sample")]
  )
  
  if (is.null(event) && is.null(targetSets)) {
    stop("Either event or targetSets is required.")
  }
  
  ## auto generate target sets
  if (is.null(targetSets)) {
    stopifnot(target.type %in% c("transcript", "gene"))
    if (verbose) cat("Infer target/control", target.type, "sets\n")
    id_column = paste0(target.type, "_id")
    targetSets <- getTargetSet(event,
                               id_column = id_column,
                               verbose = verbose)
    targetSets <- lapply(targetSets, intersect, rownames(rankings))
    if (!length(targetSets)) {
      cat("Halt: there is none SURF-inferred feature.\n")
      return(event)
    }
    controlSets <- getControlSet(event, targetSets,
                                 id_column = id_column,
                                 verbose = verbose)
    controlSets <- lapply(controlSets, intersect, rownames(rankings))
  }
  
  ## check targetSets and controlSets
  if (is.null(names(targetSets))) {
    stop("All target sets should be named.")
  }
  if (is.null(controlSets)) {
    controlSets = lapply(targetSets, as.null)
  }
  if (length(targetSets) != length(controlSets) ||
      any(names(targetSets) != names(controlSets))) {
    stop("Names of target and control sets must match.")
  }
  
  ## @AUC
  if (verbose) cat("Calculating AUC for target sets...\n")
  AUC <- calculateAUC(targetSets,
                      rankings,
                      cores = cores,
                      verbose = F, ...)
  # externalSampleData <- cbind(colData(AUC), externalSampleData)
  colData(AUC) <- externalSampleData
  
  ## non-parametric differential activity test
  if (verbose) cat("Testing for differential activity...\n")
  test.res <- mapply(
    diffAUC, targetSets, controlSets,
    MoreArgs = list(rankings = rankings,
                    sampleData = externalSampleData,
                    n.sample = n.sample,
                    cores = cores,
                    verbose = verbose, ...),
    SIMPLIFY = F
  ) %>% bind_rows()
  
  ## collect results
  res = DataFrame(id = names(targetSets),
                  size = sapply(targetSets, length),
                  set = List(targetSets),
                  control = List(controlSets),
                  test.res,
                  row.names = names(targetSets))
  res$padj <- p.adjust(res$p.value, method = "fdr")
  
  ## annotate attributes in mcols()
  mcols(res)$type = "daseq"
  conditions = levels(colData(rankings)$condition)
  mcols(res)$description = c(
    "target set identifier",
    "target set size",
    "target set",
    "full control set",
    "base AUC of target set",
    paste("target set activity (AUC) in", conditions[1]),
    paste("target set activity (AUC) in", conditions[2]),
    # paste("p-value of distinctive activity (target vs. control) in" , conditions[1]),
    # paste("p-value of distinctive activity (target vs. control) in" , conditions[2]),
    "background difference in activity (control set)",
    "statistic (parametric analogue, reference only)",
    paste0("p-value of differential activity (", conditions[1], " vs. ",conditions[2], ") contrasted to background"),
    "adjusted p-values"
  )
  metadata(res) = list(target.type = target.type,
                       n.resample = n.sample)
  
  daseqResults <- new("daseqResults", res, AUC = AUC)
  
  if (!is.null(event)) {
    event@daseqResults <- daseqResults
    event@sampleData$"External" <- externalSampleData
    metadata(event) <- metadata(event)
    metadata(event)$target.type = target.type
    metadata(event)$n.resample = n.sample
    return(event)
  } else {
    return(daseqResults)
  }
}

