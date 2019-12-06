## SURF Analysis Module 1

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## ------ construcor ------
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#' Prepare DrSeq Annotation files
#'
#' To use the "useMetaFeatures" functionality of `Rsubread`, we need the GTF input. This function is a helper function that produce the GTF files needed.
#'
#' @param anno GRanges, genome annotation
#' @param event GRanges, AS/ATS/APA event
#' @param anno.prefix prefix of annotation file
#' @param anno.format character, e.g. "gtf", as defined by `rtracklayer::export()`
#' @param remove.overlap.exon - logical, remove overlaping exons across genes (default to FALSE)
#' @param verbose logical, whether (`TRUE`) to echo progess
#' @export
prepDrseqAnno = function(anno, anno_event,
                         anno.prefix = "annotation.drseq",
                         anno.format = "gff2",
                         remove.overlap.exon = F,
                         verbose = T) {
  ## ---- "exon" flattened file ----
  # anno = import(annotation.file, annotation.format)
  if (!all(c("gene", "transcript") %in% anno$type)) anno = addGeneTx(anno)
  output.exon <- paste(anno.prefix, 'exon', anno.format, sep = ".")
  if (verbose) cat("Outputing exon flattened file to", output.exon, "...\n")

  ## merge exons by genes
  exon = anno[anno$type == 'exon']
  exon = GRangesList(split(exon, exon$gene_id))
  exon = unlist(reduce(exon))
  if (remove.overlap.exon) { ## remove overlaping exons across genes (optional)
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
  export(gr1, output.exon, anno.format)

  ## ---- "event" flattened file ----
  # anno_event = import(event.file, format = event.format)
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
  export(gr2, output.event, anno.format)
}


#' Create DrSeq dataset
#'
#' Create a `DEXSeqDataSet` object for DrSeq, so as to use the `DEXSeq` as the engine of parameter estimation.
#' It is possible to output event read counts for later usage.
#'
#' @param anno GRanges, genome annotation
#' @param anno_event GRanges, AS/ATS/APA event annotation
#' @param sampleData data.frame, describes the RNA-seq samples, contains at least two columns -- `bam` and `condition`, the row.names represent sample names
#' @param design data.frame, the design of RNA-seq experiment
#' @param remove.overlap.exon logical, remove overlaping exons across genes (default to FALSE)
#' @param create.annotation logical, whether (`TRUE`) to create a new annotation file for DrSeq. This option is useful when applying DrSeq to a large number of assays, where annotation can be shared globally once generated.
#' @param anno.prefix file names for outputting annotation files
#' @param export.count logical, whether the count matrix should be export to a file
#' @param cores integer, number of available workers, sent to `nthreads` of `featureCounts`
#' @param verbose logical, whether (`TRUE`) to echo progess
#' @param minMQS;isPairedEnd as defined in `featureCounts`, the default is customized for SURF.
#' @param ... additional parameters for `Rsubread::featureCounts`.
#' @details If you used Illumina HiSeq 2000, set `strandSpecific = 2` (reversed strand).
#' @return a `DEXSeqDataSet` object
#' @export
drseqData = function(anno, anno_event, sampleData,
                     design = ~sample + exon + condition:exon,
                     remove.overlap.exon = F,
                     create.annotation = T,
                     anno.prefix = "drseq.annotation", anno.format = "gff2",
                     export.count = F, count.prefix = "drseq.count",
                     minMQS = 10, isPairedEnd = T,
                     cores = max(1, detectCores()-2),
                     verbose = F, ...) {
  # check input
  sampleData <- as.data.frame(sampleData)
  if (is.null(sampleData$bam)) stop("sampleData must contain \"bam\" column.")
  bam.files <- sampleData$bam
  if (is.null(sampleData$condition)) stop("sampleData must contain \"condition\" column.")
  sampleData$condition <- as.factor(sampleData$condition)
  if (length(levels(sampleData$condition))!=2) stop("The condition must have two levels, for IP and input respectively.")
  if (is.null(rownames(sampleData))) {
    warning("sampleData has no row.names input, created automatically.")
    rownames(sampleData) <- paste0(
      as.character(sampleData$condition),
      unlist(lapply(rle(as.character(sampleData$condition))$lengths, seq_len)))
  }

  ## prepare annotation for exons and events
  if (create.annotation) {
    prepDrseqAnno(anno = anno, anno_event = anno_event,
                  anno.prefix = anno.prefix, anno.format = anno.format,
                  remove.overlap.exon = remove.overlap.exon,
                  verbose = verbose)
  }
  exon.file <- paste(anno.prefix, 'exon', anno.format, sep = ".")
  event.file <- paste(anno.prefix, 'event', anno.format, sep = ".")
  if (!file.exists(exon.file)) stop("Cannot find exon annotation file at ", exon.file)
  if (!file.exists(event.file)) stop("Cannot find event annotation file at ", event.file)

  ## count reads on exons ----
  count.exon = featureCounts(
    files = bam.files,
    annot.ext = exon.file, isGTFAnnotationFile = T,
    GTF.featureType = "exonic_part", GTF.attrType = "gene_id",
    allowMultiOverlap = T,
    minMQS = minMQS,
    isPairedEnd = isPairedEnd,
    nthreads = cores,
    verbose = verbose, ...)
  colnames(count.exon$counts) =
    colnames(count.exon$stat)[-1] =
    rownames(sampleData)
  # message("Total count: ", sum(count.exon$counts))

  ## count reads on events ----
  count.event = featureCounts(
    files = bam.files,
    annot.ext = event.file, isGTFAnnotationFile = T,
    GTF.featureType = "exonic_part", GTF.attrType = "event_id",
    allowMultiOverlap = T,
    minMQS = minMQS,
    isPairedEnd = isPairedEnd,
    nthreads = cores,
    verbose = verbose, ...)
  colnames(count.event$counts) =
    colnames(count.event$stat)[-1] =
    rownames(sampleData)

  ## DEXSeqDataSet
  index <- match(anno_event$event_id, rownames(count.event$counts))
  if (any(is.na(index))) {
    stop("Count results can not match the event annotation.")
  } ## just to be sure
  drd <- DEXSeqDataSet(countData = count.event$counts[index,],
                       sampleData = sampleData,
                       design = design,
                       featureID = anno_event$event_id,
                       groupID = anno_event$gene_id)

  ## replace with gene reads count
  cnt.other = count.exon$counts[rowData(drd)$groupID,] - counts(drd)[,seq_along(bam.files)]
  counts(drd)[,length(bam.files)+seq_along(bam.files)] = cnt.other
  rowMin <- rowMin(counts(drd))
  if (any(rowMin < 0)) {
    warning(sum(rowMin < 0), " event(s) prompt to invalid REU coefficients.")
    drd = drd[rowMin >= 0,]
  }

  if (export.count) {
    output.event <- paste0(count.prefix, ".txt")
    if (verbose) cat("Export count matrix to", output.event, "\n")
    write.table(counts(drd), output.event, col.names = F, quote = F)
  }

  res <- new("drseqData", anno_event, DEXSeqData = drd)
  res@params <- c(res@params, remove.overlap.exon = remove.overlap.exon)
  return(res)
}


#' DrSeq
#'
#' Fit DrSeq models using the DEXSeq engine.
#' @param drd a DEXSeqDataSet object created by `drseqData`.
#' @param cores integer, number of computing workers.
#' @param verbose logical, whether (`TRUE`) to echo progess.
#' @param ... addtional parameters for `DEXSeq`.
#' @return a `drseqResults` object.
#' @export
drseq <- function(drd,
                  cores = max(1, detectCores()-2),
                  verbose = F, ...) {
  dxd <- drd@DEXSeqData ## DEXSeqDataSet
  BPPARAM = MulticoreParam(workers = cores)
  dxdl = split(dxd, drd[mcols(dxd)$featureID, "event_name"])

  if (verbose) cat("Fitting DrSeq for", paste(names(dxdl), collapse = ", "), "...")
  t <- system.time({
    dxrl <- lapply(dxdl, DEXSeq,
                   BPPARAM = BPPARAM, quiet = !verbose, ...)
  })
  if (verbose) cat("Running time:", t[3], "sec. \n")

  ## added columns
  dxr = do.call("rbind", dxrl)
  index <- match(drd$event_id, dxr$featureID)
  if (any(is.na(index))) {
    stop("Event ID's do not match between DrSeq and annotation.")
  }
  dxr_sub <- dxr[index, c(3:10, 12)]
  dxr_sub$padj = p.adjust(dxr_sub$pvalue, method = "fdr") ## optional: re-control FDR
  names(dxr_sub)[1] <- "eventBaseMean"
  names(dxr_sub)[8] <- "logFoldChange" ## log2(condition of the 1st sample / another condition)
  mcols(dxr_sub)$type <- "DrSeq result"
  mcols(dxr_sub)["eventBaseMean", "description"] <- "mean of the counts across samples in each event"
  mcols(dxr_sub)["dispersion", "description"] <- "event dispersion estimate"
  mcols(dxr_sub)["countData", "type"] <- "RNA-seq"

  ## representation
  modelFrameBM = lapply(dxrl, slot, "modelFrameBM")
  sampleData = dxrl[[1]]@sampleData
  dispersionFunction = lapply(dxrl, slot, "dispersionFunction")
  metadata = do.call("c", lapply(dxrl, slot, "metadata"))

  drr <- new("drseqResults", cbind(drd, dxr_sub),
             modelFrameBM = DataFrameList(modelFrameBM),
             sampleData = sampleData,
             dispersionFunction = List(dispersionFunction),
             params = drd@params)
  drr@metadata <- c(drd@metadata, metadata)

  return(drr)
}

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## ------ methods ------
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#' Plot Dispersion Functions
#' @param drr a `drseqResults`` object.
#' @return a `ggplot` object.
#' @details By default, ATR type are colored with the `Paired` palette (see \url{http://colorbrewer2.org}).
#' @export
plotDispFunc <- function(drr, x.limits = c(1e-1, 1e4),
                         colrs = setNames(surf.colors, surf.events)) {
  func <- drr@dispersionFunction
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


#' Volcano plot
#'
#' Create a volcano plot for the DrSeq results, stratefied by alternative transcriptional regulation (ATR) event types.
#' A volcano plot is a scatter plot of tested units, where log2 fold change is in x-axis, and -log10(p.value) is in y-axis.
#' @param drr a `drseqResults` object.
#' @param lfc.cutoff a numeric vector of length 2, the range of log2 fold change that is consider NOT significant.
#' @param fdr.cutoff numeric, significance level of adjusted p-value.
#' @param x.limits a numeric vector of length 2, range of log2 fold change. Any values beyond this range will be projected onto the boundary.
#' @param y.limits a numeric vector of length 2, range of -log10(p.value). Any values beyond this range will be projected onto the boundary.
#' @param remove.portion.grey numeric between 0 and 1, the portion of non-significant points to be randomly remove. This is only for speeding up plotting.
#' @param remove.portion.color numeric between 0 and 1, the portion of significant points to be randomly remove. This is only for speeding up plotting.
#' @return a `ggplot` object.
#' @details By default, ATR type are colored with the `Paired` palette (see \url{http://colorbrewer2.org}).
#' @export
volcano.plot <- function(drr,
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
  dat <- data.frame(x = drr$logFoldChange,
                    y = -log10(drr$padj),
                    group = drr$event_name) %>%
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

