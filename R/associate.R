## SURF Analysis Module 2

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## ------ 1. prepare SURF data ------
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#' Screen SURF training samples
#'
#' Extract response and covariates for the analysis module 2 of SURF.
#'
#' @param event a `drseqResults` object output by by `drseq()` or a `surfData` object.
#' @param drseq.fdr FDR (BH procedure) adjusted p-value cut-off.
#' @param read.length numeric, RNA-seq read length. Default is 100 bp (e.g., Illumina TruSeq). This is used to adjust event base count, which is then used to select the representative events if replicated.
#' @param min.adjMean numeric, adjusted event base mean threshold.
#' @param filter.overlap.event logical, whether (default to `TRUE`) to select one representitive event from overlapping ones and remove the others.
#' @param verbose logical, whether (default to `FALSE`) to print out basic summary statistics.
#'
#' @return a `surfData` object, with added one columns to `event`,
#' \item{adjMean}{adjusted base mean of the event from RNA-seq data.}
#' \item{group}{group labels of ATR events, `inc` for increased REU upon RBP knock-down, and `dec` for decreased, and `ctrl` for no-changed.}
#' \item{included}{logical, indicating whether the event is included into SURF ananlysis module 2.}
#'
#' @export
screenSurfEvent = function(event,
                           drseq.fdr = .05,
                           read.length = 100,
                           min.adjMean = .05,
                           filter.overlap.event = T,
                           verbose = F) {
  if (any(c("adjMean", "group", "included") %in% colnames(event))) {
    event$adjMean = NULL
    event$group = NULL
    event$included = NULL
  }
  ## event read coverage
  eventLen = sapply(width(event$genomicData), sum)
  adjMean = event$eventBaseMean / (eventLen + (read.length-1))
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
  if (verbose) cat("Identified", length(id_evtDU), "DR events and",
                   length(id_evtEU), "ER events.\n")
  id_evt = union(id_evtDU, id_evtEU)

  ## pooling duplicated events by max{adjMean}
  ## (i) same event_name
  ## (ii) overlap more than 0%
  ## (iii) (currently not used) response and covariates are the same
  if (filter.overlap.event) {
    event1 = event[event$event_id %in% id_evt,]
    adjMean1 <- adjMean[event$event_id %in% id_evt]
    event1 = event1[order(adjMean1, decreasing = F),] ## sort by padj
    hits = findOverlaps(event1$genomicData, event1$genomicData)
    hits = hits[from(hits) > to(hits)] ## from later to former
    hits = hits[event1$event_name[from(hits)] == event1$event_name[to(hits)]] ## same event names
    if (verbose) cat(length(unique(from(hits))), "overlapping events.\n")
    id_evt = setdiff(id_evt, event1$event_id[from(hits)])
  }

  if (!length(id_evt)) stop("No usable SURF instance!")
  included = event$event_id %in% id_evt

  ## new columns
  ds = DataFrame(adjMean, group, included)
  mcols(ds)$type = c("RNA-seq", "SURF", "SURF")
  mcols(ds)$description <- c(
    "adjusted base mean of each event",
    "direction of differential regulation (upon knock-down)",
    "indicator of usable event")

  ## retrieve @shRNAsample
  if (.hasSlot(event, "sampleData")) {
    shRNAsample = event@sampleData
  } else if (.hasSlot(event, "shRNAsample")) {
    shRNAsample = event@shRNAsample
  } else {
    shRNAsample = DataFrame()
  }

  ## retrieve @eCLIPsample
  if (.hasSlot(event, "eCLIPsample")) {
    eCLIPsample = event@eCLIPsample
  } else {
    eCLIPsample = DataFrame()
  }


  ## wrap into surfData object
  ds <- new("surfData", cbind(event, ds),
            shRNAsample = shRNAsample,
            eCLIPsample = eCLIPsample,
            params = c(event@params,
                       drseq.fdr = drseq.fdr,
                       RNAseq.read.length = read.length,
                       min.adjMean = min.adjMean,
                       filter.overlap.event = filter.overlap.event))
  if (verbose) {
    cat(length(id_evt), "events included for SURF analysis.\n",
        "The distribution of AS/ATI/APA events identified:\n")
    print(table(ds[ds$included,c("event_name", "group")]))
  }
  return(ds)
}


#' Obtain AS/ATS/APA event feature signal
#'
#' Quantify signals on event features from eCLIP-seq bam files.
#'
#' @param event a `drseqResults` or `surfData` object, annotation of events
#' @param sampleData data.frame, must contain two columns -- `bam` and `condition` (for "IP" and "input", the IP should come first, for the ease of constrasting), row.names represent the sample names. `bam` is the file name of CLIP-seq bam. `condition` will be coerced to factor, whose first level will be treated as IP, and the second level as input.
#' @param signal.type character, indicate the type of feature signal wanted, support `TPM` for Transcripts Per Kilobase Million, `FPKM` for Fragments Per Kilobase Million (for paired-end reads) and Reads Per Kilobase Million (for single-end reads), and `raw.count` for raw.count read counts
#' @param FUN.aggregate function used for aggreating signals within `conditon`, default to mean.
#' @param FUN.contrast function used for contrast signals across `conditon`, default to substraction.
#' @param output.feature.signal logical, whether to output a table of feature signals
#' @param output.prefix character, prefix of output filename(s), .feature will be added, used only if `output.feature.anno` or `output.feature.signal` is TRUE
#' @param cores integer, number of available workers, sent to `nthreads` of `featureCounts`
#' @param verbose logical, whether (`TRUE`) to echo progess
#' @param minMQS;minOverlap;isPairedEnd;... parameters for `featureCounts()`. `minMQS` is default to 10, and `minOverlap` is default to 12 (25\% of the typical read length of eCLIP-seq (~50bp)), and `isPairedEnd` is default to `TRUE`. Note: If you use Illumina HiSeq 2000, set `strandSpecific = 2`.
#'
#' @return a `surfData` object.
#' @keywords signal
#' @references \url{https://www.rna-seqblog.com/rpkm-fpkm-and-tpm-clearly-explained/}
#' @export
getFeatureSignal <- function(event, sampleData,
                             signal.type = "FPKM",
                             FUN.aggregate = "mean",
                             output.signal = F, output.prefix = "feature",
                             minMQS = 10,
                             minOverlap = 12,
                             isPairedEnd = T,
                             cores = max(1, detectCores()-2),
                             verbose = F, ...) {
  ## check
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

  if (any(c("featureCount", "featureSignal") %in% colnames(event))) {
    event$featureCount = NULL
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
    verbose = verbose, ...)
  colnames(count.feature$counts) =
    colnames(count.feature$stat)[-1] =
    rownames(sampleData)

  ## eCLIP-seq read count
  featureCount <- count.feature$counts
  rownames(featureCount) <- names(unlist(event$feature, use.names = F))
  featureCount <- relist(featureCount, event$feature)

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

  ## output signal matrix
  if (output.signal) {
    output.signal <- paste0(output.prefix, ".signal.", signal.type, ".txt")
    if (verbose) cat("Save feature signals to", output.signal,"\n")
    write.table(signal, output.signal, col.names = T, quote = F)
  }

  ## aggregate by condition -- mean, then take the difference: IP - input
  if (verbose) cat("Contrasting IP sample(s) with input...\n")
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
  } else signal.IP <- apply(signal.IP, 1, FUN.aggregate)
  signal.input = signal[,sampleData$condition == sampleData$condition[nrow(sampleData)], drop = F]
  if (FUN.aggregate == "mean") {
    signal.input <- rowMeans(signal.input)
  } else signal.input <- apply(signal.input, 1, FUN.aggregate)
  signal.contrast <- signal.IP - signal.input
  names(signal.contrast) <- names(unlist(event$feature, use.names = F))
  featureSignal = relist(signal.contrast, event$feature)

  ## add attribute
  ds = DataFrame(featureCount, featureSignal)
  mcols(ds)$type <- "CLIP-seq"
  mcols(ds)$description <- c(
    "CLIP-seq read counts of each feature",
    "normalized CLIP-seq feature signals (contrasted)")

  ## prepare @eCLIPsample, a DataFrame
  eCLIPsample = DataFrame(sample = rownames(sampleData), sampleData)
  mcols(eCLIPsample)$type = "input"
  mcols(eCLIPsample)$description = ""
  eCLIPsample <- cbind(eCLIPsample, sizeFactor = f)
  mcols(eCLIPsample)["sizeFactor", "type"] = "intermediate"
  mcols(eCLIPsample)["sizeFactor", "description"] = "\"per million\" scaling factor"

  ## retrieve @shRNAsample
  if (.hasSlot(event, "sampleData")) {
    shRNAsample = event@sampleData
  } else {
    shRNAsample = DataFrame()
  }

  params = c(event@params, signal.type = signal.type, FUN.aggregate = FUN.aggregate)
  ds <- new("surfData", cbind(event, ds),
            shRNAsample = shRNAsample,
            eCLIPsample = eCLIPsample,
            params = params)
  return(ds)
}

#' SURF dataset
#'
#' Prepare SURF dataset.
#' Specifically, (i) select SURF event based on read coverage,
#' (ii) quantify feature signals using eCLIP-seq data.
#'
#' @inheritParams screenSurfEvent
#' @inheritParams getFeatureSignal
#' @param ... additional parameters for `getFeatureSignal`.
#' @return a `surfData` object.
#' @export
surfData <- function(event,
                     drseq.fdr = .05,
                     read.length = 100,
                     min.adjMean = .05,
                     filter.overlap.event = T,
                     eCLIPsample,
                     verbose = F, ...) {
  sdat <- screenSurfEvent(event,
                          drseq.fdr = drseq.fdr,
                          read.length = read.length,
                          min.adjMean = min.adjMean,
                          filter.overlap.event = filter.overlap.event,
                          verbose = verbose)
  sdat <- getFeatureSignal(sdat, sampleData = eCLIPsample,
                           verbose = verbose, ...)
  return(sdat)
}

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## ------ 2. functional association ------
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#' Learning one unit of SURF
#'
#' Train a GLM model for one RBPxAS combination. This is a helper function of testPF.
#'
#' @param data data.frame, contains training data for one RBP and one event type
#' @param min.size integer, the minimum size of "reliable" training set, default to `60`.
trainSurf = function(data, min.size = 60, trim = 0.025) {
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

#' Functional Association Testing (FAT)
#'
#' This function tests the functional association between feature signals and differential ATR.
#'
#' @param event a `surfData` obejcet.
#' @param min.size integer, the minimum size of "reliable" training set, default to `60`.
#' @param trim numeric, trimmed quantile, for head and tail `trim` portion of sample will be ignored.
#' @return a `surfResults` object
#' @export
fat <- function(event,
                min.size = 60,
                trim = 0.025,
                verbose = F, ...) {
  ## format data for SURF, group by event_name
  dat = event[event$included, c("group","featureSignal")]
  event_name = event[event$included, "event_name"]
  dat = split(dat, event_name)
  dat <- dat[!sapply(dat, is.null) & !!sapply(dat, nrow)]
  if (verbose) cat("Testing location features for", length(dat), "events:",
                   paste(names(dat), collapse = " "), "\n")

  ## after slipt, ncol() becomes the same within each event type, thus can coerce into data.frame
  dat <- lapply(dat, function(x) {
    x$featureSignal <- list_rbind(x$featureSignal)
    data.frame(group = x$group, x$featureSignal)
  })

  testing <- lapply(dat, trainSurf, min.size = min.size, trim = trim)

  res <- list_rbind(testing, save.names = "event")
  res$event <- as.factor(res$event)
  res$padj = p.adjust(res$p.value, method = "fdr")
  res <- DataFrame(res)
  mcols(res)$type = "SURF"
  mcols(res)$description = c(
    "event type/category",
    "positional feature",
    "sample size (SURF)",
    "estimated feature main effect",
    "estiamted feature standard error",
    "standardized Z value (Gaussian)",
    "p value",
    "hypothesized regulating function",
    "adjusted p value (BH)"
  )
  if (nrow(res)) {
    rownames(res) <- paste0(res$event, "-", res$feature, ":", res$functional)
  }

  new("surfResults", res,
      trainData = event[event$included, -c(15,18,19)],
      params = list(min.size = min.size, trim = trim))
}


## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## ------ 3. functional inference ------
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#' Functional association plot
#'
#' This function provides a ggplot implementation for functional association (FA Plot).
#' @param sres a `surfResults` object.
#' @param event a character vector, event type wanted. In particular, "all" means all event types.
#' @return a ggplot object.
#' @export
fa.plot <- function(sres, plot.event = "all",
                    trim = sres@params$trim,
                    fdr.cutoff = 0.05){
  if (any(plot.event == "all")) plot.event = levels(sres$event)
  levels <- outer(surf.events, surf.features, paste) %>% t %>% as.vector

  ## box plot: feature signals
  trainData <- sres@trainData
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
      left_join(as.data.frame(sres), by = c("event", "feature")) %>%
      ungroup() %>%
      mutate(event = factor(event, surf.events),
             logp = - log10(padj),
             x = factor(feature, surf.features),
             functional = ifelse(is.na(logp), "not tested", as.character(functional)) %>%
               factor(c("exclusion", "inclusion", "not tested")),
             logp = ifelse(is.na(logp), 0, logp))
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
  # ggsave("faplot.pdf", g, width = 8, height = 5)
  return(g)
}

#' Extract genomic ranges inferred by SURF
#'
#' Extract genomic ranges inferred by SURF. a targeted feature is defined by (1) event (2) group (3) inferred (4) signal
#' @param sres SURF data and results
#' @param fdr.cutoff significance cutoff for the adjusted p-values of DrSeq.
#' @param signal.cutoff threshold cut-off for the eCLIP signals, default to 20.
#' Set this to 0 if dont wnat to filter those location with low eCLIP signals of the RBP.
#' @return a GRanges object, all genomic ranges where the RBP is likely to interact and regulate some ATR events.
#' @export
inferTarget = function(sres,
                       fdr.cutoff = 0.05,
                       signal.cutoff = 20) {
  sdat = sres@trainData
  igr <- lapply(which(sres$padj < fdr.cutoff), function(i) {
    event = as.character(sres$event[i])
    feature = sres$feature[i]
    group = ifelse(sres$functional[i] == "inclusion", "decrease", "increase")
    differential <- (sdat$event_name == event) & (sdat$group == group)
    signal <- list_rbind(sdat$featureSignal[differential]) > signal.cutoff
    gr = c_granges(sdat[rownames(signal)[signal[,feature]],"feature"],
                   use.names = T, sep = ":", save.names = "event_id")
    gr = gr[sub(".*:", "",names(gr)) == sres$feature[i]]
    gr$transcript_id <- sdat[gr$event_id, "transcript_id"]
    gr$gene_id <- sdat[gr$event_id, "gene_id"]
    gr$event_name = rep(event, length(gr))
    gr$feature = rep(feature, length(gr))
    gr$functional = rep(sres$functional[i], length(gr))
    gr
  })
  igr = c_granges(igr, use.names = F)
  igr$event_name = factor(igr$event_name, surf.events)
  igr$feature = factor(igr$feature,
                       c("up3", "up2", "up1", "bd1", "bd2", "dn1", "dn2", "dn3"))
  # sum(width(reduce(igr)))
  igr
}
