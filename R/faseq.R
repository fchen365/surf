## SURF Analysis Module 2 

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## ------ constructor ------
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#' Construct FASeq Data Set
#' 
#' This function uantifies feature signals for location features using CLIP-seq data. 
#' You align CLIP-seq reads to the genome and provide FASeq with the resulting bam files. 
#' We will take care of the rest. 
#' 
#' @param event a \code{surf} object.
#' @param sampleData \code{data.frame}, must contain two columns -- `bam` and `condition` (for "IP" and "input", the IP should come first), whose \code{row.names} represent the sample names. `bam` is the file name of CLIP-seq bam. `condition` will be coerced to factor, whose first level will be treated as IP, and the second level as input.
#' @param signal.type \code{character}, indicate the type of feature signal wanted, support `TPM` for Transcripts Per Kilobase Million, `FPKM` for Fragments Per Kilobase Million (for paired-end reads) and Reads Per Kilobase Million (for single-end reads), and `raw.count` for raw.count read counts
#' @param FUN.aggregate \code{function}, used for aggreating signals within \code{conditon}, default to mean.
#' @param cores \code{integer}, number of available workers, sent to \code{nthreads} of \link{featureCounts}
#' @param verbose \code{logical}, whether (default to \code{TRUE}) to echo progess
#' @param minMQS;minOverlap;isPairedEnd;... parameters for \link{featureCounts}. \code{minMQS} is default to 10, and \code{minOverlap} is default to 12 (25\% of the typical read length of eCLIP-seq (~50bp)), and \code{isPairedEnd} is default to \code{TRUE}. 
#' @return a \code{surf} object, with (i) one column \code{featureSignal} added, (2) \code{faseqData} slot updated, and (3) \code{sampleData} slot updated.
#' @details If your sequencing platform is Illumina HiSeq 2000, set \code{strandSpecific = 2}.
#' @keywords feature signal, CLIP-seq
#' @references \url{https://www.rna-seqblog.com/rpkm-fpkm-and-tpm-clearly-explained/}
#' @export
faseqData <- function(event, sampleData, 
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
    "surf", cbind(event, ds), 
    genePartsList = event@genePartsList, 
    drseqData = event@drseqData, 
    drseqResults = event@drseqResults, 
    faseqData = faseqData, ## newly added
    faseqResults = event@faseqResults, 
    daseqResults = event@daseqResults,
    sampleData = event@sampleData
  )
  res@sampleData$"CLIP-seq" <- CLIPseqSampleData
  metadata(res) = c(
    metadata(event), 
    signal.type = signal.type, 
    FUN.aggregate = FUN.aggregate
  )
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
#' @param event a \code{faseqData} obejcet.
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
  dat = split(dat, event_name)
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
  metadata(event) <- c(metadata(event), 
                       min.size = min.size, 
                       trim = trim)
  
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
faseqInference = function(event, 
                          fdr.cutoff = 0.05, 
                          signal.cutoff = 20) {
  
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
    faseq.fdr = fdr.cutoff, 
    signal.cutoff = signal.cutoff
  )
  return(res)
}

#' DASeq
#' 
#' Perform the functional association analysis (DASeq) in a single command. 
#' This function is a wrapper that calls the necessary functions in order for DASeq.
#' 
#' @inheritParams faseqData
#' @param ... parameters for \link{Rsubread::featureCounts}.
#' @inheritParams faseqFit
#' @inheritParams faseqInference
#' @return a \code{surf} object DASeq results updated.
#' @references Chen F and Keles S. "SURF: Integrative analysis of a compendium of RNA-seq and CLIP-seq datasets highlights complex governing of alternative transcriptional regulation by RNA-binding proteins."
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
    cat("Constructing FASeq data...\n")
  event <- faseqData(
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
  event <- faseqInference(
    event, 
    fdr.cutoff = 0.05, 
    signal.cutoff = 20
  )
  
  return(event)
}


## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## ------ methods ------
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#' Functional association plot
#' 
#' This function provides a ggplot implementation for functional association (FA) plot.
#' 
#' @param object a \code{surf} object from \link{faseq} or \link{faseqFit}.
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
#' @param object a \code{surf} object from \link{faseq} or \link{faseqInference}.
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

