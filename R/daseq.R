## SURF Discovery Module 1


## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## ------ accessor ------
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#' @name getAUC
#' @rdname daseqResults-class
#' @export getAUC
setGeneric(name="getAUC",
           def=function(object) standardGeneric("getAUC"))

#' @aliases getAUC,SummarizedExperiment-method
#' @exportMethod getAUC
setMethod("getAUC",
          signature="SummarizedExperiment",
          definition = function(object) {
            if("AUC" %in% assayNames(object)) {
              SummarizedExperiment::assays(object)[["AUC"]]
            }else{
              stop("This object does not contain an AUC matrix.")
            }
          }
)

#' @rdname daseqResults-class
#' @aliases getAUC,daseqResults-method
#' @exportMethod getAUC
setMethod(
  "getAUC",
  signature = "daseqResults",
  definition = function(object) {
    getAUC(object@targetAUC)
  }
)


## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## ------ functions ------
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#' Build gene/transcript rankings for each sample
#'
#' Builds the "rankings" for each sample: expression-based ranking for all the genes/transcripts in each sample
#' The genes/transcripts with same expression value are shuffled. Therefore, genes/transcripts with expression '0' are randomly sorted at the end of the ranking.
#' These "rankings" can be seen as a new representation of the original dataset. Once they are calculated, they can be saved for future analyses.
#' @inheritParams AUCell::AUCell_buildRankings
#' @param cores integer, number of computing workers.
#' @export
getRankings <- function(exprMat, sampleData,
                        cores = max(1, detectCores()-2),
                        plotStats = F, ...) {
  ## check sample Data
  if (ncol(exprMat) != nrow(sampleData) ||
      any(colnames(exprMat) != rownames(sampleData)))
    stop("colnames(exprMat) and rownames(sampleData) must match.")
  if (is.null(sampleData$condition))
    stop("sampleData must contain a \"condition\" column indicating two groups to compare; please add.")
  sampleData$condition <- as.factor(sampleData$condition)
  if (nlevels(sampleData$condition) != 2)
    stop("The condition must have two levels.")

  rankings <- AUCell::AUCell_buildRankings(as.matrix(exprMat),
                                           plotStats = plotStats,
                                           nCores = cores, ...)
  colData(rankings) <- cbind(colData(rankings), sampleData)
  names(dimnames(assays(rankings)$ranking)) = c("genomic feature", "sample")

  daseqRanking(rankings)
}

#' Differential activity (via AUC)
#'
#' Detect differential activity using the AUC measure and RNA-seq quantification. This is a inference unit/engine of daseq.
#' This unit uses one worker only; daseq cordinates these units parallelly.
#' @inheritParams AUCell::AUCell_buildRankings
#' @inheritParams AUCell::AUCell_calcAUC
#' @param size # of times of controlSet sampling.
#' @export
diffAUC <- function(rankings,
                    targetSet, controlSet = NULL,
                    n.sample = 1000,
                    verbose = F, ...) {
  conditions = levels(colData(rankings)$condition)

  ## calculate targetSet AUC
  auc.target <- AUCell::AUCell_calcAUC(targetSet, rankings,
                               nCores = 1, verbose = verbose, ...)
  dat.obs = AUCell::getAUC(auc.target) %>% reshape2::melt(data = .) %>%
    setNames(c("txSet", "sample", "AUC")) %>%
    mutate(condition = colData(rankings)[sample, "condition"]) %>%
    summarise(stat = t.test(AUC ~ condition)$statistic,
              cond1 = median(AUC[condition == conditions[1]]),
              cond2 = median(AUC[condition == conditions[2]]),
              diff = cond1 - cond2)

  ## sample control sets
  if (is.null(controlSet)) controlSet = rownames(rankings)
  controlSetSamples = lapply(seq_len(n.sample), sample,
                            x = controlSet, size = length(targetSet))
  names(controlSetSamples) = seq_len(n.sample)

  ## calculate AUC
  auc.control <- AUCell::AUCell_calcAUC(controlSetSamples, rankings,
                                nCores = 1, verbose = verbose, ...)
  dat.null = AUCell::getAUC(auc.control) %>% data.table::data.table(.) %>%
    mutate(resampleID = rownames(auc.control)) %>%
    reshape2::melt(data = ., id.vars = c("resampleID")) %>%
    setNames(c("resampleID", "sample", "AUC")) %>%
    mutate(condition = colData(rankings)[sample, "condition"]) %>%
    group_by(resampleID) %>%
    summarize(stat = t.test(AUC ~ condition)$statistic,
              cond1 = median(AUC[condition == conditions[1]]),
              cond2 = median(AUC[condition == conditions[2]]),
              diff = cond1 - cond2)

  res <- data.frame(size = length(targetSet),
                    dat.obs$cond1, dat.obs$cond2,
                    p.cond1 = min(mean(dat.obs$cond1 > dat.null$cond1),
                                  mean(dat.obs$cond1 < dat.null$cond1)),
                    p.cond2 = min(mean(dat.obs$cond2 > dat.null$cond2),
                                  mean(dat.obs$cond2 < dat.null$cond2)),
                    stat = dat.obs$stat,
                    p.value = min(mean(dat.obs$diff > dat.null$diff),
                                  mean(dat.obs$diff < dat.null$diff)))
  names(res)[2:5] = c(conditions, paste0("p.", conditions))
  return(res)
}

#' Differential activity (via AUC)
#'
#' Detect differential activity using the AUC measure and RNA-seq quantification.
#' @return a daseqResults object, which contains a DataFrame and some testing configurations.
#' @export
daseq <- function(rankings,
                  targetSets, controlSets = NULL,
                  n.sample = 1000,
                  cores = max(1, detectCores()-2),
                  verbose = F, ...) {
  ## check sample Data
  if (is.null(colData(rankings)$condition))
    stop("colData(rankings) must contain a \"condition\" column indicating two groups to compare; please add.")
  colData(rankings)$condition <- as.factor(colData(rankings)$condition)
  if (nlevels(colData(rankings)$condition)!=2)
    stop("The condition must have two levels.")

  ## check targetSets and controlSets
  if (is.null(names(targetSets)))
    stop("All target sets should be named.")
  if (is.null(controlSets)) controlSets = lapply(targetSets, as.null)
  if (length(targetSets) != length(controlSets) ||
      any(names(targetSets) != names(controlSets)))
    stop("Names of target and control sets must match.")

  rankings <- methods::as(rankings, "aucellResults")

  ## construct target AUC slot
  AUC <- AUCell::AUCell_calcAUC(
    targetSets, rankings,
    nCores = cores, verbose = verbose, ...)
  targetAUC <- SummarizedExperiment(
    assays(AUC),
    rowData = DataFrame(
      size = sapply(targetSets, length),
      targetSet = List(targetSets)),
    colData = colData(rankings))
  names(dimnames(assays(targetAUC)$AUC)) = c("targetSet", "sample")

  registerDoParallel(cores = cores)
  dar <- foreach(targetSet = targetSets, controlSet = controlSets) %dopar% {
    diffAUC(rankings = rankings,
            targetSet = targetSet, controlSet = controlSet,
            n.sample = n.sample,
            verbose = verbose, ...)
  }
  stopImplicitCluster()
  names(dar) <- names(targetSets)
  dar <- list_rbind(dar, save.names = "targetSet")
  dar$padj <- p.adjust(dar$p.TCGA, method = "fdr")
  dar <- DataFrame(dar)
  mcols(dar)$type = "DASeq"
  conditions = levels(colData(rankings)$condition)
  mcols(dar)$description = c(
    "name of target set",
    "size of target set",
    paste0("observed activity (median AUC) in ", conditions[1]),
    paste0("observed activity (median AUC) in ", conditions[2]),
    paste("p-value of differential activity (target vs. control) in" , conditions[1]),
    paste("p-value of differential activity (target vs. control) in" , conditions[2]),
    "t statistic (for reference only)",
    paste0("p-value of differential activity (", conditions[1], " vs. ",conditions[2], ") based off the control"),
    "adjusted p-values"
  )
  new("daseqResults", dar,
      targetAUC = targetAUC,
      params = list(n.resample = n.sample))
}

#' Heatmap for DASeq results
#' @param dar a `daseqResults` object.
#' @param group an optional `factor` vector, group of `targetSet`.
#' @return a `ggplot` object.
#' @export
heatmapAUC <- function(dar, group = NULL) {
  if (is.null(group)) {
    group = rep("Target set", nrow(dar)) %>% setNames(dar$targetSet)
  } else group = setNames(group, dar$targetSet)

  mat <- getAUC(dar)
  sampleData <- colData(dar@targetAUC)
  clustered_targetSet = clusterByGroup(mat, group)
  clustered_sample <- clusterByGroup(t(mat), sampleData[colnames(mat), "condition"])
  g <- reshape2::melt(mat) %>%
    mutate(condition = sampleData[sample, "condition"],
           group = group[targetSet]) %>%
    ggplot(aes(sample, targetSet, fill = value)) +
    geom_raster() +
    scale_x_discrete(breaks = clustered_sample) +
    scale_y_discrete(breaks = clustered_targetSet) +
    scale_fill_distiller(name = "AUC", palette = "GnBu", direction = 1) +
    facet_grid(rows = vars(group), cols = vars(condition),
               scales = "free", space = "free") +
    labs(x = "Sample", y = "Target transcript set") +
    theme_bw() +
    theme(panel.spacing = unit(.2, "lines"),
          panel.border = element_blank(),
          panel.grid = element_blank(),
          axis.text = element_blank(),
          axis.line = element_blank(),
          axis.ticks = element_blank())
  return(g)
}



#' Scatter plot for DASeq results
#' @param dar a `daseqResults` object.
#' @param alpha numeric, FDR cut-off, default is 0.05.
#' @param delta numeric, lower cut-off for abs. diff. of AUC's between two conditions (default to 0.05).
#' @param min.size numeric, minimum size of the target size to color
#' @param remove.size numeric, minimum size of the target size to even display.
#' @inheritParams ggrepel::geom_label_repel
#' @return a `ggplot` object.
#' @export
scatterAUC <- function(dar,
                       alpha = 0.01, delta = 0.06, min.size = 50,
                       remove.size = 10, force = 5,
                       group = NULL, color = NULL) {
  if (is.null(group)) {
    group = rep("target set", nrow(dar))
  }
  group <- as.factor(group)
  if (length(group) != nrow(dar))
    stop("length(group) != nrow(dar).")
  if (is.null(color)) {
    if (nlevels(group) == 1)
      color <- "#E69F00"
    if (nlevels(group) == 2)
      color <- c("#56B4E9", "#E69F00")
    if (nlevels(group) >= 3)
      color <- RColorBrewer::brewer.pal(nlevels(group), "Set2")
  }
  if (nlevels(group) != length(color))
    stop("length(group) != length(color).")
  if (is.null(names(color))) names(color) = levels(group)
  if (!all(levels(group) %in% names(color)))
    stop("Missing group levels in names(color).")

  ## add "No sig." to group and color
  color = c(color, "No sig." = "grey80")

  ## filter size
  nSmallSize <- sum(dar$size < min.size)
  if (nSmallSize)
    cat("There are", nSmallSize, "small-size target set(s).\n")
  nRemove <- sum(dar$size < remove.size)
  if (nRemove) {
    cat("Remove", nRemove, "target set(s) smaller than", remove.size,".\n")
    group <- group[dar$size >= remove.size]
    dar <- dar[dar$size >= remove.size,]
  }

  cond = levels(colData(dar@targetAUC)$condition)
  names(dar)[3:4] <- c("cond1", "cond2")
  names(dar)[5:6] <- c("p.cond1", "p.cond2")

  dat <- data.frame(dar, group) %>%
    mutate(logp = -log10(padj),
           diff = cond1 - cond2,
           sig = padj < alpha & abs(diff) > delta & size >= min.size,
           color = ifelse(sig, as.character(group), "No sig.") %>%
             factor(c(levels(group), "No sig.")),
           shape = 2 * (p.cond1 < alpha) + (p.cond2 < alpha),
           shape = c("None", cond[2], cond[1], "Both")[shape + 1],
           shape = factor(shape, c("Both", cond[1], cond[2], "None")),
           label = ifelse(sig, as.character(targetSet), ""))

  ## modify if "padj == 0"
  if (any(is.infinite(dat$logp))) {
    which <- is.infinite(dat$logp)
    max <- max(dat$logp[!which])
    dat$logp[which] <- max + runif(sum(which), max = .2)
  }

  g <- dat %>%
    ggplot(aes(cond1, cond2)) +
    geom_abline(slope = 1, intercept = c(-delta, 0, delta),
                color = "grey40", linetype = c(2, 1, 2), alpha = .8) +
    geom_point(aes(color = color, shape = shape, size = size), alpha = .8) +
    scale_color_manual(name = "Signif.", values = color) +
    ggrepel::geom_label_repel(aes(label = label), size = 3, force = force,
                     label.padding = .2, box.padding = 0.4, point.padding = 0.2,
                     segment.size = 0.3, segment.alpha = .8, segment.color = 'grey50',
                     show.legend = F) +
    labs(x = paste0("median AUC (", cond[1], ")"),
         y = paste0("median AUC (", cond[2], ")"),
         shape = "DA within") +
    theme_bw() +
    theme(legend.spacing = unit(0, "lines"),
          legend.box = "horizontal",
          legend.direction = "vertical")

  return(g)
}
