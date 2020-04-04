## SURF Discovery Module 1

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## ------ constructor ------
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

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
    targeted <- sapply(isoPL[[id_column]] %in% targetSet, any)
    unlist(isoPL[[id_column]][targeted], use.names = F)
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
  names(dimnames(assays(rankings)$ranking)) = 
    c("genomic feature", "sample")
  SummarizedExperiment(assays(rankings)) 
}

#' Calculate AUC
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
  names(dimnames(assays(AUC)$AUC)) = c("set", "sample")
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
#' @references Chen F and Keles S. "SURF: Integrative analysis of a compendium of RNA-seq and CLIP-seq datasets highlights complex governing of alternative transcriptional regulation by RNA-binding proteins."
#' @export
daseq <- function(rankings, sampleData,
                  event = NULL, 
                  target.type = "transcript",
                  targetSets = NULL, 
                  controlSets = NULL, 
                  n.sample = 1000, 
                  cores = max(1, detectCores() - 2),
                  verbose = F, ...) {
  ## check rankings and sampleData
  if (ncol(rankings) != nrow(sampleData) || 
      any(colnames(rankings) != rownames(sampleData))) 
    stop("colnames(rankings) and rownames(sampleData) must match.")
  if (is.null(sampleData$condition)) 
    stop("sampleData must contain a \"condition\" column indicating two groups to compare; please add.")
  sampleData$condition <- as.factor(sampleData$condition)
  if (nlevels(sampleData$condition) != 2)
    stop("The condition must have two levels.")
  
  ## standardize sampleData 
  externalSampleData = DataFrame(
    sample = rownames(sampleData), 
    sampleData[setdiff(names(sampleData), "sample")]
  )
  
  if (is.null(event) && is.null(targetSets)) 
    stop("Either event or targetSets is required.")

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
  if (is.null(names(targetSets))) 
    stop("All target sets should be named.")
  if (is.null(controlSets)) 
    controlSets = lapply(targetSets, as.null)
  if (length(targetSets) != length(controlSets) || 
      any(names(targetSets) != names(controlSets))) 
    stop("Names of target and control sets must match.")
  
  ## @AUC 
  if (verbose) 
    cat("Calculating AUC for target sets...\n")
  AUC <- calculateAUC(targetSets,
                      rankings,
                      cores = cores,
                      verbose = F, ...)
  # externalSampleData <- cbind(colData(AUC), externalSampleData)
  colData(AUC) <- externalSampleData

  ## non-parametric differential activity test
  if (verbose) 
    cat("Testing for differential activity...\n")
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
    metadata(event) <- c(
      metadata(event), 
      target.type = target.type, 
      n.resample = n.sample
    )
    return(event)
  } else {
    return(daseqResults)
  }
}

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## ------ methods ------
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

