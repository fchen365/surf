## SURF ATR Event Parsing

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## ------ 0. main function ------
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#' Parse ATR events from genome annotation
#'
#' This function parses ATR events (including AS, ATI, and APA) from genome 
#' annotation. It also identifies location features for each event. The latter 
#' task could be computationally demanding for unfiltered (raw) genome 
#' annotation.
#'
#' @param anno.file `character`, directory to genome annotation file.
#' @param anno.format `character`, the format of the annotation file. The format 
#'     can be inferred from `anno.file` automatically, unless it is not implicit 
#'     by the extension.
#' @inheritParams getIsoPartsList
#' @inheritParams getEvent
#' @inheritParams annotateEvent
#' @inheritParams getFeature
#' @param location.feature `logical`, whether (default to `TRUE`) to add 
#'     location features for each event. This usually takes the longest time in 
#'     annotation parsing procedure.
#' @param verbose `logical`, whether (default to `TRUE`) to print out progress.
#' @return a `surf` object with ATR event annotation and updated `genePartsList` 
#'     slot.
#' @references Chen, F., & Keles, S. (2020). SURF: integrative analysis of a 
#'     compendium of RNA-seq and CLIP-seq datasets highlights complex governing 
#'     of alternative transcriptional regulation by RNA-binding proteins. 
#'     *Genome Biology*, 21(1), 1-23.
#' @export
parseEvent <- function(anno.file,
                       anno.format = tools::file_ext(anno.file),
                       cores = max(1, detectCores() - 2),
                       min.event.length = 6,
                       location.feature = TRUE,
                       depth.exon = 100,
                       depth.intron = 300,
                       remove.duplicate = TRUE,
                       verbose = FALSE) {
  ## import annotation
  if (anno.format != tools::file_ext(anno.file))
    warning("Annotation format is not implicit from file extension.")
  if (verbose)
    cat("Importing genome annotation...\n")
  anno <- import(con = anno.file, format = anno.format)
  
  if (verbose)
    cat("Extracting gene parts list... ")
  timer <- system.time({
    isoPL = getIsoPartsList(anno, depth.intron = depth.intron, cores = cores)
    isoPLas = getEvent(isoPL, cores = cores)
  })
  if (verbose)
    cat("Running time:", timer[3], "sec.\n")
  if (verbose)
    cat("Annotating AS/ATS/APA events...")
  ## (optional) start/stop codon, used to exclude some interior events 
  ## (SE/RI/A3SS/A5SS) that overlap
  anno_ss = anno[anno$type %in% c("start_codon", "stop_codon")] 
  ## staring or stoping sites
  timer <- system.time({
    anno_event_nofeat = annotateEvent(
      isoPLas,
      cores = cores,
      min.event.length = min.event.length,
      anno_ss = anno_ss,
      remove.duplicate = FALSE
    )
  })
  if (verbose)
    cat("Running time:", timer[3], "sec.\n")
  if (location.feature) {
    if (verbose)
      cat("Extracting feature windows...\n")
    timer <- system.time({
      anno_event = getFeature(
        isoPLas,
        anno_event_nofeat,
        depth.exon = depth.exon,
        depth.intron = depth.intron,
        remove.duplicate = remove.duplicate,
        cores = cores,
        verbose = verbose
      )
    })
    if (verbose)
      cat("Running time:", timer[3], "sec.\n")
  }
  
  return(anno_event)
}

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## ------ 1. generate isoform parts list ------
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#' Find multi-transciption gene
#'
#' Find genes that contain multiple transcripts.
#' @param anno `data.frame` or `GRanges`, genome annotation.
#' @param min `integer`, minimum number of transcripts.
#' @param max `integer`, maximum number of transcripts..
#' @return a `character` vector of gene identifiers that possess multiple 
#'     transcripts.
getMultiTxGene = function(anno, min = 2, max = Inf) {
  anno_tx = anno[anno$type == "transcript"]
  cnt_tx = table(anno_tx$gene_id)
  names(cnt_tx)[cnt_tx >= min & cnt_tx <= max]
}

#' Count number of `TRUE` from a logical rle
#' @param x logical rle
#' @return a integer rle
countLogicRle = function(x) {
  runs = rle(x)
  runs$values[runs$values] = seq_len(sum(runs$values))
  label = inverse.rle(runs)
}

#' Parse isoform parts list from annotation
#'
#' Parse the isoform parts list from annotation.
#'
#' @param anno a `GRanges` object of genome annotation, return by [import].
#' @param gene_id `character`, gene_id"s to analyze, default to all 
#'     multi-transcript genes.
#' @param depth.intron `integer`, extended depth into gene's flanks, default to 
#'     300 nt.
#' @param cores `integer`, number of computing cores to use.
#' @return a `list` named by gene_id, each of which contains:
#' \item{segment}{`GRanges`, genomic data of parts list.}
#' \item{label}{`integer` vector, gene model (0 for intron) and exon numbers 
#'     (coded as integer 1,2,...).}
#' \item{layout}{`lgCMatrix`, transcript structure.}
#' @export
getIsoPartsList = function(anno,
                           gene_id = getMultiTxGene(anno),
                           depth.intron = 300,
                           cores = max(1, detectCores() - 2)) {
  ## check input
  if (!is(anno, "GRanges")) {
    anno = GRanges(anno)
    warning("anno is coerced to GRanges.")
  }
  if (depth.intron < 1)
    stop("\"depth.intron\" has to be greater than 0.")
  
  registerDoParallel(cores)
  isoPartsList = foreach (g = gene_id) %dopar% {
    anno_g = anno[anno$gene_id == g]
    gene = anno_g[anno_g$type == "gene"]
    tx = anno_g[anno_g$type == "transcript"]
    exon = anno_g[anno_g$type == "exon"]
    intron = c(
      GenomicRanges::setdiff(gene, exon),
      ## exon-exon junction + exonic region of filtered isoform
      GenomicRanges::flank(gene, depth.intron),
      GenomicRanges::flank(gene, depth.intron, start = FALSE)
    ) ## add up/down-stream of genes
    strand = as.character(strand(gene))
    exon_pooled = GenomicRanges::disjoin(exon) 
    ## these are bins that DEXSeq works on (except overlapping exons)
    intron_pooled = GenomicRanges::disjoin(intron)
    
    ## valid check
    if (length(tx) == 1) {
      warning(paste("Gene", g, "has only one transcript.")) 
      return(list())
    }
    
    ## $segment - GRanges, genomic data of parts list
    segment = sort(c(exon_pooled, intron_pooled)) ## sort exonic & intronic segments
    
    ## $label - integer, exon numerator
    label = !!countOverlaps(segment, exon_pooled) ## logical rep.
    label = countLogicRle(label)
    
    ## layout - Matrix, isoform structure
    ## group exons by tx, then arrange by transcript_id
    exon_tx = S4Vectors::split(exon, exon$transcript_id)[tx$transcript_id]
    seg_tx = findOverlaps(segment, exon_tx) ## map segments to tx
    layout = Matrix::Matrix(FALSE, length(segment), length(tx), sparse = TRUE)
    colnames(layout) = tx$transcript_id
    layout[as.matrix(seg_tx)] = TRUE
    
    ## reverse everything, if "-" strand gene
    if (strand == "-") {
      segment = rev(segment)
      label = (!!rev(label)) * (max(label) + 1 - rev(label))
      layout = layout[nrow(layout):1, ]
    }
    list(segment = segment,
         label = label,
         layout = layout)
  }
  stopImplicitCluster() ## close parallel environment
  names(isoPartsList) = gene_id
  isoPartsList[!!sapply(isoPartsList, length)]
}

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## ------ 2. find ATR events ------
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#' Get exon end site from exon parts list
#'
#' This is a helper function of [getEvent].
#' @param x an element of the isoform part list.
#' @return a `logical` vector of the same length as `x`, indicating the exon end
#'   sites in the isoform part list.
getEES = function(x) {
  nEES = length(x) - match(FALSE, rev(x), nomatch = length(x) + 1) + 1
  rep(c(FALSE, TRUE), c(nEES, length(x) - nEES))
}

#' Reconstruct exon parts list to isoform parts list
#'
#' This is a helper function of [getEvent].
#' @param ex exon parts list.
#' @param gm gene model.
#' @return a `logical` vector of the same length of `gm`, indicating exonic 
#'     (`1`) or intronic (`0`) segments. 
stream = function(ex, gm) {
  for (e in seq_len(length(ex)))
    gm[gm == e] = ex[[e]]
  as.logical(gm)
}

#' Get exon start site from exon parts list
#'
#' This is a helper function of [getEvent].
#' @param x an element of the isoform part list.
#' @return a `logical` vector of the same length as `x`, indicating the exon
#'   start sites in the isoform part list.
getESS = function(x) {
  nESS = match(FALSE, x, nomatch = length(x) + 1)
  rep(c(TRUE, FALSE), c(nESS - 1, length(x) - nESS + 1))
}

#' Get AS events for one isoform (decision tree)
#' @param li boolean layout of i-th isoform.
#' @param gm gene model.
#' @return a `factor` vector of the same length as `li`, indicating the ATR 
#'     events (e.g., "SE", "A5SS", etc) or "." for no event.
getIsoAS = function(li, gm) {
  bVar = xor(li, gm) ## variable bins
  eVar = base::split(bVar, gm) ## variable exon
  if (!min(gm))
    eVar = eVar[-1]
  
  isoAS = ifelse(bVar, "", ".") ## initial
  
  ## Alternative start site 2
  nTSS = min(which(li)) ## transcript start site
  bASS = rep(c(TRUE, FALSE), c(nTSS, length(gm) - nTSS))
  isoAS[bVar & bASS] = "A5U"
  ## Alternative end site 2
  nPAS = max(which(li)) ## polyA site
  bAES = rep(c(FALSE, TRUE), c(nPAS, length(gm) - nPAS))
  isoAS[bVar & bAES] = "TAP"
  ## Skipped exon
  nSE = sapply(eVar, all)
  bSE = gm %in% which(nSE) & !bASS & !bAES
  isoAS[bVar & bSE] = "SE"
  ## Alternative start site 1 (whole exon)
  bAFE = gm %in% which(nSE) & bASS
  isoAS[bVar & bAFE] = "AFE"
  ## Alternative end site 1 (whole exon)
  bIAP = gm %in% which(nSE) & bAES
  isoAS[bVar & bIAP] = "IAP"
  ## Alternative 5" splicing (donors) site
  bEES = lapply(eVar, getEES)
  bA5SS = stream(bEES, gm) & !bSE & !bASS & !bAES
  isoAS[bVar & bA5SS] = "A5SS"
  ## Alternative 3" splicing (acceptors) site
  bESS = lapply(eVar, getESS)
  bA3SS = stream(bESS, gm) & !bSE & !bASS & !bAES
  isoAS[bVar & bA3SS] = "A3SS"
  ## Retained intron
  bRI = bVar & !bASS & !bAES & !bSE & !bA5SS & !bA3SS
  isoAS[bVar & bRI] = "RI"
  
  factor(isoAS, c(".", surf.events)) 
  ## this will be ignored if called by `apply()`
}

#' Enumerate AS events
#'
#' @details AFE, IAP, and consecutive SE events will be merged.
#' @param as a `factor` vector of AS events /isoform.
#' @param gm gene model.
#' @return a `integer` vector of the same length as `as`, indexing the order of 
#'     ATR events. 
numberAS = function(as, gm) {
  as.tmp = as.character(as[!!gm])
  runs <- rle(as.tmp)
  runs$values[runs$values != "."] = seq_len(sum(runs$values != "."))
  runs$values[runs$values == "."] = 0
  runs$values = as.numeric(runs$values)
  num <- rep(0, length(as))
  num[!!gm] = inverse.rle(runs)
  num
}

#' Enumerate AS events
#' This will not merge SE events across consecutive exon
# numberAS = function(as) {
#   ## as - factor of AS events /isoform
#   runs = rle(!(as %in% c(".", "AFE", "IAP")))
#   runs$values[runs$values] = seq_len(sum(runs$values))
#   num = inverse.rle(runs)
#   if (any(as == "AFE")) {
#     num[!!num] = num[!!num] + 1
#     num[as == "AFE"] = 1
#   }
#   num[as == "IAP"] = max(num) + 1
#   num
# }

#' Find AS events from isoform parts list
#' @param isoPL isoform parts `list`
#' @param cores `integer`, number of computing cores to use.
#' @return a `list` as input with added two components:
#' \item{asNum}{integer `matrix`, dim = dim of segment 0 means no event, and 
#'     other integers indicate event number.}
#' \item{asName}{a `list` of named `character` vector, with each element 
#'     corresponding to one isoform. name is event number, character is event 
#'     name}
getEvent = function(isoPL, cores = max(1, detectCores() - 2)) {
  registerDoParallel(cores)
  isoPLas = foreach (pl = isoPL, g = names(isoPL)) %dopar% {
    segment = pl$segment
    label = pl$label
    layout = pl$layout
    
    ## extract AS events, compress
    ## apply coerce factor columns to character
    asEvent = data.frame(apply(layout, 2, getIsoAS, gm = label)) 
    asNum = data.frame(sapply(asEvent, numberAS, gm = label))
    asName = mapply(split, asEvent, asNum, SIMPLIFY = FALSE)
    asName = lapply(asName, function(x) sapply(x, unique))
    
    ## check isolation
    asSplit = do.call("c", mapply(split, asEvent, asNum, SIMPLIFY = FALSE))
    asSplit = asSplit[!(names(asSplit) %in% c(".", "AFE", "IAP"))]
    asCheck = sapply(asSplit, function(x) length(unique(x)))
    if (any(asCheck != 1))
      message("Warning: ", g, " has adjacent AS events. Please check!")
    
    ## append
    append(pl, list(asNum = asNum, asName = asName))
  }
  stopImplicitCluster()
  setNames(isoPLas, names(isoPL))
}

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## ------ 3. annotate event (no features) ------
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#' ATR event name
#'
#' @param plas a plas object containing `asName`
#' @return a named (by event_id) vector of event type
getASName <- function(plas) {
  asName = plas$asName
  tx_id <- rep(names(asName), sapply(asName, length))
  event_name <- unlist(unname(asName), use.names = TRUE)
  event_name <- factor(event_name, c(".", surf.events))
  names(event_name) <- paste(tx_id, names(event_name), sep = "@")
  event_name <- na.omit(droplevels(event_name, exclude = "."))
  return(event_name)
}


#' Annotate AS events (generate feature bins)
#'
#' @param isoPLas isoform parts list with AS events
#' @param cores `integer`, number of computing cores to use.
#' @param min.event.length `numeric` (positive), minimum length of a valid event.
#' @param anno_ss a `GRanges` of start/stop codon, UTR, Selenocysteine, etc; 
#'     used to exclude some interior events (i.e., SE/RI/A3SS/A5SS)
#' @param remove.duplicate `logical`, whether (default to `TRUE`) to remove 
#'     identical event duplicates (by keeping one).
#' @return a `surf` object.
annotateEvent <- function(isoPLas,
                          cores = max(1, detectCores() - 2),
                          min.event.length = 6,
                          anno_ss = NULL,
                          remove.duplicate = FALSE) {
  registerDoParallel(cores)
  event.list = foreach (
    plas = isoPLas,
    g = names(isoPLas),
    .combine = "c"
    ) %dopar% {
      segment = plas$segment
      asNum = plas$asNum
      strand <- as.vector(strand(segment[1]))
      
      ## input check 1: existence of variable bins
      if (any(!sapply(asNum, max))) {
        warning(paste(
          colnames(asNum)[!sapply(asNum, max)], 
          collapse = ", "), 
          ": no (0) variable bin.")
      }
      
      ## collect event name
      event_name <- getASName(plas)
      event_name <- factor(event_name, surf.events)
      
      ## merge variable bins by events
      body = lapply(asNum, function(x) {
        ## remove 0 segments, merge segments by event
        ## note: reduce() will re-order segments by genomic coordinates.
        asSeg <- unlist(GenomicRanges::reduce(S4Vectors::split(
          segment, replace(x, x == 0, NA))))
        if (strand == "-" && length(asSeg)) {
          i <-
            unlist(aggregate(
              seq_along(asSeg),
              by = list(names(asSeg)),
              FUN = rev
            )$x)
        } else
          i <- seq_along(asSeg)
        asSeg <- asSeg[i]
        asSeg$event_part_number = unlist(lapply(
          rle(names(asSeg))$lengths, seq_len))
        return(asSeg)
      })
      
      ## construct GRangesList of events
      event <- c_granges(body, sep = "@")
      event <-
        GRangesList(S4Vectors::split(unname(event), names(event)))
      mcols(event)$event_id = event_id = names(event)
      mcols(event)$event_name = event_name[event_id]
      mcols(event)$gene_id = rep(g, length(event_id))
      mcols(event)$transcript_id = sapply(strsplit(event_id, "@"), head, 1)
      
      ## ---- clean up
      ## (1) event body length (an amino acid spans 3 bps)
      event = event[sapply(width(event), sum) >= min.event.length]
      
      ## (2) when A3SS/A5SS/RI/SE overlap with AFE/ALE, keep the latter
      if (!is.null(anno_ss)) {
        cnt = suppressWarnings(countOverlaps(
          event, anno_ss[anno_ss$gene_id == g]))
        event = event[!mcols(event)$event_name %in% 
                        c("A3SS", "A5SS", "RI", "SE") |
                        !cnt]
      }
      
      ## (3) remove duplicated event:
      ##    with the same (i) `exonic part` and (ii) `event_name`
      if (remove.duplicate) {
        hit <- findOverlaps(event, event)
        hit <- hit[from(hit) < to(hit)]
        hit <-
          hit[mcols(event)$event_name[from(hit)] == mcols(event)$event_name[to(hit)]] ## (ii)
        hit <-
          hit[as.logical(sapply(event[from(hit)] == event[to(hit)], all))] ## (i)
        if (length(hit))
          event <- event[to(hit)]
        else
          event = event[integer(0)]
      }
      
      ## this is the event annotation (w/o features) for one gene
      event
    }
  stopImplicitCluster()
  
  anno_event <- mcols(event.list)
  anno_event$genomicData <- event.list
  mcols(anno_event$genomicData) <- NULL
  
  ## add mcols()
  mcols(anno_event)$type <- "annotation"
  mcols(anno_event)["event_id", "description"] <- "event identifier"
  mcols(anno_event)["event_name", "description"] <- "event type/class"
  mcols(anno_event)["gene_id", "description"] <- "gene/group identifier"
  mcols(anno_event)["transcript_id", "description"] <- 
    "harbouring transcript identifier"
  mcols(anno_event)["genomicData", "description"] <-
    "genomic ranges of the ATR event"
  
  ## keep isoform parts list
  
  genePartsList = DataFrame(
    gene_id = names(isoPLas),
    transcript_id = List(lapply(
      lapply(isoPLas, "[[", "layout"), colnames
    )),
    segment = GRangesList(lapply(isoPLas, "[[", "segment")),
    label = List(lapply(isoPLas, "[[", "label")),
    layout = List(lapply(isoPLas, "[[", "layout")),
    row.names = names(isoPLas)
  )
  
  mcols(genePartsList)$description <- c(
    "gene identifier",
    "transcript identifiers",
    "genomic segments (parts list)",
    "exonic label (0 as intronic)",
    "transcript (isoform) layout indicators"
  )
  
  res <- new("surf",
             anno_event,
             genePartsList = genePartsList,
             sampleData = DataFrameList())
  metadata(res) = list(min.event.length = min.event.length)
  return(res)
}

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## ------ 4. add features ------
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#' Find body features
#'
#' Locate the body feature of AS events and extract the genomic bins.
#' This is a helper function to `getFeature()`.
#'
#' @param plas one element in `isoPLas`.
#' @param type `character`, either "lead" for upstream, or "lag" for downstream.
#' @param depth.exon depth extended into exon, default is 100 nt.
#' @return a `GRangesList`.
findBodyFeature = function(plas,
                           type = c("lead", "lag"),
                           depth.exon = 100) {
  segment = plas$segment
  asNum = plas$asNum
  strand <- as.vector(strand(segment[1]))
  
  ## input check 1: contain variable bins
  if (any(!sapply(asNum, max)) && type == "lead") {
    warning(paste(colnames(asNum)[!sapply(asNum, max)], 
                  collapse = ", "), 
            ": no (0) variable bin.")
  }
  
  body = lapply(asNum, function(x) {
    ## remove 0 segments, merge segments by event
    ## note: reduce() will re-order segments by genomic coordinates.
    asSeg <-
      unlist(GenomicRanges::reduce(S4Vectors::split(
        segment, replace(x, x == 0, NA))))
    if (xor(strand == "-", type == "lag")) {
      i <- cumsum(rle(names(asSeg))$lengths)
    } else {
      i <-
        data.table::shift(cumsum(rle(names(asSeg))$lengths), fill = 0) + 1
    }
    return(asSeg[i])
  })
  body <- as(body, "GRangesList")
  
  ## take flank
  bodyBin = GenomicRanges::flank(
    x = body,
    width = lapply(IRanges::width(body) * -1, pmax, -depth.exon),
    start = type == "lead"
  )
  
  as(bodyBin, "CompressedGRangesList")
}

#' Find adjacent features of given feature
#'
#' Locate the adjacent bins at the desired direction and extract the genomic 
#' bins, given some location features of AS events. This is a helper function to 
#' `getFeature()`.
#' @param plas one isoform element in `isoPLas`.
#' @param body `list`, $id_tx, $event_id = body bins indices.
#' @param type `character`, either"lead" for upstream, or "lag" for downstream.
#' @param depth.exon `integer`, extended depth into exon, default 100 nt.
#' @param depth.intron `integer`, extended depth into intron, default 300 nt.
#' @return a `GRangesList`.
findAdjacentFeature = function(plas,
                               body,
                               type = c("lead", "lag"),
                               depth.exon = 100,
                               depth.intron = 300) {
  segment = plas$segment
  label = plas$label
  layout = plas$layout
  
  adjBin = mapply(function(bd, lt) {
    ## reconstruct isoform part list
    rle <- rle(lt)
    exonic <- rle$values
    rle$values <- seq_along(rle$values)
    label.new <- inverse.rle(rle)
    mergeSeg <-
      unlist(GenomicRanges::reduce(S4Vectors::split(segment, label.new)))
    mergeSeg$exonic <- exonic
    
    ## probe into flank exon/intron
    probe = GenomicRanges::flank(bd, 1, start = type == "lead")
    adjSeg = to(findOverlaps(probe, mergeSeg))
    if (length(adjSeg) != length(bd)) {
      stop("Inproper \"body\" argument or overlapping isoform parts list \"plas\".")
    }
    
    width <-
      ifelse(mergeSeg[adjSeg]$exonic, depth.exon, depth.intron)
    if (!length(width))
      width <- 0
    tmp <- GenomicRanges::flank(x = bd,
                                width = width,
                                start = type == "lead")
    pintersect(tmp, mergeSeg[adjSeg])
  }, body, data.frame(as.matrix(layout)), SIMPLIFY = FALSE)
  
  as(adjBin, "CompressedGRangesList")
}

#' Find features on the (upstream or downstream) constitutive exons
#'
#' Locate the feature on the (upstream or downstream) constitutive exons and 
#' extract the genomic bins, given some body features of AS events. This is a 
#' helper function to `getFeature()`.
#' @param plas one isoform element in `isoPLas`.
#' @param body `list`, $id_tx$event_id = body bins indices.
#' @param type `character`, either "lead" for upstream, or "lag" for downstream.
#' @param depth.exon `integer`, extended depth into exon.
#' @param constitutive `logical`, whether to use consecutive or adjacent 
#'     neighbor exon. The default is `TRUE`.
#' @return a `GRangesList`.
findNeighborFeature = function(plas,
                               body,
                               type = c("lead", "lag"),
                               depth.exon = 100,
                               constitutive = TRUE) {
  segment = plas$segment
  layout = plas$layout
  label = plas$label
  # asNum = plas$asNum
  # asName = plas$asName
  
  strand <- as.vector(strand(segment[1]))
  nExon = max(label)
  Layout = data.frame(as.matrix(layout * label))
  
  ## exon neihbor (names: current exon index, values: neighbor exon index)
  exonIn = apply(Layout, 2, unique) ## exons included
  if (is.matrix(exonIn))
    exonIn = as.data.frame(exonIn)
  ngbExon = lapply(exonIn, function(x) {
    if (constitutive) {
      ## consecutive
      ord =  seq_len(nExon)
      if (type == "lag")
        ord = rev(ord)
      ngb = c()
      p = NA
      for (i in ord) {
        ngb[i] = p
        if (i %in% x)
          p = i
      }
      setNames(ngb, seq_len(nExon))
    } else
      ## simply adjacent
      setNames(seq_len(nExon), data.table::shift(seq_len(nExon), 1, NA, type))
  })
  
  ## exon flank (names: exon index, values: flanking bins)
  flankBins <- lapply(Layout, function(x) {
    if (constitutive) {
      exon <-
        unlist(GenomicRanges::reduce(S4Vectors::split(
          segment, replace(x, x == 0, NA)))) ## this definition is local to tx
      if (xor(strand == "-", type == "lag")) {
        i <- cumsum(rle(names(exon))$lengths)
      } else {
        i <-
          data.table::shift(cumsum(rle(names(exon))$lengths), fill = 0) + 1
      }
      exon <- exon[i]
    } else {
      exon <-
        unlist(GenomicRanges::reduce(S4Vectors::split(
          segment, replace(label, label == 0, NA)
        ))) ## this definition is global to gene
    }
    
    GenomicRanges::flank(
      exon,
      width = -pmin(IRanges::width(exon), depth.exon),
      start = type == "lag"
    )
  })
  
  ## merge
  ngbBin = mapply(function(bd, ne, fb) {
    exon <-
      unlist(GenomicRanges::reduce(S4Vectors::split(segment, replace(
        label, label == 0, NA
      )))) ## this definition is global to gene
    i.bd <-
      to(findOverlaps(bd, exon)) ## index of exon where the events fall on
    if (length(i.bd) != length(bd))
      stop("Each \"body\" should cover only one exon.")
    i.ngb <- ne[as.character(i.bd)]
    bin <- fb[as.character(na.omit(i.ngb))]
    names(bin) <- names(bd)[!is.na(i.ngb)]
    bin
  }, body, ngbExon, flankBins, SIMPLIFY = FALSE)
  
  as(ngbBin, "CompressedGRangesList")
}


#' Generate location features
#'
#' Generate location features for each ATR event.
#'
#' @param isoPLas isoform parts list with AS events.
#' @param anno_event `DataFrame`, event annotation w/o features.
#' @param depth.exon `integer`, extended depth into exon, default 50 nt.
#' @param depth.intron `integer`, extended depth into intron, default 300 nt.
#' @param cores `integer`, number of computing cores to use.
#' @param remove.duplicate `logical`, whether (default to `TRUE`) to remove 
#'     identical event duplicates (by keeping one).
#' @param verbose `logical`, whether (default to `TRUE`) to echo progress
#' @return a `surf` object with one added column called `feature`, which 
#'     contains a `GRangesList` of extracted location features.
getFeature = function(isoPLas,
                      anno_event,
                      depth.exon = 100,
                      depth.intron = 300,
                      cores = max(1, detectCores() - 2),
                      remove.duplicate = TRUE,
                      verbose = FALSE) {
  ## input check
  if (depth.intron < 1 || depth.exon < 1)
    stop("\"depth.intron\" and \"depth.exon\" must be greater than 0.")
  
  feature.names <-
    c("up3", "up2", "up1", "bd1", "bd2", "dn1", "dn2", "dn3")
  FeatName <- list(
    SE = feature.names,
    RI = feature.names[seq(3, 6)],
    A3SS = feature.names[seq(6)],
    A5SS = feature.names[seq(3, 8)],
    AFE = feature.names[seq(3, 8)],
    A5U = feature.names[seq(3, 6)],
    IAP = feature.names[seq(6)],
    TAP = feature.names[seq(3, 6)]
  )
  
  event_ids = anno_event$event_id
  
  registerDoParallel(cores)
  feature = foreach(
    g = names(isoPLas),
    plas = isoPLas,
    .combine = "c") %dopar% {
      ## collect event name
      event_name <- getASName(plas)
      event_id <- names(event_name)
      
      ## filter event
      if (!is.null(event_ids)) {
        event_id_unwanted <- base::setdiff(event_id, event_ids)
        for (i in strsplit(event_id_unwanted, "@")) {
          plas$asNum[, i[1]] = ifelse(plas$asNum[, i[1]] == i[2], 
                                      0, plas$asNum[, i[1]]) 
        }
        event_id <- base::intersect(event_id, event_ids)
      }
      
      ## local features
      bd1 = findBodyFeature(plas, "lead", depth.exon)
      bd2 = findBodyFeature(plas, "lag", depth.exon)
      up1 = findAdjacentFeature(plas, bd1, "lead", depth.exon, depth.intron)
      dn1 = findAdjacentFeature(plas, bd2, "lag", depth.exon, depth.intron)
      ## distal feature
      up3 = findNeighborFeature(plas, bd1, "lead", depth.exon)
      up2 = findAdjacentFeature(plas, up3, "lag", depth.exon, depth.intron)
      dn3 = findNeighborFeature(plas, bd2, "lag", depth.exon)
      dn2 = findAdjacentFeature(plas, dn3, "lead", depth.exon, depth.intron)
      ## collect feature bins
      Bin <- list(
        bd1 = bd1,
        bd2 = bd2,
        up1 = up1,
        dn1 = dn1,
        up3 = up3,
        up2 = up2,
        dn3 = dn3,
        dn2 = dn2
      )
      Bin <- c_granges(lapply(Bin, c_granges, sep = "@"),
                       use.names = FALSE,
                       save.names = "feature")
      Bin <- S4Vectors::split(Bin, names(Bin)) 
      ## note: the order of event_id was mixed up here
      
      ## order features (according to event type)
      featBin <- mapply(function(en, fb) {
        featname <- FeatName[[as.character(en)]]
        fb <- fb[match(featname, fb$feature)]
        names(fb) = featname
        mcols(fb) = NULL
        fb
      }, event_name[event_id], Bin[event_id], SIMPLIFY = FALSE)
      ## add event_name
      featBin <- GRangesList(featBin)
      mcols(featBin)$event_name = event_name[event_id]
      
      ## ----  remove duplicated events
      ## define with the same (i) `feature` and (ii) `event_name`
      if (remove.duplicate) {
        hit <- findOverlaps(featBin, featBin)
        hit <- hit[from(hit) < to(hit)]
        hit <- hit[mcols(featBin)$event_name[from(hit)] == 
                     mcols(featBin)$event_name[to(hit)]] 
        hit <- hit[as.logical(sapply(
          featBin[from(hit)] == featBin[to(hit)], all))] 
        if (length(hit))
          featBin <- featBin[-to(hit)]
      }
      
      ## this is the event annotation (with features) for one gene
      featBin
    }
  stopImplicitCluster()
  
  ## report # of events after dedup
  event_id_feature <- intersect(event_ids, names(feature))
  if (verbose) {
    cat("Add features to", length(event_id_feature), "events.\n")
    print(table("AS/ATI/APA Event distribution:" = anno_event$event_name))
  }
  
  ## add mcols()
  mcols(feature) = NULL
  addCols <- DataFrame(feature = feature[event_id_feature])
  mcols(addCols)$type <- "annotation"
  mcols(addCols)$description <-
    "genomic ranges of the event features"
  
  ## construct new object
  df <- as(anno_event[event_id_feature, ], "DataFrame")
  res <- new(
    "surf",
    cbind(df, addCols),
    genePartsList = anno_event@genePartsList,
    drseqData = anno_event@drseqData,
    drseqResults = anno_event@drseqResults,
    faseqData = anno_event@faseqData,
    faseqResults = anno_event@faseqResults,
    daseqResults = anno_event@daseqResults,
    sampleData = anno_event@sampleData
  )
  metadata(res) = c(
    metadata(anno_event),
    depth.exon = depth.exon,
    depth.intron = depth.intron,
    remove.duplicate.event = remove.duplicate
  )
  return(res)
}
