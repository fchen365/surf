## SURF ATR Event Parsing

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## ------ 0. main function ------
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#' Parse AS/ATS/APA events from genome annotation
#'
#' This function parses AS/ATS/APA events from genome annotation.
#' The task could be computationally demanding for unfiltered (raw) genome annotation.
#' @param annotation.file diretory to the genome annotation.
#' @param annotation.format format of annotation file, default is GTF.
#' @param cores integer, number of workers available.
#' @param save.rds logical, whether (`TRUE`) to save intermediate result to .rds file. This is useful when running SURF in batch or large scale essay.
#' @param rds.prefix output prefix, default is empty and will output hidden files.
#' @param min.event.length numeric (positive), minimum length of a valid event.
#' @param depth.exon extended depth into exon, default 50 nt.
#' @param depth.intron extended depth into intron, default 300 nt.
#' @param remove.duplicate logical, whether (`TRUE`) to remove identical event duplicates (and keep one).
#' @param verbose logical, whether (`TRUE`) to print out progress.
#' @export
parseEvent <- function(annotation.file, annotation.format = "gtf",
                       cores = max(1, detectCores() - 2),
                       save.rds = F, rds.prefix = "",
                       min.event.length = 6,
                       depth.exon = 100, depth.intron = 300,
                       remove.duplicate = T,
                       verbose = T) {
  anno = import(annotation.file, annotation.format)
  if (verbose) cat("Parsing isoform parts list... ")
  timer <- system.time({
    isoPL = getIsoPartsList(anno, cores = cores)
    isoPLas = getEvent(isoPL, cores = cores)
  })
  if (save.rds) saveRDS(isoPLas, file = paste0(rds.prefix, ".isoPLas.rds"))
  if (verbose) cat("Running time:", timer[3], "sec.\nAnnotating AS/ATS/APA events... ")
  ## (optional) start/stop codon, used to exclude some interior events (SE/RI/A3SS/A5SS) that overlap
  anno_ss = anno[anno$type %in% c("start_codon", "stop_codon")] ## staring or stoping sites
  timer <- system.time({
    anno_event_nofeat = annotateEvent(isoPLas, cores = cores,
                                      min.event.length = min.event.length,
                                      anno_ss = anno_ss,
                                      remove.duplicate = F)
  })
  if (save.rds) saveRDS(anno_event_nofeat, file = paste0(rds.prefix, ".event_nofeat.rds"))
  if (verbose) cat("Running time:", timer[3], "sec.\nExtracting feature windows... ")
  timer <- system.time({
    anno_event = getFeature(isoPLas, anno_event_nofeat,
                            depth.exon = depth.exon,
                            depth.intron = depth.intron,
                            remove.duplicate = remove.duplicate,
                            cores = cores,
                            verbose = verbose)
  })
  if (save.rds) saveRDS(anno_event, file = paste0(rds.prefix, ".event.rds"))
  if (verbose) cat("Running time:", timer[3], "sec.\nSave event annotation to", paste0(rds.prefix, ".event.rds\n"))

  return(anno_event)
}

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## ------ 1. generate isoform parts list ------
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#' Find multi-transciption gene
#'
#' Find genes that contain multiple transcripts.
#' @param anno data.frame or GRanges, annotation.
#' @param min integer, minimum number of tx"s.
#' @param max integer, maximum number of tx"s.
#' @return a character vector of gene identifiers that possess multiple transcripts.
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
  runs$values[runs$values] = 1:sum(runs$values)
  label = inverse.rle(runs)
}

#' Parse isoform parts list from annotation
#'
#' Parse the isoform parts list from annotation.
#' @param anno GRanges, annotation .
#' @param gene_id character, gene_id"s to analyze, default to all multi-transcript genes.
#' @param depth.intron integer, depth into gene"s flanks, default 300 nt.
#' @param cores integer, number of parallel worker.
#' @return a list named by gene_id, which contains:
#' \item{segment}{GRanges, genomic data of parts list.}
#' \item{label}{integer vector, gene model (0 for intron) and exon numbers (coded as integer 1,2,...).}
#' \item{layout}{lgCMatrix, transcript structure.}
#' @export
getIsoPartsList = function(anno,
                           gene_id = getMultiTxGene(anno),
                           depth.intron = 300,
                           cores = max(1, detectCores()-2)) {
  ## check input
  if (class(anno) != "GRanges") {
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
    intron = c(setdiff(gene, exon), ## exon-exon junction + exonic region of filtered isoform
               flank(gene, depth.intron),
               flank(gene, depth.intron, start = F)) ## add up/down-stream of genes
    strand = as.character(strand(gene))
    exon_pooled = disjoin(exon) ## these are bins that DEXSeq works on (except overlapping exons)
    intron_pooled = disjoin(intron)

    ## valid check
    if (length(tx) == 1) {
      warning(paste("Gene", g, "has only one transcript.")) ## ENSG00000270726.6
      return(list())
    }

    ## $segment - GRanges, genomic data of parts list
    segment = sort(c(exon_pooled, intron_pooled)) ## sort exonic & intronic segments

    ## $label - integer, exon numerator
    label = !!countOverlaps(segment, exon_pooled) ## logical rep.
    label = countLogicRle(label)

    ## layout - Matrix, isoform structure
    exon_tx = split(exon, exon$transcript_id)[tx$transcript_id] ## group exons by tx, then arrange by transcript_id
    seg_tx = findOverlaps(segment, exon_tx) ## map segments to tx
    layout = Matrix::Matrix(F, length(segment), length(tx), sparse = T)
    colnames(layout) = tx$transcript_id
    layout[as.matrix(seg_tx)] = T

    ## reverse everything, if "-" strand gene
    if (strand == "-") {
      segment = rev(segment)
      label = (!!rev(label)) * (max(label) + 1 - rev(label))
      layout = layout[nrow(layout):1,]
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
getEES = function(x) {
  nEES = length(x) - match(F, rev(x), nomatch = length(x) + 1) + 1
  rep(c(F, T), c(nEES, length(x) - nEES))
}

#' Reconstruct exon parts list to isoform parts list
#' @param ex exon parts list.
#' @param gm gene model.
stream = function(ex, gm) {
  for (e in 1:length(ex))
    gm[gm == e] = ex[[e]]
  as.logical(gm)
}

#' Get exon start site from exon parts list
getESS = function(x) {
  nESS = match(F, x, nomatch = length(x) + 1)
  rep(c(T, F), c(nESS - 1, length(x) - nESS + 1))
}

#' Get AS events for one isoform (decision tree)
#' @param li boolean layout of i-th isoform.
#' @param gm gene model.
getIsoAS = function(li, gm) {
  bVar = xor(li, gm) ## variable bins
  eVar = split(bVar, gm) ## variable exon
  if (!min(gm)) eVar = eVar[-1]

  isoAS = ifelse(bVar, "", ".") ## initial

  ## Alternative start site 2
  nTSS = min(which(li)) ## transcript start site
  bASS = rep(c(T, F), c(nTSS, length(gm) - nTSS))
  isoAS[bVar & bASS] = "A5U"
  ## Alternative end site 2
  nPAS = max(which(li)) ## polyA site
  bAES = rep(c(F, T), c(nPAS, length(gm) - nPAS))
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

  factor(isoAS,
         c(".", surf.events))
}

#' Enumerate AS events
numberAS = function(as) {
  ## as - factor of AS events /isoform
  runs = rle(!(as %in% c(".", "AFE", "IAP")))
  runs$values[runs$values] = 1:sum(runs$values)
  num = inverse.rle(runs)
  if (any(as == "AFE")) {
    num[!!num] = num[!!num] + 1
    num[as == "AFE"] = 1
  }
  num[as == "IAP"] = max(num) + 1
  num
}

#' Find AS events from isoform parts list
#' @param isoPL isoform parts list
#' @export
getEvent = function(isoPL, cores = max(1, detectCores()-2)) {
  registerDoParallel(cores)
  isoPLas = foreach (pl = isoPL, g = names(isoPL)) %dopar% {
    segment = pl$segment
    label = pl$label
    layout = pl$layout

    ## extract AS events, compress
    asEvent = data.frame(apply(layout, 2, getIsoAS, gm = label))
    asNum = data.frame(sapply(asEvent, numberAS))
    asName = mapply(split, asEvent, asNum, SIMPLIFY = F)
    asName = lapply(asName, function(x)
      sapply(x, unique))

    ## check isolation
    asSplit = do.call("c", mapply(split, asEvent, asNum, SIMPLIFY = F))
    asSplit = asSplit[!(names(asSplit) %in% c(".", "AFE", "IAP"))]
    asCheck = sapply(asSplit, function(x)
      length(unique(x)))
    if (any(asCheck != 1))
      message("Warning: ", g, " has adjacent AS events. Please check!")

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
#' get a named (by event_id) vector of event type for one gene
#' @param plas a plas object containing `asName`
getASName <- function(plas) {
  asName = plas$asName
  tx_id <- rep(names(asName), sapply(asName, length))
  event_name <- unlist(unname(asName), use.names = T)
  names(event_name) <- paste(tx_id, names(event_name), sep = "@")
  event_name <- na.omit(droplevels(event_name, exclude = "."))
  return(event_name)
}

#' annotate AS events (generate feature bins)
#"
#' @param isoPLas isoform parts list with AS events
#' @param cores # working cores
#' @param min.event.length numeric (positive), minimum length of a valid event
#' @param anno_ss GRanges of start/stop codon, UTR, Selenocysteine, etc; used to exclude some interior events (i.e., SE/RI/A3SS/A5SS)
#' @param remove.duplicate remove identical event duplicates (keep one randomly)
#' @return a `surfParsing` object
#' @export
annotateEvent <- function(isoPLas, cores = max(1, detectCores()-2),
                          min.event.length = 6,
                          anno_ss = NULL,
                          remove.duplicate = F) {
  registerDoParallel(cores)
  event.list = foreach (plas = isoPLas, g = names(isoPLas), .combine = "c") %dopar% {
    segment = plas$segment
    asNum = plas$asNum
    strand <- as.vector(strand(segment[1]))

    ## input check 1: existence of variable bins
    if (any(!sapply(asNum, max))) {
      warning(paste(colnames(asNum)[!sapply(asNum, max)], collapse = ", "), ": no (0) variable bin.")
    }

    ## collect event name
    event_name <- getASName(plas)
    event_name <- factor(event_name, surf.events)

    ## merge variable bins by events
    body = lapply(asNum, function(x) {
      ## remove 0 segments, merge segments by event
      ## note: reduce() will re-order segments by genomic coordinates.
      asSeg <- unlist(reduce(split(segment, replace(x, x == 0, NA))))
      if (strand == "-" && length(asSeg)) {
        i <- unlist(aggregate(seq_along(asSeg), by = list(names(asSeg)), FUN = rev)$x)
      } else i <- seq_along(asSeg)
      asSeg <- asSeg[i]
      asSeg$event_part_number = unlist(lapply(rle(names(asSeg))$lengths, seq_len))
      return(asSeg)
    })

    ## construct GRangesList of events
    event <- c_granges(body, sep = "@")
    event <- split(unname(event), names(event))
    mcols(event)$event_id = event_id = names(event)
    mcols(event)$event_name = event_name[event_id]
    mcols(event)$gene_id = rep(g, length(event_id))
    mcols(event)$transcript_id = sapply(strsplit(event_id, "@"), head, 1)

    ## ---- clean up ----
    ## (1) event body length (an amino acid spans 3 bps)
    event = event[sapply(width(event), sum) >= min.event.length]

    ## (2) when A3SS/A5SS/RI/SE overlap with AFE/ALE, keep the later
    if (!is.null(anno_ss)) {
      cnt = suppressWarnings(countOverlaps(event, anno_ss[anno_ss$gene_id == g]))
      event = event[!mcols(event)$event_name %in% c("A3SS","A5SS","RI","SE") | !cnt]
    }

    ## (3) remove duplicated event:
    ##    with the same (i) `exonic part` and (ii) `event_name`
    if (remove.duplicate) {
      hit <- findOverlaps(event, event)
      hit <- hit[from(hit) < to(hit)]
      hit <- hit[mcols(event)$event_name[from(hit)] == mcols(event)$event_name[to(hit)]] ## (ii)
      hit <- hit[as.logical(sapply(event[from(hit)] == event[to(hit)], all))] ## (i)
      if (length(hit)) event <- event[to(hit)] else event = event[integer(0)]
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
  mcols(anno_event)["transcript_id", "description"] <- "harbouring transcript identifier"
  mcols(anno_event)["genomicData", "description"] <- "genomic ranges of the ATR event"

  new("surfParsing",
      anno_event,
      params = list(min.event.length = min.event.length))
}

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## ------ 4. add features ------
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#' Find body features
#'
#' Locate the body feature of AS events and extract the genomic bins. This is a helper function to `getFeature()`
#"
#' @param plas one element from isoPLas.
#' @param type character, either "lead" for upstream, or "lag" for downstream
#' @param depth.exon depth extended into exon, default is 100 nt
#'
findBodyFeature = function(plas,
                           type = c("lead", "lag"),
                           depth.exon = 100) {
  segment = plas$segment
  asNum = plas$asNum
  strand <- as.vector(strand(segment[1]))

  ## input check 1: contain variable bins
  if (any(!sapply(asNum, max)) && type == "lead") {
    warning(paste(colnames(asNum)[!sapply(asNum, max)], collapse = ", "), ": no (0) variable bin.")
  }

  body = lapply(asNum, function(x) {
    ## remove 0 segments, merge segments by event
    ## note: reduce() will re-order segments by genomic coordinates.
    asSeg <- unlist(reduce(split(segment, replace(x, x == 0, NA))))
    if (xor(strand == "-", type == "lag")) {
      i <- cumsum(rle(names(asSeg))$lengths)
    } else {
      i <- data.table::shift(cumsum(rle(names(asSeg))$lengths), fill = 0) + 1
    }
    return(asSeg[i])
  })
  body <- as(body, "GRangesList")

  ## take flank
  bodyBin = flank(x = body,
                  width = lapply(IRanges::width(body) * -1, pmax, -depth.exon),
                  start = type == "lead")

  as(bodyBin, "CompressedGRangesList")
}

#' Find adjacent features of given feature
#'
#' Locate the adjacent bins at the desired direction and extract the genomic bins, given some location features of AS events. This is a helper function to `getFeature()`
#' @param plas one isoform element from isoPLas
#' @param body list, $id_tx, $event_id = body bins indices
#' @param type character, either"lead" for upstream, or "lag" for downstream
#' @param depth.exon extended depth into exon, default 100 nt
#' @param depth.intron extended depth into intron, default 300 nt
#'
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
    mergeSeg <- unlist(reduce(split(segment, label.new)))
    mergeSeg$exonic <- exonic

    ## probe into flank exon/intron
    probe = flank(bd, 1, start = type == "lead")
    adjSeg = to(findOverlaps(probe, mergeSeg))
    if (length(adjSeg) != length(bd)) {
      stop("Inproper \"body\" argument or overlapping isoform parts list \"plas\".")
    }

    width <- ifelse(mergeSeg[adjSeg]$exonic, depth.exon, depth.intron)
    if (!length(width)) width <- 0
    tmp <- flank(x = bd,
                 width = width,
                 start = type == "lead")
    pintersect(tmp, mergeSeg[adjSeg])
  }, body, data.frame(as.matrix(layout)), SIMPLIFY = F)

  as(adjBin, "CompressedGRangesList")
}

#' Find features on the (upstream or downstream) constitutive exons
#'
#' Locate the feature on the (upstream or downstream) constitutive exons and extract the genomic bins, given some body features of AS events. This is a helper function to `getFeature()`
#' @param plas one isoform element from isoPLas
#' @param body list, $id_tx$event_id = body bins indices
#' @param type character, either "lead" for upstream, or "lag" for downstream
#' @param depth.exon extended depth into exon
#' @param constitutive logical, whether to use consecutive or adjacent neighbor exon. The default is `TRUE`.
#'
findNeighborFeature = function(plas,
                               body,
                               type = c("lead", "lag"),
                               depth.exon = 100,
                               constitutive = T) {

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
      if (type == "lag") ord = rev(ord)
      ngb = c()
      p = NA
      for (i in ord) {
        ngb[i] = p
        if (i %in% x) p = i
      }
      setNames(ngb, seq_len(nExon))
    } else ## simply adjacent
      setNames(seq_len(nExon), data.table::shift(seq_len(nExon), 1, NA, type))
  })

  ## exon flank (names: exon index, values: flanking bins)
  flankBins <- lapply(Layout, function(x) {
    if (constitutive) {
      exon <- unlist(reduce(split(segment, replace(x, x == 0, NA)))) ## this definition is local to tx
      if (xor(strand == "-", type == "lag")) {
        i <- cumsum(rle(names(exon))$lengths)
      } else {
        i <- data.table::shift(cumsum(rle(names(exon))$lengths), fill = 0) + 1
      }
      exon <- exon[i]
    } else {
      exon <- unlist(reduce(split(segment, replace(label, label == 0, NA)))) ## this definition is global to gene
    }

    flank(exon,
          width = - pmin(IRanges::width(exon), depth.exon),
          start = type == "lag")
  })

  ## merge
  ngbBin = mapply(function(bd, ne, fb) {
    exon <- unlist(reduce(split(segment, replace(label, label == 0, NA)))) ## this definition is global to gene
    i.bd <- to(findOverlaps(bd, exon)) ## index of exon where the events fall on
    if (length(i.bd) != length(bd)) stop("Each \"body\" should cover only one exon.")
    i.ngb <- ne[as.character(i.bd)]
    bin <- fb[as.character(na.omit(i.ngb))]
    names(bin) <- names(bd)[!is.na(i.ngb)]
    bin
  }, body, ngbExon, flankBins, SIMPLIFY = F)

  as(ngbBin, "CompressedGRangesList")
}

#' Get event names
#'
#' Retrieve a named (by event_id) vector of event type for one gene. This is a helper function to `getFeature()`.
#' @param plas an element of isoPLas object, containing `asName`
getASName <- function(plas) {
  asName = plas$asName
  tx_id <- rep(names(asName), sapply(asName, length))
  event_name <- unlist(unname(asName), use.names = T)
  names(event_name) <- paste(tx_id, names(event_name), sep = "@")
  event_name <- na.omit(droplevels(event_name, exclude = "."))
  return(event_name)
}

#' Annotate AS events (generate feature bins)
#'
#' Annotated RBP features of AS event in `GRangesList` called `feature`
#' @param isoPLas isoform parts list with AS events
#' @param anno_event event annotation w/o features
#' @param depth.exon extended depth into exon, default 50 nt
#' @param depth.intron extended depth into intron, default 300 nt
#' @param cores # working cores
#' @param remove.duplicate remove identical event duplicates (keep one)
#' @param verbose logical, whether (`TRUE`) to echo progess
#' @return a `surfParsing` object containing the `feature` attribute.
#' @export
getFeature = function(isoPLas, anno_event,
                      depth.exon = 100,
                      depth.intron = 300,
                      cores = max(1, detectCores()-2),
                      remove.duplicate = T,
                      verbose = T) {
  ## input check
  if (depth.intron < 1 || depth.exon < 1)
    stop("\"depth.intron\" and \"depth.exon\" must be greater than 0.")

  feature.names <- c("up3", "up2", "up1", "bd1", "bd2", "dn1", "dn2", "dn3")
  FeatName <- list(SE = feature.names, RI = feature.names[3:6],
                   A3SS = feature.names[1:6], A5SS = feature.names[3:8],
                   AFE = feature.names[3:8], A5U = feature.names[3:6],
                   IAP = feature.names[1:6], TAP = feature.names[3:6])

  event_ids = anno_event$event_id

  registerDoParallel(cores)
  feature = foreach(g = names(isoPLas), plas = isoPLas, .combine = "c") %dopar% {
    ## collect event name
    event_name <- getASName(plas)
    event_id <- names(event_name)

    ## filter event
    if (!is.null(event_ids)) {
      event_id_unwanted <- setdiff(event_id, event_ids)
      for (i in strsplit(event_id_unwanted, "@")) {
        plas$asNum[,i[1]] = ifelse(plas$asNum[,i[1]] == i[2], 0, plas$asNum[,i[1]])
      }
      event_id <- intersect(event_id, event_ids)
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
    Bin <- list(bd1 = bd1, bd2 = bd2, up1 = up1, dn1 = dn1,
                up3 = up3, up2 = up2, dn3 = dn3, dn2 = dn2)
    Bin <- c_granges(lapply(Bin, c_granges, sep = "@"),
                     use.names = F, save.names = "feature")
    Bin <- split(Bin, names(Bin)) ## note: the order of event_id was mixed up here

    ## order features (according to event type)
    featBin <- mapply(function(en, fb) {
      featname <- FeatName[[as.character(en)]]
      fb <- fb[match(featname, fb$feature)]
      names(fb) = featname
      mcols(fb) = NULL
      fb
    }, event_name[event_id], Bin[event_id], SIMPLIFY = F)
    ## add event_name
    featBin <- GRangesList(featBin)
    mcols(featBin)$event_name = event_name[event_id]

    ## ----  remove duplicated events ----
    ## define with the same (i) `feature` and (ii) `event_name`
    if (remove.duplicate) {
      hit <- findOverlaps(featBin, featBin)
      hit <- hit[from(hit) < to(hit)]
      hit <- hit[mcols(featBin)$event_name[from(hit)] == mcols(featBin)$event_name[to(hit)]] ## (ii)
      hit <- hit[as.logical(sapply(featBin[from(hit)] == featBin[to(hit)], all))] ## (i)
      if (length(hit)) featBin <- featBin[-to(hit)]
    }

    ## this is the event annotation (with features) for one gene
    featBin
  }
  stopImplicitCluster()

  ## output
  event_id_feature <- intersect(event_ids, names(feature))
  if (verbose) {
    cat("Add features to", length(event_id_feature), "events.\n")
    print(table("AS/ATI/APA Event distribution:" = anno_event$event_name))
  }
  anno_event <- anno_event[event_id_feature,]

  ## add mcols()
  mcols(feature) = NULL
  addCols <- DataFrame(feature = feature[event_id_feature])
  mcols(addCols)$type <- "annotation"
  mcols(addCols)$description <- "genomic ranges of the event features"

  new("surfParsing", cbind(anno_event, addCols),
      params = c(anno_event@params,
                 depth.exon = depth.exon,
                 depth.intron = depth.intron,
                 remove.duplicate.event = remove.duplicate))
}
