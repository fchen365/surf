## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(echo = TRUE, 
                      tidy = TRUE, 
                      tidy.opts = list(comment = FALSE))
library(surf)

## ----prepare mm annotation, include=F, eval=F---------------------------------
#  library(rtracklayer)
#  library(usethis)
#  ## parse events from genome annotation (GTF)
#  ## /path/to/genome/annotation/file
#  anno_file <- "~/Downloads/gencode.vM23.primary_assembly.annotation.gtf"
#  anno_mm <- import(anno_file)
#  gene_id <- anno_mm[seqnames(anno_mm) == "chr14" &
#                       anno_mm$gene_type == "protein_coding" &
#                       anno_mm$type == "gene"]$gene_id
#  gene_id_sampled <- sample(unique(gene_id), 50)
#  anno_mm_select <- anno_mm[anno_mm$gene_id %in% gene_id_sampled]
#  
#  ## Mus musculus
#  export(anno_mm_select, "~/Downloads/gencode.vM23.primary.example.gtf")

## ----prepare hs annotation, include=F, eval=F---------------------------------
#  library(rtracklayer)
#  library(usethis)
#  ## parse events from genome annotation (GTF)
#  ## /path/to/genome/annotation/file
#  anno_file <- "~/Downloads/gencode.v32.primary_assembly.annotation.gtf"
#  anno_hs <- import(anno_file)
#  gene_id <- anno_hs[seqnames(anno_hs) == "chr16" &
#                       anno_hs$gene_type == "protein_coding" &
#                       anno_hs$type == "gene"]$gene_id
#  gene_id_sampled <- sample(unique(gene_id), 24)
#  anno_hs_select <- anno_hs[anno_hs$gene_id %in% gene_id_sampled]
#  
#  ## Homo sapiens
#  export(anno_hs_select, "~/Downloads/gencode.v32.primary.example.gtf")
#  
#  ## export gene ranges
#  genes <- reduce(anno_hs_select)
#  strand(genes) <- "*"
#  genes <- resize(genes, width(genes) + 300, fix = "center")
#  export(genes, "~/Downloads/gencode.v32.primary.region.bed")

## ----parse, eval=F------------------------------------------------------------
#  ## parse events from genome annotation (GTF)
#  event <- parseEvent("~/Downloads/gencode.v32.primary.example.gtf")

## ---- echo=F------------------------------------------------------------------
# saveRDS(event, "~/Downloads/intermediate1.rds")
event <- readRDS("~/Downloads/intermediate1.rds")

## -----------------------------------------------------------------------------
event

## -----------------------------------------------------------------------------
mcols(event)

## -----------------------------------------------------------------------------
pl <- genePartsList(event)
pl

## -----------------------------------------------------------------------------
mcols(pl)

## -----------------------------------------------------------------------------
gr0 <- GRanges(seqnames = Rle("chr16"), 
               IRanges(89710000, width = 10000))
findOverlaps(event, subject = gr0)

## -----------------------------------------------------------------------------
subsetByOverlaps(event, ranges = gr0)

## ----drseq, eval=F------------------------------------------------------------
#  rna_seq_sample <- data.frame(
#    row.names = c('sample1', 'sample2', 'sample3', 'sample4'),
#    bam = paste0("~/Downloads/",c("KD1", "KD2", "WT1", "WT2"),".bam"),
#    condition = c('knockdown', 'knockdown', 'wildtype', 'wildtype'),
#    stringsAsFactors = FALSE
#  )
#  event <- drseq(event, rna_seq_sample)

## ---- echo=F------------------------------------------------------------------
# saveRDS(event, "~/Downloads/intermediate2.rds")
event <- readRDS("~/Downloads/intermediate2.rds")

## -----------------------------------------------------------------------------
event[,7:12]

## -----------------------------------------------------------------------------
mcols(event)[7:12,]

## -----------------------------------------------------------------------------
drr <- drseqResults(event)
drr

## -----------------------------------------------------------------------------
mcols(drr)

## ---- eval=F------------------------------------------------------------------
#  clip_seq_sample = data.frame(
#    row.names = c('sample5', 'sample6', 'sample7'),
#    bam = paste0("~/Downloads/",c("IP1", "IP2", "SMI"),".bam"),
#    condition = c('IP', 'IP', 'SMI'),
#    stringsAsFactors = FALSE
#  )
#  event <- faseq(event, clip_seq_sample,
#                 min.size = 3, fdr.cutoff = 0.3, signal.cutoff = 2)

## ---- echo=F------------------------------------------------------------------
# saveRDS(event, "~/Downloads/intermediate3.rds")
event <- readRDS("~/Downloads/intermediate3.rds")

## -----------------------------------------------------------------------------
event[,13:14]

## -----------------------------------------------------------------------------
mcols(event)[13:14,]

## -----------------------------------------------------------------------------
far <- faseqResults(event)
far

## -----------------------------------------------------------------------------
mcols(far)

## ---- fig.width=8-------------------------------------------------------------
fa.plot(event)

## -----------------------------------------------------------------------------
inferredFeature(event)

## ---- eval=F------------------------------------------------------------------
#  ## rank transcripts (TPM)
#  exprMat <- readRDS('~/Downloads/TcgaTargetGtex_rsem_isoform_tpm_laml_blood_10each.rds')
#  ## sample data
#  ext_sample <- data.frame(
#    condition = rep(c('TCGA', 'GTEx'), each = 10),
#    row.names = colnames(exprMat)
#  )
#  
#  ## differential activity (transcript)
#  event <- daseq(event, getRankings(exprMat), cores = 1, ext_sample)

## ---- echo=F------------------------------------------------------------------
# saveRDS(event, "~/Downloads/intermediate4.rds")
event <- readRDS("~/Downloads/intermediate4.rds")

## -----------------------------------------------------------------------------
dar <- daseqResults(event)
dar

## -----------------------------------------------------------------------------
mcols(dar)

