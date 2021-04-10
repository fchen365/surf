## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(echo = TRUE, 
                      # tidy = TRUE, 
                      tidy.opts = list(comment = FALSE, width.cutoff = 60))
# options(width = 40)
# knitr::opts_knit$set(root.dir = '~/Downloads/surf_vignette')
# knitr::opts_chunk$set(fig.path = "figures/")
library(surf)

## ----prepare mm annotation, include=F, eval=F---------------------------------
#  library(rtracklayer)
#  library(usethis)
#  ## parse events from genome annotation (GTF)
#  ## /path/to/genome/annotation/file
#  anno_file <- "data/gencode.vM23.primary_assembly.annotation.gtf"
#  anno_mm <- import(anno_file)
#  gene_id <- anno_mm[seqnames(anno_mm) == "chr14" &
#                       anno_mm$gene_type == "protein_coding" &
#                       anno_mm$type == "gene"]$gene_id
#  gene_id_sampled <- sample(unique(gene_id), 50)
#  anno_mm_select <- anno_mm[anno_mm$gene_id %in% gene_id_sampled]
#  
#  ## Mus musculus
#  export(anno_mm_select, "data/gencode.vM23.primary.example.gtf")

## ----prepare hs annotation, include=F, eval=F---------------------------------
#  library(rtracklayer)
#  library(usethis)
#  ## parse events from genome annotation (GTF)
#  ## /path/to/genome/annotation/file
#  anno_file <- "data/gencode.v32.primary_assembly.annotation.gtf"
#  anno_hs <- import(anno_file)
#  gene_id <- anno_hs[seqnames(anno_hs) == "chr16" &
#                       anno_hs$gene_type == "protein_coding" &
#                       anno_hs$type == "gene"]$gene_id
#  gene_id_sampled <- sample(unique(gene_id), 24)
#  anno_hs_select <- anno_hs[anno_hs$gene_id %in% gene_id_sampled]
#  
#  ## Homo sapiens
#  export(anno_hs_select, "data/gencode.v32.primary.example.gtf")
#  
#  ## export gene ranges
#  genes <- reduce(anno_hs_select)
#  strand(genes) <- "*"
#  genes <- resize(genes, width(genes) + 300, fix = "center")
#  export(genes, "data/gencode.v32.primary.regSURFon.bed")

## ----pipeline, echo=FALSE, out.width="100%", fig.cap="Schematic overview of SURF pipeline. (a) In the analysis module 1, the first step parses alternative transcriptional regulation (ATR) events from genome annotation. (b) SURF quantifies relative event/exon usage (REU) using RNA-seq data following RBP knockdown (with wild-type control) and performs differential REU analysis (DrSeq). As a result, it infers regulation effects of RBPs as phenotypical labels for each differential ATR event.(c) In the analysis module 2, SURF  extracts location features for each ATR event and generates the feature signals using the complementary eCLIP-seq data. (d) Then, SURF models differential ATR labels by associating them with the feature signals and infers global positional preferences of RBPs. (e) In the discovery module, the inferred location features play a key role in downstream RBP-related discovery. The rank-based analysis of RBP transcript targets using external transcriptome datasets (e.g., TCGA and GTEx) discovers differential transcriptional activity through specific RBP and ATR event types."----
knitr::include_graphics("figures/0_pipeline.jpg")

## ---- eval=F------------------------------------------------------------------
#  library(surf)
#  
#  event <- parseEvent(anno_file)                                  # step 1
#  event <- drseq(event, rna_seq_sample)                           # step 2
#  event <- faseq(event, clip_seq_sample)                          # step 3
#  event <- daseq(event, getRankings(exprMat), ext_sample)         # step 4

## ---- eval=F------------------------------------------------------------------
#  anno_file <- "data/gencode.v32.primary.example.gtf"

## ---- eval=F------------------------------------------------------------------
#  rna_seq_sample <- data.frame(
#    row.names = c('sample1', 'sample2', 'sample3', 'sample4'),
#    bam = paste0("data/",c("KD1", "KD2", "WT1", "WT2"),".bam"),
#    condition = c('knockdown', 'knockdown', 'wildtype', 'wildtype'),
#    stringsAsFactors = FALSE
#  )

## ---- eval=F------------------------------------------------------------------
#  clip_seq_sample <- data.frame(
#    row.names = c('sample5', 'sample6', 'sample7'),
#    bam = paste0('data/', c("IP1", "IP2", "SMI"), '.bam'),
#    condition = c('IP', 'IP', 'SMI'),
#    stringsAsFactors = F
#  )

## ---- eval=F------------------------------------------------------------------
#  exprMat <- readRDS('data/TcgaTargetGtex_rsem_isoform_tpm_laml_blood_10each.rds')
#  ext_sample <- data.frame(
#    condition = rep(c('TCGA', 'GTEx'), each = 10),
#    row.names = colnames(exprMat)
#  )

## ----atr, echo=FALSE, out.width="100%", fig.cap="Illustration of eight ATR event types. The eight event types include exon skipping (SE), alternative 3' (A3SS) and 5' (A5SS) splicing, and intron retention (RI) within the AS class; alternative first exon (AFE) and alternative 5'UTR (A5U) within the ATI class; intronic (IAP) and tandem (TAP) alternative polyadenylation within the APA class. In each panel, the upper track depicts part of a gene model, and the lower track demarcates a specific ATR event in a transcript (isoform) with a white dashed box."----
knitr::include_graphics("figures/1_atr.jpg")

## ----parse, eval=F------------------------------------------------------------
#  event <- parseEvent("data/gencode.v32.primary.example.gtf")

## ---- echo=F------------------------------------------------------------------
# saveRDS(event, "results/intermediate1.rds")
event <- readRDS("results/intermediate1.rds")

## ----feature, echo=FALSE, out.width="100%", fig.cap="Illustration of location features for eight ATR event types. White boxes depict the ATR events with the event type labeled inside. Short and long curly brackets correspond to genomic regions of length 100bp and 300bp respectively."----
knitr::include_graphics("figures/2_feature.jpg")

## -----------------------------------------------------------------------------
event

## -----------------------------------------------------------------------------
mcols(event)

## -----------------------------------------------------------------------------
pl <- genePartsList(event)
pl

## -----------------------------------------------------------------------------
mcols(pl)

## ----drseq, eval=F------------------------------------------------------------
#  rna_seq_sample <- data.frame(
#    row.names = c('sample1', 'sample2', 'sample3', 'sample4'),
#    bam = paste0("data/",c("KD1", "KD2", "WT1", "WT2"),".bam"),
#    condition = c('knockdown', 'knockdown', 'wildtype', 'wildtype'),
#    stringsAsFactors = FALSE
#  )
#  event <- drseq(event, rna_seq_sample)

## ---- echo=F------------------------------------------------------------------
# saveRDS(event, "results/intermediate2.rds")
event <- readRDS("results/intermediate2.rds")

## -----------------------------------------------------------------------------
mcols(event)[7:12,]
event[,7:12]

## -----------------------------------------------------------------------------
drr <- drseqResults(event)
mcols(drr)
drr

## ---- eval=F------------------------------------------------------------------
#  clip_seq_sample = data.frame(
#    row.names = c('sample5', 'sample6', 'sample7'),
#    bam = paste0("data/",c("IP1", "IP2", "SMI"),".bam"),
#    condition = c('IP', 'IP', 'SMI'),
#    stringsAsFactors = FALSE
#  )
#  event <- faseq(event, clip_seq_sample,
#                 min.size = 3, fdr.cutoff = 0.3, signal.cutoff = 2)

## ---- echo=F------------------------------------------------------------------
# saveRDS(event, "results/intermediate3.rds")
event <- readRDS("results/intermediate3.rds")

## -----------------------------------------------------------------------------
event[,13:14]

## -----------------------------------------------------------------------------
mcols(event)[13:14,]

## -----------------------------------------------------------------------------
far <- faseqResults(event)
mcols(far)
far

## -----------------------------------------------------------------------------
inferredFeature(event)

## ---- eval=F------------------------------------------------------------------
#  exprMat <- readRDS('data/TcgaTargetGtex_rsem_isoform_tpm_laml_blood_10each.rds')
#  ext_sample <- data.frame(
#    condition = rep(c('TCGA', 'GTEx'), each = 10),
#    row.names = colnames(exprMat)
#  )
#  event <- daseq(event, getRankings(exprMat), ext_sample,
#                 cores = 1, target.type = "transcript")

## ---- echo=F------------------------------------------------------------------
# saveRDS(event, "results/intermediate4.rds")
event <- readRDS("results/intermediate4.rds")

## -----------------------------------------------------------------------------
dar <- daseqResults(event)
mcols(dar)
dar

## ----dispersion, echo=FALSE, out.width="70%", fig.cap="Estimated mean-dispersion functions from DrSeq analysis of RNA-seq datasets (ENCODE) with RBP targets CPSF6. Each line corresponds to an ATR event type from Figure 2."----
knitr::include_graphics("figures/3_disperson.jpg")

## ----volcano, echo=FALSE, out.width="110%", fig.cap="Volcano plot (-log10 transformed adjusted p-value versus log2 of fold change) of DrSeq results for CPSF6, stratified by ATR event types. Horizontal dashed line depicts FDR cut-off level of 0.01 and the vertical lines depict log2 fold change of 1 in absolute value."----
knitr::include_graphics("figures/4_volcano.jpg")

## ---- fig.width=7, fig.height=4, fig.cap="Functional association plot. The upper panels depicts the actual CLIP-seq binding signals on various location features, stratified by the differential event usage (DEU) upon the RBP knock-down (as the results of Step 2 -- DrSeq). The top strips indicate the ATR event type and the number of ATR events in each DEU group are reported in the parenthesis. The lower panels shows the p-values of the functional association test (FAT)."----
fa.plot(event, plot.event = c("AFE", "A5U", "IAP", "TAP"))

## -----------------------------------------------------------------------------
gr0 <- GRanges(seqnames = Rle("chr16"), 
               IRanges(89710000, width = 10000))
findOverlaps(event, subject = gr0)

## -----------------------------------------------------------------------------
subsetByOverlaps(event, ranges = gr0)

## -----------------------------------------------------------------------------
sessionInfo()

