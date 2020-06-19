## ------------------------------------------------------------
# what: surf vignette
# who: fan chen (fan.chen@wisc.edu)
# where: /store04/fanchen
# when: Aug 30, 2019
## ------------------------------------------------------------
## parse events from genome annotation (GTF)
## /path/to/genome/annotation/file
anno_file <- "~/Downloads/gencode.vM23.primary_assembly.annotation.gtf"
anno <- import(anno_file)
event <- parseEvent(anno_file)

## eCLIP-seq
sampleData.eclip = data.frame(
  bam = paste0("/path/to/eCLIP/bam/files"),
  condition = c('IP', 'IP', 'input'),
  stringsAsFactors = F
)
anno_event_featSignal <- featureSignal(
  anno_event,
  sampleData.eclip,
  signal.type = "TPM"
)

## shRNA-seq
sampleData.shrna = data.frame(
  bam = paste0("/path/to/shRNA/bam/files"),
  condition = rep(c('I', 'II'), each = 3),
  stringsAsFactors = F
)
drd <- drseqData(
  anno = anno,
  anno_event = anno_event,
  sampleData = sampleData.shrna,
  anno.prefix = paste(prefix, "drseq.annotation", sep = "."),
  nthreads = cores
)
drr <- drseq(drd, BPPARAM = MulticoreParam(workers = cores))

## association
sdata <- getSurfData()
sres <- trainSurf()
