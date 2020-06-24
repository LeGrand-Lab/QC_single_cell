###
# QUALITY CONTROL AND PREPROCESSING
# This workflow is valid for 10X-format raw counts.
# Third-party protocol yields matrices where GENE SYMBOLS are rownames.
# ATTENTION: Doublet detection procedures should only be applied to libraries 
# generated in the  same experimental batch.
#  *-* many thanks to Dr. L Modolo for most of this code *-*
# --
# input : data/MYEXPERMNT : barcodes.tsv.gz features.tsv.gz matrix.mtx.gz
# output : results/*.pdf and  'rdatas/MYEXPERMT_END.RData' for downstream analysis
# Joha GL 2020
##


#  ============ USER DEFINED
prloc = "~/QC_single_cell" #<<<< check working directory!!
exper="dorsowt2" # TODO change in coherence to folder input
#  ============ end user defined
exper = gsub("/",'',exper) 
listpackages <- c( "ggplot2",   "dplyr", "stringr",   "tidyverse",
                   "BSgenome",  "GenomeInfoDb", "Seurat",
                   "lubridate", # right color
                   "simpleSingleCell", # for scRNA storage & manipulation
                   "scater", # for QC control
                   "scran", # analysis pipeline
                   "uwot", # UMAP dim-red
                   "DropletUtils", #utility functions for handling single-cell (RNA-seq)
                   "AnnotationHub", # ensbl query
                   "AnnotationDbi", # ensbl query
                   "sctransform", "SingleCellExperiment", "Matrix" )

lapply(listpackages,require,character.only=TRUE)

# ================== SETTING PATHS
setwd(prloc)
resdir="results/" 
system(paste("mkdir",resdir)) # creates if not exists
system("mkdir rdatas") #creates if not exists

sink(paste0(resdir,"outputsfile.txt"), append=TRUE) 
sink(paste0(resdir,"outputsfile.txt"), append=TRUE, type="message")

# read 10X 
# ================================================================================
sce <- tryCatch({
    sce <- read10xCounts(paste0("data/",exper), type="sparse")
    return(sce)
} , error = function(e){
  print("failed DropletUtils::read10xCounts, using Seurat+SCE steps")
  matdat <- Seurat::Read10X(paste0("data/",exper))
  sce <-  SingleCellExperiment(assays=list(counts=matdat))
  rm(matdat)
  return(sce)
} , warning=function(w) {
  print("10x to sce done but check warnings")
} 
)
print("initial matrix dimensions")
dim(sce)#27998  2432
print("starting analysis")
print("loading data and adding annotations")


 
head(rowData(sce))
#DataFrame with 6 rows and 0 column, fix:
rowData(sce) <- DataFrame(
  genes_names = rownames(sce)
)
## **  DATA ANNOTATIONS **

hub_infos <- AnnotationHub() # this goes into .cache/AnnotationHub

hub_ids <- mcols(hub_infos) %>%
  data.frame() %>% 
  rownames_to_column(var = "id") %>% # we keep the rownames
  as_tibble() %>% 
  dplyr::filter(
    dataprovider %in% "Ensembl" & # Ensembl annotation
      species %in% c("Homo sapiens", "Mus musculus"), # for the two species we want
    genome %in% c("GRCh38", "GRCm38"), # on the right genome
    str_detect(title, "99"), # on the right revision
    rdataclass %in% "EnsDb",
  ) %>% 
  dplyr::select(id, species)  # id is species code in .db

# pull_ensemble id (dataset is in gene symbols instead,as we know).
# exemple: "Gsn" (symbol) --> "ENSG00000183765" (ensembl geneid) 
pull_ensembl.id <- function(id, hub_infos){
  mapIds(
    hub_infos[[id]],
    keys = str_replace(rowData(sce)$genes_names, "(.*)\\..*", "\\1"),
    keytype = "SYMBOL",   ### transform into a variable to put in function input
    column = "GENEID")
}

merge_ensembl.id <- function(id, hub_infos){
  sapply(id %>% pull(id), pull_ensembl.id, hub_infos) %>% 
    as_tibble() %>% 
    unite(col = "ensembl.id", na.rm = TRUE) %>% 
    pull(ensembl.id)
}

pull_loc <- function(id, hub_infos){
  mapIds(
    hub_infos[[id]],
    keys = str_replace(rowData(sce)$ensembl.id, "(.*)\\..*", "\\1"),
    keytype = "GENEID",   ### transform into a variable to put in function input
    column = "SEQNAME")
}

merge_loc <- function(id, hub_infos){
  sapply(id %>% pull(id), pull_loc, hub_infos) %>% 
    as_tibble() %>% 
    unite(col = "chr_pos", na.rm = TRUE) %>% 
    pull(chr_pos)
}

rowData(sce)$ensembl.id <- merge_ensembl.id(hub_ids,hub_infos)
rowData(sce)$chr_pos = merge_loc(hub_ids, hub_infos)
rowData(sce)$is_genomic <- rowData(sce)$chr_pos %in% c(as.character(1:22), "X", "Y")
rowData(sce)$species = ifelse(str_detect(rowData(sce)$ensembl.id, "^ENSMUSG"),
                   "Mus musculus", ifelse(str_detect(rowData(sce)$ensembl.id, "^ENSG"),
                   "Homo sapiens", "exoticGeneSymbol"))

rowData(sce)
print("species detected by gene symbol, before correction")
table(rowData(sce)$species)
# exoticGeneSymbol     Homo sapiens     Mus musculus 
#      1829               20            26149
# and after verification, exoticGeneSymbols belong to M musculus:
rowData(sce)$species[rowData(sce)$species=="exoticGeneSymbol"] <- "Mus musculus"

tail(rowData(sce)[rowData(sce)$species %in% "Homo sapiens",])
# WDR97       WDR97                    ENSG00000179698           8       TRUE Homo sapiens
# C2             C2 ENSG00000166278_ENSMUSG00000024371                  FALSE Homo sapiens
# C3             C3 ENSG00000125730_ENSMUSG00000024164                  FALSE Homo sapiens
# PISD         PISD                    ENSG00000241878          22       TRUE Homo sapiens
# DHRSX       DHRSX                    ENSG00000169084           X       TRUE Homo sapiens


# ==============
# QC first steps
# ==============
print("detecting contamination and inferior outliers")
colData(sce) <- DataFrame(
  n_mm_umi = colSums(counts(sce)[rowData(sce)$species %in% "Mus musculus", ]),
  n_hg_umi = colSums(counts(sce)[rowData(sce)$species %in% "Homo sapiens", ]),
  n_umi = colSums(counts(sce))
)
colData(sce)$species <- colData(sce) %>% 
  as_tibble() %>% 
  mutate(
    prop_mm = n_mm_umi / n_umi,
    prop_hg = n_hg_umi / n_umi,
    species = case_when(
      prop_mm > 0.9 ~ "Mus musculus",
      prop_hg > 0.9 ~ "Homo sapiens",
      TRUE ~ "mixed"
    )
  ) %>% pull(species)
colData(sce)

colData(sce) %>% as_tibble() %>% summary()

print("filtering out not expressed features")
rowData(sce)$expressed <- scater::nexprs(sce,byrow=TRUE)>0
per_cell <- perCellQCMetrics(sce[rowData(sce)$expressed, ], subset = list(non_genomic = !rowData(sce[rowData(sce)$expressed, ])$is_genomic))
summary(per_cell$sum) # UMI counts !
#  Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
#    567    5444    9548    9702   12948   41029
summary(per_cell$detected)
summary(per_cell$subsets_non_genomic_percent)
colData(sce) <- cbind(colData(sce),per_cell)

plotdet <- scater::plotColData(sce,x="sum",y="detected", colour_by="species")
pdf(paste0(resdir,exper,"_plotColData_sumVSdetected.pdf"))
plotdet
dev.off()
  
# inferior outliers 
colData(sce)$keep_total <- scater::isOutlier(colData(sce)$sum,type = "lower", log=TRUE)
table(colData(sce)$keep_total) # TRUE are OUTLIERS
sce <- scater::addPerFeatureQC(sce)
head(rowData(sce))
summary(rowData(sce)$detected)

outliplot <- scater::plotColData(sce, x="sum",y="detected",colour_by="keep_total")
outliplot <-outliplot + scale_fill_discrete(name="is.Outlier") + ggtitle("Cells under lower addPerFeatureQC metrics ('inferior' outliers)")
pdf(paste0(resdir,exper,"plotColData_outlier.pdf"))
outliplot
dev.off()

# save(sce, file = paste0("rdatas/",exper,".RData")) # if problems, save here and debug

# ==============
# QC continued
# ==============
print("post perfeatureQC matrix dimensions")
dim(sce)
print("colnames(colData(sce))")
colnames(colData(sce))
dim(sce) 
print("substracting empty drops, kneeplot")
knee_plot <- function(bc_rank) {
  knee_plt <- tibble(rank = bc_rank[["rank"]],
                     total = bc_rank[["total"]]) %>%
    distinct() %>%
    dplyr::filter(total > 0)
  annot <- tibble(inflection = metadata(bc_rank)[["inflection"]],
                  rank_cutoff = max(bc_rank$rank[bc_rank$total > metadata(bc_rank)[["inflection"]]]))
  p <- ggplot(knee_plt, aes(total, rank)) +
    geom_line() +
    geom_hline(aes(yintercept = rank_cutoff), data = annot, linetype = 2) +
    geom_vline(aes(xintercept = inflection), data = annot, linetype = 2) +
    scale_x_log10() +
    scale_y_log10() +
    annotation_logticks() +
    labs(y = "Rank", x = "Total UMIs")
  return(p)
}

bcrank <- DropletUtils::barcodeRanks(
  SingleCellExperiment::counts(
    sce[rowData(sce)$expressed, ]))

pdf(paste0(resdir,exper,"_knee_plot.pdf"))
knee_plot(bcrank)
dev.off()

colData(sce)$is_cell <- colData(sce)$n_umi > metadata(bcrank)$inflection
summary(colData(sce)$is_cell)  # 14 barcodes are not real cells

print("post eval empty drops, matrix dimensions")
dim(sce)

print("Finding doublets")
sce <- computeSumFactors(
  sce[rowData(sce)$expressed, sce$is_cell],
  clusters = sce$species[sce$is_cell]
)
sce <- logNormCounts(sce)

dbl_dens <- doubletCells(sce[rowData(sce)$expressed, sce$is_cell])
sce$doublet_score <- 0
sce$doublet_score[sce$is_cell] <- log10(dbl_dens + 1)

save(sce, dbl_dens, file = paste0("rdatas/",exper,"_END.RData"))

tsnepl <- plotTSNE(sce[rowData(sce)$expressed,sce$is_cell], colour_by="doublet_score")
detfeat <- scater::plotColData(sce, x="sum",y="detected",colour_by="doublet_score")
pdf(paste0(resdir,exper,"_doublets.pdf"),width=13)
tsnepl + detfeat
dev.off()

pdf(paste0(resdir,exper,"_histogram.pdf"))
qplot(sce$doublet_score, geom="histogram")
hist(colData(sce)$doublet_score)
dev.off()
print("END")

print("FINAL (post doublets detection) matrix dimensions")
dim(sce)
sink()
sink(type="message")
# END
# ================================================================================
## NOTES:
# note that exemple query:
#   mapIds(org.Hs.eg.db, keys=MYVECTOR, column="SYMBOL", keytype="ENTREZID")
# is the same as : mapIds(hub_infos[["AH78783"]], keys= ...)
#   because "AH78783" is the id accession for H sapiens database

