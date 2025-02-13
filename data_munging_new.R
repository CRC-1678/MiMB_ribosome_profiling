##### R setup #####

library(parallel)

##### Global variables #####

IN_PATH <- "./data/input/"
OUT_PATH <- "./data/output/"

##### Data preparation for RiboSeq project #####
##### Load input #####
# Data was retrieved from ensembl.org/biomart on 2023/01/24
# with following query:
# Dataset: Human genes (GRCh38.p13)
# Filters: Transcript type: protein_coding; GENCODE basic annotation: Only;
# Chromosome/scaffold: 1 - 22 , MT , X , Y

pcoding_exons_GRCh38p13 <- read.table(
  paste0(IN_PATH, "240730_biomaRt_GENCODEbasic_pcoding_exons_GRCh38p13.rds"),
  header = T, sep = "\t", quote = "\"")


# Select only genes with annotated 5' UTR and 3' UTR present:

parallel::stopCluster(clust)
clust <- parallel::makeCluster(8)
parallel::clusterExport(clust, varlist = "pcoding_exons_GRCh38p13")

pcoding_hasUTR <- parSapply(clust, unique(pcoding_exons_GRCh38p13$ensembl_transcript_id),
                         FUN = function(transID) 
                           return(any(pcoding_exons_GRCh38p13[
                             pcoding_exons_GRCh38p13$ensembl_transcript_id == transID
                           ,"3_utr_start"]) &
                             any(pcoding_exons_GRCh38p13[
                               pcoding_exons_GRCh38p13$ensembl_transcript_id == transID
                               ,"5_utr_start"])))

parallel::stopCluster(clust)
rm(clust)

pcoding_hasUTR_trans <- unique(pcoding_exons_GRCh38p13$ensembl_transcript_id)[pcoding_hasUTR]

pcoding_exons_GRCh38p13_wUTR <- pcoding_exons_GRCh38p13[
  sapply(pcoding_exons_GRCh38p13$ensembl_transcript_id,
  FUN = function(transID) is.element(transID, pcoding_hasUTR_trans)),]

rm(pcoding_exons_GRCh38p13)

# select the transcripts which encode the longest protein isoform:

pcoding_longestCDS <- sapply(unique(pcoding_exons_GRCh38p13_wUTR$ensembl_gene_id),
                         FUN = function(geneID) 
                           return(max(pcoding_exons_GRCh38p13_wUTR[
                             pcoding_exons_GRCh38p13_wUTR$ensembl_gene_id == geneID
                             ,"cds_length"])))

pcoding_isLongest <- lapply(unique(pcoding_exons_GRCh38p13_wUTR$ensembl_gene_id),
                            FUN = function(geneID)
                           return(which(pcoding_exons_GRCh38p13_wUTR$ensembl_gene_id == geneID &
                                          pcoding_exons_GRCh38p13_wUTR$cds_length ==
                                          pcoding_longestCDS[geneID])))

pcoding_longestOnly <- pcoding_exons_GRCh38p13_wUTR[unlist(pcoding_isLongest),]

rm(pcoding_exons_GRCh38p13_wUTR)

identical(length(unique(pcoding_longestOnly$ensembl_gene_id)),
          length(unique(pcoding_longestOnly$ensembl_transcript_id)))

# There are genes with several transcripts remaining. Hence, now pick
# the longest transcript:

parallel::stopCluster(clust)
clust <- parallel::makeCluster(8)
parallel::clusterExport(clust, "pcoding_longestOnly")

pcoding_Longest_trans <- parSapply(clust,
  unique(pcoding_longestOnly$ensembl_gene_id),
  FUN = function(geneID) {
    geneRows <- pcoding_longestOnly$ensembl_gene_id == geneID
    geneTrans <- unique(pcoding_longestOnly[geneRows,"ensembl_transcript_id"])
    if(length(geneTrans) == 1) return(geneTrans)
    else if(length(geneTrans) > 1) {
      geneTrans_length <- sapply(geneTrans, FUN = function(TransID){
        return(as.numeric(pcoding_longestOnly[
          geneRows & pcoding_longestOnly$ensembl_transcript_id == TransID,
          "transcript_length"][1]))
      })
      return(geneTrans[which.max(geneTrans_length)])
    }
    else return(NA)
})

pcoding_longestOnly <- pcoding_longestOnly[
  is.element(el = pcoding_longestOnly$ensembl_transcript_id,pcoding_Longest_trans),]

identical(length(unique(pcoding_longestOnly$ensembl_gene_id)),
length(unique(pcoding_longestOnly$ensembl_transcript_id)))

saveRDS(pcoding_longestOnly, file = "./data/input/240730_biomaRt_GENCODEbasic_pcoding_longestOnly_GRCh38p13.rds")
pcoding_longestOnly <- readRDS(file = "./data/input/240730_biomaRt_GENCODEbasic_pcoding_longestOnly_GRCh38p13.rds")

# Export transcript list for sequence query:

saveRDS(paste(unique(pcoding_longestOnly$ensembl_gene_id),
              unique(pcoding_longestOnly$ensembl_transcript_id),
              sep = ";"),
        file = "./data/output/pcoding_Longest_transcrIDs.rds")

##### Export Length file for all transcript entries #####

parallel::stopCluster(clust)
clust <- parallel::makeCluster(8)
parallel::clusterExport(clust, "pcoding_longestOnly")

GENE_LENGTH <- t(parSapply(clust, unique(pcoding_longestOnly$ensembl_gene_id), FUN = function(geneID) {
  allLines <- grep(geneID,pcoding_longestOnly$ensembl_gene_id)
  firstLine <- allLines[1]
  lastLine <- allLines[length(allLines)]
  chromosome <- pcoding_longestOnly[firstLine,"chromosome_name"]
  strand <- pcoding_longestOnly[firstLine,"strand"]
  transID <- pcoding_longestOnly[firstLine,"ensembl_transcript_id"]
  transStart <- pcoding_longestOnly[firstLine,"transcript_start"]
  transEnd <- pcoding_longestOnly[firstLine,"transcript_end"]
  CDS_length <- pcoding_longestOnly[firstLine,"cds_length"]
  X5_UTR <- sum(apply(pcoding_longestOnly[allLines,c("5_utr_start","5_utr_end")],
                      1, FUN=function(X5UTRpos) {
                        return(X5UTRpos[2]-X5UTRpos[1]+1)
                        }), na.rm = T)
  X3_UTR <- sum(apply(pcoding_longestOnly[allLines,c("3_utr_start","3_utr_end")],
                      1, FUN=function(X5UTRpos) {
                        return(X5UTRpos[2]-X5UTRpos[1]+1)
                      }), na.rm = T)
  full_trans <- sum(c(CDS_length,X5_UTR,X3_UTR))
  return(setNames(c(paste(c(geneID,transID),collapse = ";"),chromosome,strand,transStart,transEnd,X5_UTR,CDS_length,X3_UTR,full_trans),
                  c("Gene;Transcript","chromosome","strand","transcript_start","transcript_end","5_UTR","CDS","3_UTR","FullTranscript")))
  }))

GENE_LENGTH[,"Gene;Transcript"] <- stringi::stri_replace(GENE_LENGTH[,"Gene;Transcript"],replacement = ";",fixed = ".")

write.table(GENE_LENGTH,
            file = "./data/output/pcoding_Longest_Lengths.tab",
            quote = F, col.names = T, row.names = F, sep = "\t")
# NB: RefSeq IDs will be added by biomaRt_retrieval.R

##### Re-format for use in samtools coverage #####

parallel::stopCluster(clust)
clust <- parallel::makeCluster(24)


pcoding_Longest_Exons_BED <- t(parallel::parApply(clust,pcoding_longestOnly,1,FUN = function(exonEntry) {
  exonBED <- setNames(vector(mode = "character", length = 6),
                      c("chr","EXONstart","EXONend","Gene.Transcript","dot","strand"))
  exonBED["strand"] <- c("-","+")[as.numeric(as.numeric(exonEntry[["strand"]]) == 1)+1]
  exonBED["chr"] <- exonEntry[["chromosome_name"]]
  exonBED["Gene.Transcript"] <- paste(exonEntry[c("ensembl_gene_id","ensembl_transcript_id")], collapse = ";")
  exonBED["dot"] <- "."
  exonBED["EXONstart"] <- as.integer(exonEntry[["exon_chrom_start"]])-1 # For BED format, reduce position by -1
  exonBED["EXONend"] <- as.integer(exonEntry[["exon_chrom_end"]])
  return(exonBED)
}))

write.table(pcoding_Longest_Exons_BED,
            file = "./data/output/pcoding_Longest_Exons_BED.bed",
            quote = F, col.names = F, row.names = F, sep = "\t")

# For usage with .bam files that use "chr..." as chromosome names:

write.table(cbind(chr = paste0("chr",pcoding_Longest_Exons_BED[,1]),
                  pcoding_Longest_Exons_BED[,-1]),
            file = "./data/output/pcoding_Longest_Exons_BED_chr.bed",
            quote = F, col.names = F, row.names = F, sep = "\t")


pcoding_Longest_Exons_CDSonly <- t(parallel::parApply(clust,pcoding_longestOnly,1,FUN = function(exonEntry) {
  exonBED <- setNames(vector(mode = "character", length = 6),
                      c("chr","CDSstart","CDSend","Gene.Transcript","dot","strand"))
  if(is.na(exonEntry[["cds_start"]])) return(rep(NA,6))
  exonBED["strand"] <- c("-","+")[as.numeric(as.numeric(exonEntry[["strand"]]) == 1)+1]
  exonBED["chr"] <- exonEntry[["chromosome_name"]]
  exonBED["Gene.Transcript"] <- paste(exonEntry[c("ensembl_gene_id","ensembl_transcript_id")], collapse = ";")
  exonBED["dot"] <- "."
  X5UTR_start <- as.numeric(exonEntry[["5_utr_start"]])
  X5UTR_end <- as.numeric(exonEntry[["5_utr_end"]])
  X3UTR_start <- as.numeric(exonEntry[["3_utr_start"]])
  X3UTR_end <- as.numeric(exonEntry[["3_utr_end"]])
  if(exonBED["strand"] == "+") {
    if(!is.na(X5UTR_end)) {
      print(paste(c(X5UTR_end, exonEntry[["cds_start"]]), collapse = ";"))
      if(as.numeric(exonEntry[["cds_start"]]) != 1) print(paste0("For gene ", exonEntry[["ensembl_gene_id"]],
      " 1st exon does not contain START!"))
      exonBED["CDSstart"] <- X5UTR_end + 1
    }
    else exonBED["CDSstart"] <- exonEntry[["exon_chrom_start"]]
    if(!is.na(X3UTR_start)) {
      exonBED["CDSend"] <- X3UTR_start - 1
    }
    else exonBED["CDSend"] <- exonEntry[["exon_chrom_end"]]
  }
  else if (exonBED["strand"] == "-") {
    if(!is.na(X5UTR_start)) {
      if(as.numeric(exonEntry[["cds_start"]]) != 1) print(paste0("For gene ", exonEntry[["ensembl_gene_id"]],
                                                     " 1st exon does not contain START!"))
      exonBED["CDSend"] <- X5UTR_start-1
    }
    else exonBED["CDSend"] <- exonEntry[["exon_chrom_end"]]
    if(!is.na(X3UTR_start)) {
      exonBED["CDSstart"] <- X3UTR_end + 1
    }
    else exonBED["CDSstart"] <- exonEntry[["exon_chrom_start"]]
  }
  return(exonBED)
}))

pcoding_Longest_Exons_CDSonly_BED <- as.data.frame(pcoding_Longest_Exons_CDSonly)
pcoding_Longest_Exons_CDSonly_BED <- pcoding_Longest_Exons_CDSonly_BED[!is.na(pcoding_Longest_Exons_CDSonly_BED[,1]),]
colnames(pcoding_Longest_Exons_CDSonly_BED) <- c("chr","CDSstart","CDSend","Gene.Transcript","dot","strand")
pcoding_Longest_Exons_CDSonly_BED[,"CDSstart"] <- as.integer(pcoding_Longest_Exons_CDSonly_BED[,"CDSstart"])-1 # For BED format, reduce position by -1
pcoding_Longest_Exons_CDSonly_BED[,"CDSend"] <- as.integer(pcoding_Longest_Exons_CDSonly_BED[,"CDSend"])

write.table(pcoding_Longest_Exons_CDSonly_BED,
            file = "./data/output/pcoding_Longest_Exons_CDSonly.bed",
            quote = F, col.names = F, row.names = F, sep = "\t")

# For usage with .bam files that use "chr..." as chromosome names:

write.table(cbind(chr = paste0("chr",pcoding_Longest_Exons_CDSonly_BED[,1]),
      pcoding_Longest_Exons_CDSonly_BED[,-1]),
      file = "./data/output/pcoding_Longest_Exons_CDSonly_chr.bed",
      quote = F, col.names = F, row.names = F, sep = "\t")


##### Compare stable genes depending on counting on CDS or entire transcript #####

stable_genes_EXON <- read.table("/mnt/d/data_repository/RiboSeq_SG_LS/Annotation/highest_expressed_on_complete/4_stable_genes_sorted.bed", sep = "\t", header = F)

stable_genes_CDS <- read.table("/mnt/d/data_repository/RiboSeq_SG_LS/Annotation/highest_expressed_on_CDS/4_stable_genes_sorted.bed", sep = "\t", header = F)

stable_genes_both <- read.table("/mnt/d/data_repository/RiboSeq_SG_LS/Annotation/4_stable_genes_both.bed", sep = "\t", header = F)

setdiff(unique(stable_genes_EXON$V4),unique(stable_genes_CDS$V4))
setdiff(unique(stable_genes_EXON$V4),unique(stable_genes_both$V4))
setdiff(unique(stable_genes_CDS$V4),unique(stable_genes_both$V4))

# There are 435 genes with stably expressed transcripts but not CDS in sample 0 uM H2O2 + Harr
# including cathepsin L -> interesting for mTOR regulation

stable_genes_EXON_25uM <- read.table("/mnt/d/data_repository/RiboSeq_SG_LS/Annotation/5_stable_genes.bed", sep = "\t", header = F)
stable_genes_CDS_25uM <- read.table("/mnt/d/data_repository/RiboSeq_SG_LS/Annotation/highest_expressed_on_CDS/5_stable_genes.bed", sep = "\t", header = F)

saveRDS(setdiff(unique(stable_genes_EXON_25uM$V4),unique(stable_genes_CDS_25uM$V4)), "../RiboSeq_SG_LS_pipeline1/data/other_input/NTEonlyReadProteins.rds")

# 431 for 25 uM H2O2 + Harr

##### Workspace cleanup for saving #####

# List all objects in the global environment
all_objects <- ls()

# Select objects that start with "genes_"
genes_objects <- all_objects[grep("^genes_", all_objects)]

# Define a threshold size (change as needed)
size_threshold <- 100000000

# Filter objects based on both name and size
selected_objects <- all_objects[sapply(all_objects, function(obj) object.size(get(obj)) > size_threshold)]

# Remove selected objects
rm(list = selected_objects)


