##### R setup #####
library(biomaRt)
library(stringi)
library(Biostrings)

##### biomaRt setup #####
# # Having problems with SSL certificate verification. solved by doing:
# options(download.file.method="curl", download.file.extra="-k -L")
# backup_curl <- curl::curl_options("verify")

BMhost <- "https://oct2022.archive.ensembl.org"

ensembl_archive <- useMart("ensembl", host = BMhost, verbose = TRUE)
ensembl_current <- useMart("ensembl", verbose = TRUE)

# Choose ensembl mart depending on application:
ensembl <- ensembl_archive
# ensembl <- ensembl_current

ensembl_Hs = useDataset("hsapiens_gene_ensembl", ensembl)
ensembl_Hs_filters <- biomaRt::listFilters(ensembl_Hs,what = c("name","description","fullDescription"))

ensembl_Hs_filters[grep("protein", ensembl_Hs_filters$description),]
ensembl_Hs_filters[grep("GENCODE", ensembl_Hs_filters$description),]
ensembl_Hs_filters[grep("Canonical", ensembl_Hs_filters$description),]

filterSet <- list(with_protein_id = TRUE,
                  transcript_gencode_basic = TRUE, 
#                  transcript_is_canonical = TRUE, 
                  chromosome_name = c(1:22,"MT","X","Y"))

ensembl_pages <- attributePages(ensembl_Hs)

listAttributes(ensembl_Hs, page = "feature_page")#[grep("ID",listAttributes(ensembl_Hs, page = "feature_page")),]

listAttributes(ensembl_Hs, page = "structure")#[grep("transcript",listAttributes(ensembl_Hs, page = "structure")),]

listAttributes(ensembl_Hs, page = "sequences")

attributeSet <- c("ensembl_gene_id","ensembl_transcript_id","ensembl_exon_id",
                  "chromosome_name","strand",
                  "transcript_start","transcript_end","transcript_length",
                  "exon_chrom_start","exon_chrom_end",
                  "5_utr_start","5_utr_end",
                  "cds_start","cds_end",
                  "3_utr_start","3_utr_end",
                  "cds_length")

##### Helper functions #####

orderSeqByStrand <- function(SeqDf) {
  if(SeqDf[1,]$strand == "1"){
    SeqDf <- SeqDf[order(SeqDf$exon_chrom_start,SeqDf$exon_chrom_end),] # Order the gene sequence
  } 
  else if(SeqDf[1,]$strand == "-1") {
    SeqDf <- SeqDf[order(-SeqDf$exon_chrom_start,-SeqDf$exon_chrom_end),]
  }
  return(SeqDf)
}
##### Queries #####

getCanonicalTrIDs <- getBM(mart = ensembl_Hs, 
      attributes = c("ensembl_gene_id","ensembl_transcript_id"), 
      filters = c(setNames(TRUE,"transcript_is_canonical"),filterSet))

pcoding_exons_GRCh38p13 <- getBM(mart = ensembl_Hs,       
                                 attributes = attributeSet, 
                                 filters = filterSet)

saveRDS(pcoding_exons_GRCh38p13, "./data/input/240730_biomaRt_GENCODEbasic_pcoding_exons_GRCh38p13.rds")

# Process pcoding_exons_GRCh38.p13 using script "data_munging_new.R" to
# extract transcripts encoding the longest proteins with GENCODE basic annotation.
##### Import selected transcript set  #####

selectedGeneTrans <- readRDS(file = "./data/output/pcoding_Longest_transcrIDs.rds")
selectedGIDTrID <- sapply(strsplit(selectedGeneTrans, split = ";"),FUN = function(x) return(x))
GENE_LENGTH <- read.table("./data/output/pcoding_Longest_Lengths.tab",header = T, sep = "\t")


##### Query corresponding refSeq IDs #####
# Retrieve refseq mRNA IDs and transcript length to identify correct match for multi-matches:

selectedTrRefSeqID <- getBM(mart = ensembl_Hs,       
                           attributes = c("ensembl_gene_id","ensembl_transcript_id","refseq_mrna", "transcript_length"), 
                           filters = list(ensembl_transcript_id = selectedGIDTrID[2,]))

sum(duplicated(selectedTrRefSeqID[,"ensembl_transcript_id"]))
# Many multi-matches to refSeq IDs (matching between ensembl and NCBI database is similarity-based)
sum(selectedTrRefSeqID[,"refseq_mrna"] == "")
# 1304 transcripts without matching refseqID on biomaRt catalog
saveRDS(selectedTrRefSeqID[selectedTrRefSeqID[,"refseq_mrna"] == "","ensembl_transcript_id"],
        file = "./data/output/selectedTrs_no_refSeqID_found.rds")
# Remove:
selectedTrRefSeqID <- selectedTrRefSeqID[!(selectedTrRefSeqID[,"refseq_mrna"] == ""),]

# Identify correct match using refSeq catalog:
refSeq_catalog_human <- read.table("/mnt/d/NGS_files/genomes/NCBI_refseq_225/filtered_refseq_lengths.tsv", sep = "\t")
colnames(refSeq_catalog_human) <- c("taxonomy_ID","species_name","accession.version","refseq_release","refseq_status","length")
refSeq_catalog_human[,"accession_stripped"] <- sapply(refSeq_catalog_human$accession.version, FUN = function(x) stri_split(x, fixed = ".")[[1]][1])
anyDuplicated(refSeq_catalog_human[,"accession_stripped"])
# Stripped accessions are unique.

# Look for transcripts whose matched refSeq ID is not part of the catalogue:

saveRDS(selectedTrRefSeqID[!(selectedTrRefSeqID[,"refseq_mrna"] %in% refSeq_catalog_human$accession_stripped),"ensembl_transcript_id"],
        file = "./data/output/selectedTrs_refSeqID_noLength.rds")

# Remove:
selectedTrRefSeqID <- selectedTrRefSeqID[selectedTrRefSeqID[,"refseq_mrna"] %in% refSeq_catalog_human$accession_stripped,]

# Process duplicate results for ensembl transcript IDs:
selectedTrRefSeqID <- cbind(selectedTrRefSeqID,
                            refSeq_ID_lengthmatched =rep(NA,nrow(selectedTrRefSeqID)),
                            refSeq_ID_length =rep(NA,nrow(selectedTrRefSeqID)),
                            refSeq_match_unique = rep(TRUE,nrow(selectedTrRefSeqID)))

for(ensembl_id in selectedTrRefSeqID[!(selectedTrRefSeqID[,"refseq_mrna"] == ""),"ensembl_transcript_id"]) {
  #print(paste("Looking up", ensembl_id))
  our_length <- GENE_LENGTH[grep(ensembl_id, GENE_LENGTH$Gene.Transcript),"FullTranscript"]
  matched_refseqs <- selectedTrRefSeqID[selectedTrRefSeqID[,"ensembl_transcript_id"] == ensembl_id,"refseq_mrna"]
  if(!any(matched_refseqs %in% refSeq_catalog_human$accession_stripped)) {
    selectedTrRefSeqID[selectedTrRefSeqID[,"ensembl_transcript_id"] == ensembl_id,"refSeq_ID_lengthmatched"] <- "none"
    next
  }
  # Only transcripts with matched refSeq IDs are processed... the others keep NA for all columns.
  length_matched  <- which(refSeq_catalog_human[refSeq_catalog_human$accession_stripped %in% matched_refseqs,"length"] == our_length) 
  if(!any(length_matched)) {
    selectedTrRefSeqID[selectedTrRefSeqID[,"ensembl_transcript_id"] == ensembl_id,"refSeq_ID_lengthmatched"] <- "none"
    selectedTrRefSeqID[selectedTrRefSeqID[,"ensembl_transcript_id"] == ensembl_id,"refSeq_ID_length"] <-
      refSeq_catalog_human[refSeq_catalog_human$accession_stripped == matched_refseqs[1],"length"]
    if(length(matched_refseqs) > 1) 
      selectedTrRefSeqID[selectedTrRefSeqID[,"ensembl_transcript_id"] == ensembl_id,"refSeq_match_unique"] <- FALSE
    next
  }
  if(length(length_matched) > 1) {
      # print(paste0("More than one matching refseq with correct length for ", ensembl_id))
      length_matched <- length_matched[1]
      selectedTrRefSeqID[selectedTrRefSeqID[,"ensembl_transcript_id"] == ensembl_id,"refSeq_match_unique"] <- FALSE
  }
  refseq_matched <- refSeq_catalog_human[refSeq_catalog_human$accession_stripped %in% matched_refseqs,"accession_stripped"][length_matched]
  selectedTrRefSeqID[selectedTrRefSeqID[,"ensembl_transcript_id"] == ensembl_id,"refSeq_ID_lengthmatched"] <- refseq_matched
  selectedTrRefSeqID[selectedTrRefSeqID[,"ensembl_transcript_id"] == ensembl_id,"refSeq_ID_length"] <- 
    refSeq_catalog_human[refSeq_catalog_human$accession_stripped == refseq_matched,"length"]
}


length(unique(selectedTrRefSeqID[selectedTrRefSeqID[,"refSeq_ID_lengthmatched"] == "none","ensembl_transcript_id"]))
# 4041 lines/ 2693 cases with refSeq ID match but without length match:

saveRDS(unique(selectedTrRefSeqID[selectedTrRefSeqID[,"refSeq_ID_lengthmatched"] == "none","ensembl_transcript_id"]),
        file = "./data/output/selectedTrs_refSeqID_wrongLength.rds")

# Remove:
selectedTrRefSeqID <- selectedTrRefSeqID[!selectedTrRefSeqID[,"refSeq_ID_lengthmatched"] == "none",]

# Remove 6128 wrong multi-matches:
sum(!(selectedTrRefSeqID[,"refseq_mrna"] == selectedTrRefSeqID[,"refSeq_ID_lengthmatched"]) &
        selectedTrRefSeqID[,"refSeq_match_unique"])

# Also remove 134 multi-matches for which several refSeq IDs have the correct length:
sum(!selectedTrRefSeqID[,"refSeq_match_unique"])

selectedTrRefSeqID <- selectedTrRefSeqID[selectedTrRefSeqID[,"refseq_mrna"] == selectedTrRefSeqID[,"refSeq_ID_lengthmatched"],]

# 13973 matches (of 17974 transcripts in total) remaining.
# Attach the refSeq IDs that showed a length match to the GENE_LENGTH information:

GENE_LENGTH[,"refSeq_match"] <- selectedTrRefSeqID[match(GENE_LENGTH$Gene.Transcript, apply(selectedTrRefSeqID[,c(1,2)],1,FUN=function(x) paste(x, collapse = ";"))),"refSeq_ID_lengthmatched"]

write.table(GENE_LENGTH,
            file = "./data/output/pcoding_Longest_Lengths.tab",
            quote = F, col.names = T, row.names = F, sep = "\t")

##### Retrieve exon information for selected transcripts if 5' UTR spans several exons #####

pcoding_exons_GRCh38p13_reduced <- pcoding_exons_GRCh38p13[pcoding_exons_GRCh38p13[,"ensembl_transcript_id"] %in% selectedGIDTrID[2,],]
# Transcripts in which the CDS does not start on Exon 1-x have NA entries in column "cds_start" for exon 1-x.
# For these cases, retrieve genomic coordinates for all exons with cds_start == NA and for the first exon that
# has the CDS (cds_start == 1).

pcoding_multiexon_5UTR <- pcoding_exons_GRCh38p13_reduced[is.na(pcoding_exons_GRCh38p13_reduced$cds_start) |
                                                            pcoding_exons_GRCh38p13_reduced$cds_start == 1,
                                                          c("ensembl_transcript_id","exon_chrom_start","exon_chrom_end")]
# Remove single-line entries (each transcript included due to cds_start == 1):

pcoding_multiexon_5UTR <- pcoding_multiexon_5UTR[!pcoding_multiexon_5UTR$ensembl_transcript_id %in% 
                                                   names(table(pcoding_multiexon_5UTR$ensembl_transcript_id)[table(pcoding_multiexon_5UTR$ensembl_transcript_id) == 1]),]

# Collapse the exon spans for each transcript ID with multiexon-5UTR:

pcoding_multiexon_5UTR_collapse <- setNames(sapply(unique(pcoding_multiexon_5UTR$ensembl_transcript_id), FUN = function(transID){
  paste(apply(pcoding_multiexon_5UTR[pcoding_multiexon_5UTR$ensembl_transcript_id == transID,c("exon_chrom_start","exon_chrom_end")],1, 
        FUN = function(transID_exon) paste(transID_exon, collapse = "-")),collapse = "|")
}), unique(pcoding_multiexon_5UTR$ensembl_transcript_id))

# Attach this info to GENE_LENGTH file:

GENE_LENGTH[,"multiexon_5UTR_spans"] <- rep(NA,nrow(GENE_LENGTH))

for(transID in names(pcoding_multiexon_5UTR_collapse)) {
  GENE_LENGTH[grep(transID,GENE_LENGTH$Gene.Transcript),"multiexon_5UTR_spans"] <- pcoding_multiexon_5UTR_collapse[transID]
}

write.table(GENE_LENGTH,
            file = "./data/output/pcoding_Longest_Lengths.tab",
            quote = F, col.names = T, row.names = F, sep = "\t")


##### Compare selected transcripts to canonical transcript IDs #####

match(selectedGIDTrID[1,],getCanonicalTrIDs[,1])

CanTrIDmatched <- getCanonicalTrIDs[match(selectedGIDTrID[1,],getCanonicalTrIDs[,1]),2]

selectedGIDTrID <- rbind(selectedGIDTrID, selectedGIDTrID[2,] == CanTrIDmatched)
selectedGIDTrID[3,is.na(selectedGIDTrID[3,])] <- "missing"

pie(table(selectedGIDTrID[3,]))

# For the selected transcripts, get the CDS start on canonical and non-canonical:

selectedTrIDstart <- getBM(mart = ensembl_Hs,       
                          attributes = c("ensembl_gene_id","ensembl_transcript_id","genomic_coding_start","strand"), 
                          filters = list(ensembl_transcript_id = selectedGIDTrID[2,]))
  
selectedTrIDstart_uniq <- sapply(unique(selectedTrIDstart[,"ensembl_gene_id"]), FUN = function(geneID){
  return(min(selectedTrIDstart[selectedTrIDstart$ensembl_gene_id == geneID,"genomic_coding_start"], na.rm = T))
})

canonicalTrIDstart <- getBM(mart = ensembl_Hs, 
                            attributes = c("ensembl_gene_id","ensembl_transcript_id","genomic_coding_start","strand"), 
                           filters = c(setNames(TRUE,"transcript_is_canonical"),filterSet))

canonicalTrIDstart_uniq <- sapply(unique(canonicalTrIDstart[,"ensembl_gene_id"]), FUN = function(geneID){
  return(min(canonicalTrIDstart[canonicalTrIDstart$ensembl_gene_id == geneID,"genomic_coding_start"], na.rm = T))
})

# Compare Start sites:

sum(is.na(selectedTrIDstart_uniq))
sum(is.na(canonicalTrIDstart_uniq))

selectedVsCanStart <- selectedTrIDstart_uniq - canonicalTrIDstart_uniq[match(names(selectedTrIDstart_uniq),names(canonicalTrIDstart_uniq))]

sum(is.na(selectedVsCanStart))

pie(table(selectedVsCanStart == 0))

# Export list of genes for which the chosen transcript has a different CDS start annotation:

write.csv(names(selectedTrIDstart_uniq)[c(which(is.na(selectedVsCanStart)),which(!selectedVsCanStart == 0))],
          file = "./data/output/selectedTrs_With_diffSTART.csv")


##### Retrieveexon sequences for selected transcripts #####
# Use downloaded FASTA file (ensembl release 107 2022-10-24):

allTransF_path <- "./data/input/Homo_sapiens.GRCh38.cdna.all.fa"

allTrans_seq <- readDNAStringSet(allTransF_path)

allTrans_simpleNames <- sapply(names(allTrans_seq), FUN = function(fheader) {
  strsplit(fheader, split = "\\.")[[1]][1]
})

any(duplicated(allTrans_simpleNames))

selectedTrans_seq <- allTrans_seq[allTrans_simpleNames %in% 
                                    sapply(selectedGeneTrans, FUN = function(GeneTrans) {
                                      strsplit(GeneTrans, split = ";")[[1]][2]
                                    })]


GeneTransFastaHead <- sapply(names(selectedTrans_seq), FUN = function(fheader) {
  trans.Ver <- strsplit(fheader, split = " ")[[1]][1]
  gene.Ver <- strsplit(fheader, split = " ")[[1]][4]
  transNoVer <- strsplit(trans.Ver, split = "\\.")[[1]][1]
  geneNoVer <- strsplit(substring(gene.Ver,6), split = "\\.")[[1]][1]
  return(paste0(geneNoVer,";",transNoVer))
})

names(selectedTrans_seq) <- GeneTransFastaHead

setdiff(selectedGeneTrans,GeneTransFastaHead)
setdiff(GeneTransFastaHead,selectedGeneTrans)

# A few transcripts were assigned to different genes (4) and some (8) were not present.

# Remove the entries corresponding to 4 wrongly assigned transcripts from the sequence set:

which(!(names(selectedTrans_seq) %in% selectedGeneTrans))
selectedTrans_seq <- selectedTrans_seq[-which(!(names(selectedTrans_seq) %in% selectedGeneTrans))]

# For all failed cases, retrieve the transcript sequence from biomaRt:

selectFailed <- sapply(setdiff(selectedGeneTrans,names(selectedTrans_seq)), FUN = function(GeneTrans) {
     TransFASTA <- paste(orderSeqByStrand(getBM(mart = ensembl_Hs, 
                               attributes = c("gene_exon","strand","exon_chrom_start","exon_chrom_end"), 
                               filters = list(ensembl_gene_id = strsplit(GeneTrans, split = ";")[[1]][1],
                                              ensembl_transcript_id = strsplit(GeneTrans, split = ";")[[1]][2])))[,1],
                         collapse = "")
     return(TransFASTA)
})

# Append and sort by the selectedGeneTrans vector:

selectedTrans_seq_sort <- c(as.character(selectedTrans_seq),selectFailed)[selectedGeneTrans]

identical(names(selectedTrans_seq_sort),selectedGeneTrans)

write.table(file = "./data/output/pcoding_Longest_Exons.fasta",
            paste0(paste0(names(selectedTrans_seq_sort), "\t"),selectedTrans_seq_sort),
            quote = F, row.names = F, col.names = F, sep = "\t")


# Querying biomaRt as below is very slow.
# group_size <- ceiling(length(selectedGeneTrans) / 24)
# group_indices <- rep(1:24, each = group_size, length.out = length(selectedGeneTrans))
# 
# # Split the vector into 24 equally sized groups
# selectedGeneTrans_split <- split(selectedGeneTrans,group_indices)
# 
# parallel::stopCluster(clust)
# clust <- parallel::makeCluster(24)
# parallel::clusterExport(clust,c("selectedGeneTrans_split","orderSeqByStrand","ensembl_Hs"))
# 
# parSapply(clust, 1:24, FUN = function(TransGroup) {
#   require(biomaRt)
#   ensembl = useMart("ensembl", host = BMhost, verbose = TRUE)
#   ensembl_Hs = useDataset("hsapiens_gene_ensembl", ensembl)
#   FastaFile <- file(paste0("./data/output/pcoding_Longest_Exons_",TransGroup,".fasta"),open = "w")
#   sapply(selectedGeneTrans_split[[TransGroup]], FUN = function(GeneTrans) {
#     TransFASTA <- paste(orderSeqByStrand(getBM(mart = ensembl_Hs, 
#                               attributes = c("gene_exon","strand","exon_chrom_start","exon_chrom_end"), 
#                               filters = list(ensembl_gene_id = strsplit(GeneTrans, split = ";")[[1]][1],
#                                              ensembl_transcript_id = strsplit(GeneTrans, split = ";")[[1]][2])))[,1],
#                         collapse = "")
#     writeLines(paste0(GeneTrans, "\t", TransFASTA), FastaFile)
#     return(NULL)
#   })
#   return(NULL)
# })



