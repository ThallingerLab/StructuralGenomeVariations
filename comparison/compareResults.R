
library(stringr)
library(seqinr)

#BiocManager::install("vcfR")
#BiocManager::install("VariantAnnotation")
#BiocManager::install("StructuralVariantAnnotation")

#library(vcfR)
library(VariantAnnotation)
library(StructuralVariantAnnotation)

loadTruthGR <- function(bed_folder){
  truth_files <- list.files(bed_folder, recursive = T, pattern = "*_summary.bed")
  
  TruthSet <- list()
  
  for(file in truth_files){
   base <- unlist(strsplit(file, "/"))[1]
   
  }
  
  
}

dir = "/Data/Analyses/2022/202205_SV-SIM/accuracy_testing/"
setwd(dir)

delly_vcf <- readVcf("svs/MSv1_f100_l150_m350_s105/DEL-1/delly_50/delly.vcf")
delly_gr <- breakpointRanges(delly_vcf)

gridss_vcf <- readVcf("svs/MSv1_f100_l150_m350_s105/DEL-1/gridss_50/svs.vcf")
gridss_gr <- c(breakpointRanges(gridss_vcf), breakendRanges(gridss_vcf))


manta_vcf <- readVcf("svs/MSv1_f100_l150_m350_s105/DEL-1/manta_50/results/variants/candidateSV.vcf.gz")
manta_gr <- breakpointRanges(manta_vcf)

lumpy_vcf <- readVcf("svs/MSv1_f100_l150_m350_s105/DEL-1/lumpy_50/lumpy_50-smoove.genotyped.vcf.gz")
lumpy_gr <- breakpointRanges(lumpy_vcf)

softsv_vcf <- readVcf("svs/MSv1_f100_l150_m350_s105/DEL-1/softsv_50/svsoft_parsed.vcf")
softsv_gr <- breakpointRanges(softsv_vcf)


manta_gr@elementMetadata

info(header(vcf)) = unique(as(rbind(as.data.frame(info(header(vcf))), data.frame(
  row.names=c("SIMPLE_TYPE"),
  Number=c("1"),
  Type=c("String"),
  Description=c("Simple event type annotation based purely on breakend position and orientation."))), "DataFrame"))

gr <- breakpointRanges(vcf)
svtype <- simpleEventType(gr)

info(vcf)$SIMPLE_TYPE <- "BE"
info(vcf[gr$sourceId])$SIMPLE_TYPE <- svtype
info(vcf[gr$sourceId])$SVLEN <- gr$svLen
ranges(vcf[gr$sourceId]) <- ranges(gr)
strand(vcf[gr$sourceId]) <- strand(gr)
