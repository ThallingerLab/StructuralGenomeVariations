library(stringr)
library(seqinr)
library(VariantAnnotation)
library(dplyr)
# devtools::install("/home/veronika/Projects/StructuralVariantAnnotation/")

rootdir <- ifelse(as.character(Sys.info())[1] == "Windows", "H://Analyses/", "/Data/Analyses/")
setwd(rootdir)

toolsdir <- ifelse(as.character(Sys.info())[1] == "Windows", "H://Programs/", "~/Projects/")

source(paste0(toolsdir,"StructuralVariantAnnotation/R/util.R"))
source(paste0(toolsdir,"StructuralGenomeVariations/evaluation/sv_benchmark.R"))
#adapted from https://github.com/PapenfussLab/sv_benchmark/blob/master/R/sv_benchmark.R


# findBreakpointOverlaps()

# used only in na12878.R
# .distance <- function(r1, r2) {
#   return(data.frame(
#     min=pmax(0, pmax(start(r1), start(r2)) - pmin(end(r1), end(r2))),
#     mean=abs(mean(start(r1),end(r1))-mean(start(r2),end(r2))),
#     max=pmax(end(r2) - start(r1), end(r1) - start(r2))
#     ))
# }

################################################################################

all_callers <- c("bdmax",
                 "breseq",
                 "delly",
                 "dysgu",
                 "gridss",
                 "lumpy",
                 "manta",
                 "pindel",
                 "softsv",
                 "svaba",
                 "wham")

vcf_patterns <- c("bdmax\\.vcf$",
                  "output\\.vcf$",
                  "delly\\.vcf$",
                  "dysgu\\.vcf$",
                  "svs\\.vcf$",
                  "smoove\\.genotyped\\.vcf\\.gz$",
                  "diploidSV\\.vcf\\.gz$",
                  "\\.sorted.*?\\.vcf$",
                  "softsv_parsed\\.vcf$",
                  "svaba\\.sv\\.vcf$",
                  "wham\\.vcf$")

names(vcf_patterns) <- all_callers

# manta alternative: candidateSV.vcf.gz (in results/variants/)

################################################################################
#### SIMULATIONS

bedFolder <- "2022/202205_SV-SIM/accuracy_testing/20230223_bed_all/"
vcfFolder <- "2022/202205_SV-SIM/accuracy_testing/20230223_svs_all/MSv1_f100_l150_m350_s105/"

truthSet <- loadTruthGR(paste0(rootdir,bedFolder))

VcfCallMetadata <- lapply(vcf_patterns, function(x) ParseMetadata(paste0(rootdir,vcfFolder),x))

vcfCallSummary <- bind_rows(VcfCallMetadata, .id = "tool")

vcfCallSubset <- vcfCallSummary[vcfCallSummary$readlength == 150 & 
                                  vcfCallSummary$fragmentlength == 350 &
                                  vcfCallSummary$fl_sd == 105 &
                                  vcfCallSummary$coverage == 100 &
                                  vcfCallSummary$machine == "MSv1", ]
SubResults <- list()

if(nrow(vcfCallSubset) > 0){
  
  print(paste0("Analysing Directory: ", sub(rootdir,"", dirname(vcfCallSubset$abs_file[1])), ", has: ", nrow(vcfCallSubset), " entries."))
  
  for(ind in 1:nrow(vcfCallSubset)){
    vcf <- readVcf(vcfCallSubset$abs_file[ind])
    sample_id <- vcfCallSubset$sample[ind]
    # vcf_id <- paste("gridss",paste(metadata[ind,2:7], collapse = ":"), sep = ":")
    vcf_id <- sub(rootdir,"", vcfCallSubset$abs_file[ind])
    maxgap = max(as.numeric(vcfCallSubset$readlength[ind]),
                 (as.numeric(vcfCallSubset$fragmentlength)-as.numeric(vcfCallSubset$readlength[ind])*2))
    
    test_gr <- c(breakpointRanges(vcf, inferMissingBreakends=TRUE),breakendRanges(vcf))
    
    SubResults[[vcf_id]] <- ScoreVariantsFromTruthVCF(test_gr, truthgr = truthSet[[sample_id]], 
                                                     includeFiltered = T, sizemargin = NULL, 
                                                     maxgap = maxgap ,
                                                     ignore.strand = F,
                                                     id = vcf_id)
  }
}

names(SubResults) <- sub(rootdir,"", vcfCallSubset$abs_file)

vcfCallSubset$tp_called_bp <- unlist(lapply(SubResults, function(res){
  sum(res$calls$tp)
}))

vcfCallSubset$tp_called_sourceId <- unlist(lapply(SubResults, function(res){
  sum((res$calls %>%                               # Summary by group using dplyr
           group_by(sourceId) %>% 
           summarize(tp = any(tp)))$tp)
}))

vcfCallSubset$tp_truth_bp <- unlist(lapply(SubResults, function(res){
  sum(res$truth$tp)
}))

vcfCallSubset$tp_truth_sourceId <- unlist(lapply(SubResults, function(res){
  sum((res$truth %>%                               # Summary by group using dplyr
         group_by(sourceId) %>% 
         summarize(tp = any(tp)))$tp)
}))

vcfCallSubset$fp_bp <- unlist(lapply(SubResults, function(res){
  sum(res$calls$fp)
}))

vcfCallSubset$fp_sourceId <- unlist(lapply(SubResults, function(res){
  sum((res$calls %>%                               # Summary by group using dplyr
         group_by(sourceId) %>% 
         summarize(fp = all(fp)))$fp)
}))

vcfCallSubset$fn_bp <- unlist(lapply(SubResults, function(res){
  sum(res$truth$fn)
}))

vcfCallSubset$fp_sourceId <- unlist(lapply(SubResults, function(res){
  sum((res$truth %>%                               # Summary by group using dplyr
         group_by(sourceId) %>% 
         summarize(fn = all(fn)))$fn)
}))

vcfCallSubset$precision_bp <- unlist(lapply(SubResults, function(res){
  sum(res$calls$tp)/nrow(res$calls)
}))

vcfCallSubset$recall_bp <- unlist(lapply(SubResults, function(res){
  sum(res$truth$tp)/nrow(res$truth)
}))

vcfCallSubset$precision_sourceId <- unlist(lapply(SubResults, function(res){
  so <- res$calls %>%                               # Summary by group using dplyr
         group_by(sourceId) %>% 
         summarize(tp = any(tp))
  
  return(sum(so$tp)/nrow(so))
}))

vcfCallSubset$recall_sourceId <- unlist(lapply(SubResults, function(res){
  so <- res$truth %>%                               # Summary by group using dplyr
    group_by(sourceId) %>% 
    summarize(tp = any(tp))
  
  return(sum(so$tp)/nrow(so))
}))

vcfCallSubset[vcfCallSubset$precision_bp == 1,]

vcfCallSubset[vcfCallSubset$precision_bp != vcfCallSubset$precision_sourceId,]


measures <- vcfCallSubset %>%                               # Summary by group using dplyr
  group_by(tool) %>% 
  summarize(rec = sum(tp_truth_bp)/(sum(tp_truth_bp)+sum(fn_bp)),
            pre = sum(tp_called_bp)/(sum(tp_called_bp)+sum(fp_bp)))

measures$f1 <- (measures$rec*measures$pre)/(measures$rec+measures$pre)

RecallPerID <- lapply(unique(vcfCallSubset$sample), function(sample_id){
  vcfIds <- sub(rootdir,"",vcfCallSubset$abs_file[vcfCallSubset$sample == sample_id])
  
  tpCall <- as.data.frame(bind_cols(lapply(SubResults[vcfIds], function(res){
    res$truth$tp
  })))
  
  rownames(tpCall) <- rownames(SubResults[vcfIds][[1]]$truth)
  
  tpCall$sumTP <- rowSums(tpCall)
  
  return(tpCall)
})

names(RecallPerID) <- unique(vcfCallSubset$sample)

notCalled <- lapply(RecallPerID, function(x){
  x$ID <- rownames(x)
  return(x[x$sumTP == 0,c("ID","sumTP")])
})

notCalled

test_sample <- "SUBKO-9"
#test_sample <- "INSP-7"

sub_test_vcf <- paste0("2022/202205_SV-SIM/accuracy_testing/20230223_svs_all/HSXn_f100_l150_m350_s35/",test_sample,"/gridss_100/svs.vcf")
#sub_test_vcf <- paste0("2022/202205_SV-SIM/accuracy_testing/20230223_svs_all/HSXn_f100_l150_m350_s35/",test_sample,"/dysgu_100/dysgu.vcf")
#sub_test_vcf <- paste0("2022/202205_SV-SIM/accuracy_testing/20230223_svs_all/HSXn_f100_l150_m350_s35/",test_sample,"/manta_100/results/variants/diploidSV.vcf.gz")

test_vcf <- readVcf(sub_test_vcf)
test_gr <- c(breakpointRanges(test_vcf, inferMissingBreakends=TRUE),breakendRanges(test_vcf))

res <- ScoreVariantsFromTruthVCF(test_gr, truthgr = truthSet[[test_sample]], maxgap = 10, ignore.strand = T, id = vcf_id)

findBreakpointOverlaps(test_gr, truthSet[[test_sample]], maxgap=10, ignore.strand=FALSE, sizemargin=0.25)

query <- test_gr
subject <- truthSet[[test_sample]]

callgr <- test_gr
truthgr <- truthSet[[test_sample]]

#clearCache()

# AllResults <- lapply(VcfCallMetadata, function(metadata){
#   res <- list()
#   for(ind in 1:nrow(metadata)){
#     vcf <- readVcf(metadata$abs_file[ind])
#     sample_id <- metadata$sample[ind]
#     vcf_id <- paste("gridss",paste(metadata[ind,2:7], collapse = ":"), sep = ":")
#     
#     test_gr <- c(breakpointRanges(vcf, inferMissingBreakends=TRUE),breakendRanges(vcf))
#     
#     res <- c(res, ScoreVariantsFromTruthVCF(test_gr, truthgr = truthSet[[sample_id]], maxgap = 10, ignore.strand = T, id = vcf_id))
#   }
#   return(res)
# })


#Formulas
# R = TP/(TP + FN) ... (TP + FN) == All True variants
# P = TP/(TP + FP) ... (TP + FP) == All called variants
# F1 = 2*P*R/(P + R)

vcf <- readVcf(VcfCallMetadata$gridss$abs_file[ind])
sample_id <- VcfCallMetadata$gridss$sample[ind]
vcf_id <- paste("gridss",paste(VcfCallMetadata$gridss[ind,2:7], collapse = ":"), sep = ":")

test_gr <- c(breakpointRanges(vcf, inferMissingBreakends=TRUE),breakendRanges(vcf))

test_scores <- ScoreVariantsFromTruthVCF(test_gr,truthgr = truthSet[[sample_id]], maxgap = 100, ignore.strand = T, id = vcf_id)

sum(test_scores$calls$tp+test_scores$calls$duptp)
sum(test_scores$truth$tp)


hits <- as.data.frame(findBreakpointOverlaps(test_gr,  truthSet[[sample_id]], maxgap=100, ignore.strand=T, sizemargin=0.25))

test_gr[hits$queryHits]

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
