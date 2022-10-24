library(stringr)
library(seqinr)
library(VariantAnnotation)
library(StructuralVariantAnnotation)

rootdir <- ifelse(as.character(Sys.info())[1] == "Windows", "H://Analyses/", "/Data/Analyses/")
setwd(rootdir)

toolsdir <- ifelse(as.character(Sys.info())[1] == "Windows", "H://Programs/", "~/Projects/StructuralGenomeVariations/comparison/")

source(paste0(toolsdir,"sv_benchmark.R"))
#adapted from https://github.com/PapenfussLab/sv_benchmark/blob/master/R/sv_benchmark.R


findBreakpointOverlaps()

# used only in na12878.R
.distance <- function(r1, r2) {
  return(data.frame(
    min=pmax(0, pmax(start(r1), start(r2)) - pmin(end(r1), end(r2))),
    max=pmax(end(r2) - start(r1), end(r1) - start(r2))))
}

################################################################################



all_callers <- c("bdmax",
                # "breseq",
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
                  # "",
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

truth_set <- loadTruthGR(paste0(rootdir,"2022/202205_SV-SIM/bed"))

VcfCallMetadata <- lapply(vcf_patterns, function(x) ParseMetadata(paste0(rootdir,"2022/202205_SV-SIM/accuracy_testing/svs/HSXn_f100_l150_m550_s165/"),x))

lapply(VcfCallMetadata, length)

vcf <- readVcf(VcfCallMetadata$abs_file[715])
sample_id <- VcfCallMetadata$sample[715]

test_gr <- c(breakpointRanges(vcf, inferMissingBreakends=TRUE),breakendRanges(vcf))

test_scores <- ScoreVariantsFromTruthVCF(test_gr,truthgr = truth_set[[sample_id]], maxgap = 100, ignore.strand = T, id = sample_id)

hits <- as.data.frame(findBreakpointOverlaps(test_gr,  truth_set[[sample_id]], maxgap=10, ignore.strand=T, sizemargin=0.25))

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
