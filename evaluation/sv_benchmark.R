# library(GenomicRanges)
library(devtools)
library(StructuralVariantAnnotation) #install_github("d-cameron/StructuralVariantAnnotation")
library(testthat)
library(tidyverse)
library(R.cache)

#' Replaces the NA values in a with corresponding values in b
#' @param a,b objects to be tested or coerced.
#' @return The altered object.
'%na%' <- function(a, b) {
	if (is.null(a) || length(a) == 0) return(b)
	if (is.null(b) || length(b) == 0) return(a)
	return(ifelse(is.na(a), b, a))
}

#' Uses b if a is NULL
#' @param a,b objects to be tested or coerced.
#' @return An un-null object.
'%null%' <- function(a, b) {
	if (is.null(a)) return(b)
	return (a)
}

toolsdir <- ifelse(as.character(Sys.info())[1] == "Windows", "W:/Projects/StructuralGenomeVariations/comparison/", "~/Projects/StructuralGenomeVariations/comparison/")

loadTruthGR <- function(directory){
  truth_files <- list.files(directory, recursive = T, pattern = "*_summary.vcf")
  
  TruthSet <- list()
  options(warn=1)
  for(file in truth_files){
    print(paste("Reading file:", file))
    base <- unlist(strsplit(file, "/"))[1]
    vcf <- readVcf(paste(directory,file, sep = "/"))
    TruthSet[[base]] <- c(breakpointRanges(vcf, inferMissingBreakends=TRUE),breakendRanges(vcf, inferMissingBreakends=TRUE))
  }
  options(warn=0)
  return(TruthSet) 
}

GetId <- function(filenames) {
  cf <- as.character(filenames)
  if (length(cf) == 0) {
    return(character(0))
  } else if (all(is.na(cf))) {
    return(cf)
  } else {
    return(str_match(basename(as.character(filenames)), "([^.]+)(\\..*)*$")[,2])
  }
}

ParseMetadata <- function(directory, pattern) {
  key <- list(directory, pattern)
  vcf_table <- loadCache(key=key, dirs=".Rcache/ParseMetadata")
  
  if (is.null(vcf_table)) {
    vcf_table <- .ParseMetadata(directory, pattern)
    saveCache(vcf_table, key=key, dirs=".Rcache/ParseMetadata")
  }
  
  return(vcf_table)
}

.ParseMetadata <- function(directory, pat) {
  
  write("Loading metadata", stderr())
  
  if(endsWith(directory,"/")) directory <- substring(directory,1, nchar(directory)-1)
  
  simSettings <- strsplit(basename(directory),"_")[[1]][1:5]
  simSettings[2:5] <- substr(simSettings[2:5],2,4)
  
  names(simSettings) <- c("machine","coverage","readlength", "fragmentlength","fl_sd")
  
  vcf_table <- data.frame(matrix(ncol = 8, nrow = 0))

  cur_files <- list.files(directory, recursive = T, pattern = pat)
  
  cur_table <- data.frame(abs_file = paste0(directory,"/",cur_files),
                          sample = sub("/.*","",cur_files),
                          machine = simSettings["machine"],
                          readlength = simSettings["readlength"],
                          fragmentlength = simSettings["fragmentlength"],
                          fl_sd = simSettings["fl_sd"],
                          coverage = str_match(cur_files,"/.*_([0-9]*)/")[,2],
                          row.names = NULL)
  
  vcf_table <- rbind(vcf_table, cur_table)

  return(vcf_table)
}

GetId <- function(filenames) {
	cf <- as.character(filenames)
	if (length(cf) == 0) {
		return(character(0))
	} else if (all(is.na(cf))) {
		return(cf)
	} else {
		return(str_match(basename(as.character(filenames)), "([^.]+)(\\..*)*$")[,2])
	}
}
# Load metadata into a dataframe
LoadMetadata <- function(directory) {
	write("Loading metadata", stderr())
	filenames <- list.files(directory, pattern="*.metadata$", full.names=TRUE)
	zeroSizeFiles = file.info(filenames)$size == 0
	if (any(zeroSizeFiles)) {
		warning(paste("Skipping files", filenames[zeroSizeFiles], "as they have 0 size."))
		filenames <- filenames[!zeroSizeFiles]
	}
	#metadata <- foreach (filename=filenames, .export=c("GetId"), .combine=rbind) %dopar% {
	metadata <- lapply(filenames, function(filename) {
		md <- read.csv(filename, header=FALSE, sep="=", quote = "\"'", col.names=c("CX", "V"))
		md$File <- filename
		md$Id <- GetId(filename)
		md
	})
	metadata <- do.call(rbind, metadata)
	metadata <- data.frame(lapply(metadata, as.character), stringsAsFactors=FALSE)
	if (nrow(metadata) == 0) {
		return(data.frame(Id=c()))
	}
	metadata <- spread(metadata, CX, V)
	# convert data from older format
	if (!is.null(metadata$CX_ALIGNER_SOFTCLIP)) {
		metadata$CX_ALIGNER_MODE <- metadata$CX_ALIGNER_MODE %na% ifelse(metadata$CX_ALIGNER_SOFTCLIP == 1, "local", "global")
	}
	# transform known numeric data to expected type
	if (!is.null(metadata$CX_READ_FRAGMENT_LENGTH)) {
		metadata$CX_READ_FRAGMENT_LENGTH <- as.numeric(as.character(metadata$CX_READ_FRAGMENT_LENGTH))
	}
	if (!is.null(metadata$CX_READ_LENGTH)) {
		metadata$CX_READ_LENGTH <- as.numeric(as.character(metadata$CX_READ_LENGTH))
	}
	if (!is.null(metadata$CX_READ_DEPTH)) {
		metadata$CX_READ_DEPTH <- as.numeric(as.character(metadata$CX_READ_DEPTH))
	}
	if (is.null(metadata$CX_MULTIMAPPING_LOCATIONS)) {
		metadata$CX_MULTIMAPPING_LOCATIONS <- NA_integer_
	}
	rownames(metadata) <- metadata$Id
	write(paste(nrow(metadata), "metadata files loaded"), stderr())
	return(metadata)
}
# infer proxy quality scores for ROC purposes based on strength of support
withqual <- function(vcf, caller) {
	if (is.null(VariantAnnotation::fixed(vcf)$QUAL)) {
		fixed(vcf)$QUAL <- NA_real_
	}
	if (any(is.na(VariantAnnotation::fixed(vcf)$QUAL))) {
		if (!is.na(caller) && !is.null(caller)) {
			caller <- str_extract(caller, "^[^/]+") # strip version
			# use total read support as a qual proxy
			if (caller %in% c("delly")) {
				altqual <- ifelse(is.na(info(vcf)$PE), 0, info(vcf)$PE) + ifelse(is.na(info(vcf)$SR), 0, info(vcf)$SR)
			} else if (caller %in% c("crest")) {
				altqual <- ifelse(is.na(info(vcf)$right_softclipped_read_count), 0, info(vcf)$right_softclipped_read_count) + ifelse(is.na(info(vcf)$left_softclipped_read_count), 0, info(vcf)$left_softclipped_read_count)
			} else if (caller %in% c("pindel")) {
				altqual <- geno(vcf)$AD[,1,2]
			} else if (caller %in% c("lumpy")) {
				altqual <- unlist(info(vcf)$SU)
			} else if (caller %in% c("cortex")) {
				# Assuming 1 is ref, 2 is alt
				altqual <- geno(vcf)$COV[,1,2]
				altqual[is.na(altqual)] <- 0
				print(str_c("The number of Cortex calls with altqual 0 is ", sum(is.na(altqual))))
			} else if (caller %in% c("manta")) {
			    altqual <- 0
			}
			VariantAnnotation::fixed(vcf)$QUAL <- ifelse(is.na(VariantAnnotation::fixed(vcf)$QUAL), altqual, VariantAnnotation::fixed(vcf)$QUAL)
		} else {
			# use a placeholder QUAL for truth sets
			VariantAnnotation::fixed(vcf)$QUAL <- 1
		}
	}
	if (any(is.na(VariantAnnotation::fixed(vcf)$QUAL))) {
		#if (is.null(caller) && is.na(caller)) {
		warning(paste("Missing QUAL scores for", caller))
		stop(paste("Missing QUAL scores for", caller))
	}
	return(vcf)
}

StripCallerVersion <- function(caller, gridssfirst = FALSE) {
	if (length(caller) == 0) return(caller)
	caller <- paste0(str_extract(caller, "^([^/]+)"), str_match(caller, "^([^/]+)\\/[^/]+(/[^/]+)?")[,3] %na% "") %na% caller
	if (gridssfirst && any(caller=="gridss")) {
		caller <- relevel(factor(caller), "gridss")
	}
	return(caller)
}
PrettyVariants <- function(x) {
	x[x=="hetDEL"] <- "Deletion"
	x[x=="hetINS"] <- "Insertion"
	x[x=="hetDUP"] <- "Tandem Duplication"
	x[x=="hetINV"] <- "Inversion"
	x[x=="hetBP"] <- "Breakpoint"
	x[x=="hetBP_SINE"] <- "SINE/ALU Breakpoint"
	return(x)
}
#' Loads a minimal structural variant GRanges from the VCF
LoadMinimalSVs <- function(filename, caller, transform, nominalPosition) {
	key <- list(filename, caller, transform, nominalPosition)
	gr <- loadCache(key=key, dirs=".Rcache/LoadMinimalSVs")
	if (is.null(gr)) {
		gr <- .LoadMinimalSVs(filename, caller, transform, nominalPosition)
		saveCache(gr, key=key, dirs=".Rcache/LoadMinimalSVs")
	}
	if (length(gr) > 0) {
		seqlevelsStyle(gr) <- "UCSC"
	}
	return(gr)
}
.LoadMinimalSVs <- function(filename, caller, transform, nominalPosition) {
	vcf <- readVcf(filename, "hg19")
	if (!is.null(transform)) {
		vcf <- transform(vcf)
	}
	vcf <- withqual(vcf, caller)
	gr <- breakpointRanges(vcf, nominalPosition)
	gr$paramRangeID <- NULL
	gr$REF <- NULL
	gr$ALT <- NULL
	#gr$svtype <- NULL
	#gr$svLen <- NULL
	gr$insSeq <- NULL
	#gr$insLen <- NULL
	return(gr)
}
#vcf <- readVcf("C:/dev/sv_benchmark/data.aligner/5afa7ffdf2cc32602476526d5b477c5c.vcf", "hg19")
#' Loads structural variant GRanges from the VCFs in the given directory
LoadMinimalSVFromVCF <- function(directory, pattern="*.vcf$", metadata=NULL, existingList=NULL, transform=NULL, nominalPosition=FALSE) {
	write("Loading VCFs", stderr())
	filenames <- list.files(directory, pattern=pattern, full.names=TRUE)
	zeroSizeFiles = file.info(filenames)$size == 0
	if (any(zeroSizeFiles)) {
		write(paste("Skipping file", filenames[zeroSizeFiles], "due to 0 size.\n"))
		warning(paste("Skipping files", paste(filenames[zeroSizeFiles]), "due to 0 size.\n"))
		filenames <- filenames[!zeroSizeFiles]
	}
	# exclude already loaded VCFs
	filenames <- filenames[!(GetId(filenames) %in% names(existingList))]
	# exclude VCFs without metadata
	filenames <- filenames[GetId(filenames) %in% metadata$Id]
	# only load VCFS that have metadata
	#if (!is.null(metadata)) {
	#	filenames <- filenames[GetId(filenames) %in% metadata$Id]
	#}
	#vcfs <- foreach (filename=filenames, .packages="VariantAnnotation") %dopar% { # Parallel load of VCFs
	grlist <- lapply(filenames, function(filename) {
		write(paste0("Loading ", filename), stderr())
		caller <- NULL
		if (!is.null(metadata)) {
			caller <- metadata$CX_CALLER[metadata$Id == GetId(filename)]
		}
		gr <- LoadMinimalSVs(filename, caller, transform, nominalPosition)
		gr$Id <- rep(GetId(filename), length(gr))
		return (gr)
	})
	names(grlist) <- GetId(filenames)
	grlist[sapply(grlist, is.null)] <- NULL # Remove NULL VCFs list
	write(paste("Loaded", length(grlist), "VCFs"), stderr())
	return(c(existingList, grlist))
}
.interval_distance <- function(s1, e1, s2, e2) {
	return (ifelse(s2 >= s1 & s2 <= e1, 0,
					ifelse(s1 >= s2 & s1 <= e2, 0,
					ifelse(s1 < s2, s2 - e1, s1 - e2))))
}
# used only in na12878.R
# .distance <- function(r1, r2) {
#   return(data.frame(
#     min=pmax(0, pmax(start(r1), start(r2)) - pmin(end(r1), end(r2))),
#     max=pmax(end(r2) - start(r1), end(r1) - start(r2))
#   ))
# }

.distance <- function(r1, r2) {
  return(data.frame(
    min=pmax(0, pmax(start(r1), start(r2)) - pmin(end(r1), end(r2))),
    mean=abs((start(r1)+end(r1))/2-(start(r2)+end(r2))/2),
    max=pmax(end(r2) - start(r1), end(r1) - start(r2))
  ))
}

# used only in na12878.R
findMatchingBreakpoints <- function(query, subject, maxgap=0L, ignore.strand=FALSE, sizemargin=0.25, restrictMarginToSizeMultiple=0.5) {
  hits <- as.data.frame(findBreakpointOverlaps(query, subject, maxgap=maxgap, ignore.strand=ignore.strand))
  # take into account confidence intervals when calculating event size
  callwidth <- .distance(query, partner(query))
  truthwidth <- .distance(subject, partner(subject))
  callsize <- callwidth + (query$insLen %na% 0)
  truthsize <- truthwidth + (subject$insLen %na% 0)
  hits$sizeerror <- .distance(IRanges(start=callsize$min[hits$queryHits], end=callsize$max[hits$queryHits]),
                              IRanges(start=truthsize$min[hits$subjectHits], end=truthsize$max[hits$subjectHits]))$min
  # event sizes must be within sizemargin
  if (!is.null(sizemargin)) {
    hits <- hits[hits$sizeerror - 2 < sizemargin * pmin(callsize$max[hits$queryHits], truthsize$max[hits$subjectHits]),]
  }
  # further restrict breakpoint positions for small events
  hits$localbperror <- .distance(query[hits$queryHits], subject[hits$subjectHits])$mean
  hits$remotebperror <- .distance(partner(query)[hits$queryHits], partner(subject)[hits$subjectHits])$mean
  if (!is.null(restrictMarginToSizeMultiple)) {
    allowablePositionError <- (pmin(callsize$max[hits$queryHits], truthsize$max[hits$subjectHits]) * restrictMarginToSizeMultiple + 2)
    hits <- hits[hits$localbperror <= allowablePositionError & hits$remotebperror <= allowablePositionError, ]
  }
  return(hits)
}

#' Finds pairs of query breakpoints that span a matching subject breakpoint
findSpanningBreakpoints <- function(query, subject, maxgap=0L, ignore.strand=FALSE, sizemargin=0.25, restrictMarginToSizeMultiple=0.5, maxSpanningFragmentSize,
																		matchDirection=TRUE) {
	# must be larger than min fragment size (if not, small indel + large event get matched to large event)
	query <- query[is.na(query$svLen) | abs(query$svLen) >= maxSpanningFragmentSize]
	# <A1-------A2>	 <B1-----B2>
	#						 |-d-|
	# find A2-B1 matches
	# d must be < maxSpanningFragmentSize
	spanningHits <- as.data.frame(findOverlaps(query, query, maxgap=maxSpanningFragmentSize, ignore.strand=TRUE))
	# removes self-intersections and duplication due to picking up both sides
	spanningHits <- spanningHits[spanningHits$subjectHits < spanningHits$queryHits,]
	a2 <- query[spanningHits$subjectHits]
	b1 <- query[spanningHits$queryHits]
	spanningHits$fragmentSize <- abs((start(a2) + end(a2)) / 2 - (start(b1) + end(b1)) / 2)
	# directions must be appropriate for a small fragment
	if (matchDirection) {
		spanningHits <- spanningHits[ifelse(paste0(strand(a2), strand(b1)) == "+-", start(a2) < end(b1),
																 ifelse(paste0(strand(b1), strand(a2)) == "+-", start(b1) < end(a2), FALSE)),]
	}
	# generate A1-B2 breakpoint gr
	a1 <- partner(query)[spanningHits$subjectHits]
	b2 <- partner(query)[spanningHits$queryHits]
	spanninggr <- c(a1, b2)
	spanninggr$partner <- c(names(b2), names(a1))
	spanninggr$fragmentSize <- spanningHits$fragmentSize
	# spanning against subject
	hits <- findMatchingBreakpoints(spanninggr, subject, maxgap, ignore.strand, sizemargin, restrictMarginToSizeMultiple)
	hits$localBreakend <- names(spanninggr)[hits$queryHits]
	hits$remoteBreakend <- names(partner(spanninggr))[hits$queryHits]
	hits$fragmentSize <- spanninggr$fragmentSize[hits$queryHits]
	hits$queryHits <- NULL
	return(hits)
}

LoadAndScoreCalls <- function(metadata, truthdata, includeFiltered=FALSE, maxgap, ignore.strand, sizemargin=0.25, id=NULL, requiredHits=1){
  
  scores <- list()
  
  for(i in 1:nrow(metadata)){
    vcf <- readVcf(metadata$abs_file[i])
    gr <- c(breakpointRanges(vcf, inferMissingBreakends=TRUE),breakendRanges(vcf))
# 
    scores <- c(scores, ScoreVariantsFromTruthVCF(callgr = gr,
                                             truthgr = truthdata[[metadata$sample[i]]],
                                             maxgap = maxgap,
                                             ignore.strand = ignore.strand,
                                             sizemargin = sizemargin,
                                             id = metadata$abs_file[i],
                                             requiredHits = requiredHits ))
  }
  
  result <- list(
    calls=bind_rows(lapply(scores, function(x) x$calls)),
    truth=bind_rows(lapply(scores, function(x) x$truth)))
  if (length(result$truth) == 0) {
    result$truth <- NULL
  }
  return(result)  
}


#' @param keytruth unique identifier of truthgr if it is not the 'natural' truth to compare to
ScoreVariantsFromTruthVCF <- function(callgr, truthgr, includeFiltered=FALSE, maxgap, ignore.strand, sizemargin=0.25, id=NULL, requiredHits=1, keytruth=NULL, keycalls=NULL) {
	
  if (length(callgr) == 0) {
		return(.ScoreVariantsFromTruthVCF(callgr, truthgr, includeFiltered, maxgap, ignore.strand, sizemargin, id %null% NA_character_, requiredHits))
  }
  
	id <- id %null% callgr$Id[1]
	key <- list(includeFiltered, maxgap, ignore.strand, sizemargin, id, requiredHits, keytruth, keycalls)
	result <- loadCache(key=key, dirs=".Rcache/ScoreVariantsFromTruthVCF")
	
	if (is.null(result)) {
		write(paste0("ScoreVariantsFromTruth ", id), stderr())
		result <- .ScoreVariantsFromTruthVCF(callgr, truthgr, includeFiltered, maxgap, ignore.strand, sizemargin, id, requiredHits)
		saveCache(result, key=key, dirs=".Rcache/ScoreVariantsFromTruthVCF")
	} else {
		cat(".", file=stderr())
	}
	
	return(result)
}

#' @param requiredHits number of matching truth hits before being considered a true positive.
#' A value greater than 1 is useful when comparing directly to split long reads
#' and long read indels
.ScoreVariantsFromTruthVCF <- function(callgr, truthgr, includeFiltered, maxgap, ignore.strand, sizemargin, id, requiredHits) {
	if (!includeFiltered) {
		callgr <- callgr[callgr$FILTER %in% c("PASS", "."),]
	}
	if (is.null(truthgr)) {
		stop(paste("Missing truth ", truthid, " for ", id))
	}
	if (!all(callgr$partner %in% names(callgr))) {
		warning(paste(sum(!(callgr$partner %in% names(callgr))), "breakends missing partners. Ignoring."))
		callgr <- callgr[callgr$partner %in% names(callgr)]
	}
	if (is.null(callgr$ihomlen)) {
		callgr$ihomlen <- rep(NA_integer_, length(callgr))
	}
	if (is.null(truthgr$ihomlen)) {
		truthgr$ihomlen <- rep(NA_integer_, length(truthgr))
	}
	hits <- as.data.frame(findBreakpointOverlaps(callgr, truthgr, maxgap=maxgap, ignore.strand=ignore.strand, sizemargin=sizemargin))

	hits$QUAL <- callgr$QUAL[hits$queryHits]
	hits <- hits[order(-hits$QUAL),] # sort by qual so the highest QUAL writes last when doing hit assignments on subjectHits or queryHits
	hits$dup <- duplicated(hits$subjectHits)
	hitcount <- hits %>%
		group_by(queryHits) %>%
		summarise(besthits=sum(!dup), allhits=n())

	calldf <- as.data.frame(callgr) %>%
		dplyr::select(QUAL, svLen, insLen, sourceId, HOMLEN, ihomlen) %>%
		mutate(
			Id=id,
			tp=rep(FALSE, nrow(.)),
			duptp=rep(FALSE, nrow(.)),
			fp=rep(FALSE, nrow(.)),
			fn=rep(FALSE, nrow(.)),
			sizeerror=rep(NA, nrow(.)),
			localbperror=NA,
			remotebperror=NA,
			includeFiltered=rep(includeFiltered, nrow(.)),
			maxgap=rep(maxgap, nrow(.)),
			ignore.strand=rep(ignore.strand, nrow(.)))
	
	calldf$tp[hitcount$queryHits[hitcount$besthits >= requiredHits]] <- TRUE
	calldf$duptp[hitcount$queryHits[hitcount$allhits >= requiredHits]] <- TRUE
	calldf$duptp <- calldf$duptp & !calldf$tp
	calldf$fp <- !calldf$tp
	calldf$localbperror[hits$queryHits] <- hits$localbperror
	calldf$remotebperror[hits$queryHits] <- hits$remotebperror
	#calldf$sizeerror[hits$queryHits] <- hits$sizeerror
	calldf$simpleEvent <- simpleEventType(callgr)
	calldf$repeatClass <- callgr$repeatClass
	calldf$breakendId <- names(callgr)

	truthdf <- NULL
	if (requiredHits == 1) {
		truthdf <- as.data.frame(truthgr) %>%
			dplyr::select(svLen, insLen, sourceId, HOMLEN, ihomlen) %>%
			mutate(Id=id, QUAL=0, tp=FALSE, fp=FALSE, fn=FALSE, localbperror=NA, remotebperror=NA) %>%
			mutate(
				includeFiltered=includeFiltered,
				maxgap=maxgap,
				ignore.strand=ignore.strand)
		truthdf$tp[hits$subjectHits] <- TRUE
		truthdf$fn <- !truthdf$tp
		truthdf$QUAL[hits$subjectHits] <- hits$QUAL
		truthdf$localbperror[hits$subjectHits] <- hits$localbperror
		truthdf$remotebperror[hits$subjectHits] <- hits$remotebperror
		#truthdf$sizeerror[hits$subjectHits] <- hits$sizeerror
		truthdf$simpleEvent <- simpleEventType(truthgr)
		truthdf$repeatClass <- truthgr$repeatClass
		# using caller-defined homology length
		truthdf$HOMLEN[hits$subjectHits] <- callgr[hits$queryHits]$HOMLEN
		truthdf$breakendId <- names(truthgr)
	}
	return(list(calls=calldf, truth=truthdf))
}

simpleEventType <- function(gr) {

  gr$et <- ifelse(seqnames(gr) != seqnames(partner(gr)), "BP",
          ifelse(gr$insLen >= abs(gr$svLen) * 0.7, "INS",
           ifelse(strand(gr) == strand(partner(gr)), "INV",
            ifelse(xor(start(gr) < start(partner(gr)), strand(gr) == "-"), "DEL",
             "DUP"))))

  etp <- partner(gr)$et
  return(ifelse(gr$et < etp, gr$et, etp))
}

ScoreVariantsFromTruth <- function(vcfs, metadata, includeFiltered=FALSE, maxgap, ignore.strand, sizemargin=0.25, requiredHits=1, truthgr=NULL, keytruth=NULL, keycalls=NULL) {
	ids <- metadata$Id[!is.na(metadata$CX_CALLER) & metadata$Id %in% names(vcfs)]
	scores <- lapply(ids, function(id) {
		callgr <- vcfs[[id]]
		if (is.null(truthgr)) {
			truthid <- GetId((metadata %>% filter(Id==id))$CX_REFERENCE_VCF)
			truthgr <- vcfs[[truthid]]
		}
		if (is.null(truthgr)) {
			stop("Missing truth for ", id)
		}
		if (length(callgr) == 0) {
			return(list(calls=NULL))
		}
		return(ScoreVariantsFromTruthVCF(callgr, truthgr, includeFiltered, maxgap, ignore.strand, sizemargin, id, requiredHits=requiredHits, keytruth=keytruth, keycalls=keycalls))
	})
	result <- list(
		calls=bind_rows(lapply(scores, function(x) x$calls)),
		truth=bind_rows(lapply(scores, function(x) x$truth)))
	if (length(result$truth) == 0) {
		result$truth <- NULL
	}
	return(result)
}

#' subsets the breakpoints to only include breakpoints in which both breakends
#' occur within the specified bed regions
subsetbed <- function(gr, bed, maxgap) {
	gr <- subsetByOverlaps(gr, bed, maxgap=maxgap, ignore.strand=TRUE)
	gr <- gr[gr$partner %in% names(gr)]
	return(gr)
}
#' imports SVs in BEDPE breakpoint format
import.sv.bedpe <- function(file, placeholderName="bedpe") {
	df <- read.table(file, col.names=c("chr1", "start1", "end1", "chr2", "start2", "end2", "id", "score", "strand1", "strand2", "info"), stringsAsFactors=FALSE)
	# ensure row names are unique
	row.names(df) <- ifelse(duplicated(df$id), paste0(placeholderName, seq_along(df$chr1)), df$id)
	gro <- GRanges(seqnames=df$chr1, strand=df$strand1, ranges=IRanges(df$start1, df$end1), id=df$id, score=df$score, info=df$info)
	grh <- GRanges(seqnames=df$chr2, strand=df$strand2, ranges=IRanges(df$start2, df$end2), id=df$id, score=df$score, info=df$info)
	names(gro) <- paste0(row.names(df), "_1")
	names(grh) <- paste0(row.names(df), "_2")
	gro$partner <- names(grh)
	grh$partner <- names(gro)
	return(c(gro, grh))
}

# find . -name '*.fa.out' -exec tail -n +4 {} \;  > merged.fa.out
# list.files(paste0(referenceLocation, "/UCSC/repeatmasker/"), pattern="*.fa.out", recursive=TRUE, full.names=TRUE),
import.repeatmasker.fa.out <- function(repeatmasker.fa.out) {
	rmdt <- read_table2(repeatmasker.fa.out, col_names=FALSE, skip=3)
	grrm <- GRanges(
	  seqnames=rmdt$X5,
	  ranges=IRanges(start=rmdt$X6 + 1, end=rmdt$X7),
	  repeatType=rmdt$X10,
	  repeatClass=rmdt$X11)
	grrm$repeatClass <- str_replace(str_replace(grrm$repeatClass, "[?]", ""), "/.*", "")
	return(grrm)
}

findBreakpointOverlaps <- function(query, subject, maxgap=-1L, minoverlap=0L, ignore.strand=FALSE, sizemargin=NULL, restrictMarginToSizeMultiple=NULL) {
  # .assertValidBreakpointGRanges(query)
  # .assertValidBreakpointGRanges(subject)
  pquery = partner(query, selfPartnerSingleBreakends = T)
  squery = partner(subject, selfPartnerSingleBreakends = T)
  localhits = findOverlaps(query, subject, maxgap=maxgap, minoverlap=minoverlap, type="any", select="all", ignore.strand=ignore.strand)
  remotehits = findOverlaps(pquery, squery, maxgap=maxgap, minoverlap=minoverlap, type="any", select="all", ignore.strand=ignore.strand)
  ## duplicated() version:
  #hits = Hits(c(S4Vectors::queryHits(localhits), S4Vectors::queryHits(remotehits)), c(S4Vectors::subjectHits(localhits), S4Vectors::subjectHits(remotehits)), nLnode=nLnode(localhits), nRnode=nRnode(localhits), sort.by.query=TRUE)
  #hits = hits[duplicated(hits)]
  
  ## intersect() version:
  hits = BiocGenerics::intersect(localhits, remotehits)
  
  ## dplyr() version:
  #hits <- dplyr::bind_rows(
  #	as.data.frame(localhits, row.names=NULL),
  #	as.data.frame(remotehits, row.names=NULL))
  #hits = hits %>% dplyr::arrange(queryHits, subjectHits) %>%
  #	dplyr::filter(!is.na(dplyr::lead(.$queryHits)) & !is.na(dplyr::lead(.$subjectHits)) & dplyr::lead(.$queryHits) == .$queryHits & dplyr::lead(.$subjectHits) == .$subjectHits)
  
  ## dplyr() exploiting the sorted nature of the findOverlaps():
  #hits = Hits(c(S4Vectors::queryHits(localhits), S4Vectors::queryHits(remotehits)), c(S4Vectors::subjectHits(localhits), S4Vectors::subjectHits(remotehits)), nLnode=nLnode(localhits), nRnode=nRnode(localhits), sort.by.query=TRUE)
  #queryLead  = dplyr::lead(S4Vectors::queryHits(hits))
  #querySubject  = dplyr::lead(S4Vectors::queryHits(hits))
  #hits = hits[
  #	!is.na(queryLead) &d
  #	!is.na(querySubject) &
  #	queryLead == S4Vectors::queryHits(hits) &
  #	querySubject == S4Vectors::subjectHits(hits)]
  if (!is.null(sizemargin) && !is.na(sizemargin)) {
    # take into account confidence intervals when calculating event size
    callwidth <- .distance(query, pquery)
    truthwidth <- .distance(subject, squery)
    callsize <- callwidth + .replaceNa(query$insLen, 0)
    truthsize <- truthwidth + .replaceNa(subject$insLen, 0)
    sizeerror <- .distance(
      IRanges::IRanges(start=callsize$min[S4Vectors::queryHits(hits)], end=callsize$max[S4Vectors::queryHits(hits)]),
      IRanges::IRanges(start=truthsize$min[S4Vectors::subjectHits(hits)], end=truthsize$max[S4Vectors::subjectHits(hits)])
    )$min
    # event sizes must be within sizemargin
    hits <- hits[sizeerror - 1 < sizemargin * pmin(callsize$max[S4Vectors::queryHits(hits)], truthsize$max[S4Vectors::subjectHits(hits)]),]
    
    # further restrict breakpoint positions for small events
    if (!is.null(restrictMarginToSizeMultiple)) {
      localbperror <- .distance(query[S4Vectors::queryHits(hits)], subject[S4Vectors::subjectHits(hits)])$min
      remotebperror <- .distance(pquery[S4Vectors::queryHits(hits)], squery[S4Vectors::subjectHits(hits)])$min
      
      allowablePositionError <- (pmin(callsize$max[S4Vectors::queryHits(hits)], truthsize$max[S4Vectors::subjectHits(hits)]) * restrictMarginToSizeMultiple + 1)
      hits <- hits[localbperror <= allowablePositionError & remotebperror <= allowablePositionError, ]
    }
  }
  
  hits_frame <- as.data.frame(hits)
  hits_frame$localbperror <- .distance(query[S4Vectors::queryHits(hits)], subject[S4Vectors::subjectHits(hits)])$mean
  hits_frame$remotebperror <- .distance(pquery[S4Vectors::queryHits(hits)], squery[S4Vectors::subjectHits(hits)])$mean
  
  return(hits_frame)
}
# 
# findBreakEndOverlaps <- function (query, subject, maxgap = -1L, minoverlap = 0L, ignore.strand = FALSE, 
#                                     sizemargin = NULL, restrictMarginToSizeMultiple = NULL) 
# {
#   .assertValidBreakpointGRanges(subject)
#   squery = partner(subject)
#   
#   localhits = findOverlaps(query, subject, maxgap = maxgap, 
#                            minoverlap = minoverlap, type = "any", select = "all", 
#                            ignore.strand = ignore.strand)
#   
#   remotehits = findOverlaps(query, squery, maxgap = maxgap, 
#                             minoverlap = minoverlap, type = "any", select = "all", 
#                             ignore.strand = ignore.strand)
#   
#   hits = BiocGenerics::intersect(localhits, remotehits)
#   
#   
#   if (!is.null(sizemargin) && !is.na(sizemargin)) {
#     localbperror <- .distance(query[S4Vectors::queryHits(hits)], 
#                               subject[S4Vectors::subjectHits(hits)])$min
#   }
#   
#   return(hits)
# }