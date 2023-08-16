setwd("/Data/Analyses/2023/202308_CloveBiotech/")


################################################################################
library(stringr)

annotCompareRes <- function(compareTab){
  if("vcffile" %in% colnames(compareTab)){
    
    list <- lapply(str_split(compareTab$vcffile,"/"), function(x) {
      n <- which(str_detect(x,".*_f.*_l.*_m.*_s.*"))
      
      machine <- str_match(x[n],"(.*)_f")[2]
      len <- str_match(x[n],"_l([0-9]*)_")[2]
      ins <- str_match(x[n],"_m([0-9]*)_")[2]
      insstd <- str_match(x[n],"_s([0-9]*)")[2]
      
      maintype <- x[n+1]
      
      toolCov <- strsplit(x[n+2],"_")[[1]]
      
      fildered <- ifelse(str_detect(x[n+3],"pass"),"pass","raw")
      
      namedVec <- c(machine,maintype,len, ins, insstd, toolCov, fildered)
      names(namedVec) <- c("machine","type","readLen","insLen","insStd","tool","coverage","pass")
      
      return(namedVec)
      
    })
    return(list)
  } else {
    stop("Error: No column called \'vcffile\'")
  }
} 

################################################################################

compareHeader <- c("bedfiled","vcffile","ntrue","ncalled","TP","FN","FP")

summaries <- c("20221120_bed/summary-clovebiotech.tsv","20220819_bed/summary-clovebiotech.tsv",
               "20221120_bed/summary-clove.tsv","20220819_bed/summary-clove.tsv")#,"20230310_bed/summary-clovebiotech.tsv")

compareRes <- do.call(rbind, lapply(summaries, function(x){
  tab <- read.delim2(x, header = F, stringsAsFactors = F)
  colnames(tab) <- c("bedfiled","vcffile","ntrue","ncalled","TP","FN","FP")
  annot <- annotCompareRes(tab)
  fullTab <- cbind(tab, as.data.frame(do.call(rbind, annot)))
  return( fullTab )
}))

#############################################

compareRes <- transform(compareRes, ntrue = as.numeric(ntrue), 
                        ncalled = as.numeric(ncalled),
                        TP = as.numeric(TP),
                        FP = as.numeric(FP),
                        FN = as.numeric(FN),
                        readLen = as.numeric(readLen),
                      #  insLen = as.numeric(insLen),
                      #  insStd = as.numeric(insStd),
                        coverage = as.numeric(coverage),
                        pass = relevel(as.factor(pass),"pass")
                        )

comparSummary <- aggregate(compareRes[c("ntrue","ncalled","TP","FN","FP")],
          by=compareRes[c("tool", "machine", "readLen","insLen","insStd","coverage","pass")],
          FUN = sum)

comparSummary <- comparSummary[order(comparSummary$coverage, comparSummary$insLen, comparSummary$tool, comparSummary$pass),]

comparSummary$group <- paste0(comparSummary$tool, comparSummary$pass)

groupColors <- c("#f768a1", "#c51b8a","#67a9cf", "#1c9099")
names(groupColors)<- unique(comparSummary$group)

comparSummary$tpr <- (comparSummary$TP*100/comparSummary$ntrue)

library(ggplot2)

miseq_tp <- ggplot(comparSummary[comparSummary$readLen == "150",], aes(fill = group, x = insLen, y = TP)) +
    geom_col(position=position_dodge(c(0.8)), width = 0.7)  +
    geom_text(aes(label = TP), position = position_dodge(c(0.8)), color = "white", hjust = 1.5) +
    scale_fill_manual("legend", values = groupColors)+
    facet_wrap(~coverage) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.title = element_text(size = 20)) +
  xlab("insert size") + ylab("true positives") +
  ylim(0,max(comparSummary$ntrue))

svg("~/Nextcloud/Publications/2023_CloveBiotech/CloveBiotech-150_TP.svg")

miseq_tp + coord_flip()

dev.off()

miseq_fp <- ggplot(comparSummary[comparSummary$readLen == "150",], aes(fill = group, x = insLen, y = FP)) +
  geom_col(position=position_dodge(0.8), width = 0.7)  +
  geom_text(aes(label = FP), position = position_dodge(.8), color = "white", hjust = 1.5) +
  scale_fill_manual("legend", values = groupColors)+
  facet_wrap(~coverage) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.title = element_text(size = 20)) +
  xlab("insert size") +  ylab("false positives") +
  ylim(0,max(comparSummary$FP)+5)

# 
# hiseq_n <- ggplot(compar_summary[compar_summary$machine == "HSXn" & compar_summary$group != "exact", ], aes(fill = group, x = insert_size, y = called_additional)) +
#   geom_col(position=position_dodge(0.75), width = 0.7)  +
#   geom_text(aes(label = called_additional), position = position_dodge(.7), color = "black", hjust = -0.15) +
#   scale_fill_manual("legend", values = c("raw" = "#67a9cf", "pass" = "#1c9099", "exact" = "#016c59"))+
#   facet_wrap(~cov) +
#   theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.title = element_text(size = 20)) +
#   xlab("insert size") +  ylab("false positives") +
#   ylim(0,max(compar_summary$called_additional)+5)

svg("~/Nextcloud/Publications/2023_CloveBiotech/CloveBiotech-150_FP.svg")

miseq_fp + coord_flip()

dev.off()

miseq_tp <- ggplot(comparSummary[comparSummary$readLen == "250",], aes(fill = group, x = insLen, y = TP)) +
  geom_col(position=position_dodge(c(0.8)), width = 0.7)  +
  geom_text(aes(label = TP), position = position_dodge(c(0.8)), color = "white", hjust = 1.5) +
  scale_fill_manual("legend", values = groupColors)+
  facet_wrap(~coverage) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.title = element_text(size = 20)) +
  xlab("insert size") + ylab("true positives") +
  ylim(0,max(comparSummary$ntrue))

svg("~/Nextcloud/Publications/2023_CloveBiotech/CloveBiotech-250_TP.svg")

miseq_tp + coord_flip()

dev.off()

miseq_fp <- ggplot(comparSummary[comparSummary$readLen == "250",], aes(fill = group, x = insLen, y = FP)) +
  geom_col(position=position_dodge(0.8), width = 0.7)  +
  geom_text(aes(label = FP), position = position_dodge(.8), color = "white", hjust = 1.5) +
  scale_fill_manual("legend", values = groupColors)+
  facet_wrap(~coverage) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.title = element_text(size = 20)) +
  xlab("insert size") +  ylab("false positives") +
  ylim(0,max(comparSummary$FP)+5)

# 
# hiseq_n <- ggplot(compar_summary[compar_summary$machine == "HSXn" & compar_summary$group != "exact", ], aes(fill = group, x = insert_size, y = called_additional)) +
#   geom_col(position=position_dodge(0.75), width = 0.7)  +
#   geom_text(aes(label = called_additional), position = position_dodge(.7), color = "black", hjust = -0.15) +
#   scale_fill_manual("legend", values = c("raw" = "#67a9cf", "pass" = "#1c9099", "exact" = "#016c59"))+
#   facet_wrap(~cov) +
#   theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.title = element_text(size = 20)) +
#   xlab("insert size") +  ylab("false positives") +
#   ylim(0,max(compar_summary$called_additional)+5)

svg("~/Nextcloud/Publications/2023_CloveBiotech/CloveBiotech-250_FP.svg")

miseq_fp + coord_flip()

dev.off()

comparSummary$recall <- comparSummary$TP/comparSummary$ntrue
comparSummary$precion <- comparSummary$TP/comparSummary$ncalled
comparSummary$fscore <- 2*(comparSummary$recall*comparSummary$precion)/(comparSummary$recall+comparSummary$precion)

comparSummary[order(-comparSummary$precion),]
comparSummary[order(-comparSummary$recall),]
comparSummary[order(-comparSummary$fscore),]


comparSummaryCBR <- comparSummary[comparSummary$pass == "raw" & comparSummary$tool == "clovebiotech",]
cor(as.numeric(comparSummaryCBR$insLen),comparSummaryCBR$FP)
cor(as.numeric(comparSummaryCBR$coverage),comparSummaryCBR$FP)



################################################################################

comparType <- aggregate(compareRes[c("ntrue","ncalled","TP","FN","FP")],
                           by=compareRes[c("type","tool","pass")],
                           FUN = sum)

comparTypeCB <- comparType[comparType$tool == "clovebiotech" & comparType$pass == "raw",]

comparTypeCB[order(-comparTypeCB$FN),] 

################################################################################

summariesPlasmids <- c("20230310_bed/summary-clove.tsv","20230310_bed/summary-clovebiotech.tsv")

compareResPlasmids <- do.call(rbind, lapply(summariesPlasmids, function(x){
  tab <- read.delim2(x, header = F, stringsAsFactors = F)
  colnames(tab) <- c("bedfiled","vcffile","ntrue","ncalled","TP","FN","FP")
  annot <- annotCompareRes(tab)
  fullTab <- cbind(tab, as.data.frame(do.call(rbind, annot)))
  return( fullTab )
}))


comparSummaryPlasmids <- aggregate(compareResPlasmids[c("ntrue","ncalled","TP","FN","FP")],
                           by=compareResPlasmids[c("tool", "machine", "readLen","insLen","insStd","coverage","pass")],
                           FUN = sum)

comparSummaryPlasmids <- comparSummaryPlasmids[order(comparSummaryPlasmids$coverage,
                                                     comparSummaryPlasmids$insLen, comparSummaryPlasmids$tool, comparSummaryPlasmids$pass),]

comparSummaryPlasmids$group <- paste0(comparSummaryPlasmids$tool, comparSummaryPlasmids$pass, "plasmids")

comparSummaryPlasmids$plasmid <- T

#######################

compareResNoPlasmids <- compareRes[sub("-.*","",compareRes$type) %in% c("INSP","SUBP"), ]

comparSummaryNoPlasmids <- aggregate(compareResNoPlasmids[c("ntrue","ncalled","TP","FN","FP")],
                                   by=compareResNoPlasmids[c("tool", "machine", "readLen","insLen","insStd","coverage","pass")],
                                   FUN = sum)


comparSummaryNoPlasmids <- comparSummaryNoPlasmids[order(comparSummaryNoPlasmids$coverage,
                                                     comparSummaryNoPlasmids$insLen, comparSummaryNoPlasmids$tool, comparSummaryNoPlasmids$pass),]

comparSummaryNoPlasmids$group <- paste0(comparSummaryNoPlasmids$tool, comparSummaryNoPlasmids$pass)
comparSummaryNoPlasmids$plasmid <- F

plasmidSummaries <- rbind(comparSummaryPlasmids,comparSummaryNoPlasmids)
plasmidSummaries <- plasmidSummaries[plasmidSummaries$pass == "raw",]

plasmidSummaries <- plasmidSummaries[order(plasmidSummaries$tool),]

plasmidSummaries$recall <- plasmidSummaries$TP/plasmidSummaries$ntrue
plasmidSummaries$precion <- plasmidSummaries$TP/plasmidSummaries$ncalled

########################################

groupColorsPlasmids <- c("#f768a1", "#c51b8a","#67a9cf", "#1c9099")
names(groupColorsPlasmids)<- unique(plasmidSummaries$group)

plasmis_tp <- ggplot(plasmidSummaries[plasmidSummaries$readLen == "150",], aes(fill = group, x = insLen, y = TP)) +
  geom_col(position=position_dodge(c(0.8)), width = 0.7)  +
  geom_text(aes(label = TP), position = position_dodge(c(0.8)), color = "white", hjust = 1.5) +
  scale_fill_manual("legend", values = groupColorsPlasmids)+
  facet_wrap(~coverage) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.title = element_text(size = 20)) +
  xlab("insert size") + ylab("true positives") +
  ylim(0,max(plasmidSummaries$ntrue))

svg("~/Nextcloud/Publications/2023_CloveBiotech/CloveBiotech-150_TP_Plasmids.svg")

plasmis_tp + coord_flip()

t.test(x = plasmidSummaries$precion[plasmidSummaries$tool == "clovebiotech" & plasmidSummaries$plasmid], 
       y = plasmidSummaries$precion[plasmidSummaries$tool == "clovebiotech" & !plasmidSummaries$plasmid], alternative = "greater")

t.test(x = plasmidSummaries$recall[plasmidSummaries$tool == "clovebiotech" & plasmidSummaries$plasmid], 
       y = plasmidSummaries$recall[plasmidSummaries$tool == "clovebiotech" & !plasmidSummaries$plasmid], alternative = "greater")

wilcox.test(recall~plasmid, data = plasmidSummaries[plasmidSummaries$tool == "clovebiotech",] , exact = FALSE, correct = FALSE, conf.int = FALSE)
wilcox.test(precion~plasmid, data = plasmidSummaries[plasmidSummaries$tool == "clovebiotech",] , exact = FALSE, correct = FALSE, conf.int = FALSE)

dev.off()
