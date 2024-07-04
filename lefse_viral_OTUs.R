
###LEfSe###

library(tidyverse)
library(data.table)

## Read in vOTU summary counts table with checkV data
viral_summary_checkV <- read.delim("vOTUs_filtered_summary_table_VIRUSES_checkv_counts.tsv", header=T, sep='\t')

# Check if contigIG contains vOTU string and subset # 922756 --> 660899
votu.table <- viral_summary_checkV[viral_summary_checkV$contig_ID  %like% "vOTU", ]
# Subset table by vOTUs with counts >1 # 660889 --> 7915 (non dereplicated)
votu.table.sub <- votu.table[which(votu.table$checkv_viral_genes > 1),]

votu.sub.tpm <- votu.table[, -c(2:88,102:114)]

votu.sub.tpm.dup <- votu.sub.tpm %>%
  filter(!duplicated(.)) %>%
  column_to_rownames("contig_ID")

votu.sub.tpm.dup.sort = votu.sub.tpm.dup[order(rowSums(votu.sub.tpm.dup), decreasing = TRUE),]
write.table(votu.sub.tpm.dup.sort, "viral.tpm.summary.table.sorted.tsv", sep = '\t')

## Check if there are any vOTUs present in all samples
nrow(votu.sub.tpm.dup.sort)
allOK = apply(votu.sub.tpm.dup.sort[1:6699,] >= 1, 1, all) 
data2 = votu.sub.tpm.dup.sort[allOK, ]
common.votus = rownames(data2)

common.votu.data = votu.table %>% 
  filter(contig_ID == common.votus[1] | contig_ID == common.votus[2] | contig_ID == common.votus[3])
write.csv(common.votu.data, "common.votus.csv")

## Subset vOTU table by those with overall counts >0
votu.sub.tpm.dup.nonzero = votu.sub.tpm.dup[which(rowSums(votu.sub.tpm.dup) > 0),]
sample.info <- data.frame(colnames(votu.sub.tpm.dup.nonzero))
names(sample.info) <- c('SampleID')
sample.info$status <- ifelse(sample.info$SampleID %in% c("counts_S12.wts_TPM","counts_S13.wts_TPM","counts_S8.wts_TPM","counts_S9.wts_TPM"), 'Healthy', 'Affected')

## LEfSe analysis
library(SummarizedExperiment)

se.df <- SummarizedExperiment(assays =  list(votu.sub.tpm.dup.nonzero),
                              colData = sample.info)

library(lefser)
library(hrbrthemes)

lfse <- lefser(se.df, kruskal.threshold = 0.05, wilcox.threshold = 0.05, groupCol = "status") 

lefse.plot <- lefserPlot(lfse) + theme_classic() + 
  scale_fill_manual(labels = c("Affected","Healthy"), values = c("#4b5e1e", "#7D9D33"), name = "Cloacitis status") + 
  labs(x = "Viral OTUs") +
  ggtitle("Differentially abundant viral OTUs")

pdf("lefse_viral_OTUs.pdf", 
    width = 5, height = 6)
lefse.plot
dev.off()

ggsave("lefse_viral_OTUs.png", lefse.plot, 
       bg = "white", dpi = 300, width = 5, height = 6)

viral.summary.sub <- filter(votu.table, votu.table$contig_ID %in% lfse$Names)
write.table(viral.summary.sub, "viral_summary_lefse_subset.tsv", sep = '\t')


##Subsetting
votu.sub.tpm.tp = data.frame(t(votu.sub.tpm.dup))
votu.sub.tpm.tp$Status = ifelse(rownames(votu.sub.tpm.tp) %in% c("counts_S12.wts_TPM","counts_S13.wts_TPM","counts_S8.wts_TPM","counts_S9.wts_TPM"), 'Healthy', 'Affected')

# vOTUs found in healthy birds
votu.healthy = data.frame(votu.sub.tpm.tp[which(votu.sub.tpm.tp$Status=='Healthy'),])
votu.healthy$Status = NULL
votu.sub.tpm.tp.healthy.zero = votu.sub.tpm.tp[,which(colSums(votu.healthy) < 1)]
# vOTUs not present in any healthy birds
votus.not.healthy = names(votu.sub.tpm.tp.healthy.zero)

# vOTUs found in affected birds (and only affected birds)
votu.affected = data.frame(votu.sub.tpm.tp[which(votu.sub.tpm.tp$Status=='Affected'),])
votu.affected$Status = NULL 
votu.affected.not.healthy.votus = data.frame(votu.affected[,which(colnames(votu.affected) %in% votus.not.healthy)])
votu.affected.not.healthy.votus.nonzero = votu.affected.not.healthy.votus[, which(colSums(votu.affected.not.healthy.votus) > 0)]

votu.affected.nonzero.sort <- data.frame(votu.affected.not.healthy.votus.nonzero[,order(colSums(votu.affected.not.healthy.votus.nonzero), decreasing = TRUE)])

write.table(votu.affected.nonzero.sort, "votu.affected.not.healthy.nonzero.sort.tsv", sep = '\t')

viral.summary.affected.sub <- filter(votu.table, votu.table$contig_ID %in% colnames(votu.affected.nonzero.sort))
write.table(viral.summary.affected.sub, "viral.summary.affected.sub.tsv", sep = '\t')
