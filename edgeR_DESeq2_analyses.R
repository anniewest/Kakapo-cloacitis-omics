
###Differentially expressed genes across WTS samples###

library(tidyverse)

## Read in wts read counts table
wts.summary.df = read.delim("wts tables/wts_summary_count_table.tsv", sep = '\t', header = T)
names(wts.summary.df)

wts.summary.df.tmm = wts.summary.df[,-c(2:42)]

wts.summary.df.tmm = wts.summary.df.tmm %>%
  `rownames<-`(.$geneID)

colnames(wts.summary.df.tmm) = str_replace(colnames(wts.summary.df.tmm),'.wts_TMM','')

# edgeR significantly DA genes 
DE.genes = c("co_assembly_NODE_88297_length_213_cov_2.006329_1","co_assembly_NODE_9027_length_685_cov_5.193651_1",
             "co_assembly_NODE_601_length_5923_cov_6.513804_7","D2_assembly_NODE_1132_length_600_cov_3.662385_1",
             "D4_assembly_NODE_3516_length_318_cov_1.771863_1")

## DESeq2

library("DESeq2"); packageVersion("DESeq2")

wts.counts.sub = wts.summary.df[,c(1,4:16)] %>% 
  column_to_rownames(var = "geneID")

sample.info <- data.frame(colnames(wts.counts.sub))
names(sample.info) <- c('SampleID')
sample.info$status <- ifelse(sample.info$SampleID %in% c("S12.wts_featurecounts_counts","S13.wts_featurecounts_counts","S8.wts_featurecounts_counts","S9.wts_featurecounts_counts"), 'Healthy', 'Diseased')

# Run DESeq2
dds <- DESeqDataSetFromMatrix(countData = wts.counts.sub,
                              colData = sample.info,
                              design= ~ status)
#dds$status = relevel(dds$status, ref = "Healthy")
dds <- DESeq(dds)
resultsNames(dds) # lists the coefficients
res <- results(dds)
summary(res)
sum(res$padj < 0.1, na.rm=TRUE)
resOrdered <- res[order(res$pvalue),]

# Investigate output
plotCounts(dds, gene=which.min(res$padj), intgroup="status")

res05 <- results(dds, alpha=0.05)
summary(res05)
sum(res05$padj < 0.05, na.rm=TRUE)

res05.df = data.frame(res05) %>% 
  subset(padj < 0.05)

# Check for edgeR genes in DESeq2 output
any(DE.genes %in% rownames(res05.df))
res05.df.subset <- res05.df[rownames(res05.df) %in% DE.genes, ]

## Append annotations
wgs.annot = read.delim("wgs tables/dereplicated_contigs.annotation_summary.aa", sep = '\t', header = T) %>% 
  filter(Query.Gene %in% rownames(res05.df)) %>% 
  as.data.frame()
wgs.annot.subset = wgs.annot %>% select(Query.Gene, Description.2)
names(wgs.annot.subset) = c("geneID","KEGG_desc")

res05.df.annot = res05.df %>%
  rownames_to_column(var = "geneID") %>% 
  left_join(wgs.annot.subset, by = "geneID") %>% 
  mutate(DifferentiallyAbundant = levels(factor(sample.info[,"status"]))[as.numeric(res05.df$log2FoldChange>0)+1])

write.csv(res05.df.annot,"cloacitis_status_de2seq_gene_comparison.csv")

## Hok/gef genes
deseq2.hok.genes = res05.df.annot %>% 
  filter(grepl("K18919|K18920", KEGG_desc)) %>% 
  pull(geneID)

edgeR.hok.genes = c("D2_assembly_NODE_1132_length_600_cov_3.662385_1",
             "D4_assembly_NODE_3516_length_318_cov_1.771863_1")
print(edgeR.hok.genes %in% deseq2.hok.genes)

hok.wgs.genes = read.delim("wgs tables/dereplicated_contigs.annotation_summary.aa", sep = '\t', header = T) %>% 
  filter(str_detect(Description.2, "K18919|K18920")) %>%
  as.data.frame() %>%
  select(Query.Gene,Description.2) 

names(hok.wgs.genes) = c("geneID","KEGG_desc")

hok.genes.df = wts.summary.df.tmm[rownames(wts.summary.df.tmm) %in% hok.wgs.genes$geneID,]

hok.genes.df.plot = pivot_longer(hok.genes.df, 
                                cols = -geneID,  
                                names_to = "Sample",  
                                values_to = "TMM") %>%
  left_join(hok.wgs.genes) %>% 
  separate(KEGG_desc, c("K0_number", "gene_name"), ": ")

hok.genes.df.plot$status <- ifelse(hok.genes.df.plot$Sample %in% c("S12","S13","S8","S9"), 'Healthy', 'Diseased')
hok.genes.df.plot$Sample <- factor(hok.genes.df.plot$Sample, level = c("S1","S2","S3","S4","S5","S6","S7","S8","S9","S10",
                                                                     "S11","S12","S13"))

hok.genes.df.plot.edit <- hok.genes.df.plot %>%
  group_by(gene_name, Sample, status) %>%
  dplyr::summarise(sum_TMM = sum(TMM))

library(ggplot2)
library(hrbrthemes)
library(ggsci)

hok.gene.plot = ggplot(hok.genes.df.plot.edit, aes(Sample, log10(sum_TMM), fill=gene_name)) +
  geom_bar(aes(), stat="identity", position="stack", width = 1, alpha = 0.75)  +  guides(fill=guide_legend(title="Gene", nrow = 2)) +
  theme_ipsum() + theme(legend.position = "bottom", 
                        axis.line.x=element_line(color="black",linewidth=1.0,linetype=1),
                        axis.line.y=element_line(color="black",linewidth=1.0,linetype=1),
                        panel.grid.major = element_blank(),
                        panel.grid.minor = element_blank(),
                        strip.background = element_rect(color="black", fill=NA, linewidth=1, linetype="solid"),
                        strip.text = element_text(size = 10),
                        axis.text.x = element_text(size = 10),
                        axis.text.y = element_text(size = 10)) +
  labs(x="\nK\u101k\u101p\u14d transcriptome sample\n", y="log10 TMM transcript abundance") +  
  facet_grid(.~status, scales = "free_x", space='free') + 
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
  scale_fill_futurama()
hok.gene.plot
ggsave("hok.de.genes.plot.png", hok.gene.plot, 
       dpi = 300, bg = "white", width = 9, height = 5)

pdf("hok.de.genes.plot.pdf", width=10,height=5)
hok.gene.plot
dev.off()
