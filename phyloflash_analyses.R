
###phyloFlash analyses###

library(data.table)
library(dplyr)
library(stringr)

## Read in wts phyloFlash table
phyloflash.df = fread("WTS_phyloFlash_compare.ntu_table.tsv", header=F, sep = "\t") %>%
  data.frame()
  
phyloflash.df$V1 <- str_replace(phyloflash.df$V1,"\"Bacteria", "Bacteria") #str_replace_all(phyloflash.df$V1,';','; ') %>% 
phyloflash.df$V2 <- str_replace(phyloflash.df$V2,'_R1_','') 
head(phyloflash.df)

# Filter Eukaryotes
choanozoa.df <- phyloflash.df %>% 
  filter(grepl("choanozoa", V1, ignore.case = T))

platnts.df <- phyloflash.df %>% 
  filter(grepl("streptophyta", V1, ignore.case = T))

phyloflash.df.noeuk <- phyloflash.df %>% 
  filter(!grepl("Metazoa", V1, ignore.case = T)) %>% 
  filter(!grepl("Streptophyta", V1, ignore.case = T))

phyloflash.df.wide = phyloflash.df.noeuk %>%
  reshape(idvar = "V1", timevar = "V2", direction = "wide") %>%
  `rownames<-`(.$V1) %>%
  select(-V1)

colnames(phyloflash.df.wide) <- str_replace(colnames(phyloflash.df.wide),'V3.','')
head(phyloflash.df.wide)

# Change NA to 0
phyloflash.df.wide[is.na(phyloflash.df.wide)] <- 0

## Create phyloseq object
library(phyloseq)

sample.info <- data.frame(colnames(phyloflash.df.wide)) %>% 
  arrange(.[[1]]) %>%
  rename(SampleID = 1)

sample.info$status <- ifelse(sample.info$SampleID %in% c("Kakapo_S12","Kakapo_S13","Kakapo_S8","Kakapo_S9"), 'Healthy', 'Diseased')
sample.info$name <- c("F-Bravo","S-Uri","S-Hikoi","S-Mukeke","S-Atareta","F-Alice","F-Taeatanga",
                      "F-Cyndy","F-Merv","F-Bella","F-Sinbad","F-Nora","F-Scratch")
sample.info$type <- ifelse(sample.info$SampleID %in% c("Kakapo_S10","Kakapo_S11","Kakapo_S12","Kakapo_S13"), 'Swab', 'Faecal')
rownames(sample.info) = sample.info$SampleID

taxa = phyloflash.df.wide %>%
  rownames() %>%
  as.character() %>%
  strsplit(split = ";") %>%
  data.frame() %>%
  t() %>%
  `colnames<-`(c("Kingdom","Phylum","Order","Class","Family","Genus","Species")) %>%
  `rownames<-`(rownames(phyloflash.df.wide)) %>%
  as.matrix()

ps = phyloseq(otu_table(phyloflash.df.wide, taxa_are_rows = T),
              sample_data(sample.info),
              tax_table(taxa))


rank_names(ps)
ps_tax_table <- as.data.frame(tax_table(ps))

str(ps_tax_table$Phylum)
sapply(ps_tax_table, n_distinct)
levels(factor(ps_tax_table$Kingdom))
count(ps_tax_table, Kingdom)

# Filter chloroplasts, mitochondria, unassigned taxa
ps0 <- subset_taxa(ps, !Phylum %in% c("(Bacteria)","(Eukaryota)","uncultured eukaryote","uncultured")) 
ps0

table(is.na(as.data.frame(tax_table(ps))$Class))
ps0 = subset_taxa(ps0, !Class %in% c("Chloroplast","Rickettsiales"))
ps0
rank_names(ps0)
table(tax_table(ps0)[,"Phylum"], exclude = NULL)

# Filter low-abundance taxa
minTotRelAbun = 1e-5

x = taxa_sums(ps0)
keepTaxa = (x / sum(x)) > minTotRelAbun
ps1 = prune_taxa(keepTaxa, ps0)
ps1

## SRS normalisation
library(SRS)
srs.otu = data.frame(otu_table(ps1))
read.counts = data.frame(colSums(srs.otu)) # Check reads per sample

SRS_output <- SRS(srs.otu, Cmin = 90000) 
rownames(SRS_output) = row.names(srs.otu)
SRS_output = SRS_output[order(rowSums(SRS_output),decreasing = T),]

ps.SRS = phyloseq(otu_table(SRS_output, taxa_are_rows = T),
                  sample_data(sample.info),
                  tax_table(taxa))

## Plot taxa
library(hrbrthemes)
library(ggplot2)
library(plyr)
library(extrafont)

taxaplot_theme <-  theme_ipsum() + theme(plot.title = element_text(size = 35),
                                         axis.text.x = element_text(size = 20),
                                         axis.ticks.x = element_blank(),
                                         axis.text.y=element_text(size=28),
                                         axis.title.y=element_text(size=35),
                                         axis.title.x = element_text(size = 35),
                                         strip.text.x = element_text(size = 26),
                                         panel.grid.major = element_blank(),
                                         panel.grid.minor = element_blank(),
                                         panel.background = element_blank(),
                                         legend.box.background = element_rect(),
                                         legend.box.margin = margin(5, 5, 5, 5),
                                         legend.position="bottom",
                                         legend.title = element_text(face = "bold", size = 35),
                                         legend.text = element_text(size = 30))


# Create 'Others' group for taxa <0.3% mean abundance
others <- transform_sample_counts(ps.SRS, function(x) x/sum(x))
glom_g <- tax_glom(others, taxrank = 'Genus')
df_g <- psmelt(glom_g)
df_g$Genus <- as.character(df_g$Genus)
means_g <- ddply(df_g, ~Genus, function(x) c(mean=mean(x$Abundance)))
remainder_g <- means_g[means_g$mean <= 0.003,]$Genus 
df_g[df_g$Genus %in% remainder_g,]$Genus <- 'Other genera <0.3%' 
print(levels(factor(df_g$Genus)))

tax.colours = c("#FF0F39","#FF99AB","#402751","#a593a5","#666092","#f4de58","#0b5e65","#139eaa","#8fde5d","#0b3165","#c09473","#ecded5","#F6C141",
                    "#e0a10b","#fbff86","#ffd1d5","#429058","#d1d5ff","#CAE0AB","#E8601C","#F9ECA0","#ffb570","#7BAFDE","#4f7190","#787878","#9f9f9f")

genus_taxa_plotbar <- ggplot(data=df_g, aes(x=name, y=Abundance,  
                                             fill = factor(Genus, 
                                                           c("Coccidia","Giardiinae","Fungi","(Actinobacteria)",
                                                             "(Actinomycetaceae)","(Clostridia)",
                                                             "(Corynebacteriaceae)","(Enterobacterales)","(Enterobacteriaceae)",
                                                             "(Firmicutes)","(Lachnospiraceae)","(Pasteurellaceae)","Atopobium","Bacteroides",
                                                             "Clostridium sensu stricto 1","Dickeya","Enterobacter","Enterococcus",
                                                             "Escherichia-Shigella","Klebsiella","NK4A214 group",
                                                             "Salmonella","Streptococcus","Tyzzerella","Varibaculum","Other genera <0.3%")))) +
  geom_bar(aes(), stat="identity", position="fill", width = 1)  + 
  guides(fill=guide_legend(title="Genera", nrow = 6)) + taxaplot_theme + 
  #theme(axis.title.y = element_blank(), axis.text.y = element_blank()) + 
  labs(x="\nK\u101k\u101p\u14d transcriptome sample\n", y="Relative transcript abundance") +  
  facet_grid(.~status, scales = "free_x", space='free') +
  scale_fill_manual(values = tax.colours)
genus_taxa_plotbar
ggsave("phyloflash_genus_taxaplot_Weukaryotes_SRS_JUNE2024.png", genus_taxa_plotbar, width = 25, height = 15, bg = "white", dpi = 300)
ggsave("phyloflash_genus_taxaplot_Weukaryotes_SRS_JUNE2024.pdf", genus_taxa_plotbar, width = 25, height = 15, bg = "white", dpi = 300)

## Alpha diversity

rich_ps <- estimate_richness(ps.SRS, measures = c("InvSimpson","Shannon","Observed"))
richF <- list(rich_ps, sample_data(ps.SRS)$status, sample_data(ps.SRS)$type)
names(richF) <- c("alpha_diversity", "status", "type")

shapiro.test(richF$alpha_diversity$InvSimpson) 
shapiro.test(richF$alpha_diversity$Shannon)
shapiro.test(richF$alpha_diversity$Observed) 

rich.df = data.frame(richF) %>% 
  mutate(sample = str_replace(rownames(.), "Kakapo_", ""))
names(rich.df) = c("Observed","Shannon","InvSimpson","status","type","sample")

library(lme4)
obs_model_null <- glmer(Observed ~ (1|type), data = rich.df, family = poisson(link=log))
obs_model <- glmer(Observed ~ status + (1|type), data = rich.df, family = poisson(link=log)) 
summary(obs_model) 
anova(obs_model, obs_model_null)
confint(obs_model)


rich.df.lmm <- rich.df
rich.df.lmm[,c(2,3)] <- scale(rich.df.lmm[,c(2,3)])

shannon_model_null <- lmer(Shannon ~ (1|type), data = rich.df.lmm)
shannon_model <- lmer(Shannon ~ status + (1|type), data = rich.df.lmm) 
summary(shannon_model) 
anova(shannon_model, shannon_model_null)
confint(shannon_model)

invsimp_model_null <- lmer(InvSimpson ~ (1|type), data = rich.df.lmm)
invsimp_model <- lmer(InvSimpson ~ status + (1|type), data = rich.df.lmm) 
summary(invsimp_model) 
anova(invsimp_model, invsimp_model_null)
confint(invsimp_model)

## Differential taxa (via expression) analysis
library(DESeq2)
library(microbiome)

my_theme <- theme_ipsum() + 
  theme(axis.text.y = element_text(face="italic", size = 15),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(color="black",linewidth=0.7,linetype=1, fill = "NA"), 
        panel.background = element_blank(),
        axis.title.y = element_text(size = 25),
        axis.title.x = element_text(size=20),
        legend.box.background = element_rect(),
        legend.box.margin = margin(5, 5, 5, 5),
        legend.title = element_text(size = 25),
        legend.text = element_text(size = 20), 
        legend.position = "top") 


ps1.deseq = ps1
meta.st <- meta(ps1.deseq)
meta.st$status <- as.factor(meta.st$status)
diagdds_sta = phyloseq_to_deseq2(ps1.deseq, ~ status)
gm_mean = function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}
geoMeans = apply(counts(diagdds_sta), 1, gm_mean)
diagdds_sta = estimateSizeFactors(diagdds_sta, geoMeans = geoMeans)
dds_st = DESeq(diagdds_sta, test="Wald", fitType="local")

otu.ab1 <- abundances(ps1.deseq)
res1 = results(dds_st, cooksCutoff = FALSE)
res_tax1 = cbind(as.data.frame(res1), as.matrix(rownames(otu.ab1)), OTU = rownames(res1))
res_tax1 = cbind(as(res_tax1, "data.frame"), as(tax_table(ps1.deseq)[rownames(res_tax1), ], "matrix"))
res_tax_sig1 = subset(res_tax1, padj < 0.05 & 0 < abs(log2FoldChange))
res_tax1$Significant <- ifelse(rownames(res_tax1) %in% rownames(res_tax_sig1) , "Yes", "No")
res_tax1$Significant[is.na(res_tax1$Significant)] <- "No"
sig_res1 <- res_tax1[rownames(res_tax_sig1),"OTU"]
res_table1 <- data.frame(res_tax_sig1$baseMean , res_tax_sig1$log2FoldChange,res_tax_sig1$padj)
row.names(res_table1) <- rownames(res_tax_sig1)

data_to_write1 <-res_tax_sig1[,c("baseMean","log2FoldChange","pvalue","padj","Phylum", "Class", "Order", "Family", "Genus","Species","OTU")]
data_to_write1 = data_to_write1 %>%
  mutate(Genus = str_replace(Genus,"uncultured","(Pasteurellaceae)")) %>%
  mutate(Genus = str_replace(Genus,"(uncultured bacterium)","(Coriobacteriia)"))

data_to_write1$DifferentiallyAbundant <-levels(meta.st[,"status"])[as.numeric(data_to_write1$log2FoldChange>0)+1]

# Total numer of differentially present OTUs
nrow(data_to_write1) 
length(unique(data_to_write1$Genus)) 
write.csv(data_to_write1,"cloacitis_status_deseq_comparison_0.05.csv")

status.deqseq.plot<- ggplot(data_to_write1, aes(log2FoldChange, Genus)) + 
  geom_point(aes(color = DifferentiallyAbundant), shape = 21, size = 3, stroke = 2) + 
  scale_color_manual(values= c("#4b5e1e", "#7D9D33"), name = "Cloacitis status") + my_theme +
  labs(y = "Differentially active bacterial genera") +  geom_vline(xintercept = 0) 

status.deqseq.plot
ggsave("cloacitis_wts_genera_deseq_plot_0.05.png", status.deqseq.plot, width = 10, height = 10, bg = "white", dpi = 300)
ggsave("cloacitis_wts_genera_deseq_plot_0.05.pdf", status.deqseq.plot, width = 12, height = 12, bg = "white", dpi = 300)


############################################################################################################################################################

##WGS##

## Read in wgs phyloFlash table
phyloflash.wgs.df = fread("WGS_phyloflash_compare.ntu_table.tsv", header=F, sep = "\t") %>%
  data.frame()
head(phyloflash.wgs.df)

# Filter Eukaryotes
choanozoa.wgs.df <- phyloflash.wgs.df %>% 
  filter(grepl("choanozoa", V1, ignore.case = T))

platnts.wgs.df <- phyloflash.wgs.df %>% 
  filter(grepl("streptophyta", V1, ignore.case = T))

phyloflash.wgs.df.noeuk <- phyloflash.wgs.df %>% 
  filter(!grepl("Metazoa", V1, ignore.case = T)) %>% 
  filter(!grepl("Streptophyta", V1, ignore.case = T))

phyloflash.wgs.df.wide = phyloflash.wgs.df.noeuk %>%
  reshape(idvar = "V1", timevar = "V2", direction = "wide") %>%
  `rownames<-`(.$V1) %>%
  select(-V1)

colnames(phyloflash.wgs.df.wide) <- str_replace(colnames(phyloflash.wgs.df.wide),'V3.','')
head(phyloflash.wgs.df.wide)

#change NA to 0
phyloflash.wgs.df.wide[is.na(phyloflash.wgs.df.wide)] <- 0
read.counts.wgs = data.frame(colSums(phyloflash.wgs.df.wide))

## Create a phyloseq object
sample.info.wgs <- data.frame(colnames(phyloflash.wgs.df.wide))
names(sample.info.wgs) <- c('SampleID')
sample.info.wgs$Run <- ifelse(sample.info.wgs$SampleID %in% c("Kakapo_WGS_D7","Kakapo_WGS_D8","Kakapo_WGS_D9"), 'New', 'Old')
sample.info.wgs <- sample.info.wgs[with(sample.info.wgs, order(SampleID)), ]
sample.info.wgs$name <- c("Bonus","Rakiura","Ellie","Tia","Waa","Tutoko","Merv","Taeatanga","Bravo")
sample.info.wgs$age <- ifelse(sample.info.wgs$SampleID %in% c("Kakapo_WGS_D4","Kakapo_WGS_D5","Kakapo_WGS_D6"), 'Chick', 'Adult')

rownames(sample.info.wgs) = sample.info.wgs$SampleID

taxa.wgs = phyloflash.wgs.df.wide %>%
  rownames() %>%
  as.character() %>%
  strsplit(split = ";") %>%
  data.frame() %>%
  t() %>%
  `colnames<-`(c("Kingdom","Phylum","Order","Class","Family","Genus","Species")) %>%
  `rownames<-`(rownames(phyloflash.wgs.df.wide)) %>%
  as.matrix()

ps.wgs = phyloseq(otu_table(phyloflash.wgs.df.wide, taxa_are_rows = T),
              sample_data(sample.info.wgs),
              tax_table(taxa.wgs))
ps.wgs

rank_names(ps.wgs)
ps_wgs_tax_table <- as.data.frame(tax_table(ps.wgs))

str(ps_wgs_tax_table$Phylum)
sapply(ps_wgs_tax_table, n_distinct)
levels(factor(ps_wgs_tax_table$Kingdom))

# Filter chloroplasts, mitochondria, unclassified taxa
ps0.wgs <- subset_taxa(ps.wgs, !Phylum %in% c("(Bacteria)","(Eukaryota)","uncultured eukaryote","uncultured"))

phyla2Filter = c("Chloroplast","Rickettsiales")
ps0.wgs = subset_taxa(ps0.wgs, !Class %in% phyla2Filter)
ps0.wgs
rank_names(ps0.wgs)
table(tax_table(ps0.wgs)[,"Phylum"], exclude = NULL)

# Filter low-abundance taxa
x = taxa_sums(ps0.wgs)
keepTaxa = (x / sum(x)) > minTotRelAbun  
ps1.wgs = prune_taxa(keepTaxa, ps0.wgs)
ps1.wgs

## SRS normalisation
srs.otu.wgs = data.frame(otu_table(ps1.wgs))
read.counts = data.frame(colSums(srs.otu.wgs)) # Check reads per sample

SRS_output.wgs <- SRS(srs.otu.wgs, Cmin = 490)
rownames(SRS_output.wgs) = row.names(srs.otu.wgs)
SRS_output.wgs = SRS_output.wgs[order(rowSums(SRS_output.wgs),decreasing = T),]

ps.SRS.wgs = phyloseq(otu_table(SRS_output.wgs, taxa_are_rows = T),
                  sample_data(sample.info.wgs),
                  tax_table(taxa.wgs))
ps.SRS.wgs

## Plot taxa
others.wgs <- transform_sample_counts(ps.SRS.wgs, function(x) x/sum(x))
glom_g.wgs <- tax_glom(others.wgs, taxrank = 'Genus')
df_g.wgs <- psmelt(glom_g.wgs)
df_g.wgs$Genus <- as.character(df_g.wgs$Genus)
means_g.wgs <- ddply(df_g.wgs, ~Genus, function(x) c(mean=mean(x$Abundance)))
remainder_g.wgs <- means_g.wgs[means_g.wgs$mean <= 0.003,]$Genus 
df_g.wgs[df_g.wgs$Genus %in% remainder_g.wgs,]$Genus <- 'Other species <0.3%' 
print(levels(factor(df_g.wgs$Genus)))

tax.colours.wgs = c("#4646ff","#139eaa","#8fde5d","#d6e272", "#ff858f","#fbff86","#ffd1d5","#429058",
                    "#CAE0AB","#666092","#E8601C","#cc0000","#ffb570","#f4de58","#546e93","#6666e0","#7BAFDE","#9a4d76","#9f9f9f")

genus_taxa_plotbar.wgs <- ggplot(data=df_g.wgs, aes(x=name, y=Abundance,  
                                            fill = factor(Genus, 
                                                          c("(Bacilli)","(Enterobacterales)","(Enterobacteriaceae)",
                                                            "Bacillus","Caenimonas","Clostridium sensu stricto 1",
                                                            "Dickeya","Enterobacter","Escherichia-Shigella","Hymenobacter",
                                                            "Klebsiella","Nostoc PCC-73102","Salmonella","Serratia","Sphingomonas",
                                                            "Spirosoma","Streptococcus","Vibrio","Other species <0.3%")))) +
  geom_bar(aes(), stat="identity", position="fill", width = 1)  + 
  guides(fill=guide_legend(title="Genera", nrow = 5)) + taxaplot_theme + 
  labs(x="\nK\u101k\u101p\u14d metagenome sample\n", y="Relative sequence abundance") +  
  facet_grid(.~Run, scales = "free_x", space='free') +
  scale_fill_manual(values = tax.colours.wgs)
genus_taxa_plotbar.wgs
ggsave("phyloflash_wgs_genus_taxaplot_Weukaryotes_SRS.png", genus_taxa_plotbar.wgs, width = 25, height = 15, dpi = 300)
ggsave("phyloflash_wgs_genus_taxaplot_Weukaryotes_SRS.pdf", genus_taxa_plotbar.wgs, width = 25, height = 15, dpi = 300)

## Alpha diversity

rich_ps.wgs <- estimate_richness(ps.SRS.wgs, measures = c("InvSimpson","Shannon","Observed"))
richF.wgs <- list(rich_ps.wgs, sample_data(ps.SRS.wgs)$Run, sample_data(ps.SRS.wgs)$age, sample_data(ps.SRS.wgs)$name)
names(richF.wgs) <- c("alpha_diversity", "Run", "age","name")

shapiro.test(richF.wgs$alpha_diversity$InvSimpson) 
shapiro.test(richF.wgs$alpha_diversity$Shannon)
shapiro.test(richF.wgs$alpha_diversity$Observed) 

rich.df.wgs = data.frame(richF.wgs)
names(rich.df.wgs) = c("Observed","Shannon","InvSimpson","Run","age","name")

obs_model_null.wgs <- glmer(Observed ~ (1|age), data = rich.df.wgs, family = poisson(link=log))
obs_model.wgs <- glmer(Observed ~ Run + (1|age), data = rich.df.wgs, family = poisson(link=log)) 
summary(obs_model.wgs) 
anova(obs_model.wgs, obs_model_null.wgs)
confint(obs_model.wgs)


rich.df.lmm.wgs <- rich.df.wgs
rich.df.lmm.wgs[,c(2,3)] <- scale(rich.df.lmm.wgs[,c(2,3)])

shannon_model_null <- lmer(Shannon ~ (1|age), data = rich.df.lmm.wgs)
shannon_model <- lmer(Shannon ~ Run + (1|age), data = rich.df.lmm.wgs) 
summary(shannon_model) 
anova(shannon_model, shannon_model_null)
confint(shannon_model)

invsimp_model_null <- lmer(InvSimpson ~ (1|age), data = rich.df.lmm.wgs)
invsimp_model <- lmer(InvSimpson ~ Run + (1|age), data = rich.df.lmm.wgs) 
summary(invsimp_model) 
anova(invsimp_model, invsimp_model_null)
confint(invsimp_model)

obs.wgs <- ggplot(rich.df.wgs, aes(x=name, y=Observed, fill=Run))
obs.wgs.plot <- obs.wgs + geom_bar(aes(), stat="identity", position="stack", width = 1) + #geom_boxplot(size=1.0, aes(fill=Run)) + 
  myadptheme + theme(plot.title = element_text(size = 35)) +
  labs(y = "Observed species", x = "\nSequencing run") +
  ggtitle("glmm p < 0.001***") +
  facet_grid(~Run, scales = "free_x", space = "free") +
  scale_fill_manual(values = c("#7D9D33", "#CED38C"))
obs.wgs.plot
ggsave("adp_obs_wgs.png", obs.wgs.plot, width = 12, height = 10, bg = "white", dpi = 300)

## Differential abundance analysis
library(DESeq2)
library(microbiome)

ps1.wgs.bac = subset_taxa(ps1.wgs, Kingdom=="Bacteria")

meta.wgs <- meta(ps1.wgs.bac)
meta.wgs$Run <- as.factor(meta.wgs$Run)
diagdds_wgs = phyloseq_to_deseq2(ps1.wgs.bac, ~ Run)
gm_mean = function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}
geoMeans.wgs = apply(counts(diagdds_wgs), 1, gm_mean)
diagdds_wgs = estimateSizeFactors(diagdds_wgs, geoMeans = geoMeans.wgs)
dds_wgs = DESeq(diagdds_wgs, test="Wald", fitType="local")

otu.ab1.wgs <- abundances(ps1.wgs.bac)
res1.wgs = results(dds_wgs, cooksCutoff = FALSE)
res_tax1.wgs = cbind(as.data.frame(res1.wgs), as.matrix(rownames(otu.ab1.wgs)), OTU = rownames(res1.wgs))
res_tax1.wgs = cbind(as(res_tax1.wgs, "data.frame"), as(tax_table(ps1.wgs.bac)[rownames(res_tax1.wgs), ], "matrix"))
res_tax_sig1.wgs = subset(res_tax1.wgs, padj < 0.01 & 0 < abs(log2FoldChange))
res_tax1.wgs$Significant <- ifelse(rownames(res_tax1.wgs) %in% rownames(res_tax_sig1.wgs) , "Yes", "No")
res_tax1.wgs$Significant[is.na(res_tax1.wgs$Significant)] <- "No"
sig_res1.wgs <- res_tax1.wgs[rownames(res_tax_sig1.wgs),"OTU"]
res_table1.wgs <- data.frame(res_tax_sig1.wgs$baseMean , res_tax_sig1.wgs$log2FoldChange,res_tax_sig1.wgs$padj)
row.names(res_table1.wgs) <- rownames(res_tax_sig1.wgs)

data_to_write1.wgs <-res_tax_sig1.wgs[,c("baseMean","log2FoldChange","pvalue","padj","Phylum", "Class", "Order", "Family", "Genus","Species","OTU")]
data_to_write1.wgs$DifferentiallyAbundant <-levels(meta.wgs[,"Run"])[as.numeric(data_to_write1.wgs$log2FoldChange>0)+1]

## Total number of differentially abundant OTUs 
nrow(data_to_write1.wgs) 
length(unique(data_to_write1.wgs$Genus)) 
write.csv(data_to_write1.wgs,"wgs_cloacitis_status_deseq_comparison.csv")

wgs.deqseq.plot<- ggplot(data_to_write1.wgs, aes(log2FoldChange, Genus)) + 
  geom_point(aes(color = DifferentiallyAbundant), shape = 21, size = 3, stroke = 2) + 
  scale_color_manual(values= c("#7BAFDE", "#7D9D33"), name = "Sequencing run") + my_theme +
  labs(y = "Differentially abundant genera") +  geom_vline(xintercept = 0) 

wgs.deqseq.plot
ggsave("cloacitis_wgs_bac_deseq_plot.png", wgs.deqseq.plot, width = 10, height = 8, bg = "white", dpi = 300)
