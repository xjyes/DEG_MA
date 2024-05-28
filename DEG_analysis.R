setwd("D:\\STUDY\\szbl\\MA_HA_CO_vs._AL\\code")

# Load all data as data
data <- read.csv("D:\\STUDY\\szbl\\MA_HA_CO_vs._AL\\data\\all_compare_MA_Co_vs_Al.csv")
colnames(data)[1] = "gene_id"
head(data)

# Delete all data with no expression
data = na.omit(data)
save(data, file = "data_MA.Rdata")

# Extract expression matrix
exp = data[,-c(8:27)]
rownames(exp) = exp[,1]
exp = exp[,-1]
head(exp)

# Change rownames to gene names
tmp = by(exp,
         data$gene_name,
         function(x) rownames(x)[which.max(rowMeans(x))])
probes = as.character(tmp)
dim(exp)
exp = exp[rownames(exp) %in% probes, ]
dim(exp)
rownames(exp) = data[match(rownames(exp),data$gene_id),19]
head(exp)
colnames(exp)[1:6] = c("Co_1","Co_2","Co_3","Al_1","Al_2","Al_3")
exp = exp[,c(4,5,6,1,2,3)]
head(exp)
save(exp, file = "exp_MA.Rdata")

# Prepare matrices for DEG analysis
rm(list = ls())
options(stringsAsFactors = F)
load("exp_HA.Rdata")
library(limma)
group <- factor(rep(c("Al","Co"),each=3))
design <- model.matrix(~0+group)
colnames(design)=levels(group)
rownames(design)=colnames(exp)
contrast.matrix <- makeContrasts(paste0(c("Co","Al"),collapse = "-"),levels = design)
print(design)
print(contrast.matrix)
library(edgeR)
y <- DGEList(counts = exp)
head(y)

# Filter genes with low counts
keep <- rowSums(cpm(y)>1) >= 2
y <- y[keep, , keep.lib.sizes=FALSE]
head(y)
y <- calcNormFactors(y)
head(y)
logCPM <- cpm(y, log=TRUE, prior.count=3)
head(logCPM)
save(logCPM, file = "logCPM_MA.Rdata")

# Perform PCA to check grouping
library(ggfortify)
df = as.data.frame((t(exp)))
df$group = group
autoplot(prcomp(df[,1:(ncol(df)-1)]), data = df, colour = 'group',
         label = TRUE, label.size = 3
         ) 

# DEG analysis
fit <- lmFit(logCPM, design)
fit2 <- contrasts.fit(fit,contrast.matrix)
fit2 <- eBayes(fit2, trend=TRUE)
tempOutput <- topTable(fit2,coef=1,n=Inf)
nrDEG = na.omit(tempOutput)
head(nrDEG)
save(nrDEG, file = "nrDEG_MA.Rdata")


# Plot heatmap
library(pheatmap)
rm(list = ls())
options(stringsAsFactors = F)
load(file = "nrDEG_HA.Rdata")
load(file = "exp_HA.Rdata")
choose_gene1 = rownames(nrDEG)
choose_matrix1 = exp[choose_gene1,]
matrix_scale1 = t(scale(t(choose_matrix1)))
annotation_col = data.frame(sampletype = c(rep("Al",3),rep("Co",3)))
row.names(annotation_col) <- colnames(choose_matrix1)
pheatmap(matrix_scale1, scale = "row", annotation_col = annotation_col,
         cluster_rows = T, cluster_cols = T,show_rownames = F,border_color = NA)

choose_gene2 = c(head(rownames(nrDEG)[which(nrDEG$logFC < 0, arr.ind = T)],20),
                 head(rownames(nrDEG)[which(nrDEG$logFC > 0, arr.ind = T)],20))
choose_matrix2 = exp[choose_gene2,]
matrix_scale2 = t(scale(t(choose_matrix2)))
pheatmap(matrix_scale2, scale = "row",annotation_col = annotation_col,
         cluster_cols = F ,border_color = NA)

choose_gene3 = head(rownames(nrDEG),40)
choose_matrix3 = exp[choose_gene3,]
matrix_scale3 = t(scale(t(choose_matrix3)))
pheatmap(matrix_scale3, scale ="row", annotation_col = annotation_col,border_color = NA)

# Export the top 20 UP and DOWN regulated genes
load("data_HA.Rdata")
top20 = data[match(rownames(choose_matrix2),data$gene_name),]
top20$label[1:20] =  "Al"
top20$label[21:40] = "Co"
save(top20,file = "top20_HA.Rdata")
write.csv(top20, file = "top20_HA_UPandDOWN.csv",row.names = F)

# Plot Volcano map
library(ggplot2)
library(ggrepel)
rm(list = ls())
options(stringsAsFactors = F)
load(file = "nrDEG_MA.Rdata")
plot(nrDEG$logFC,-log10(nrDEG$P.Value))
DEG = nrDEG
logFC_cutoff <- with(DEG, mean(abs(logFC),na.rm = TRUE) + 2*sd(abs(logFC),na.rm = TRUE))
DEG$change = as.factor(ifelse(DEG$P.Value < 0.05 & abs(DEG$logFC) > logFC_cutoff,
                              ifelse(DEG$logFC >logFC_cutoff, 'UP', 'DOWN'),'NOT')
)
this_tile <- paste0('Cutoff for logFC is ',round(logFC_cutoff,3),
                    '\nThe number of up gene is ',nrow(DEG[DEG$change =='UP',]) ,
                    '\nThe number of down gene is ',nrow(DEG[DEG$change =='DOWN',]),
                    '\nThe number of not classified gene is ', nrow(DEG[DEG$change == 'NOT',]))

this_tile
head(DEG)
save(DEG, file = "DEG_MA.Rdata")

g = ggplot(data=DEG, aes(x=logFC, y=-log10(P.Value), color=change)) +
  geom_point(alpha=0.4, size=1.75) +
  theme_set(theme_set(theme_bw(base_size=20)))+
  xlab("log2 fold change") + ylab("-log10 p-value") +
  ggtitle( this_tile  ) + theme(plot.title = element_text(size=10,hjust = 0.5))+
  scale_colour_manual(values = c('blue','black','red'))
print(g)



# Plot multiple DEG
library(reshape2)
library(stringr)
rm(list = ls())
options(stringsAsFactors = F)
load("top20_MA.Rdata")
load("logCPM_MA.Rdata")
load("exp_MA.Rdata")
ttop = top20
ttop[,2:7] = logCPM[match(ttop$gene_name,rownames(logCPM)),1:6]
ttop = ttop[,-c(1,8:18,20:27)]
colnames(ttop)[1:6] = c("Al_1","Al_2","Al_3","Co_1","Co_2","Co_3")
ttop20 = melt(ttop[1:20,])
colnames(ttop20)[2:4] = c("group","sample","value")
ttop20$group = ifelse(str_detect(ttop20$sample,"Co")==TRUE,"Co","Al")
ttop20$gene_name <- factor(ttop20$gene_name, levels = ttop20$gene_name[1:20])

ggplot(ttop20) +
  geom_point(aes(x = gene_name, y = value, color = group, shape = group),stat = "identity") +
  xlab("Genes") +
  ylab("Normalized expression level") +
  ggtitle("Top 20 significant DE genes in Al") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(plot.title = element_text(hjust = 0.5))

ttop40 = melt(ttop[21:40,],value.names=rownames(ttop))
colnames(ttop40)[2:4] = c("group","sample","value")
ttop40$group = ifelse(str_detect(ttop40$sample,"Co")==TRUE,"Co","Al")
ttop40$gene_name <- factor(ttop40$gene_name, levels = ttop40$gene_name[1:20])

ggplot(ttop40) +
  geom_point(aes(x = gene_name, y = value, color = group, shape = group),stat = "identity") +
  #scale_y_log10() +
  xlab("Genes") +
  ylab("Normalized expression level") +
  ggtitle("Top 20 significant DE genes in Co") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(plot.title = element_text(hjust = 0.5))


# # Plot Venn graph of common genes
# library(VennDiagram)
# rm(list = ls())
# options(stringsAsFactors = F)
# load("exp_HA.Rdata")
# exp_HA = exp
# load("exp_MA.Rdata")
# exp_MA = exp
# venn_list <- list(group1 = rownames(exp_HA),group2 = rownames(exp_MA))
# venn.diagram(venn_list, fill = c("red","blue"), filename = "Venn common genes.png", 
#              imagetype = 'png', category = c("HA","MA"),alpha =0.5, cat.col = rep('black', 2),
#              col ='black', cex = 1.5, fontfamily = 'serif',
#               ext.pos = 15, ext.line.lty = "dashed", ext.length = 0.8, ext.dist = -0.1,
#              cat.dist = 0.05, cat.just = list(c(-1,-1),c(1,-1)),cat.cex = 1.5, cat.fontfamily = 'serif')
# inter <- get.venn.partitions(venn_list)
# for (i in 1:nrow(inter)) inter[i,'values'] <- paste(inter[[i,'..values..']], collapse = ', ')
# write.csv(inter[,-c(3,4)], 'venn_common_gene.csv', row.names = F, quote = F)



# Enrichment analysis: data preparation
# add the symbol column
rm(list = ls())
options(stringsAsFactors = F)
load("nrDEG_HA.Rdata")
library(dplyr)
deg = nrDEG
deg = mutate(deg,symbol = rownames(deg))
head(deg)

# add the change column: up and down regulation
logFC_t <- with(deg, mean(abs(logFC),na.rm = TRUE) + 2*sd(abs(logFC),na.rm = TRUE))
change = ifelse(deg$logFC > logFC_t,'up',
                ifelse(deg$logFC < -logFC_t,'down','stable'))
deg <- mutate(deg, change)
head(deg)
table(deg$change)

# add the ENTREZID column: map gene symbols to ENTREZID
library(ggplot2)
library(R.utils)
library(clusterProfiler)
library(org.Mm.eg.db) # Load mouse gene annotation
library(org.Hs.eg.db) # Load human gene annotation
R.utils::setOption("clusterProfiler.download.method","auto")
s2e <- bitr(unique(deg$symbol),fromType = "SYMBOL",
            toType = c("ENTREZID"),
            OrgDb = org.Hs.eg.db)
head(s2e)
head(deg)
deg <- inner_join(deg,s2e,by=c("symbol" = "SYMBOL"))

head(deg)
save(deg, file = "deg_HA.Rdata")


# Enrichment analysis #

# KEGG pathway analysis
rm(list = ls())
options(stringsAsFactors = F)
library(R.utils)
library(clusterProfiler)
library(ggplot2)
R.utils::setOption("clusterProfiler.download.method","auto")
load("deg_MA.Rdata")
gene_up = deg[deg$change == 'up','ENTREZID']
gene_down = deg[deg$change == 'down', 'ENTREZID']
gene_diff = c(gene_up, gene_down)
gene_all = deg[,'ENTREZID']
kk.up <- enrichKEGG(gene = gene_up,
                    organism = 'mmu',
                    universe = gene_all,
                    pvalueCutoff = 0.9,
                    qvalueCutoff = 0.9)
head(kk.up)[,1:6]
dim(kk.up)
kk.down <- enrichKEGG(gene = gene_down,
                      organism = 'mmu',
                      universe = gene_all,
                      pvalueCutoff = 0.9,
                      qvalueCutoff = 0.9)
head(kk.down)[,1:6]
dim(kk.down)
kk.diff <- enrichKEGG(gene = gene_diff,
                      organism = 'mmu',
                      pvalueCutoff = 0.05)
head(kk.diff)[,1:6]
class(kk.diff)
# Change entrez id to gene symbol
kk.up =DOSE::setReadable(kk.up, OrgDb='org.Mm.eg.db',keyType='ENTREZID')
kk.down =DOSE::setReadable(kk.down, OrgDb='org.Mm.eg.db',keyType='ENTREZID')
kk.diff =DOSE::setReadable(kk.diff, OrgDb='org.Mm.eg.db',keyType='ENTREZID')

# extract enrichment result data frame
kegg_diff_dt <- kk.diff@result
# select by pvalue
down_kegg <- kk.down@result %>%
  filter(pvalue<0.05) %>%
  mutate(group=-1)
up_kegg <- kk.up@result %>%
  filter(pvalue<0.05) %>%
  mutate(group=1)


# Visualization
kegg_plot <- function(up_kegg, down_kegg){
  dat=rbind(up_kegg,down_kegg)
  colnames(dat)
  dat$pvalue = -log10(dat$pvalue)
  dat$pvalue = dat$pvalue*dat$group
  
  dat = dat[order(dat$pvalue, decreasing = F),]
  
  g_kegg <- ggplot(dat, aes(x=reorder(Description,order(pvalue,decreasing = F)),y = pvalue, fill = group)) +
    geom_bar(stat="identity", width = 0.8) +
    scale_fill_gradient(low="blue",high="red",guide = FALSE) +
    scale_x_discrete(name = "Pathway names") +
    scale_y_continuous(name = "log10P-value")+
    coord_flip() + theme_bw() + theme(plot.title = element_text(hjust = 0.5))+
    ggtitle("Pathway Enrichment")
}

g_kegg <- kegg_plot(up_kegg, down_kegg)
g_kegg


# Get enriched genes and the logFC of KEGG
library(dplyr)
library(tidyr)
library(org.Mm.eg.db) # Load mouse gene annotation
library(org.Hs.eg.db) # Load human gene annotation

# up regulated gene
gene_ID_up <- data.frame("gene_id"=up_kegg$geneID)
gene_ID_up <- gene_ID_up %>% separate(gene_id, sep = "/", into = as.character(1:max(up_kegg$Count)))
gene_ID_up["logFC"] = deg[match(gene_ID_up[,1],deg$symbol),"logFC"]
gene_ID_up["logFC"] = substr(gene_ID_up$logFC,1,4)
gene_ID_up = tidyr::unite(gene_ID_up,"symbol_logFC",as.character(1),logFC,sep = "_")

for (i in 2:max(up_kegg$Count)){
  gene_ID_up["logFC"] = deg[match(gene_ID_up[,2],deg$symbol),"logFC"]
  gene_ID_up["logFC"] = substr(gene_ID_up$logFC,1,4)
  gene_ID_up = tidyr::unite(gene_ID_up,"temp",2,logFC,sep = "_",na.rm = T)
  gene_ID_up = tidyr::unite(gene_ID_up,"symbol_logFC",symbol_logFC,temp,sep = "/", na.rm = T)
}
gene_ID_up["symbol_logFC"] = gsub("/*$","",gene_ID_up$symbol_logFC)
up_kegg["symbol_logFC"] = gene_ID_up$symbol_logFC


# down regulated gene
gene_ID_down <- data.frame("gene_id"=down_kegg$geneID)
gene_ID_down <- gene_ID_down %>% separate(gene_id, sep = "/", into = as.character(1:max(down_kegg$Count)))
gene_ID_down["logFC"] = deg[match(gene_ID_down[,1],deg$symbol),"logFC"]
gene_ID_down["logFC"] = substr(gene_ID_down$logFC,1,4)
gene_ID_down = tidyr::unite(gene_ID_down,"symbol_logFC",as.character(1),logFC,sep = "_")

for (i in 2:max(down_kegg$Count)){
  gene_ID_down["logFC"] = deg[match(gene_ID_down[,2],deg$symbol),"logFC"]
  gene_ID_down["logFC"] = substr(gene_ID_down$logFC,1,4)
  gene_ID_down = tidyr::unite(gene_ID_down,"temp",2,logFC,sep = "_",na.rm = T)
  gene_ID_down = tidyr::unite(gene_ID_down,"symbol_logFC",symbol_logFC,temp,sep = "/", na.rm = T)
}
gene_ID_down["symbol_logFC"] = gsub("/*$","",gene_ID_down$symbol_logFC)
down_kegg["symbol_logFC"] = gene_ID_down$symbol_logFC

# Generate the kegg enrichment analysis result
total_kegg = rbind(up_kegg,down_kegg)

total_kegg$group[which(total_kegg$group == '-1')] <- "DOWN"
total_kegg$group[which(total_kegg$group == "1")] <- "UP"
write.table(total_kegg,"kegg_MA.csv",row.names=FALSE,col.names=TRUE,sep=",")

# GO database analysis #
rm(list = ls())
library(ggplot2)
library(clusterProfiler)
library(org.Mm.eg.db) #load mouse gene annotation
library(org.Hs.eg.db) # load human gene annotation
options(stringsAsFactors = F)
load("deg_HA.Rdata")

gene_up = deg[deg$change == 'up', 'ENTREZID']
gene_down = deg[deg$change == 'down', 'ENTREZID']
gene_diff = c(gene_up,gene_down)
head(deg)

# Cellular component
ego_CC <- enrichGO(gene = gene_diff,
                   OrgDb = org.Hs.eg.db,
                   ont = "CC",
                   pAdjustMethod = "BH",
                   minGSSize = 1,
                   pvalueCutoff = 0.01,
                   qvalueCutoff = 0.01,
                   readable = TRUE)



# Biological process
ego_BP <- enrichGO(gene = gene_diff,
                   OrgDb = org.Hs.eg.db,
                   ont = "BP",
                   pAdjustMethod = "BH",
                   minGSSize = 1,
                   pvalueCutoff = 0.01,
                   qvalueCutoff = 0.01,
                   readable = TRUE)


# Molecular function
ego_MF <- enrichGO(gene = gene_diff,
                   OrgDb = org.Hs.eg.db,
                   ont = "MF",
                   pAdjustMethod = "BH",
                   minGSSize = 1,
                   pvalueCutoff = 0.01,
                   qvalueCutoff = 0.01,
                   readable = TRUE)


# All
ego_ALL <- enrichGO(gene = gene_diff,
                   OrgDb = org.Hs.eg.db,
                   ont = "ALL",
                   pAdjustMethod = "BH",
                   minGSSize = 1,
                   pvalueCutoff = 0.01,
                   qvalueCutoff = 0.01,
                   readable = TRUE)

# Plot bar
barplot(ego_CC, showCategory=20, font.size = 8, title = "EnrichmentGO_CC")
barplot(ego_BP, showCategory=20, font.size = 8, title = "EnrichmentGO_BP")
barplot(ego_MF, showCategory=20, font.size = 8, title = "EnrichmentGO_MF")
barplot(ego_ALL, font.size = 8, title = "EnrichmentGO_ALL",split = "ONTOLOGY") + 
  facet_grid(ONTOLOGY~.,scale="free")
# Plot dot
dotplot(ego_CC, showCategory=20, font.size = 8, title="EnrichmentGO_CC")
dotplot(ego_BP, showCategory=20, font.size = 8, title="EnrichmentGO_BP")
dotplot(ego_MF, showCategory=20, font.size = 8, title="EnrichmentGO_MF")
dotplot(ego_ALL, showCategory=8, font.size = 8, title = "EnrichmentGO_ALL",split = "ONTOLOGY") + 
  facet_grid(ONTOLOGY~.,scale="free")

# Get enriched genes and logFC in GO
library(dplyr)
library(tidyr)
GO = ego_ALL@result
GO_ID <- data.frame("gene_id"=GO$geneID)
GO_ID <- GO_ID %>% separate(gene_id, sep = "/", into = as.character(1:max(GO$Count)))
GO_ID["logFC"] = deg[match(GO_ID[,1],deg$symbol),"logFC"]
GO_ID["logFC"] = substr(GO_ID$logFC,1,4)
GO_ID = tidyr::unite(GO_ID,"symbol_logFC",as.character(1),logFC,sep = "_")

for (i in 2:max(GO$Count)){
  GO_ID["logFC"] = deg[match(GO_ID[,2],deg$symbol),"logFC"]
  GO_ID["logFC"] = substr(GO_ID$logFC,1,4)
  GO_ID = tidyr::unite(GO_ID,"temp",2,logFC,sep = "_",na.rm = T)
  GO_ID = tidyr::unite(GO_ID,"symbol_logFC",symbol_logFC,temp,sep = "/", na.rm = T)
}
GO_ID["symbol_logFC"] = gsub("/*$","",GO_ID$symbol_logFC)
GO["symbol_logFC"] = GO_ID$symbol_logFC

write.table(GO,"GO_HA.csv",row.names=FALSE,col.names=TRUE,sep=",")











