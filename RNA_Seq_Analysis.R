# TFM MADOBIS
#Alumno: Sergio Bustamante

#Analisis Rna-seq
#BiocManager::install("DESeq2")

library(DESeq2)
library(dplyr)

#Cargar datos
pheno.data <- read.csv("pheno_data2.csv")
pheno.data

gene.count.matrix <- read.table(file = "UNnormalized_counts.csv",header = T,sep = ";")
colnames(gene.count.matrix) = c("gene.ID","GFP1","GFP2","GFP3","S1","S2","S3")
#gene.count.matrix <- gene.count.matrix[,-1]

genesid = gene.count.matrix$gene.ID
head(genesid)

gene.count.matrix = dplyr::select(gene.count.matrix, c("GFP1","GFP2","GFP3","S1","S2","S3"))
rownames(gene.count.matrix) = genesid
head(gene.count.matrix)


#### observación de los datos
boxplot(log2(gene.count.matrix), outline=F,col=c(rep("green",3),rep("blue",3)),
        ylab="log2(Gene Expression)",cex.lab=1.5)



#PCA and dendogram
library(FactoMineR)
library(factoextra)


pca.gene.expression <- data.frame(colnames(gene.count.matrix),
                                  t(gene.count.matrix))
colnames(pca.gene.expression)[1] <- "Sample"
head(pca.gene.expression)

res.pca <- PCA(pca.gene.expression, graph = F,scale.unit = TRUE,quali.sup = 1 )
res.hcpc <- HCPC(res.pca, graph=T,nb.clust = 2)   
fviz_dend(res.hcpc,k=2,
          cex = 0.75,                       # Label size
          palette = "jco",               # Color palette see ?ggpubr::ggpar
          rect = TRUE, rect_fill = TRUE, # Add rectangle around groups
          rect_border = "jco",           # Rectangle color
          type="rectangle",
          labels_track_height = 1400      # Augment the room for labels
)

fviz_pca_ind(res.pca, col.ind = c("gf1","gf2","gf3","s1","s2","s3"), 
             pointsize=2, pointshape=21,fill="black",
             repel = TRUE, 
             #addEllipses = TRUE,ellipse.type = "confidence",
             legend.title="Samples",
             title="Principal Component Analysis GFP vs S",
             show_legend=TRUE,show_guide=TRUE)


#Deseq2
dds <- DESeqDataSetFromMatrix(countData=gene.count.matrix, colData=pheno.data, design = ~ genotype)
dds

dds <- DESeq(dds) 
res <- results(dds) 
head(res)
summary(res)


#Cuantificacion de genes expresados diferencialmente
log.fold.change <- res$log2FoldChange
q.value <- res$padj
names(log.fold.change) <- genesid
names(q.value) <- genesid

activated.genes.deseq2 <- genesid[log.fold.change > 2 & q.value < 0.05]
activated.genes.deseq2 <- activated.genes.deseq2[!is.na(activated.genes.deseq2)]

repressed.genes.deseq2 <- genesid[log.fold.change < - 2 & q.value < 0.05]
repressed.genes.deseq2 <- repressed.genes.deseq2[!is.na(repressed.genes.deseq2)]

length(activated.genes.deseq2)
length(repressed.genes.deseq2)

###############


## Volcano plot
log.q.val <- -log10(q.value)
plot(log.fold.change,log.q.val,pch=19,col="grey",cex=0.8,
     xlim=c(-12,12),ylim = c(-5,300), 
     xlab="log2(Fold-chage)",ylab="-log10(q-value)",cex.lab=1.5,main = "GFP vs s")

points(x = log.fold.change[activated.genes.deseq2],
       y = log.q.val[activated.genes.deseq2],col="red",cex=0.8,pch=19)
points(x = log.fold.change[repressed.genes.deseq2],
       y = log.q.val[repressed.genes.deseq2],col="blue",cex=0.8,pch=19)




#Enriquecimiento funcional
#BiocManager::install("clusterProfiler")
library(clusterProfiler)
#BiocManager::install("org.Mm.eg.db")
#Libreria para organismo mus musculus
library(org.Mm.eg.db)

#Activated genes
#ontology BP
activated.atha.enrich.go_BP <- enrichGO(gene          = activated.genes.deseq2,
                                     OrgDb         = org.Mm.eg.db,
                                     ont           = "BP",
                                     pAdjustMethod = "BH",
                                     pvalueCutoff  = 0.05,
                                     readable      = FALSE,
                                     keyType = "ENSEMBL")

barplot(activated.atha.enrich.go_BP,showCategory = 22)
dotplot(activated.atha.enrich.go,showCategory = 20)


#ontology CC
activated.atha.enrich.go_CC <- enrichGO(gene          = activated.genes.deseq2,
                                     OrgDb         = org.Mm.eg.db,
                                     ont           = "CC",
                                     pAdjustMethod = "BH",
                                     pvalueCutoff  = 0.05,
                                     readable      = FALSE,
                                     keyType = "ENSEMBL")

barplot(activated.atha.enrich.go_CC,showCategory = 22)

#ontology MF
activated.atha.enrich.go_MF <- enrichGO(gene          = activated.genes.deseq2,
                                     OrgDb         = org.Mm.eg.db,
                                     ont           = "MF",
                                     pAdjustMethod = "BH",
                                     pvalueCutoff  = 0.05,
                                     readable      = FALSE,
                                     keyType = "ENSEMBL")

barplot(activated.atha.enrich.go_MF,showCategory = 22)





#Repressed genes
#ontology BP
repressed.atha.enrich.go_BP <- enrichGO(gene          = repressed.genes.deseq2,
                                     OrgDb         = org.Mm.eg.db,
                                     ont           = "BP",
                                     pAdjustMethod = "BH",
                                     pvalueCutoff  = 0.05,
                                     readable      = FALSE,
                                     keyType = "ENSEMBL")

barplot(repressed.atha.enrich.go_BP,showCategory = 22)
dotplot(repressed.atha.enrich.go_BP,showCategory = 20)

#ontology CC
repressed.atha.enrich.go_CC <- enrichGO(gene          = repressed.genes.deseq2,
                                        OrgDb         = org.Mm.eg.db,
                                        ont           = "CC",
                                        pAdjustMethod = "BH",
                                        pvalueCutoff  = 0.05,
                                        readable      = FALSE,
                                        keyType = "ENSEMBL")

barplot(repressed.atha.enrich.go_CC,showCategory = 22)

#ontology MF
repressed.atha.enrich.go_MF <- enrichGO(gene          = repressed.genes.deseq2,
                                        OrgDb         = org.Mm.eg.db,
                                        ont           = "MF",
                                        pAdjustMethod = "BH",
                                        pvalueCutoff  = 0.05,
                                        readable      = FALSE,
                                        keyType = "ENSEMBL")

barplot(repressed.atha.enrich.go_MF,showCategory = 22)

































