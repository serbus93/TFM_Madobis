# TFM MADOBIS
#Alumno: Sergio Bustamante

#Analisis Proteomica


library(readxl)
library(tidyverse)
library(dplyr)
library(tidyr)
#importing data
data = read_excel("Protein_intensity.xlsx")
#En este caso, se utilizan datos brutos de entrada, normalizados previamente

# Remove the missing values
data_tidy <- data %>% drop_na()
# Number of proteins in original data
data %>% summarise(Number_of_proteins = n())
# Number of proteins without missing values
data_tidy %>% summarise(Number_of_proteins = n())

#exploring data
glimpse(data_tidy)
str(data_tidy)
summary(data_tidy)
boxplot(log2(data_tidy[,3:8]),outline = F,col=c(rep("green",3),rep("blue",3,)),
        ylab="log2(Protein intensity)",cex.lab=1.5)


#PCA
library(FactoMineR)
library(factoextra)

data_pca = data_tidy
rownames(data_pca) = protid
data_pca = data_pca[,3:8]
head(data_pca)

pca_proteins <- data.frame(colnames(data_pca),
                                  t(data_pca))
colnames(pca_proteins)[1] <- "Sample"
head(pca_proteins)

res.pca <- PCA(pca_proteins, graph = F,scale.unit = TRUE,quali.sup = 1 )
fviz_pca_ind(res.pca, col.ind = c("gf1","gf2","gf3","s1","s2","s3"), 
             pointsize=2, pointshape=21,fill="black",
             repel = TRUE, 
             #addEllipses = TRUE,ellipse.type = "confidence",
             legend.title="Samples",
             title="Principal Component Analysis GFP vs S",
             show_legend=TRUE,show_guide=TRUE)



#Hypothesis testing with the t-test
# T-test function for multiple experiments
t_test <- function(dt,grp1,grp2){
  # Subset control group and convert to numeric
  x <- dt[grp1] %>% unlist %>% as.numeric()
  # Subset treatment group and convert to numeric
  y <- dt[grp2] %>% unlist %>% as.numeric()
  # Perform t-test using the mean of x and y
  result <- t.test(x, y)
  # Extract p-values from the results
  p_vals <- tibble(p_val = result$p.value)
  # Return p-values
  return(p_vals)
} 



# Apply t-test function to data using plyr adply
#  .margins = 1, slice by rows, .fun = t_test plus t_test arguements
dat_pvals <- plyr::adply(data_tidy,.margins = 1, .fun = t_test, 
                         grp1 = c(3:5), grp2 = c(6:8)) %>% as.tibble()

# Perform t-test on first protein
# Columnas 3 a 5: control-GFP
# Columnas 6 a 8: S-GFP
t.test(as.numeric(data_tidy[1,3:5]),
       as.numeric(data_tidy[1,6:8]))$p.value


# Plot histogram
dat_pvals %>% 
  ggplot(aes(p_val)) + 
  geom_histogram(binwidth = 0.05, 
                 boundary = 0.5, 
                 fill = "darkblue",
                 colour = "white") +
  xlab("p-value") +
  ylab("Frequency") +
  theme_minimal()


#Calculate the Fold Change (FC)
# Select columns and log data
dat_log <- dat_pvals %>% 
  select(-c(Protein,gene.names,p_val)) %>% 
  log2()

# Bind columns to create transformed data frame
dat_combine <- bind_cols(dat_pvals[,c(1:2)], dat_log, dat_pvals[,9]) 

dat_fc <- dat_combine %>% 
  group_by(Protein) %>% 
  mutate(mean_GFP = mean(c(GFP.1,
                           GFP.3,
                           GFP.4)),
         mean_S= mean(c(S.1,
                                S.2,
                                S.3)),
         log_fc = mean_S - mean_GFP,
         log_pval = -1*log10(p_val))

# Final transformed data
dat_tf <- dat_fc %>% dplyr::select(Protein,
                            gene.names,
                            log_fc, log_pval)

head(dat_tf)

#Visualising the transformed data
# Plot a histogram to look at the distribution.
dat_tf %>%
  ggplot(aes(log_fc)) + 
  geom_histogram(binwidth = 0.5,
                 boundary = 0.5,
                 fill = "darkblue",
                 colour = "white") +
  xlab("log2 fold change") +
  ylab("Frequency") +
  theme_minimal()


#Cuantificacion de proteinas expresadas diferencialmente
protid = dat_tf$Protein
log.fold.change <- dat_tf$log_fc
log.p.value <- dat_tf$log_pval
names(log.fold.change) <- protid
names(log.p.value) <- protid

#logFC = 1 equals to FC = 2
#logFC = 0.5 equals to +-45% expression
#-log Pvalue = 1.3 equals to Pvalue = 0.05
#Comparison GFP vs S. Activated genes = genes activated in GFP (control) and 
#less expressed in condition S
activated.proteins <- protid[log.fold.change  > 0.5 & log.p.value > 1.3]
activated.proteins <- activated.proteins[!is.na(activated.proteins)]

repressed.proteins <- protid[log.fold.change < - 0.5 & log.p.value > 1.3]
repressed.proteins <- repressed.proteins[!is.na(repressed.proteins)]

length(activated.proteins)
length(repressed.proteins)



#Volcano plot
dat_tf %>% ggplot(aes(log_fc,log_pval)) + geom_point()

dat_tf %>%
  # Add a threshold for significant observations
  #threshold log_pval = 0.6 = pvalue 0.05
  
  mutate(threshold = if_else(log_fc >= 0.5 & log_pval >= 1.3 |
                               log_fc <= -0.5 & log_pval >= 1.3,"A", "B")) %>%
  # Plot with points coloured according to the threshold
  ggplot(aes(log_fc,log_pval, colour = threshold)) +
  geom_point(alpha = 0.5) + # Alpha sets the transparency of the points
  # Add dotted lines to indicate the threshold, semi-transparent
  geom_hline(yintercept = 1.3, linetype = 2, alpha = 0.5) + 
  geom_vline(xintercept = 0.5, linetype = 2, alpha = 0.5) +
  geom_vline(xintercept = -0.5, linetype = 2, alpha = 0.5) +
  # Set the colour of the points
  scale_colour_manual(values = c("A"= "red", "B"= "black")) +
  xlab("log2 fold change") + ylab("-log10 p-value") + # Relabel the axes
  theme_minimal() + # Set the theme
  theme(legend.position="none") # Hide the legend

#Extract significant proteins
sig_proteins = dat_tf %>%
  # Filter for significant observations
  filter(log_pval >= 1.3 & (log_fc >= 0.5 | log_fc <= -0.5)) %>% 
  # Ungroup the data
  ungroup() %>% 
  # Select columns of interest
  dplyr::select(Protein,gene.names,log_fc,log_pval)

head(sig_proteins)
write.table(sig_proteins, file = "sig_proteins.tsv", sep = "\t",row.names = F,quote = F)

#Heatmap
# Keep the same p-val cut-off, And SAME LOG FC.
dat_filt <- dat_fc %>%
  filter(log_pval >= 1.3 & (log_fc >= 0.5 | log_fc <= -0.5))

# Convert to matrix data frame
dat_matrix <- as.matrix.data.frame(dat_filt[,3:8]) 
# Name the rows with protein ids
row.names(dat_matrix) <- dat_filt$Protein
# Transpose and scale the data to a mean of zero and sd of one
dat_scaled <- scale(t(dat_matrix)) %>% t()


# Transpose the matrix to calculate distance between experiments, row-wise
d1 <- dat_scaled %>% t() %>%
  dist(.,method = "euclidean", diag = FALSE, upper = FALSE)
# Calculate the distance between proteins row-wise 
d2 <- dat_scaled %>%
  dist(.,method = "euclidean", diag = FALSE, upper = FALSE)

# Show the values for d1
round(d1,2)


# Clustering distance between experiments using Ward linkage
c1 <- hclust(d1, method = "ward.D2", members = NULL)
# Clustering distance between proteins using Ward linkage
c2 <- hclust(d2, method = "ward.D2", members = NULL)

# Check clustering by plotting dendrograms
par(mfrow=c(2,1),cex=0.5) # Make 2 rows, 1 col plot frame and shrink labels
plot(c1); plot(c2) # Plot both cluster dendrograms



# Set colours for heatmap, 20 increments
my_palette <- colorRampPalette(c("blue","white","red"))(n = 20)

# Plot heatmap with heatmap.2
par(cex.main=0.75) # Shrink title fonts on plot
dat_scaled %>% 
  # Rename the comlums
  magrittr::set_colnames(c("GFP 1", "GFP 3", "GFP 4",
                           "S 1", "S 2", "S 3")) %>% 
  # Plot heatmap
  gplots::heatmap.2(.,                     
                    Colv=as.dendrogram(c1),     # Experiments clusters in cols
                    Rowv=as.dendrogram(c2),     # Protein clusters in rows
                    revC=TRUE,                  # Flip plot to match pheatmap
                    density.info="histogram",   # Plot histogram of data and colour key
                    trace="none",               # Turn of trace lines from heat map
                    col = my_palette,           # Use my colour scheme
                    cexRow=0.6,cexCol=0.75)     # Amend row and column label fonts


#Enriquecimiento funcional
#BiocManager::install("clusterProfiler")
library(clusterProfiler)
#BiocManager::install("org.Mm.eg.db")
library(org.Mm.eg.db)

#Activated proteins
#ontology BP
activated.atha.enrich.go_BP <- enrichGO(gene          = activated.proteins,
                                        OrgDb         = org.Mm.eg.db,
                                        ont           = "BP",
                                        pAdjustMethod = "BH",
                                        pvalueCutoff  = 0.05,
                                        readable      = FALSE,
                                        keyType = "UNIPROT")

barplot(activated.atha.enrich.go_BP,showCategory = 20)

#Activated genes
#ontology CC
activated.atha.enrich.go_CC <- enrichGO(gene          = activated.proteins,
                                        OrgDb         = org.Mm.eg.db,
                                        ont           = "CC",
                                        pAdjustMethod = "BH",
                                        pvalueCutoff  = 0.05,
                                        readable      = FALSE,
                                        keyType = "UNIPROT")

barplot(activated.atha.enrich.go_CC,showCategory = 20)
#ontology MF
activated.atha.enrich.go_MF <- enrichGO(gene          = activated.proteins,
                                        OrgDb         = org.Mm.eg.db,
                                        ont           = "MF",
                                        pAdjustMethod = "BH",
                                        pvalueCutoff  = 0.05,
                                        readable      = FALSE,
                                        keyType = "UNIPROT")

barplot(activated.atha.enrich.go_MF,showCategory = 15)



#Repressed proteins
#ontology BP
repressed.atha.enrich.go_BP <- enrichGO(gene          = repressed.proteins,
                                        OrgDb         = org.Mm.eg.db,
                                        ont           = "BP",
                                        pAdjustMethod = "BH",
                                        pvalueCutoff  = 0.05,
                                        readable      = FALSE,
                                        keyType = "UNIPROT")

barplot(repressed.atha.enrich.go_BP,showCategory = 15)

#ontology CC
repressed.atha.enrich.go_CC <- enrichGO(gene          = repressed.proteins,
                                        OrgDb         = org.Mm.eg.db,
                                        ont           = "CC",
                                        pAdjustMethod = "BH",
                                        pvalueCutoff  = 0.05,
                                        readable      = FALSE,
                                        keyType = "UNIPROT")

barplot(repressed.atha.enrich.go_CC,showCategory = 15)

#ontology MF
repressed.atha.enrich.go_MF <- enrichGO(gene          = repressed.proteins,
                                        OrgDb         = org.Mm.eg.db,
                                        ont           = "MF",
                                        pAdjustMethod = "BH",
                                        pvalueCutoff  = 0.05,
                                        readable      = FALSE,
                                        keyType = "UNIPROT")

barplot(repressed.atha.enrich.go_MF,showCategory = 15)



