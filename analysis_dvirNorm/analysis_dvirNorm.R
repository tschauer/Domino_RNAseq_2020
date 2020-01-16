


rm(list = ls())


setwd(dirname(rstudioapi::getActiveDocumentContext()$path))



########################################################  libraries ######################################################## 


library(GenomicFeatures)
library(AnnotationDbi)
library(org.Dm.eg.db)
library(GO.db)
library(RColorBrewer)
library(DESeq2)
library(affy)
library(pheatmap)
library(scales)
library(genefilter)
library(ggplot2)
library(sva)
library(grid)
library(GenomicFeatures)
library(gdata)
library(gridExtra)
library(GenomicAlignments)
library(rtracklayer)



source("../functions/functions_DESeq.R")


cbPalette <- c("#999999", "#0072B2", "#CC79A7", "#009E73", "#E69F00", "#D55E00", "#56B4E9", "#F0E442")




############################################################################################################################# 
############################################################################################################################# 
############################################################################################################################# 




############################################################################################################################# 
############################################################################################################################# 
############################################################################################################################# 

########################################################  annotation ######################################################## 



# dir.gtf <- "../../genome/"
# gtffile <- file.path(dir.gtf,"dmel-all-r6.17.gtf")
# txdb <- makeTxDbFromGFF(gtffile, format="gtf")
# 
# saveDb(txdb, file = "txdb.sqlite")


txdb <- loadDb("txdb.sqlite")




############################################################################################################################# 
############################################################################################################################# 
############################################################################################################################# 

########################################################  SampleTable  ######################################################



SampleTable <- read.xls("../SampleTable.xlsx", sheet = 1, header = TRUE)


SampleTable$Conditions <- factor(gsub("_r.*", "", SampleTable$Name))
SampleTable$Conditions <- relevel(SampleTable$Conditions, ref = "Control")



########################################################  count data  ########################################################


myCount_Dir <- "../counts/"

myCount_Files <- list.files(path = myCount_Dir, pattern = ".*dmel.*.out.tab")
myCount_File_Path <- file.path(myCount_Dir, myCount_Files)



dmel.counts <-  makeCountTable(count_files = myCount_Files, 
                             count_file_path = myCount_File_Path,
                             stranded = TRUE)


dmel.counts_genes <- dmel.counts[grep("FBgn", rownames(dmel.counts)),]



write.csv(dmel.counts_genes, "dmel.counts_genes.csv")




######################################################## 




myCount_Files <- list.files(path = myCount_Dir, pattern = ".*dvir.*.out.tab")
myCount_File_Path <- file.path(myCount_Dir, myCount_Files)



dvir.counts <-  makeCountTable(count_files = myCount_Files, 
                               count_file_path = myCount_File_Path,
                               stranded = TRUE)


dvir.counts_genes <- dvir.counts[grep("FBgn", rownames(dvir.counts)),]



write.csv(dvir.counts_genes, "dvir.counts_genes.csv")




############################################################################################################################# 
############################################################################################################################# 
############################################################################################################################# 



############################################################################################################################# 
############################################################################################################################# 
############################################################################################################################# 













############################################################################################################################# 
############################################################################################################################# 
############################################################################################################################# 



############################################################################################################################# 
############################################################################################################################# 
############################################################################################################################# 


############################################            Setup DESeq           ###############################################




dds.dmel <- setupDDS(CountTableName = "dmel.counts_genes", n_samples_for_filtering = 6)

dds.dvir <- setupDDS(CountTableName = "dvir.counts_genes", n_samples_for_filtering = 6)



my_conditions <- colData(dds.dmel)$Sample





######################################################## 




pdf("DESeq_sizefactors.pdf", width = 8, height = 8, useDingbats = FALSE)

par(oma=c(5,5,5,5), mar=c(5,5,2,2), mgp=c(2.5,1,0), bg=NA,
    cex.axis = 1.2, cex.main = 1.5, cex.lab=1.33, pch=19, cex=1)

par(mfrow=c(2,2))


plot(sizeFactors(dds.dmel), sizeFactors(dds.dvir),
     xlab = "dmel", ylab = "dvir",
     xlim = c(0.6,1.4), ylim = c(0.6,1.4),
     pch = 19, col = cbPalette[my_conditions],
     main = "Norm. Factors")

abline(coef = c(0,1), col="grey32", lty=2)

plot(0,0, type="n", xlab="", ylab="", xaxt="n", yaxt="n", bty="n")

legend("center", legend = levels(my_conditions), 
       pch =19, col = cbPalette[seq_along(levels(my_conditions))])

dev.off()




######################################################## 





sizeFactors(dds.dmel) <- sizeFactors(dds.dvir)


dds.dmel <- DESeq(dds.dmel)



######################################################## 




############################################################################################################################# 
############################################################################################################################# 
############################################################################################################################# 



############################################################################################################################# 
############################################################################################################################# 
############################################################################################################################# 


############################################        Batch Correction          ###############################################



batchVar <- colData(dds.dmel)$Batch

modcombat <- model.matrix(~Sample, data = colData(dds.dmel))

log2_counts <- ComBat(dat = log2(counts(dds.dmel, normalized = TRUE)+1), batch = batchVar, mod = modcombat, par.prior = TRUE, prior.plots = FALSE)





######################################################## 



pdf("DESeq_PCA.pdf", width = 8, height = 8, useDingbats = FALSE)

par(oma=c(5,5,5,5), mar=c(5,5,2,2), mgp=c(2.5,1,0), 
    cex.axis = 1.2, cex.main = 1.5, cex.lab=1.33, pch=19, cex=1)

par(mfrow=c(2,2))


plottingPCA(my_data = log2_counts,
            color_palette = cbPalette,
            conditions = my_conditions,
            quantiles = c(0,1),
            point_size = 1,
            show_labels = FALSE)



plottingPCA(my_data = log2_counts, 
            xcomp = 1, ycomp = 3,
            color_palette = cbPalette,
            conditions = my_conditions,
            quantiles = c(0,1),
            point_size = 1,
            show_labels = FALSE)

plot(0,0, type="n", xlab="", ylab="", xaxt="n", yaxt="n", bty="n")

legend("center", legend = levels(my_conditions), 
       pch = 19, col = cbPalette[seq_along(levels(my_conditions))])





dev.off()


######################################################## 



############################################################################################################################# 
############################################################################################################################# 
############################################################################################################################# 



############################################################################################################################# 
############################################################################################################################# 
############################################################################################################################# 


############################################        Favorite Genes          ###############################################




pdf("DESeq_example_genes.pdf", width = 8, height = 8, useDingbats = FALSE)

par(oma=c(5,5,5,5), mar=c(5,5,2,2), mgp=c(2.5,1,0), 
    cex.axis = 1.2, cex.main = 1.5, cex.lab=1.33, pch=19, cex=1)

par(mfrow=c(2,2))



my_favorite_genes <- c("dom", "His2Av", "Tip60")


for(my_favorite_gene in my_favorite_genes){
      
      my_gene_id <- mapIds(org.Dm.eg.db, my_favorite_gene, keytype="SYMBOL", column = "FLYBASE", multiVals="first")
      
      log2_counts_favorite <- log2_counts[rownames(log2_counts)==my_gene_id,]
      log2_counts_favorite_control <- mean(log2_counts_favorite[names(log2_counts_favorite) %in% SampleTable$SampleID[SampleTable$Conditions == "Control"]])
      log2_counts_favorite_norma <- log2_counts_favorite - log2_counts_favorite_control
      
      
      plotDots(my_data = log2_counts_favorite_norma, 
               my_title = my_favorite_gene, 
               color_palette = cbPalette, 
               color_groups = my_conditions, 
               conditions = my_conditions, 
               point_size = 0.7,
               x_label = levels(my_conditions),
               ylims = c(-6,6),
               xlims = c(0.5,length(levels(my_conditions))+0.5))
      
      axis(2, at =0, labels = "relative log2 counts", line = 1.5, lwd = 0)
      
}

plot(0,0, type="n", xlab="", ylab="", xaxt="n", yaxt="n", bty="n")

legend("center", legend = levels(my_conditions), 
       pch = 19, col = cbPalette[seq_along(levels(my_conditions))])



dev.off()



############################################################################################################################# 
############################################################################################################################# 
############################################################################################################################# 



############################################################################################################################# 
############################################################################################################################# 
############################################################################################################################# 


############################################      Tests and MA plots          ###############################################






pdf("DESeq_MAplots.pdf", width = 8, height = 8, useDingbats = FALSE)

par(oma=c(5,5,5,5), mar=c(5,5,2,2), mgp=c(2.5,1,0), 
    cex.axis = 1.2, cex.main = 1.5, cex.lab=1.33, pch=19, cex=1)



par(mfrow=c(2,2))



my_contrast <- c("Sample", "DomA", "Control")

summarizeResults(dds = dds.dmel, 
                 contrast = my_contrast, 
                 favorite_gene = NULL,
                 pval_cutoff = 0.01,
                 lfc_cutoff = 0,
                 ylims = c(-8,8),
                 xlims = c(1, 10^6), 
                 plotMA = TRUE)


my_contrast <- c("Sample", "Tip60", "Control")

summarizeResults(dds = dds.dmel, 
                 contrast = my_contrast, 
                 favorite_gene = NULL,
                 pval_cutoff = 0.01,
                 lfc_cutoff = 0,
                 ylims = c(-8,8),
                 xlims = c(1, 10^6), 
                 plotMA = TRUE)








my_contrast <- c("Sample", "DomB", "Control")

summarizeResults(dds = dds.dmel, 
                 contrast = my_contrast, 
                 favorite_gene = NULL,
                 pval_cutoff = 0.01,
                 lfc_cutoff = 0,
                 ylims = c(-8,8),
                 xlims = c(1, 10^6), 
                 plotMA = TRUE)



my_contrast <- c("Sample", "H2Av", "Control")

summarizeResults(dds = dds.dmel, 
                 contrast = my_contrast, 
                 favorite_gene = NULL,
                 pval_cutoff = 0.01,
                 lfc_cutoff = 0,
                 ylims = c(-8,8),
                 xlims = c(1, 10^6), 
                 plotMA = TRUE)


my_contrast <- c("Sample", "DomA", "DomB")

summarizeResults(dds = dds.dmel, 
                 contrast = my_contrast, 
                 favorite_gene = NULL,
                 pval_cutoff = 0.01,
                 lfc_cutoff = 0,
                 ylims = c(-8,8),
                 xlims = c(1, 10^6), 
                 plotMA = TRUE)




dev.off()








############################################################################################################################# 
############################################################################################################################# 
############################################################################################################################# 



############################################################################################################################# 
############################################################################################################################# 
############################################################################################################################# 

############################################           log2FC plot            ###############################################




######################################################## 




pdf("DESeq_log2FC.pdf", width = 8, height = 8, useDingbats = FALSE)

par(oma=c(5,5,5,5), mar=c(5,5,2,2), mgp=c(2.5,1,0), 
    cex.axis = 1.2, cex.main = 1.5, cex.lab=1.33, pch=19, cex=1)



par(mfrow=c(2,2))



plotLog2FC(dds = dds.dmel, 
           contrast1 = c("Sample", "H2Av", "Control"), 
           contrast2 = c("Sample", "DomA", "Control"), 
           favorite_genes1 = "", 
           favorite_genes2 = "", 
           selection_name1 = "", 
           selection_color1 = rgb(0.8,0,0,1),
           selection_color2 = rgb(1,1,1,1),
           selection_name2 = "",
           text_label1 = TRUE,
           text_label2 = FALSE,
           pval_cutoff = 0.01, 
           lfc_cutoff = 0, 
           lims = c(-8,8)
)


plotLog2FC(dds = dds.dmel, 
           contrast1 = c("Sample", "H2Av", "Control"), 
           contrast2 = c("Sample", "DomB", "Control"), 
           favorite_genes1 = "", 
           favorite_genes2 = "", 
           selection_name1 = "", 
           selection_color1 = rgb(0.8,0,0,1),
           selection_color2 = rgb(1,1,1,1),
           selection_name2 = "",
           text_label1 = TRUE,
           text_label2 = FALSE,
           pval_cutoff = 0.01, 
           lfc_cutoff = 0, 
           lims = c(-8,8)
)


plotLog2FC(dds = dds.dmel, 
           contrast1 = c("Sample", "Tip60", "Control"), 
           contrast2 = c("Sample", "DomA", "Control"), 
           favorite_genes1 = "", 
           favorite_genes2 = "", 
           selection_name1 = "", 
           selection_color1 = rgb(0.8,0,0,1),
           selection_color2 = rgb(1,1,1,1),
           selection_name2 = "",
           text_label1 = TRUE,
           text_label2 = FALSE,
           pval_cutoff = 0.01, 
           lfc_cutoff = 0, 
           lims = c(-8,8)
)


plotLog2FC(dds = dds.dmel, 
           contrast1 = c("Sample", "Tip60", "Control"), 
           contrast2 = c("Sample", "DomB", "Control"), 
           favorite_genes1 = "", 
           favorite_genes2 = "", 
           selection_name1 = "", 
           selection_color1 = rgb(0.8,0,0,1),
           selection_color2 = rgb(1,1,1,1),
           selection_name2 = "",
           text_label1 = TRUE,
           text_label2 = FALSE,
           pval_cutoff = 0.01, 
           lfc_cutoff = 0, 
           lims = c(-8,8)
)



plotLog2FC(dds = dds.dmel, 
           contrast1 = c("Sample", "DomA", "Control"), 
           contrast2 = c("Sample", "DomB", "Control"), 
           favorite_genes1 = "", 
           favorite_genes2 = "", 
           selection_name1 = "", 
           selection_color1 = rgb(0.8,0,0,1),
           selection_color2 = rgb(1,1,1,1),
           selection_name2 = "",
           text_label1 = TRUE,
           text_label2 = FALSE,
           pval_cutoff = 0.01, 
           lfc_cutoff = 0, 
           lims = c(-8,8)
)



dev.off()

######################################################## 







############################################################################################################################# 
############################################################################################################################# 
############################################################################################################################# 



############################################################################################################################# 
############################################################################################################################# 
############################################################################################################################# 





