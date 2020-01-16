

library(GenomicFeatures)
library(AnnotationDbi)
library(org.Dm.eg.db)
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
library(lattice)
library(magick)



######################################################## make Count Table ######################################################## 



makeCountTable <- function(count_files, count_file_path, stranded = FALSE){
  
      if(stranded){
        cidx <- 4
      } else {
        cidx <- 2
      }
  
      for(i in seq_along(count_files)){
        
        tmp <- read.table(count_file_path[i])
        
        if(i == 1){
          my_counts <- tmp[,cidx] 
        } else {
          my_counts <- cbind(my_counts, tmp[,cidx])
        }
      }
      
      rownames(my_counts) <- tmp[,1]
      colnames(my_counts) <- gsub("_[G,A,T,C].*","", count_files)
      
      return(my_counts)
}





######################################################## make Kmer Table ######################################################## 



makeKmerTable <- function(kmer_files, kmer_file_path, sample_names){
      
      for(i in seq_along(kmer_files)){
            
            if(i == 1){
                  
                  tmp1 <- read.table(kmer_file_path[i])      
                  colnames(tmp1) <- c("id", sample_names[i])
                  
            } else {
                  
                  tmp2 <- read.table(kmer_file_path[i])
                  colnames(tmp2) <- c("id", sample_names[i])
                  
                  
                  tmp1 <- merge(data.table(tmp1), 
                                data.table(tmp2), by = "id")
            }
      }
      
      tmp1 <- as.data.frame(tmp1)
      rownames(tmp1) <- tmp1$id
      tmp1 <- tmp1[,-1]
      
      return(tmp1)
      
}







######################################################## Summarize Results and MAplot ######################################################## 



# magickPoints <- function(x, y, ...) {
#       w <- convertWidth(unit(1, "npc"), "in", valueOnly=TRUE)
#       h <- convertHeight(unit(1, "npc"), "in", valueOnly=TRUE)
#       cvp <- current.viewport()
#       dev <- dev.cur()
#       raster <- image_graph(width=w*72, height=h*72)
#       pushViewport(viewport(xscale=cvp$xscale, yscale=cvp$yscale))
#       panelDefault(x, y, ...)
#       dev.off()
#       dev.set(dev)
#       grid.raster(raster)
# }


summarizeResults <- function(dds, contrast, favorite_gene, 
                             pval_cutoff, lfc_cutoff = 0,
                             ylims = c(-5,5), xlims = c(1, 10^6),
                             plotMA = TRUE){
  
    res <- results(dds, 
                   contrast = contrast, 
                   lfcThreshold = lfc_cutoff, independentFiltering = FALSE
                   #, alpha = pval_cutoff
                   )
    
    #res <- lfcShrink(dds, contrast = contrast, res = res)
    
    
    res$chr <- select(txdb, rownames(res), "TXCHROM", keytype="GENEID", multiVals="first")$TXCHROM
    res$symbol <- mapIds(org.Dm.eg.db, rownames(res), "SYMBOL", keytype="FLYBASE", multiVals="first")
    
    
    res$padj[is.na(res$padj)] <- 1
    res$log2FoldChange[is.na(res$log2FoldChange)] <- 0
    res.sign <- subset(res, padj < pval_cutoff)
    
    
    my_main <- paste(contrast[2:3], collapse = " - ")
    my_main <- gsub("plus", "+", my_main)
    my_main <- gsub("minus", "-", my_main)
    my_main <- gsub("\\_", "", my_main)
    my_main <- gsub("\\.", " ", my_main)
    

    # p1 <- xyplot(res$log2FoldChange ~ log10(res$baseMean+1), panel = magickPoints,
    #      main = my_main, xlab = "mean counts", ylab = "log2 fold change",
    #      col=rgb(0,0,0,0.1),  ylim=ylims, xlim=xlims, pch=19, cex = 0.25)
    # 
    # update(p1, panel = function(...) {
    #       panel.xyplot(...)
    #       xyplot(res.sign$log2FoldChange ~ log10(res.sign$baseMean+1), #panel = magickPoints,
    #              #main = my_main, xlab = "mean counts", ylab = "log2 fold change",
    #              col=rgb(0.8,0,0,0.5),  ylim=ylims, xlim=xlims, pch=19, cex = 0.25)
    # })

    if(plotMA){
          plot(res$baseMean, res$log2FoldChange, log="x",
              xlab = "log10 mean counts", ylab = "log2 fold change", xaxt="n",
              col=rgb(0,0,0,0.1),  ylim=ylims, xlim=xlims, pch=19, cex = 0.25)
          
          axis(side = 1, at = c(1,10^2,10^4,10^6), labels = log10(c(1,10^2,10^4,10^6)))
          #axis(side = 1, at = c(100), labels = c(2))
          
          
          points(res.sign$baseMean, res.sign$log2FoldChange, col=rgb(0.8,0,0,0.5), pch=19, cex = 0.25)
          
          mtext(text = my_main, side = 3, line = 0.5, adj = 0.5, font=2, cex = 1.2)

          for(i in seq_along(favorite_gene)){
                
                my_favorite <- grep(paste0("^",favorite_gene[i], "$"), res$symbol)
                
                if(res$padj[my_favorite] < pval_cutoff){
                      
                      points(res$baseMean[my_favorite],
                             res$log2FoldChange[my_favorite],
                             col=rgb(0.9,0.6,0,1), pch=19)
                      
                      text(res$baseMean[my_favorite],
                           res$log2FoldChange[my_favorite],
                           labels = res$symbol[my_favorite], adj = c(0,-0.5),
                           col=rgb(0.9,0.6,0,1))
                } else {
                      
                      points(res$baseMean[my_favorite],
                             res$log2FoldChange[my_favorite],
                             col="grey", pch=19)
                      
                      text(res$baseMean[my_favorite],
                           res$log2FoldChange[my_favorite],
                           labels = res$symbol[my_favorite], adj = c(0,-0.5),
                           col="grey")
                }
          }
          
          
          
          abline(h=0, col="grey32")
          
          legend("bottomright", legend =  c(paste("up", sum(res.sign$log2FoldChange > 0), sep = " - "),
                                            paste("down", sum(res.sign$log2FoldChange < 0), sep = " - ")),
                 bg = NA, cex = 0.8, border = NA, bty = "n")
          
    }
    

    
    my_main <- gsub(" ","", my_main)
    
    write.table(res[order(res$padj),], file = paste("res", my_main, "txt", sep="."), quote = F, sep = "\t", row.names = T, col.names = NA)
    #write.csv(res.sign[order(res.sign$padj),], file = paste("res.sign", my_main, "csv", sep="."))
    
  
}





######################################################## Summarize Kmer Results and MAplot ######################################################## 



summarizeKmerResults <- function(dds, contrast, favorite_gene, pval_cutoff, lfc_cutoff = 0, ylims = c(-5,5)){
      
      res <- results(dds, 
                     contrast = contrast, 
                     lfcThreshold = lfc_cutoff, independentFiltering = FALSE
                     #, alpha = pval_cutoff
      )
      
      #res <- lfcShrink(dds, contrast = contrast, res = res)
      
      
      res$symbol <- rownames(res)
      
      
      res$padj[is.na(res$padj)] <- 1
      res$log2FoldChange[is.na(res$log2FoldChange)] <- 0
      res.sign <- subset(res, padj < pval_cutoff)
      
      
      my_main <- paste(contrast[2:3], collapse = " - ")
      my_main <- gsub("plus", "+", my_main)
      my_main <- gsub("minus", "-", my_main)
      
      
      plot(res$baseMean, res$log2FoldChange, log="x", 
           main = my_main, xlab = "mean expression", ylab = "log2 fold change",
           col=rgb(0,0,0,0.1),  ylim=ylims, pch=19)
      points(res.sign$baseMean, res.sign$log2FoldChange, col=rgb(0.8,0,0,0.5), pch=19, cex = 0.5)
      
      
      for(i in seq_along(favorite_gene)){
            
            my_favorite <- grep(favorite_gene[i], res$symbol)
            
            if(res$padj[my_favorite] < pval_cutoff){
                  text(res$baseMean[my_favorite], 
                       res$log2FoldChange[my_favorite],
                       labels = res$symbol[my_favorite], adj = c(0,-0.5), 
                       col=rgb(0.8,0,0,1))
            } else {
                  text(res$baseMean[my_favorite], 
                       res$log2FoldChange[my_favorite],
                       labels = res$symbol[my_favorite], adj = c(0,-0.5), 
                       col="grey")
            }
      }
      
      
      
      abline(h=0, col="grey32")
      
      legend("bottomright", legend =  c(paste("up", sum(res.sign$log2FoldChange > 0), sep = " - "),
                                        paste("down", sum(res.sign$log2FoldChange < 0), sep = " - ")) )
      
      
      my_main <- gsub(" ","", my_main)
      
      write.csv(res[order(res$padj),], file = paste("res", my_main, "csv", sep="."))
      #write.csv(res.sign[order(res.sign$padj),], file = paste("res.sign", my_main, "csv", sep="."))
      
      
}









######################################################## log2FC vs log2FC ######################################################## 





plotLog2FC <- function(dds, 
                       contrast1,
                       contrast2,
                       favorite_genes1,
                       favorite_genes2,
                       pval_cutoff,
                       lfc_cutoff = 0,
                       lims = c(-5,5),
                       selection_color1 = rgb(0,0,0.9,1),
                       selection_color2 = rgb(0.8,0,0,1),
                       selection_name1 = "",
                       selection_name2 = "",
                       text_label1 = FALSE,
                       text_label2 = FALSE,
                       show_legend = FALSE,
                       show_corr = TRUE){
      
      
      ######
      
      res1 <- results(dds, 
                     contrast = contrast1, 
                     lfcThreshold = lfc_cutoff, independentFiltering = FALSE
                     #, alpha = pval_cutoff
      )
      
      #res1 <- lfcShrink(dds, contrast = contrast1, res1 = res1,)
      
      res1$chr <- select(txdb, rownames(res1), "TXCHROM", keytype="GENEID", multiVals="first")$TXCHROM
      res1$symbol <- mapIds(org.Dm.eg.db, rownames(res1), "SYMBOL", keytype="FLYBASE", multiVals="first")
      
      res1$padj[is.na(res1$padj)] <- 1
      res1$log2FoldChange[is.na(res1$log2FoldChange)] <- 0
      res1.sign <- subset(res1, padj < pval_cutoff)
      
      
      my_label1 <- paste("log2FC", paste(contrast1[2:3], collapse = " - "))
      my_label1 <- gsub("plus", "+", my_label1)
      my_label1 <- gsub("minus", "-", my_label1)
      my_label1 <- gsub("\\.", " ", my_label1)
      my_label1 <- gsub("\\_", "", my_label1)
      
      #######
      
      res2 <- results(dds, 
                      contrast = contrast2, 
                      lfcThreshold = lfc_cutoff, independentFiltering = FALSE
                      #, alpha = pval_cutoff
      )
      
      #res2 <- lfcShrink(dds, contrast = contrast2, res2 = res2,)
      
      res2$chr <- select(txdb, rownames(res2), "TXCHROM", keytype="GENEID", multiVals="first")$TXCHROM
      res2$symbol <- mapIds(org.Dm.eg.db, rownames(res2), "SYMBOL", keytype="FLYBASE", multiVals="first")
      
      res2$padj[is.na(res2$padj)] <- 1
      res2$log2FoldChange[is.na(res2$log2FoldChange)] <- 0
      res2.sign <- subset(res2, padj < pval_cutoff)
      
      
      my_label2 <- paste("log2FC", paste(contrast2[2:3], collapse = " - "))
      my_label2 <- gsub("plus", "+", my_label2)
      my_label2 <- gsub("minus", "-", my_label2)
      my_label2 <- gsub("\\.", " ", my_label2)
      my_label2 <- gsub("\\_", "", my_label2)
      
      #####
      
      stopifnot(
            identical(rownames(res1), rownames(res2))
      )
      
      #####
      
      if(show_corr){
            my_main <- paste("Spearman´s r = ", round(cor(res1$log2FoldChange,res2$log2FoldChange, method = "spearman"), 2))
      } else {
            my_main <- ""
      }
      
      plot(res1$log2FoldChange,
           res2$log2FoldChange,
           main = my_main, 
           xlab = my_label1, 
           ylab = my_label2,
           col=rgb(0,0,0,0.1),  ylim=lims, xlim=lims, pch=19, cex = 0.25)
      
      
      # for(i in seq_along(res1$symbol)){
      #       
      #       my_favorite <- i
      #       
      #       stopifnot(
      #             identical(res1$symbol[my_favorite], res2$symbol[my_favorite])
      #       )
      #       
      #       if(res1$padj[my_favorite] < pval_cutoff & res2$padj[my_favorite] < pval_cutoff){
      #             points(res1$log2FoldChange[my_favorite], 
      #                  res2$log2FoldChange[my_favorite],
      #                  col=rgb(0.8,0,0,0.5), pch=19)
      #       } else if(res1$padj[my_favorite] < pval_cutoff | res2$padj[my_favorite] < pval_cutoff){
      #             points(res1$log2FoldChange[my_favorite], 
      #                  res2$log2FoldChange[my_favorite],
      #                  col=rgb(0.8,0.5,0,0.5), pch=19)
      #       } else {
      #             points(res1$log2FoldChange[my_favorite], 
      #                  res2$log2FoldChange[my_favorite],
      #                  col=rgb(0,0,0,0), pch=19)
      #       }
      # }
      # 
            
      
      for(i in seq_along(favorite_genes1)){
            
            my_favorite <- favorite_genes1[i] == res1$symbol
            
            if(length(my_favorite) == 0){
                  next()
            }
            
            stopifnot(
                  identical(res1$symbol[my_favorite], res2$symbol[my_favorite])
            )
            
            
            #if(res1$padj[my_favorite] < pval_cutoff & res2$padj[my_favorite] < pval_cutoff){
                  
                  points(res1$log2FoldChange[my_favorite], 
                       res2$log2FoldChange[my_favorite],
                       col=selection_color1, pch=19, cex = 0.5)

                  if(text_label1){
                        text(res1$log2FoldChange[my_favorite],
                             res2$log2FoldChange[my_favorite],
                             labels = res1$symbol[my_favorite], adj = c(0,-0.5),
                             col=selection_color1)
                  }                  
                  
            # } else if(res1$padj[my_favorite] < pval_cutoff | res2$padj[my_favorite] < pval_cutoff){
            #       
            #       points(res1$log2FoldChange[my_favorite], 
            #            res2$log2FoldChange[my_favorite],
            #            col=rgb(0.8,0.5,0,1), pch =19)
            #       
            #       text(res1$log2FoldChange[my_favorite], 
            #            res2$log2FoldChange[my_favorite],
            #            labels = res1$symbol[my_favorite], adj = c(0,1), 
            #            col=rgb(0.8,0.5,0,1))
            # } else {
            #       text(res1$log2FoldChange[my_favorite], 
            #            res2$log2FoldChange[my_favorite],
            #            labels = res1$symbol[my_favorite], adj = c(0,1), 
            #            col="grey")
            # }
      }
      
      
      
      
      for(i in seq_along(favorite_genes2)){
            
            my_favorite <- favorite_genes2[i] == res1$symbol
            
            if(length(my_favorite) == 0){
                  next()
            }
            
            stopifnot(
                  identical(res1$symbol[my_favorite], res2$symbol[my_favorite])
            )
            
            
            #if(res1$padj[my_favorite] < pval_cutoff & res2$padj[my_favorite] < pval_cutoff){
            
            points(res1$log2FoldChange[my_favorite], 
                   res2$log2FoldChange[my_favorite],
                   col=selection_color2, pch=19, cex = 0.5)
            
            if(text_label2){
                  text(res1$log2FoldChange[my_favorite],
                       res2$log2FoldChange[my_favorite],
                       labels = res1$symbol[my_favorite], adj = c(0,-0.5),
                       col=selection_color2)
            }                  
            
            # } else if(res1$padj[my_favorite] < pval_cutoff | res2$padj[my_favorite] < pval_cutoff){
            #       
            #       points(res1$log2FoldChange[my_favorite], 
            #            res2$log2FoldChange[my_favorite],
            #            col=rgb(0.8,0.5,0,1), pch =19)
            #       
            #       text(res1$log2FoldChange[my_favorite], 
            #            res2$log2FoldChange[my_favorite],
            #            labels = res1$symbol[my_favorite], adj = c(0,1), 
            #            col=rgb(0.8,0.5,0,1))
            # } else {
            #       text(res1$log2FoldChange[my_favorite], 
            #            res2$log2FoldChange[my_favorite],
            #            labels = res1$symbol[my_favorite], adj = c(0,1), 
            #            col="grey")
            # }
      }
      
      abline(h=0, v=0, col="grey32")
      abline(coef = c(0,1), col="grey32", lty=2)
      
      
      if(show_legend){
            legend("topleft", legend = c(selection_name1, selection_name2), bg = "white",
                   #legend =  c("sign. in both", "sign. in one"),
                   col = c(selection_color1, selection_color2), pch = 19, cex=0.8)      
      }
            
      
      
      
      
      

}




######################################################## Kmer log2FC vs log2FC ######################################################## 





plotKmerLog2FC <- function(dds, 
                       contrast1,
                       contrast2,
                       favorite_gene,
                       pval_cutoff,
                       lfc_cutoff = 0,
                       lims = c(-5,5),
                       selection_name = "",
                       text_label = FALSE){
      
      
      ######
      
      res1 <- results(dds, 
                      contrast = contrast1, 
                      lfcThreshold = lfc_cutoff, independentFiltering = FALSE
                      #, alpha = pval_cutoff
      )
      
      #res1 <- lfcShrink(dds, contrast = contrast1, res1 = res1,)
      
      res1$symbol <- rownames(res1)
      
      res1$padj[is.na(res1$padj)] <- 1
      res1$log2FoldChange[is.na(res1$log2FoldChange)] <- 0
      res1.sign <- subset(res1, padj < pval_cutoff)
      
      
      my_label1 <- paste("log2FC", paste(contrast1[2:3], collapse = " - "))
      my_label1 <- gsub("plus", "+", my_label1)
      my_label1 <- gsub("minus", "-", my_label1)
      my_label1 <- gsub("\\.", " ", my_label1)
      my_label1 <- gsub("\\_", "", my_label1)
      
      #######
      
      res2 <- results(dds, 
                      contrast = contrast2, 
                      lfcThreshold = lfc_cutoff, independentFiltering = FALSE
                      #, alpha = pval_cutoff
      )
      
      #res2 <- lfcShrink(dds, contrast = contrast2, res2 = res2,)
      
      res2$symbol <-rownames(res2)
      
      res2$padj[is.na(res2$padj)] <- 1
      res2$log2FoldChange[is.na(res2$log2FoldChange)] <- 0
      res2.sign <- subset(res2, padj < pval_cutoff)
      
      
      my_label2 <- paste("log2FC", paste(contrast2[2:3], collapse = " - "))
      my_label2 <- gsub("plus", "+", my_label2)
      my_label2 <- gsub("minus", "-", my_label2)
      my_label2 <- gsub("\\.", " ", my_label2)
      my_label2 <- gsub("\\_", "", my_label2)
      
      
      #####
      
      stopifnot(
            identical(rownames(res1), rownames(res2))
      )
      
      #####
      
      plot(res1$log2FoldChange,
           res2$log2FoldChange,
           main = "Comparison", 
           xlab = my_label1, 
           ylab = my_label2,
           col=rgb(0,0,0,0.1),  ylim=lims, xlim=lims, pch=19, cex = 0.5)
      
      
      # for(i in seq_along(res1$symbol)){
      #       
      #       my_favorite <- i
      #       
      #       stopifnot(
      #             identical(res1$symbol[my_favorite], res2$symbol[my_favorite])
      #       )
      #       
      #       if(res1$padj[my_favorite] < pval_cutoff & res2$padj[my_favorite] < pval_cutoff){
      #             points(res1$log2FoldChange[my_favorite], 
      #                  res2$log2FoldChange[my_favorite],
      #                  col=rgb(0.8,0,0,0.5), pch=19)
      #       } else if(res1$padj[my_favorite] < pval_cutoff | res2$padj[my_favorite] < pval_cutoff){
      #             points(res1$log2FoldChange[my_favorite], 
      #                  res2$log2FoldChange[my_favorite],
      #                  col=rgb(0.8,0.5,0,0.5), pch=19)
      #       } else {
      #             points(res1$log2FoldChange[my_favorite], 
      #                  res2$log2FoldChange[my_favorite],
      #                  col=rgb(0,0,0,0), pch=19)
      #       }
      # }
      # 
      
      
      for(i in seq_along(favorite_gene)){
            
            my_favorite <- favorite_gene[i] == res1$symbol
            
            if(sum(my_favorite) == 0){
                  next()
            }
            
            stopifnot(
                  identical(res1$symbol[my_favorite], res2$symbol[my_favorite])
            )
            
            
            #if(res1$padj[my_favorite] < pval_cutoff & res2$padj[my_favorite] < pval_cutoff){
            
            points(res1$log2FoldChange[my_favorite], 
                   res2$log2FoldChange[my_favorite],
                   col=rgb(0.8,0,0,1), pch=19, cex = 0.75)
            
            if(text_label){
                  text(res1$log2FoldChange[my_favorite],
                       res2$log2FoldChange[my_favorite],
                       labels = res1$symbol[my_favorite], adj = c(0,-0.5),
                       col=rgb(0.8,0,0,1))
            }                  
            
            # } else if(res1$padj[my_favorite] < pval_cutoff | res2$padj[my_favorite] < pval_cutoff){
            #       
            #       points(res1$log2FoldChange[my_favorite], 
            #            res2$log2FoldChange[my_favorite],
            #            col=rgb(0.8,0.5,0,1), pch =19)
            #       
            #       text(res1$log2FoldChange[my_favorite], 
            #            res2$log2FoldChange[my_favorite],
            #            labels = res1$symbol[my_favorite], adj = c(0,1), 
            #            col=rgb(0.8,0.5,0,1))
            # } else {
            #       text(res1$log2FoldChange[my_favorite], 
            #            res2$log2FoldChange[my_favorite],
            #            labels = res1$symbol[my_favorite], adj = c(0,1), 
            #            col="grey")
            # }
      }
      
      abline(h=0, v=0, col="grey32")
      
      legend("topleft", legend = selection_name, bg = "white",
             #legend =  c("sign. in both", "sign. in one"),
             col = c(rgb(0.8,0,0,1)), pch = 19)
      
      
      
      
      
}





######################################################## log2FC vs oligo frequency ######################################################## 






plotOligoFreq <- function(dds, 
                       contrast1,
                       freq_ranges,
                       oligo = "TTTT",
                       favorite_genes1,
                       favorite_genes2,
                       pval_cutoff,
                       lfc_cutoff = 0,
                       lims = c(-5,5),
                       selection_color1 = rgb(0,0,0.9,1),
                       selection_color2 = rgb(0.8,0,0,1),
                       selection_name1 = "",
                       selection_name2 = "",
                       text_label1 = FALSE,
                       text_label2 = FALSE,
                       show_cor = FALSE){
      
      
      ######
      
      res1 <- results(dds, 
                      contrast = contrast1, 
                      lfcThreshold = lfc_cutoff, independentFiltering = FALSE
                      #, alpha = pval_cutoff
      )
      
      #res1 <- lfcShrink(dds, contrast = contrast1, res1 = res1,)
      
      #res1$chr <- select(txdb, rownames(res1), "TXCHROM", keytype="GENEID", multiVals="first")$TXCHROM
      res1$symbol <- mapIds(org.Dm.eg.db, rownames(res1), "SYMBOL", keytype="FLYBASE", multiVals="first")
      res1 <- res1[!(is.na(res1$symbol)),]
      
      res1$padj[is.na(res1$padj)] <- 1
      res1$log2FoldChange[is.na(res1$log2FoldChange)] <- 0
      res1.sign <- subset(res1, padj < pval_cutoff)
      
      my_label1 <- paste("log2FC", paste(contrast1[2:3], collapse = " - "))
      my_label1 <- gsub("plus", "+", my_label1)
      my_label1 <- gsub("minus", "-", my_label1)
      my_label1 <- gsub("\\.", " ", my_label1)
      my_label1 <- gsub("\\_", "", my_label1)
      
      #######
      
      
      res_merged <- merge(as.data.frame(res1), 
                          as.data.frame(freq_ranges), by.x = "row.names", by.y = "gene_id")
      
      #######
      
      plot(res_merged$log2FoldChange,
           res_merged[oligo][,1],
           main = "",
           xlim = lims,
           xlab = my_label1, 
           ylab = paste("Frequeny", oligo),
           col=rgb(0,0,0,0.1), pch=19, cex = 0.25)
      
      if(show_cor){
            title(paste("cor =", round(cor( res_merged$log2FoldChange,res_merged[oligo][,1], method = "pearson"),2)), line = 0.5, cex.main=1 )
            
      }
      
      #####
      
      
      
      
      for(i in seq_along(favorite_genes1)){
            
            my_favorite <- favorite_genes1[i] == res_merged$symbol
            
            if(sum(my_favorite) == 0){
                  next()
            }

            points(res_merged$log2FoldChange[my_favorite],
                   res_merged[oligo][,1][my_favorite],
                   col=selection_color1, pch=19, cex = 0.5)
            
            if(text_label1){
                  text(res_merged$log2FoldChange[my_favorite],
                       res_merged[oligo][,1][my_favorite],
                       labels = res_merged$symbol[my_favorite], adj = c(0,-0.5),
                       col=selection_color1)
            }                  
            

      }
      
      for(i in seq_along(favorite_genes2)){
            
            my_favorite <- favorite_genes2[i] == res_merged$symbol
            
            if(sum(my_favorite) == 0){
                  next()
            }
            
            points(res_merged$log2FoldChange[my_favorite],
                   res_merged[oligo][,1][my_favorite],
                   col=selection_color2, pch=19, cex = 0.5)
            
            if(text_label2){
                  text(res_merged$log2FoldChange[my_favorite],
                       res_merged[oligo][,1][my_favorite],
                       labels = res_merged$symbol[my_favorite], adj = c(0,-0.5),
                       col=selection_color2)
            }                  
            
            
      }
      abline(v=0, col="grey32")
      
      legend("topleft", legend = c(selection_name1, selection_name2), bg = "white",
             #legend =  c("sign. in both", "sign. in one"),
             col = c(selection_color1, selection_color2), pch = 19, cex=0.8)
      
      
      # if(oligo == "TTTT"){
      #       
      #       my_main <- paste(contrast1[2:3], collapse = "-")
      #       my_main <- gsub("plus", "+", my_main)
      #       my_main <- gsub("minus", "-", my_main)
      #       
      #       write.csv(res_merged[order(res_merged$TTTT, decreasing = TRUE),], file = paste("res", my_main, "freqs", "csv", sep="."))
      # 
      # }
  
}







######################################################## Heatmap ######################################################## 






createHeatmap <- function(my_res_name, data_counts, enriched = TRUE, n_genes = 200, favorite_genes = NULL, horizontal = FALSE){

      
      my_res <- get(my_res_name)
      
      if(nrow(my_res) == 1){next()}
      
      my_res <- my_res[order(my_res$padj),]

      if(enriched){
            my_res <- my_res[my_res$log2FoldChange > 0,]
      } else {
            my_res <- my_res[my_res$log2FoldChange < 0,]
      }
            
      
      if(!(is.null(favorite_genes))){
            my_res <- my_res[rownames(my_res) %in% favorite_genes,]
      }
      
      if(nrow(my_res) > n_genes){
            my_res <- my_res[1:n_genes,]
      }
      
      ############################
      
      mat <- data_counts[rownames(data_counts) %in% rownames(my_res),]
      mat <- mat - rowMeans(mat)
      my_gene_id <- rownames(mat)
      rownames(mat) <- mapIds(org.Dm.eg.db, rownames(mat), "SYMBOL", keytype="FLYBASE", multiVals="first")
      rownames(mat)[is.na(rownames(mat))] <-  my_gene_id[is.na(rownames(mat))] 
      
      ############################      
      
      my_mat_size <- nrow(mat)
      my_ht_size <- (1/(my_mat_size/100))*4.5
      if(my_ht_size > 12){my_ht_size <- 12}
      
      ############################
      
      if(horizontal){
            mat <- t(mat)
      }
      
      hm = pheatmap(mat,
                     cellwidth = my_ht_size, cellheight = my_ht_size, fontsize = 8,
                     fontsize_row = my_ht_size-1,  fontsize_col = my_ht_size-1, 
                     cluster_cols = TRUE, cluster_rows = TRUE,
                     border_color = NA, silent = TRUE, 
                     color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdBu")))(100),
                     breaks = seq(-1*max(abs(mat)), max(abs(mat)), length.out = 101))
      
      return(hm)
}








############################################################################################################################
############################################################################################################################
############################################################################################################################


#####################################################      PCA      ######################################################## 






plottingPCA <- function(my_data, 
                        xcomp = 1,
                        ycomp = 2,
                        color_palette,
                        conditions,
                        quantiles = c(0,1),
                        show_labels = TRUE,
                        point_size = 1.5,
                        my_limits = c(-100,100)){
      
      rv <- rowVars(my_data)
      
      selection <- (rv >  quantile(rv, quantiles[1])  & rv < quantile(rv, quantiles[2]))
      
      
      pca <- prcomp(t(my_data[selection, ]), scale. = TRUE)
      
      percentVar <- round(pca$sdev^2/sum(pca$sdev^2)*100,1)[1:10]
      
      
      
      plot(pca$x[, xcomp], pca$x[, ycomp]*-1, 
           col = color_palette[conditions], 
           pch=16, cex = point_size,
           xlab = paste("PC",xcomp," (", percentVar[xcomp], "%)", sep=""),
           ylab = paste("PC",ycomp," (", percentVar[ycomp], "%)", sep=""),
           xlim= my_limits, ylim=my_limits)
      
      title(main = "PCA", line = 0.5)
      
      
      if(show_labels){
            text(pca$x[, xcomp], pca$x[, ycomp]*-1, labels = rownames(pca$x), 
                 adj = -0.5, col = "gray32", cex=0.5)      
      }
      
      
      
}



############################################################################################################################
############################################################################################################################
############################################################################################################################

###################################################       dot plots      ################################################### 



plotDots <- function(my_data, 
                     my_title,
                     color_palette,
                     color_groups,
                     conditions,
                     point_size = 1.5,
                     ylims = c(0, 15),
                     xlims = NULL,
                     x_label = NULL){
      
      
      
      grouped_data = data.frame(exp = as.numeric(my_data),
                                group =  as.numeric(conditions))
      
      plot(grouped_data$group + runif(length(grouped_data$group), -0.25, 0.25),
           grouped_data$exp,
           pch=19, cex = point_size, 
           ylim =  ylims, xlim = xlims,
           xaxt = "n", xlab = "", ylab = "",
           col =  color_palette[color_groups])
      
      title(main = my_title, line = 0.5)
      
      if(!(is.null(x_label))){
            axis(side = 1, at = seq_along(unique(grouped_data$group)), labels = x_label, las=2)
            
      } else {
            axis(side = 1, labels = FALSE)
            
      }
      
      
}








############################################################################################################################
############################################################################################################################
############################################################################################################################

###################################################       Heatmaps       ################################################### 


plotHeatmap <- function(my_mat,
                        my_row_order,
                        my_col_order,
                        min_value = 0,
                        max_value = 15,
                        my_title,
                        my_color_palette,
                        show_xaxis = FALSE,
                        useRaster = TRUE,
                        cex.main = 1.25){
      
      my_mat[my_mat < min_value] <- min_value
      my_mat[my_mat > max_value] <- max_value
      
      image(t(my_mat[my_row_order, my_col_order]), 
            #main = my_title, cex.main = cex.main,
            col = my_color_palette, 
            breaks =  seq(min_value, max_value, length.out = 101),
            axes=FALSE, useRaster = useRaster)
      
      title(main = my_title, line = 1.5, cex.main = cex.main)
      
      if(show_xaxis){
            axis(side = 1, at = c(0,0.5,1), lwd = 0, lwd.ticks = 1, las=1,tck = -0.15, cex.axis = 0.8,
                 labels = c( paste(round(as.integer(colnames(my_mat))[1] / 1000), "kb", sep=""), "",
                             paste("+",round(as.integer(colnames(my_mat))[ncol(my_mat)] / 1000), "kb", sep="")
                 ))
      }
      
      
}





plotHeatmapKey <- function(my_mat,
                           min_value = 0,
                           max_value = 15,
                           my_title,
                           my_color_palette,
                           show_yaxis = FALSE,
                           show_xaxis = FALSE,
                           useRaster = TRUE,
                           cex.main = 1.25,
                           cex.axis = 0.6){
      
      my_mat[my_mat < min_value] <- min_value
      my_mat[my_mat > max_value] <- max_value
      
      image((my_mat), 
            main = my_title, cex.main = cex.main,
            col = my_color_palette, 
            breaks =  seq(min_value, max_value, length.out = 101),
            axes=FALSE, useRaster = useRaster)
      
      if(show_yaxis){
            axis(side = 4, at = c(0,1), lwd = 0, lwd.ticks = 1, las=1, cex.axis = cex.axis, tck = -0.25,
                 labels = round(c(min_value,max_value)))
      }
      
      if(show_xaxis){
            axis(side = 1, at = c(0,1), lwd = 0, lwd.ticks = 1, las=1, cex.axis = cex.axis, tck = -0.25,
                 labels = round(c(min_value,max_value)))
      }
      
      
}



############################################################################################################################
############################################################################################################################
############################################################################################################################






############################################################################################################################
############################################################################################################################
############################################################################################################################







setupDDS <- function(SampleTableName = "SampleTable",
                     CountTableName = "dmel.counts_genes",
                     SampleIdName = "SampleID",
                     ConditionsName = "Conditions",
                     BatchName = "Replicate",
                     n_samples_for_filtering = 3
                     ){
      
      SampleTable <- get(SampleTableName)
      
      my_counts_genes <- get(CountTableName)
      
      stopifnot(
            identical(colnames(my_counts_genes), as.character(SampleTable[,SampleIdName]))
      )
      
      
      
      ######################################################## 
      
      
      
      filter <- apply(my_counts_genes, 1, function(x) length(x[x>1]) >= n_samples_for_filtering)
      
      my_counts_filtered <- my_counts_genes[filter,]
      
      
      
      ######################################################## 
      
      
      my_colData <- DataFrame(Sample = SampleTable[,ConditionsName],
                              Batch = factor(SampleTable[,BatchName]))
      
      rownames(my_colData) <- SampleTable[,SampleIdName]
      
      
      ######################################################## 
      
      dds <- DESeqDataSetFromMatrix(countData = my_counts_filtered, 
                                     colData = my_colData, 
                                     design = ~Batch+Sample)
      
      
      
      dds <- DESeq(dds)
      
      return(dds)
}












############################################################################################################################
############################################################################################################################
############################################################################################################################






