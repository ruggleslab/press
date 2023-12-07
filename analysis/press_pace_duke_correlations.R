# HEADER --------------------------------------------------
# Author: Matthew Muller
# 
# Date: 2022-12-02
# 
# Script Name: Press Pace and Duke Longitudinal Groups Correlations
# 
# Notes:
#  This script makes a correlation plot for the genetic data between pace and duke

# SET WORKING DIRECTORY -----------------------------------
wd <- FALSE
if (wd != FALSE) {setwd(wd)}

# LOAD LIBRARIES ------------------------------------------
packages <- c(
  "tidyverse",
  "ggplot2",
  "devtools",
  "DESeq2",
  "reshape2",
  "scales",
  "corrplot",
  "cowplot"
)

for (pkg in packages) {
  paste0(pkg)[[1]]
  library(pkg, character.only = T, quietly = T)
}

source_url('https://raw.githubusercontent.com/mattmuller0/scripts/main/Rtools/general_functions.R')

###########################################################################
#
#                       Correlation Analysis
#
###########################################################################
par(mfrow=c(1,2)) 

# Load pace and duke, then subset based on press
pace_dge <- read.csv('data/hyper_v_hypo_deseqoutput.csv')
duke_dge <- read.csv('output/duke_dge/duke_longitudinal_dge_output.csv')
press_genes <- read.csv('data/clean/press_genes.csv', header = F)$V1

pace_dge <- pace_dge[pace_dge$X %in% press_genes,]
row.names(pace_dge) <- pace_dge$X
duke_dge <- duke_dge[duke_dge$X %in% press_genes,]
row.names(duke_dge) <- duke_dge$X
duke_dge <- add_missing_rows(duke_dge, pace_dge)

# Corrs
cor.test(pace_dge$log2FoldChange, duke_dge$log2FoldChange)


plot(pace_dge$log2FoldChange, duke_dge$log2FoldChange, pch=19, col=2, cex=0.75, ylim=c(-1.5,1),
     main = 'Pace and Duke DGE Correlations \nDuke Longitudinal Group1') + 
  abline(h=0, v=0, lty=3) + abline(lm(pace_dge$log2FoldChange ~ duke_dge$log2FoldChange), col = "red")



# Load pace and duke, then subset based on press
pace_dge <- read.csv('data/hyper_v_hypo_deseqoutput.csv')
duke_dge <- read.csv('data/duke_validation_run3/dge_analysis/comp_group1__hyper_v_nothyper_AGCONTROL_deseqout.csv')
press_genes <- read.csv('data/clean/press_genes.csv', header = F)$V1

pace_dge <- pace_dge[pace_dge$X %in% press_genes,]
row.names(pace_dge) <- pace_dge$X
duke_dge <- duke_dge[duke_dge$X %in% press_genes,]
row.names(duke_dge) <- duke_dge$X
duke_dge <- add_missing_rows(duke_dge, pace_dge)

# Corrs
cor.test(pace_dge$log2FoldChange, duke_dge$log2FoldChange)


plot(pace_dge$log2FoldChange, duke_dge$log2FoldChange, pch=19, col=2, cex=0.75, ylim=c(-1.5,1),
     main = 'Pace and Duke DGE Correlations \nDuke Group 1') + 
  abline(h=0, v=0, lty=3) + abline(lm(pace_dge$log2FoldChange ~ duke_dge$log2FoldChange), col = "red")
