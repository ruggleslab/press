################################
# Name: MacIntosh Cornwell
# Email: mcornwell1957@gmail.com
################################
## Master Script for platelet pace data

## Library load
packages = c("ggplot2", "reshape2", "DESeq2", "grid", "gridExtra", "scales", "ggrepel", "tools", "psych", "Hmisc", "devtools", "VennDiagram")
pkgs <- lapply(packages, function(x) suppressMessages(require(x, character.only=T,quietly=T))) # nolint

# print any packages that failed to load
print(packages[!sapply(pkgs, isTRUE)])

source("src/tosh_rnaseq_scripts/deseq_functions.R")
source("src/tosh_rnaseq_scripts/overlap_finder_function.R")
source("src/tosh_rnaseq_scripts/mgc_plotting_functions.R")
source("src/tosh_rnaseq_scripts/geneset_analysis_functions.R")

## Outfolder
runnumber = 18
hypervalue = 60
hypovalue = 40
outfilepathmaster <- paste0("output/hyper_geneset_creation/run", runnumber, "_hyper", hypervalue, "_hypo", hypovalue, "_AGRCONTROL", "/")
dir.create(outfilepathmaster, showWarnings = FALSE, recursive = TRUE)

relative_path <- "/Users/muller/Library/CloudStorage/GoogleDrive-mm12865@nyu.edu/My Drive/RugglesLab/projects/platelet-pace"

## Save and load data
# save.image(file = paste0(outfilepathmaster, "/platelet_hyper_testing.RData"))
# load(file = paste0(outfilepathmaster, "/platelet_hyper_testing.RData"))

# make some notes on the run
notes <- "This is the final run using press_1"

# save notes
write(notes, file = paste0(outfilepathmaster, "run_notes.txt"))

# --------------------------------- Read in files ---------------------------------

# hypersamples <- c("PACE_202","PACE_148","PACE_248","PACE_076","PACE_095","PACE_212","PACE_005","PACE_109","PACE_167","PACE_282","PACE_157","PACE_178")
# nothypersamples <- c("PACE_024","PACE_114","PACE_136","PACE_080","PACE_065","PACE_019","PACE_007","PACE_186","PACE_143","PACE_263","PACE_241","PACE_014","PACE_022","PACE_113","PACE_185","PACE_246","PACE_094","PACE_147","PACE_184","PACE_230","PACE_199","PACE_188","PACE_088","PACE_092","PACE_176","PACE_071","PACE_106","PACE_161","PACE_093","PACE_125","PACE_227","PACE_269","PACE_029","PACE_235","PACE_165","PACE_289","PACE_104","PACE_155","PACE_032","PACE_020","PACE_036","PACE_067","PACE_121","PACE_236","PACE_128","PACE_284","PACE_051","PACE_129","PACE_300","PACE_073","PACE_096","PACE_013","PACE_059","PACE_069","PACE_146")

# Countfile
normcounttable_file <- "output/run6_rmoutliers2_agesexcontrol_withpaired_20201007/rna_processing/normcounttab.txt"
normcounttable <- read.table(normcounttable_file, sep = "\t", header = TRUE, row.names = 1, check.names = FALSE)
rawcounttable_file <- "output/run6_rmoutliers2_agesexcontrol_withpaired_20201007/rna_processing/filtrawcounttab.txt"
rawcounttable <- read.table(rawcounttable_file, sep = "\t", header = TRUE, row.names = 1)

## Add minimum filtering:
# minreadfilter <- 100
# normcounttable <- normcounttable[rowMeans(normcounttable) >= minreadfilter,]
# rawcounttable <- rawcounttable[rowMeans(rawcounttable) >= minreadfilter,]

## File loading 
metafile = "data/metadata_7_20210810_expanded.txt"
metatable = read.table(metafile, sep = "\t", header=TRUE, stringsAsFactors = FALSE, na.strings = c(NA, "NA", ""),check.names = FALSE)
rownames(metatable) <- metatable[,"FINAL_PACE_ID"]

# --------------------------------- Define hyper and normo PACE samples and define PRESS sets (mgc_custom) ---------------------------------

## Select our hyper hypo groups
epi_metatable <- metatable[,c("FINAL_PACE_ID", "epi_04um_300s_n")]
epi_metatable[,"hypercohort"] <- NA
epi_metatable[epi_metatable[,"epi_04um_300s_n"] >= hypervalue & !is.na(epi_metatable[,"epi_04um_300s_n"]),"hypercohort"] <- "hyper"
epi_metatable[epi_metatable[,"epi_04um_300s_n"] < hypovalue & !is.na(epi_metatable[,"epi_04um_300s_n"]),"hypercohort"] <- "nothyper"
rownames(epi_metatable) <- epi_metatable[,"FINAL_PACE_ID"]

hypersamples <- epi_metatable[epi_metatable[,"hypercohort"] == "hyper" & !is.na(epi_metatable[,"hypercohort"]),"FINAL_PACE_ID"]
nothypersamples <- epi_metatable[epi_metatable[,"hypercohort"] == "nothyper" & !is.na(epi_metatable[,"hypercohort"]),"FINAL_PACE_ID"]
hypersamples_rna <- hypersamples[hypersamples %in% colnames(normcounttable)]
nothypersamples_rna <- nothypersamples[nothypersamples %in% colnames(normcounttable)]

## Adding additional criteria to ONLY TAKE PATIENTS ON APS
APsamples <- metatable[metatable[,"antiplatelet_therapy"] %in% "1:Yes","FINAL_PACE_ID"]
hypersamples <- hypersamples_rna[hypersamples_rna %in% APsamples]
nothypersamples <- nothypersamples_rna[nothypersamples_rna %in% APsamples]

## Write this out so we can track it.
epi_metatable_writeout <- merge(epi_metatable, metatable[,"antiplatelet_therapy",drop=FALSE], 
                                by.x = "FINAL_PACE_ID", by.y = "row.names", all.x = TRUE)
rownames(epi_metatable_writeout) <- epi_metatable_writeout[,"FINAL_PACE_ID"]
epi_metatable_writeout[,"hypercohort_inrnaseq"] <- epi_metatable_writeout[,"hypercohort"]
epi_metatable_writeout[!rownames(epi_metatable_writeout) %in% colnames(normcounttable), "hypercohort_inrnaseq"] <- NA
epi_metatable_writeout[,"hypercohort_inrnaseq_AP"] <- epi_metatable_writeout[,"hypercohort_inrnaseq"]
epi_metatable_writeout[!rownames(epi_metatable_writeout) %in% APsamples, "hypercohort_inrnaseq_AP"] <- NA
write.table(epi_metatable_writeout, paste0(outfilepathmaster, "hypercohort_metatable.csv"), sep = ",", col.names = NA, row.names = TRUE)

## PLot histogram of samples by metric - only those on AP though
pout <- plot_histogram(epi_metatable[rownames(epi_metatable) %in% APsamples,"epi_04um_300s_n",drop=FALSE], limity = c(0,5), binparam = 94,
                       labsparam = list(title = "Platelet Aggregation in Response to 0.4uM Epinephrine at 300s", x = "Value", y = "Number of PACE Patients"))
pout <- pout + geom_vline(xintercept = c(40,60), linetype = 2, color = "red")
pout <- pout + scale_y_continuous(breaks = seq(0,10,2))
pdf(paste0(outfilepathmaster, "metric_APonly_histogram.pdf"))
print(pout)
junk <- dev.off()

## Lets do deseq, t.test, wilcox to see what genes come out from each test
# 1 t.test and wilcox test
testout1 <- apply(normcounttable, 1, function(x) {
  group1 <- unlist(x[hypersamples])
  group2 <- unlist(x[nothypersamples])
  ttestout <- t.test(group1, group2)
  out1 <- c(ttestout[["estimate"]][1], ttestout[["estimate"]][2], ttestout[["p.value"]], log2(ttestout[["estimate"]][1]/ttestout[["estimate"]][2]))
  out2 <- c(wilcox.test(group1, group2)$p.value, cohens_d(group1, group2))
  c(out1, out2)
})
rownames(testout1) <- c("hypermean", "nothypermean", "ttest_pval", "log2_hyper_v_nothyper", "wilcox_pval", "cohens")
testout1 <- t(testout1)

# 2 DESeq
## Run all of our DESeq analyses based on the various comps I put in
compcols <- data.frame(comp_hyper__hyper_v_nothyper = c(rep(1, length(hypersamples)), rep(0, length(nothypersamples))), row.names = c(hypersamples, nothypersamples))
# Adding in age and sex control - because i think it gives the best results.............
controlcols <- metatable[c(hypersamples, nothypersamples),c("sex1", "age", "race1")]
controlcols[,"age_cat"] <- ifelse(controlcols[,"age"] >= 74, "74_98",
                           ifelse(controlcols[,"age"] >= 61, "61_73",
                           ifelse(controlcols[,"age"] >= 28, "28_60", NA)))
deseqout <- DESeq_analysis(compcols = compcols[c(hypersamples, nothypersamples),,drop=FALSE], controlcols = controlcols[,c("sex1", "age_cat", "race1")],
                           rawcounttab = rawcounttable[,c(hypersamples, nothypersamples)])
# deseqout <- DESeq_analysis(compcols = compcols[c(hypersamples, nothypersamples),,drop=FALSE], controlcols = NULL,
#                            rawcounttab = rawcounttable[,c(hypersamples, nothypersamples)])
testout2 <- deseqout[[1]]
write.table(testout2, paste0(outfilepathmaster, "hyper_v_hypo_deseqoutput.csv"), sep = ",", col.names = NA, row.names = TRUE)

## Combine all results
fulltestout1 <- merge(testout1, testout2[,c("pvalue", "log2FoldChange")], by = "row.names")
colnames(fulltestout1)[c(8,9)] <- c("deseq_pval", "deseq_log2")
# apply(fulltestout[,grepl("_pval", colnames(fulltestout))], 2, function(x) sum(x < 0.05))

## Now, lets see which genes ALSO correlate with Epi values
epiPACEIDs <- epi_metatable[!is.na(epi_metatable[,"epi_04um_300s_n"]), "FINAL_PACE_ID"]
epiPACEIDs <- epiPACEIDs[epiPACEIDs %in% c(hypersamples, nothypersamples)]
# epiPACEIDs <- epiPACEIDs[epiPACEIDs %in% intersect(APsamples, colnames(normcounttable))]
cortestout <- corr.test(epi_metatable[epiPACEIDs, "epi_04um_300s_n",drop=FALSE], t(normcounttable[,epiPACEIDs]), method = "spearman", adjust = "none")
testout3 <- merge(t(cortestout$p), t(cortestout$r), by = "row.names")
colnames(testout3) <- c("Row.names", "epi_04_um_300s_pval", "epi_04um_300s_spearman")

## Merge this to the full table
fulltestout2 <- merge(fulltestout1, testout3, by = "Row.names")

# Lets do a correlation test for all AP samples
cortestout <- corr.test(epi_metatable[
    intersect(APsamples, colnames(normcounttable)), "epi_04um_300s_n",drop=FALSE], 
    t(normcounttable[,intersect(APsamples, colnames(normcounttable))]), 
    method = "spearman", adjust = "none")
epiPACEIDs <- epiPACEIDs[epiPACEIDs %in% intersect(APsamples, colnames(normcounttable))]
testout4 <- merge(t(cortestout$p), t(cortestout$r), by = "row.names")
colnames(testout4) <- c("Row.names", "epi_04_um_300s_HYPER_pval", "epi_04um_300s_HYPER_spearman")

# Merge this additional test now too
fulltestout3 <- merge(fulltestout2, testout4, by = "Row.names")

# So now whats the overlap of genes here:
fulltestoutfinal <- fulltestout3
rownames(fulltestoutfinal) <- fulltestoutfinal[,"Row.names"]
overlapinlist <- list(
                      # ttest_sig = fulltestoutfinal[fulltestoutfinal[,"ttest_pval"] < 0.05,"Row.names"], 
                      wilcox_sig = fulltestoutfinal[fulltestoutfinal[,"wilcox_pval"] < 0.05,"Row.names"],
                      deseq_sig = fulltestoutfinal[fulltestoutfinal[,"deseq_pval"] < 0.05 & !is.na(fulltestoutfinal[,"deseq_pval"]),"Row.names"],
                      # epicorr_sig = fulltestoutfinal[fulltestoutfinal[,"epi_04_um_300s_pval"] < 0.05,"Row.names"],
                      epicorrALL_sig = fulltestoutfinal[fulltestoutfinal[,"epi_04_um_300s_HYPER_pval"] < 0.05,"Row.names"]
                      )

overlapout <- overlap_finder(overlapinlist)
overlaptable <- overlapout$overlaptable
vennplot <- overlapout$vennplot
overlapgrouptab <- overlapout$overlapgrouptab
colnames(overlapgrouptab)[1] <- "Row.names"
overlapsummary <- overlapout$overlapsummary

write.table(overlaptable, paste0(outfilepathmaster, "testoverlap_overlaptable.csv"), sep = ",", col.names = TRUE, row.names = FALSE)
write.table(overlapgrouptab, paste0(outfilepathmaster, "testoverlap_overlapgrouptable.csv"), sep = ",", col.names = TRUE, row.names = FALSE)
write.table(overlapsummary, paste0(outfilepathmaster, "testoverlap_overlapsummary.csv"), sep = ",", col.names = TRUE, row.names = TRUE)
pdf(paste0(outfilepathmaster, "testoverlap_overlapvenn.pdf"))
grid.draw(vennplot)
junk <- dev.off()

## Get the up and down for the genes that passed our tests of means
# meanstestgenes <- overlapgrouptab[overlapgrouptab[,2] %in% c("ttest_sig;wilcox_sig;deseq_sig;epicorr_sig", "ttest_sig;wilcox_sig;deseq_sig"), 1]
meanstestgenes <- overlapgrouptab[overlapgrouptab[,2] %in% c("wilcox_sig;deseq_sig;epicorr_sig", "ttest_sig;wilcox_sig;deseq_sig", "wilcox_sig;deseq_sig;epicorr_sig;epicorrALL_sig", 'wilcox_sig;deseq_sig;epicorrALL_sig'), 1]
table(sign(testout2[meanstestgenes,"log2FoldChange"]))

## 330 is a pretty good number - now lets go back and make sure these all go the same direction (so all pos or all neg)
# GOI <- overlapgrouptab[overlapgrouptab[,"overlapgroup"] == "ttest_sig;wilcox_sig;deseq_sig;epicorr_sig", "Row.names"]
# GOI <- overlapgrouptab[overlapgrouptab[,"overlapgroup"] == "wilcox_sig;deseq_sig;epicorr_sig;epicorrALL_sig", "Row.names"]
GOI <- meanstestgenes
GOIuptable <- fulltestoutfinal[sign(fulltestoutfinal[,"log2_hyper_v_nothyper"]) == 1 & sign(fulltestoutfinal[,"epi_04um_300s_spearman"]) == 1 & rownames(fulltestoutfinal) %in% GOI, ]
GOIup <- rownames(GOIuptable)
GOIdowntable <- fulltestoutfinal[sign(fulltestoutfinal[,"log2_hyper_v_nothyper"]) == -1 & sign(fulltestoutfinal[,"epi_04um_300s_spearman"]) == -1 & rownames(fulltestoutfinal) %in% GOI, ]
GOIdown <- rownames(GOIdowntable)

## Lets write out some things now
# writeouttable <- merge(fulltestoutfinal, overlapgrouptab, by.x = "Row.names", by.y = "features", all.x = TRUE)
writeouttable = Reduce(function(dtf1, dtf2)
  merge(dtf1, dtf2, by = "Row.names", all.x = TRUE, sort = FALSE), 
  list(
    fulltestoutfinal, 
    overlapgrouptab, 
    data.frame(GOIup = rep("GOIup", length(GOIup)), Row.names = GOIup),
    data.frame(GOIdown = rep("GOIdown", length(GOIdown)), Row.names = GOIdown)
  ))

write.table(writeouttable, paste0(outfilepathmaster, "full_hypertest_stattable.csv"), sep = ",", col.names = TRUE, row.names = FALSE)
write.table(data.frame(writeouttable[!is.na(writeouttable[,"GOIup"]),"Row.names"]), paste0(outfilepathmaster, "custom_mgc_hyper_up.txt"), sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)
write.table(data.frame(writeouttable[!is.na(writeouttable[,"GOIdown"]),"Row.names"]), paste0(outfilepathmaster, "custom_mgc_hyper_down.txt"), quote = FALSE, sep = "\t", col.names = FALSE, row.names = FALSE)
write.table(data.frame(GOI), paste0(outfilepathmaster, "custom_mgc_all_GOI.txt"), quote = FALSE, sep = "\t", col.names = FALSE, row.names = FALSE)


# Writing out another table here - where we have our samples on one side, and epi level and genes on the other
feature_outtable <- merge(epi_metatable[c(hypersamples, nothypersamples), c("epi_04um_300s_n", "hypercohort")],
                          t(normcounttable[GOI, c(hypersamples, nothypersamples)]), by = "row.names")
rownames(feature_outtable) <- feature_outtable[,"Row.names"]
feature_outtable <- feature_outtable[,!grepl("Row.names", colnames(feature_outtable))]
write.table(feature_outtable, paste0(outfilepathmaster, "hyper_feature_outtable.csv"), sep = ",", col.names  = NA, row.names = TRUE)

# --------------------------------- Viz and test our new PRESS sets ---------------------------------

## Some additional illustrations of this:
dir.create(paste0(outfilepathmaster, "figures/"), showWarnings = FALSE, recursive = TRUE)
# 1 scatterplot of deseq p stat vs correlation p stat
# plottable1 <- cbind(fulltestoutfinal[,c("deseq_pval", "epi_04_um_300s_pval", "deseq_log2")], -log10(fulltestoutfinal[,c("deseq_pval", "epi_04_um_300s_pval")]))
# colnames(plottable1)[c(4,5)] <- c("deseq_pstat", "epi_pstat")
# plottable1[,"signif"] <- "notsig"
# plottable1[GOIup,"signif"] <- "sig_UP"
# plottable1[GOIdown,"signif"] <- "sig_DN"
# pout1 <- scatter_plotter(plottable1[,c("deseq_pstat", "epi_pstat")], colorvar = plottable1[,"signif",drop=FALSE], 
#                          labsparam = list(title = "deseq vs corr pstat", x = "deseq p stat", y = "corr p stat"), plotstats = FALSE)
# pout1 <- pout1 + scale_color_manual(breaks = c("notsig", "sig_UP", "sig_DN"), values = c("grey", "darkred", "darkblue"))
# pout1 <- pout1 + geom_hline(yintercept = -log10(0.05), linetype = 2) + geom_vline(xintercept = -log10(0.05), linetype = 2)
# pout1

## Need to get the pvalue of rval cutoff - so taking the avg of the two values will do a very accurate estimate (not exact, thats fine - purely viz)
rvalcutoff <- mean(min(fulltestoutfinal[fulltestoutfinal[,"epi_04_um_300s_HYPER_pval"] < 0.05 & 
                                        sign(fulltestoutfinal[,"epi_04um_300s_HYPER_spearman"]) == 1,"epi_04um_300s_spearman"]),
                   abs(max(fulltestoutfinal[fulltestoutfinal[,"epi_04_um_300s_HYPER_pval"] < 0.05 & 
                                            sign(fulltestoutfinal[,"epi_04um_300s_HYPER_spearman"]) == -1,"epi_04um_300s_HYPER_spearman"]))
)

rvalcutoff <- 0.1858201
stattype = "pvalue"
pvalcutoff = 0.05
log2fccutoff = rvalcutoff
deseqanalysislabel = "deseq v corr"
plottab = fulltestoutfinal[c("epi_04um_300s_HYPER_spearman", "deseq_pval")]

## PLOT CODE
plottab[is.na(plottab[,2]),2] <- 1
plottab[is.na(plottab[,1]),1] <- 0

plottab$color = "notsig"
plottab[plottab[,2] < pvalcutoff & plottab[,1,drop=FALSE] > 0, "color"] <- "sigpvalup" # pink
plottab[plottab[,1,drop=FALSE] > log2fccutoff, "color"] <- "siglog2fcup" # red
plottab[plottab[,1,drop=FALSE] > log2fccutoff & plottab[,2] < pvalcutoff, "color"] <- "sigbothup" # darkred
plottab[plottab[,2] < pvalcutoff & plottab[,1,drop=FALSE] < 0, "color"] <- "sigpvaldown" # light blue
plottab[plottab[,1,drop=FALSE] < -log2fccutoff, "color"] <- "siglog2fcdown" # blue
plottab[plottab[,1,drop=FALSE] < -log2fccutoff & plottab[,2] < pvalcutoff, "color"] <- "sigbothdown" # dark blue

plottab[,2] = -log10(plottab[,2])
plottabin = as.data.frame(plottab)
# great, pink, red, darkred, light blue, blue, dark blue
# colorparam = c("grey", "#ff9999", "#ff0000", "#990000", "#b2b2ff", "#1919ff", "#0000b2")
# colorparam = c("grey", "#e5af73", "#8a181a", "#176533", "#e5af73", "#8a181a", "#176533")
colorparam = c("grey", "#925e9f", "#3d5488", "#42b53f", "#925e9f", "#3d5488", "#42b53f")
names(colorparam) = c("notsig","sigpvalup", "siglog2fcup", "sigbothup", "sigpvaldown", "siglog2fcdown", "sigbothdown" )

#e74b35 # red
#3d5488 # blue
#f29b7e # tan
#925e9f # purple
#42b53f # green

## MANUALLY LABELING GENES
pout <- ggplot(plottabin, aes(x = plottabin[,1], y = plottabin[,2], color = plottabin[,3]))
pout <- pout + geom_point(size = 1) + scale_color_manual(values = colorparam)
pout <- pout + labs(x = "R Value", y = "-Log10(Adj.P.Value)", color = "significance")
# labeledgenes = c(
# )
# labeltab = plottabin[rownames(plottabin) %in% labeledgenes,]
# pout <- pout + geom_label_repel(
#   data = labeltab,
#   aes(x = labeltab[,1], y=labeltab[,2], label = rownames(labeltab), color=labeltab[,3]),
#   size=5, segment.size = 0.2, fontface="bold", segment.color = "black",
#   box.padding = unit(0.5, "lines"), point.padding = unit(0.2, "lines"), max.overlaps = 1000
# )
pout <- pout + geom_hline(yintercept = -log10(0.05), linetype = 2) + geom_vline(xintercept = c(-rvalcutoff, rvalcutoff), linetype = 2)
pout <- pout + theme_pubr(base_size = 16, legend = "none")
pdf(paste0(outfilepathmaster, "figures/", "deseqpval_v_corr_volcano.pdf"), width = 5, height = 4, useDingbats = FALSE)
print(pout)
junk <- dev.off()

## Histogram of the individal values for deseq
deseq_histotable <- -log10(fulltestoutfinal[,"deseq_pval",drop=FALSE])
deseq_histotable[,"pvalsig"] <- ifelse(deseq_histotable[,"deseq_pval"] > -log10(0.05), "#e5af73", "grey")
labsparam = list(title = "DESeq pval for hyper v hypo", x = "Pval", y = "value", fill = "significance")
pout <- ggplot(deseq_histotable, aes(x = deseq_histotable[,"deseq_pval"], fill = deseq_histotable[,"pvalsig"]))
pout <- pout + scale_fill_manual(breaks = c("#e5af73", "grey"), values = c("#e5af73", "grey"))
pout <- pout + scale_x_continuous(breaks = c(0,3,6,9,12))
pout <- pout + geom_histogram(color = "black", bins = 100, alpha = 1, position="stack")
pout <- pout + labs(title = labsparam$title, x = labsparam$x, y = labsparam$y, fill = labsparam$fill )
pout <- pout + theme_bw()
pdf(paste0(outfilepathmaster, "figures/", "deseq_pval_histogram.pdf"))
print(pout)
junk <- dev.off()

## Histogram of the individal values for epi corr
corr_histotable <- fulltestoutfinal[,"epi_04um_300s_spearman",drop=FALSE]
corr_histotable[,"rvalsig"] <- ifelse(abs(corr_histotable[,"epi_04um_300s_spearman"]) > rvalcutoff, "#8a181a", "grey")
labsparam = list(title = "DESeq pval for hyper v hypo", x = "Pval", y = "value", fill = "significance")
pout <- ggplot(corr_histotable, aes(x = corr_histotable[,"epi_04um_300s_spearman"], fill = corr_histotable[,"rvalsig"]))
pout <- pout + scale_fill_manual(breaks = c("#8a181a", "grey"), values = c("#8a181a", "grey"))
# pout <- pout + scale_x_continuous(breaks = c(0,3,6,9,12))
pout <- pout + geom_histogram(color = "black", bins = 100, alpha = 1, position="stack")
pout <- pout + labs(title = labsparam$title, x = labsparam$x, y = labsparam$y, fill = labsparam$fill )
pout <- pout + theme_bw()
pdf(paste0(outfilepathmaster, "figures/", "corr_rval_histogram.pdf"))
print(pout)
junk <- dev.off()





## PCA?
pcadata = t(normcounttable[GOI,c(hypersamples, nothypersamples)])
colorvar = epi_metatable[c(hypersamples, nothypersamples),"hypercohort",drop=FALSE]
labsparam = list(title = "PCA", x = "PC1", y = "PC2", color="hypercohort")
pcaplotout <- pca_plotter(pcadata = pcadata, colorvar = colorvar, scalecolor=FALSE, separatelegend = FALSE,
                            labelpoints = FALSE, labsparam = labsparam, returnoutliers = FALSE)

pdf(paste0(outfilepathmaster, "figures/", "GOI_PCA.pdf"))
print(pcaplotout[[1]] + stat_ellipse())
junk <- dev.off()
pcamodeldata <- pcaplotout[["pcdata"]]

## What if i did a simple RF using the PCs, and then returned the variable importance for prediction....?
# Create modeling input data
library("caret")
outcomelabel <- "hypercohort"
outcomelevels <- c("hyper", "nothyper")
SOI <- c(hypersamples, nothypersamples)
featurelabels <- colnames(pcamodeldata[[5]])
outcometable <- epi_metatable[SOI,"hypercohort",drop=FALSE]
featuretable <- pcamodeldata[[5]][SOI,]
modeling_intable <- cbind(outcometable, featuretable)
colnames(modeling_intable)[1] <- outcomelabel
modeling_intable[,outcomelabel] <- factor(modeling_intable[,outcomelabel], levels = outcomelevels)

# create data partition and prepare training scheme
trainindex <- createDataPartition(modeling_intable[,outcomelabel], p = 1, list = FALSE)
training <- modeling_intable[trainindex,]
testing <- modeling_intable[-trainindex,]
modellabel = "rf"
trcontrol <- trainControl(method="repeatedcv", number=10, repeats=5, 
                          classProbs = TRUE, savePredictions = TRUE)

modfit <- train(x=training[,featurelabels], y = training[,outcomelabel], method = modellabel, 
                trControl=trcontrol, preProc = c("center", "scale"), metric = "ROC")

PCimportance <- varImp(modfit, scale=FALSE)$importance
PCimportance <- PCimportance[order(PCimportance[,1], decreasing = TRUE),,drop=FALSE]

# top PCs are 1 and 9
pcadata = t(normcounttable[GOI,c(hypersamples, nothypersamples)])
colorvar = epi_metatable[c(hypersamples, nothypersamples),"hypercohort",drop=FALSE]
labsparam = list(title = "PCA", x = "PC1", y = "PC20", color="hypercohort")
pcaplotout <- pca_plotter(pcadata = pcadata, colorvar = colorvar, scalecolor=FALSE, separatelegend = FALSE,
                          labelpoints = FALSE, labsparam = labsparam, returnoutliers = FALSE, custompcs = c("PC1", "PC20"))

pdf(paste0(outfilepathmaster, "figures/", "GOI_PCA_PC1andPC20.pdf"))
print(pcaplotout[[1]] + stat_ellipse())
junk <- dev.off()







## What about a standard heatmap (annotated) for the hyper hypo and our GOI
# Trying other genes for our HM viz
## Using our full GOI
# hm_GOI <- GOI
## Using our DESeq values
deseq_results_table <- deseqout[[1]]
# hm_GOI <- rownames(deseq_results_table[deseq_results_table[,"padj"] < 0.05 & !is.na(deseq_results_table[,"padj"]) &
#                                        abs(deseq_results_table[,"log2FoldChange"]) > 0 ,])
# hm_GOI <- overlapgrouptab[str_count(overlapgrouptab[,"overlapgroup"],";") > 0,"Row.names"]
# hm_GOI <- overlapgrouptab[grepl("deseq_sig", overlapgrouptab[,"overlapgroup"]) &
#                           grepl("wilcox_sig", overlapgrouptab[,"overlapgroup"]),"Row.names"]
hm_GOI <- overlapgrouptab[grepl("deseq_sig", overlapgrouptab[,"overlapgroup"]) & 
                          grepl("epicorr_sig", overlapgrouptab[,"overlapgroup"]), "Row.names"]
hm_GOI <- GOI # not sure why we'd want to use anything else
length(hm_GOI)



plottable2 <- normcounttable[hm_GOI,c(hypersamples, nothypersamples)]
hmmetatable <- na.omit(epi_metatable[c(hypersamples, nothypersamples),c("hypercohort", "epi_04um_300s_n")])
colannotation <- annotationlist_builder(hmmetatable)
# hmout <- create_heatmap(counttab = plottable2, scale_data = TRUE, colmetatable = hmmetatable, colannotationlist = colannotation, colclusterparam = FALSE, rowclusterparam = TRUE, separate_legend = FALSE)

sampleannottable <- hmmetatable
sampleannotatiaonlist <- colannotation
hmplottab <- plottable2
colsplitparam <- sampleannottable[,"hypercohort"]
# rowsplitparam

## Top annotation
temp1 <- vector("list", length(sampleannotatiaonlist))
names(temp1) = names(sampleannotatiaonlist)
annotlegendlist = lapply(temp1, function(x) x[[1]] =
                           list(title_gp=gpar(fontsize=5, fontface="bold"), labels_gp=gpar(fontsize=4)))
## Define a param that will go through each annotation - and keep a legend if its continuous or has less than 10 discrete terms, otherwise hide the legend
showlegendparam = unname(unlist(lapply(sampleannotatiaonlist, function(x) {
  numterms = tryCatch(length(na.omit(unique(x))), error=function(e) NULL)
  is.null(numterms) || numterms <= 10})))
hatop = HeatmapAnnotation(df = sampleannottable,
                          col = sampleannotatiaonlist,
                          na_col = "white",
                          show_annotation_name = TRUE,
                          annotation_name_gp = gpar(fontsize = 8, fontface="bold"),
                          annotation_name_side = "left",
                          simple_anno_size = unit(min(60/length(sampleannotatiaonlist), 5),"mm"),
                          show_legend = TRUE,
                          annotation_legend_param = annotlegendlist)


# ## Side annotation
# ## Define parameters for each of the labels on the annotation bars
# temp1 <- vector("list", length(geneannotatiaonlist))
# names(temp1) = names(geneannotatiaonlist)
# annotlegendlist = lapply(temp1, function(x)
#   x[[1]] = list(title_gp=gpar(fontsize=8, fontface="bold"), labels_gp=gpar(fontsize=8)))

# ## Define a param that will go through each annotation - and keep a legend if its continuous or has less than 10 discrete terms, otherwise hide the legend
# showlegendparam = unname(unlist(lapply(geneannotatiaonlist, function(x) {
#   numterms = tryCatch(length(na.omit(unique(x))), error=function(e) NULL)
#   is.null(numterms) || numterms <= 10})))
# 
# ## Look for any empty annotations - fill them with white, and later, make sure to hide their legend
# emptyannots = names(sapply(geneannotatiaonlist, length)[sapply(geneannotatiaonlist, length)==0])
# if (length(emptyannots) > 0){
#   for (i in 1:length(emptyannots)) {
#     temp1 = "white"
#     names(temp1) = emptyannots[i]
#     geneannotatiaonlist[[emptyannots[i]]] = temp1
#   }
#   showlegendparam[which(names(geneannotatiaonlist) %in% emptyannots)] = FALSE
# }
# ## Add param that will bolden the side annotation bars if it is <100, and omit the grid lines if more
# if (nrow(geneannottable) < 100) {
#   sideannotation_linebold_param <- gpar(fontsize = 0.5)} else {sideannotation_linebold_param <- NULL}
# haside = rowAnnotation(df = geneannottable,
#                        col = geneannotatiaonlist,
#                        na_col = "white",
#                        # gp = gpar(fontsize = 0.01),
#                        gp = sideannotation_linebold_param,
#                        show_annotation_name=TRUE,
#                        annotation_name_gp = gpar(fontsize = 8, fontface="bold"),
#                        annotation_name_side = "top",
#                        simple_anno_size = unit(min(60/length(geneannotatiaonlist), 5),"mm"),
#                        show_legend = TRUE,
#                        annotation_legend_param = annotlegendlist)


hmplottab = as.matrix(t(apply(hmplottab, 1, function(x) zscore(x))))
heatmapcolorparam <- colorRamp2(breaks = c(2, 0, -2), 
                                c("red", "white", "blue"))
summary_outhm <- Heatmap(as.matrix(hmplottab),
                         col = heatmapcolorparam, row_title = "Genes", column_title = "Samples",
                         # column_split = colsplitparam,
                         top_annotation = hatop,
                         cluster_columns = TRUE, cluster_rows = TRUE,

                         show_column_names = TRUE, column_names_gp = gpar(fontsize = 6),
                         show_row_names = FALSE, row_names_side = "left", row_names_gp = gpar(fontsize=6),
                         
                         heatmap_legend_param = list(
                           title = "Zscore",
                           title_gp = gpar(fontsize = 8, fontface = "bold")),
                         height = unit(min((nrow(hmplottab)/2), 12),"cm"),
                         width = unit(min(ncol(hmplottab), 18),"cm")
                         
)
# pdf(paste0(outfilepathmaster, "figures/", "GOI_heatmap_clustered_split.pdf"), 15, 10, useDingbats = FALSE)
# pdf(paste0(outfilepathmaster, "figures/", "GOI_heatmap_clustered.pdf"), 15, 10, useDingbats = FALSE)
# pdf(paste0(outfilepathmaster, "figures/", "GOI_heatmap_sorted.pdf"), 15, 10, useDingbats = FALSE)
pdf(paste0(outfilepathmaster, "figures/", "custGOI_deseqANDcorr_heatmap_clustered.pdf"), 15, 10, useDingbats = FALSE)
draw(summary_outhm)
junk <- dev.off()





pcadata = t(normcounttable[hm_GOI,c(hypersamples, nothypersamples)])
colorvar = epi_metatable[c(hypersamples, nothypersamples),"hypercohort",drop=FALSE]
labsparam = list(title = "PCA", x = "PC1", y = "PC2", color="hypercohort")
pcaplotout <- pca_plotter(pcadata = pcadata, colorvar = colorvar, scalecolor=FALSE, separatelegend = FALSE,
                          labelpoints = FALSE, labsparam = labsparam, returnoutliers = FALSE)

pdf(paste0(outfilepathmaster, "figures/", "custGOI_PCA.pdf"))
print(pcaplotout[[1]] + stat_ellipse())
junk <- dev.off()





## Another figure comparing our GOI to existing GO terms for platelets
customgeneset <- data.frame(gs_name = c(rep("custom_mgc_hyper_up", length(GOIup)), 
                                        rep("custom_mgc_hyper_down", length(GOIdown))), 
                            gene_symbol = c(GOIup, GOIdown))
## Platelet GO term annotation table
speciesparam <- "Homo sapiens"
m_t2g_cat <- as.data.frame(msigdbr(species = speciesparam, category = c("C5"))[,c("gs_name", "gene_symbol")])
plateletgotermtable <- unique(m_t2g_cat[grepl("PLATELET", m_t2g_cat[,1]),])
plateletgotermtable <- rbind(plateletgotermtable[!grepl("PLATELET_DERIVED_GROWTH_FACTOR|HP_", plateletgotermtable[,1]),], customgeneset)
plateletgotermtable[,"gene"] <- plateletgotermtable[,"gene_symbol"]
plateletgoterm_annottable <- dcast(plateletgotermtable,gene_symbol ~ gs_name, value.var = "gene")
rownames(plateletgoterm_annottable) <- plateletgoterm_annottable[,"gene_symbol"]
plateletgoterm_annottable <- plateletgoterm_annottable[,!grepl("gene_symbol", colnames(plateletgoterm_annottable))]

## upset plot viz
upset_plot_inlist <- apply(plateletgoterm_annottable, 2, function(x) unname(x[!is.na(x)]))
pout <- create_upset_plot(upset_indata = upset_plot_inlist)
## This was too much - so I think instead, I should just do a comparison of in ANY platelet pathway
speciesparam <- "Homo sapiens"
m_t2g_cat <- as.data.frame(msigdbr(species = speciesparam, category = c("C5"))[,c("gs_name", "gene_symbol")])
plateletgotermtable <- unique(m_t2g_cat[grepl("PLATELET", m_t2g_cat[,1]),])
plateletgotermtable <- plateletgotermtable[!grepl("PLATELET_DERIVED_GROWTH_FACTOR|HP_", plateletgotermtable[,1]),]
overlapinlist <- list(GOIup = GOIup, GOIdown = GOIdown, anyplateletgoterm = unique(plateletgotermtable[,"gene_symbol"]))
overlapout <- overlap_finder(overlapinlist)
write.table(overlapout$overlaptable, paste0(outfilepathmaster, "figures/", "GOI_plateletgoterm_overlaptable.csv"),
            sep = ",", col.names = TRUE, row.names = FALSE)
write.table(overlapout$overlapgrouptab, paste0(outfilepathmaster, "figures/", "GOI_plateletgoterm_overlapgrouptab.csv"),
            sep = ",", col.names = TRUE, row.names = FALSE)
write.table(overlapout$overlapsummary, paste0(outfilepathmaster, "figures/", "GOI_plateletgoterm_overlapsummary.csv"),
            sep = ",", col.names = TRUE, row.names = TRUE)
pdf(paste0(outfilepathmaster, "figures/", "GOI_plateletgoterm_overlap_venn.pdf"))
grid.draw(overlapout$vennplot)
junk <- dev.off()










## Whats the GSEA results from our hyper hypo deseq
outfilepathGSEA = paste0(outfilepathmaster, "gsea/hyper_gseacheck/")
dir.create(outfilepathGSEA, recursive = TRUE, showWarnings = FALSE)
## {"Homo sapiens", "Mus musculus"}
speciesparam = "Homo sapiens"
pstatparam = "pvalue"
numpathways_plotparam <- 10

## Seeding to get rid of the change in results between runs
seedparam = 1118065690
## GSEA run through for HALLMARK terms
geneset_analysis_out_HALL = geneset_analysis(testout2, pvalcutoffparam = 1, genesetparam = c("H"), speciesparam = speciesparam, seedparam = seedparam)
write.table(geneset_analysis_out_HALL, file = paste0(outfilepathGSEA, "hyper_v_nothyperdeseq", "_gsea_HALL.csv"), 
            sep = ",", row.names = TRUE, col.names = NA)
geneset_analysis_out_GO = geneset_analysis(testout2, pvalcutoffparam = 1, genesetparam = c("C5"), speciesparam = speciesparam, seedparam = seedparam)
write.table(geneset_analysis_out_GO, file = paste0(outfilepathGSEA, "hyper_v_nothyperdeseq", "_gsea_GO.csv"), 
            sep = ",", row.names = TRUE, col.names = NA)

## GSEA Hyper
statcutoffparamlist = c("stattype" = "pvalue", "pstatcutoff" = 0.05, "log2fccutoff" = 0)
hypergeo_genetest_out_HALL = hypergeo_genetest(testout2,
                                               statcutoffparam = statcutoffparamlist, 
                                               genesetparam = c("H"), speciesparam = speciesparam)
write.table(hypergeo_genetest_out_HALL$enricherUPout, 
            file = paste0(outfilepathGSEA, "hyper_v_nothyperdeseq", "_hypergeo_UP_HALL.csv"), 
            sep = ",", row.names = TRUE, col.names = NA)
write.table(hypergeo_genetest_out_HALL$enricherDOWNout, 
            file = paste0(outfilepathGSEA, "hyper_v_nothyperdeseq", "_hypergeo_DOWN_HALL.csv"), 
            sep = ",", row.names = TRUE, col.names = NA)

hypergeo_genetest_out_GO = hypergeo_genetest(testout2,
                                             statcutoffparam = statcutoffparamlist, 
                                             genesetparam = c("C5"), speciesparam = speciesparam)
write.table(hypergeo_genetest_out_GO$enricherUPout, 
            file = paste0(outfilepathGSEA, "hyper_v_nothyperdeseq", "_hypergeo_UP_GO.csv"), 
            sep = ",", row.names = TRUE, col.names = NA)
write.table(hypergeo_genetest_out_GO$enricherDOWNout, 
            file = paste0(outfilepathGSEA, "hyper_v_nothyperdeseq", "_hypergeo_DOWN_GO.csv"), 
            sep = ",", row.names = TRUE, col.names = NA)

## Recalc the p.adjust values for upGSEA
# platelet_gogsea <- hypergeo_genetest_out_GO$enricherUPout
# platelet_gogsea <- platelet_gogsea[grepl("PLATELET_", platelet_gogsea[,"ID"]),]
# platelet_gogsea <- platelet_gogsea[!grepl("PLATELET_DERIVED_", platelet_gogsea[,"ID"]),]
# platelet_gogsea[,"new_padj"] <- p.adjust(platelet_gogsea[,"pvalue"], method = "fdr")


## Shoot - all tests before were using the DESeq output, what about just selecting our GOI...........
hypergeo_GOIup_out_GO = hypergeo_genetest(data.frame(GOIup),
                                             statcutoffparam = statcutoffparamlist, 
                                             genesetparam = c("C5"), speciesparam = speciesparam)
write.table(hypergeo_GOIup_out_GO$enricherUPout, 
            file = paste0(outfilepathGSEA, "hyper_v_nothyper", "_GOIup_hypergeo_GO.csv"), 
            sep = ",", row.names = TRUE, col.names = NA)

hypergeo_GOIdown_out_GO = hypergeo_genetest(data.frame(GOIdown),
                                          statcutoffparam = statcutoffparamlist, 
                                          genesetparam = c("C5"), speciesparam = speciesparam)
write.table(hypergeo_GOIdown_out_GO$enricherUPout, 
            file = paste0(outfilepathGSEA, "hyper_v_nothyper", "_GOIdown_hypergeo_GO.csv"), 
            sep = ",", row.names = TRUE, col.names = NA)

hypergeo_GOIup_out_H = hypergeo_genetest(data.frame(GOIup),
                                          statcutoffparam = statcutoffparamlist, 
                                          genesetparam = c("H"), speciesparam = speciesparam)
write.table(hypergeo_GOIup_out_H$enricherUPout, 
            file = paste0(outfilepathGSEA, "hyper_v_nothyper", "_GOIup_hypergeo_H.csv"), 
            sep = ",", row.names = TRUE, col.names = NA)

hypergeo_GOIdown_out_H = hypergeo_genetest(data.frame(GOIdown),
                                          statcutoffparam = statcutoffparamlist, 
                                          genesetparam = c("H"), speciesparam = speciesparam)
write.table(hypergeo_GOIdown_out_H$enricherUPout, 
            file = paste0(outfilepathGSEA, "hyper_v_nothyper", "_GOIdown_hypergeo_H.csv"), 
            sep = ",", row.names = TRUE, col.names = NA)

## What if we enrich for all of the genes... not something i would normally do? but lets see - NOPE, JUST NOISIER
hypergeo_GOIall_out_GO = hypergeo_genetest(data.frame(GOI),
                                            statcutoffparam = statcutoffparamlist, 
                                            genesetparam = c("C5"), speciesparam = speciesparam)
write.table(hypergeo_GOIall_out_GO$enricherUPout, 
            file = paste0(outfilepathGSEA, "hyper_v_nothyper", "_GOIall_hypergeo_GO.csv"), 
            sep = ",", row.names = TRUE, col.names = NA)

hypergeo_GOIall_out_H = hypergeo_genetest(data.frame(GOI),
                                         statcutoffparam = statcutoffparamlist, 
                                         genesetparam = c("H"), speciesparam = speciesparam)
write.table(hypergeo_GOIall_out_H$enricherUPout, 
            file = paste0(outfilepathGSEA, "hyper_v_nothyper", "_GOIall_hypergeo_H.csv"), 
            sep = ",", row.names = TRUE, col.names = NA)


## Recalc the p.adjust values for upGSEA
platelet_gogsea <- hypergeo_GOIup_out_GO$enricherUPout
platelet_gogsea <- platelet_gogsea[grepl("PLATELET_", platelet_gogsea[,"ID"]),]
platelet_gogsea <- platelet_gogsea[!grepl("PLATELET_DERIVED_", platelet_gogsea[,"ID"]),]
platelet_gogsea[,"new_padj"] <- p.adjust(platelet_gogsea[,"pvalue"], method = "fdr")



## The best one is the following:
# hypergeo_GOIup_out_GO - so lets doa barplot with some chosen pathways
POI <- c("GOCC_SECRETORY_GRANULE", "GOBP_MYELOID_LEUKOCYTE_MEDIATED_IMMUNITY", "GOCC_FICOLIN_1_RICH_GRANULE", 
         "GOCC_ANCHORING_JUNCTION", "GOBP_CELL_ACTIVATION_INVOLVED_IN_IMMUNE_RESPONSE", "GOBP_EXOCYTOSIS",
         "GOBP_RHO_PROTEIN_SIGNAL_TRANSDUCTION", 
         "GOBP_MYELOID_LEUKOCYTE_ACTIVATION", 
         "GOMF_ENZYME_BINDING", 
         "GOMF_CYTOSKELETAL_PROTEIN_BINDING", 
         "GOBP_LEUKOCYTE_MEDIATED_IMMUNITY", 
         "GOBP_AUTOPHAGY_OF_MITOCHONDRION", 
         "GOCC_ACTIN_CYTOSKELETON", 
         "GOBP_IMMUNE_EFFECTOR_PROCESS", 
         "GOBP_WOUND_HEALING")
gseaplotin <- data.frame(hypergeo_GOIup_out_GO$enricherUPout)[POI,c("ID", "p.adjust")]
gseaplotin[,"pstat"] <- -log10(gseaplotin[,"p.adjust"])
gseaplotin[,"ID"] <- apply(data.frame(tolower(gsub("GOBP |GOCC |GOMF ", "", gsub("_", " ", gseaplotin[,"ID"])))), 1, simpleCap)
gseaplotin <- gseaplotin[order(gseaplotin[,"pstat"], decreasing = TRUE),]
pout <- ggplot(gseaplotin, mapping = aes(x = gseaplotin[,"ID"], y = gseaplotin[,"pstat"]))
pout <- pout + geom_bar(stat = "identity")
# pout <- pout + scale_x_discrete(limits=rev(gseaplotin[,"ID"]))
pout <- pout + coord_flip()
pout <- pout + labs(x="Geneset", y = "-log10(p.adjust)", title = "hypergeo_GOIup_out_GO_POI")
pout <- pout + theme_pubr()

pdf(paste0(outfilepathGSEA, "hypergeo_GOIup_out_GO_POI.pdf"), width = 10, height = 5)
print(pout)
junk <- dev.off()

## Also do down GOI - with a NOMINAL P VALUE
POI <- c("GOBP_PROTEIN_MODIFICATION_BY_SMALL_PROTEIN_CONJUGATION_OR_REMOVAL", "GOBP_REGULATION_OF_RNA_SPLICING",
         "GOBP_REGULATION_OF_TRANSCRIPTION_INITIATION_FROM_RNA_POLYMERASE_II_PROMOTER", "GOMF_CATALYTIC_ACTIVITY_ACTING_ON_RNA",
         "GOBP_RNA_POLYMERASE_II_PREINITIATION_COMPLEX_ASSEMBLY", "GOBP_PROTEIN_TRANSPORT_ALONG_MICROTUBULE",
         "GOBP_CELLULAR_IRON_ION_HOMEOSTASIS", "GOMF_OXYGEN_CARRIER_ACTIVITY", "GOBP_RNA_SPLICING", "GOMF_ENZYME_BINDING")
gseaplotin <- data.frame(hypergeo_GOIdown_out_GO$enricherUPout)[POI,c("ID", "pvalue")]
gseaplotin[,"pstat"] <- -log10(gseaplotin[,"pvalue"])
gseaplotin[,"ID"] <- apply(data.frame(tolower(gsub("GOBP |GOCC |GOMF ", "", gsub("_", " ", gseaplotin[,"ID"])))), 1, simpleCap)
gseaplotin <- gseaplotin[order(gseaplotin[,"pstat"], decreasing = TRUE),]
pout <- ggplot(gseaplotin, mapping = aes(x = gseaplotin[,"ID"], y = gseaplotin[,"pstat"]))
pout <- pout + geom_bar(stat = "identity")
pout <- pout + scale_x_discrete(limits=rev(gseaplotin[,"ID"]))
pout <- pout + coord_flip()
pout <- pout + labs(x="Geneset", y = "-log10(pvalue)", title = "hypergeo_GOIdown_out_GO_POI")
pout <- pout + theme_pubr()

pdf(paste0(outfilepathGSEA, "hypergeo_GOIdown_out_GO_POI.pdf"), width = 10, height = 5)
print(pout)
junk <- dev.off()



# --------------------------------- Validation across all of our data sets---------------------------------
## Validation!
# -          I think the VERY OLD sequencing we did many years ago where we used 12 very HYPER versus 12 very HYPO
# o   Of note, we chose 6 HYPER males vs 6 HYPO males, and 6 HYPER females vs 6 HYPO females (sex was diff in hypr and hypo so you my want to consider)
# o   While there was a small issue with this dataset as it related to sex – I think the HYPER and HYPO were robust (and this was in controls)
# o   The nice think about this dataset was that people came back several times so we were confident in their HYPER and HYPO category
# -          When we did the PACE Platelet RNA, we also did THR and EMD – both of those datasets have Epi 0.4 300s
# o   We would need ot figure out who was on/off antiplatelet therapy – as that may be important
# o   I think looking at PACE vs control (THR/EMD) you would expect to see platelet related genes to be diff expressed between these populations
# -          HARP – I would not look at PLT aggregation as use of antiplatelet therapy is robust and some used P2Y12i and others not – but I would definitely look at MI vs control (at the very least would compare MI-CAD vs controls)
# -          On Friday, we can talk to Deepak about the Duke dataset….. but I think we should have enough to start with
# 

## Validation barplot
validation_barplot <- function(gseabarplot_intable, numberofterms = 10, titleparam = NULL, testtype = "hyper") {
  
  if (testtype == "hyper") {
    gseaplotin <- data.frame(gseabarplot_intable)[order(abs(gseabarplot_intable[,"p.adjust"]), decreasing = FALSE)[1:numberofterms],c("ID", "p.adjust")]
    ## if the mgc one isnt there - then add it on the bottom
    if (sum(grepl("mgc", rownames(gseaplotin))) == 0) {
      gseaplotin <- rbind(gseaplotin, gseabarplot_intable[grepl("mgc", gseabarplot_intable[,"ID"]), c("ID", "p.adjust")])
    }
  } else {
    ## Want the number of terms /2, and then that many top and bottom gene sets.
    splitnumberofterms <- numberofterms/2
    gseaplotin <- gseabarplot_intable[order(gseabarplot_intable[,"NES"], decreasing = TRUE), c("ID", "p.adjust", "NES")][c(1:splitnumberofterms, (nrow(gseabarplot_intable)-(splitnumberofterms-1)):nrow(gseabarplot_intable)),]
    ## if the mgc one isnt there - then add it on the bottom
    if (sum(grepl("mgc", rownames(gseaplotin))) != 2) {
      gseaplotin <- unique(rbind(gseaplotin, gseabarplot_intable[grepl("mgc", gseabarplot_intable[,"ID"]), c("ID", "p.adjust", "NES")]))
      gseaplotin <- gseaplotin[order(gseaplotin[,"NES"], decreasing = TRUE),]
    }
  }

  ## Perform transforming and data cleaning
  gseaplotin[,"pstat"] <- -log10(gseaplotin[,"p.adjust"])
  gseaplotin[,"ID"] <- apply(data.frame(tolower(gsub("GOBP |GOCC |GOMF ", "", gsub("_", " ", gseaplotin[,"ID"])))), 1, simpleCap)
  
  #pout <- ggplot(gseaplotin, mapping = aes(x = gseaplotin[,1], y = gseaplotin[,2], fill = gseaplotin[,3]))
  # pout <- ggplot(gseaplotin, mapping = aes(x = str_wrap(gseaplotin[,"ID"], 30), y = gseaplotin[,"pstat"]))
  if (testtype == "hyper") {
    pout <- ggplot(gseaplotin, mapping = aes(x = gseaplotin[,"ID"], y = gseaplotin[,"pstat"]))
  }
  if(testtype == "fullgsea") {
    pout <- ggplot(gseaplotin, mapping = aes(x = gseaplotin[,"ID"], y = gseaplotin[,"NES"], fill = gseaplotin[,"pstat"]))
    pout <- pout + scale_fill_viridis(option = "magma")
  }
  pout <- pout + geom_bar(stat = "identity", color = "black")
  pout <- pout + scale_x_discrete(limits=rev(gseaplotin[,"ID"]))
  pout <- pout + coord_flip() + labs(x="Geneset", y = "-log10(p.adjust)", title = titleparam) + theme_pubr()
  pout
  return(pout)
}

validation_gsea_function <- function(indeseqtable, validationlabel, validationoutfilepath) {
  statcutoffparamlist <- c("stattype" = "pvalue", "pstatcutoff" = 0.05, "log2fccutoff" = 0)
  customgeneset <- data.frame(gs_name = c(rep("custom_mgc_hyper_up", length(GOIup_validation)), 
                                          rep("custom_mgc_hyper_down", length(GOIdown_validation))), 
                              gene_symbol = c(GOIup_validation, GOIdown_validation))
  custom_hypergseatest = hypergeo_genetest(DEseqtable = indeseqtable,
                                           statcutoffparam = statcutoffparamlist, 
                                           genesetparam = c("C5"), speciesparam = speciesparam, customgeneset = customgeneset)
  write.table(custom_hypergseatest$enricherUPout, paste0(validationoutfilepath, validationlabel, "_GO_w_custom_UP.csv"), 
              sep = ",", col.names = TRUE, row.names = FALSE)
  pout1 <- validation_barplot(custom_hypergseatest$enricherUPout, numberofterms = 10, titleparam = "MI vs Control")
  pdf(paste0(validationoutfilepath, validationlabel, "_GO_w_custom_UP.pdf"), width = 10, height = 5)
  print(pout1)
  junk <- dev.off()
  write.table(custom_hypergseatest$enricherDOWNout, paste0(validationoutfilepath, validationlabel, "_GO_w_custom_DN.csv"), 
              sep = ",", col.names = TRUE, row.names = FALSE)
  pout2 <- validation_barplot(gseabarplot_intable = custom_hypergseatest$enricherDOWNout, numberofterms = 10, titleparam = "MI vs Control")
  pdf(paste0(validationoutfilepath, validationlabel, "_GO_w_custom_DN.pdf"), width = 10, height = 5)
  print(pout2)
  junk <- dev.off()
  
  ## Doing this with preranked full GSEA to including our custom set (yikes...)
  custom_fullgseatest = geneset_analysis(indeseqtable, 
                                         pvalcutoffparam = 1, genesetparam = c("C5"), speciesparam = speciesparam,
                                         seedparam = seedparam, customgeneset = customgeneset)
  
  write.table(custom_fullgseatest, file = paste0(validationoutfilepath, validationlabel, "_GO_w_custom_fullgsea.csv"), 
              sep = ",", row.names = TRUE, col.names = NA)
  pout3 <- validation_barplot(gseabarplot_intable = custom_fullgseatest, numberofterms = 10, titleparam = "MI vs Control", testtype = "fullgsea")
  pdf(paste0(validationoutfilepath, validationlabel, "_GO_w_custom_fullgsea.pdf"), width = 10, height = 5)
  print(pout3)
  junk <- dev.off()
  
  if ("custom_mgc_hyper_up" %in% custom_fullgseatest[,1]) {
      gseaplot_UP <- gsea_custom_randomwalkplot(gseaobject = custom_fullgseatest, geneSetID = "custom_mgc_hyper_up", addNESvalue = TRUE)
  } else { gseaplot_UP <- NULL }
  if ("custom_mgc_hyper_down" %in% custom_fullgseatest[,1]) {
      gseaplot_DN <- gsea_custom_randomwalkplot(gseaobject = custom_fullgseatest, geneSetID = "custom_mgc_hyper_down", addNESvalue = TRUE)
  } else { gseaplot_DN <- NULL }
  # gseaplot_UP <- gseaplot(custom_fullgseatest, geneSetID = "custom_mgc_hyper_up", title = "custom_mgc_hyper_up")
  # gseaplot_DN <- gseaplot(custom_fullgseatest, geneSetID = "custom_mgc_hyper_down", title = "custom_mgc_hyper_down")
  pdf(width = 11, height = 8.5, paste0(validationoutfilepath, validationlabel, "_GO_w_custom_gseawalkplots.pdf"))
  print(gseaplot_UP)
  print(gseaplot_DN)
  junk <- dev.off()
  
  return(list(hypergsea_out = custom_hypergseatest, fullgsea_out = custom_fullgseatest))
}


write_out_intermediate_files <- function(rawcounttable, metatable, compcols, controlcols, outfilepath, labelparam) {
    dir.create(outfilepath, showWarnings = FALSE, recursive = TRUE)
    write.table(rawcounttable, paste0(outfilepath, labelparam, "__rawcounttable.csv"), sep = ",", col.names = NA, row.names = TRUE)
    write.table(metatable, paste0(outfilepath, labelparam, "__metatable.csv"), sep = ",", col.names = NA, row.names = TRUE)
    write.table(compcols, paste0(outfilepath, labelparam, "__compcols.csv"), sep = ",", col.names = NA, row.names = TRUE)
    write.table(controlcols, paste0(outfilepath, labelparam, "__controlcols.csv"), sep = ",", col.names = NA, row.names = TRUE)
}


## Set valiation outfilepath
validationoutfilepath <- paste0(outfilepathmaster, "validation/")
dir.create(validationoutfilepath, showWarnings = FALSE, recursive = TRUE)

## Set universe of GO term genes
speciesparam <- "Homo sapiens"
GOtermuniversegenes_table <- as.data.frame(msigdbr(species = speciesparam, category = c("C5"))[,c("gs_name", "gene_symbol")])
GOtermuniversegenes <- unique(GOtermuniversegenes_table[,"gene_symbol"])

## Just for validation purposes - to make sure we arent grossly overbiasing, let me ONLY select the genes that are in our universe already
GOIup_validation <- GOIup[GOIup %in% GOtermuniversegenes]
GOIdown_validation <- GOIdown[GOIdown %in% GOtermuniversegenes]




# --------------------------------- PACE vs THR - the original dataset, going back and testing ---------------------------------

## 1 - PACE vs THR or EMD I think would be a good place to start
dir.create(paste0(outfilepathmaster, "validation/", "pace_pacethr/"), showWarnings = FALSE, recursive = TRUE)




# HAVE to RERUN WITH AGR (AGE, GENDER, RACE) CONTROLS
pace_pacethr_rawcount_file <- "output/run7_rmoutliers2_agesexcontrol_withpaired_20220422/rna_processing/filtrawcounttab.txt"
pace_pacethr_rawcount_table <- read.table(pace_pacethr_rawcount_file, sep = "\t", header = TRUE, row.names = 1, check.names = FALSE)
pace_pacethr_meta_file <- "output/run7_rmoutliers2_agesexcontrol_withpaired_20220422/rna_processing/metatable_in.csv"
pace_pacethr_metatable <- read.table(pace_pacethr_meta_file, sep =",", header = TRUE, row.names = 1, comment.char = "", quote = "")
SOI <- rownames(pace_pacethr_metatable)

compcols <- na.omit(pace_pacethr_metatable[SOI,"comp_cohort__pace_v_thr",drop=FALSE])
controlcols <- pace_pacethr_metatable[SOI,c("sex1","age_cat", "race")]
pace_pacethr_deseq <- DESeq_analysis(compcols = compcols[SOI,,drop=FALSE], controlcols = controlcols[SOI,],
                           rawcounttab = pace_pacethr_rawcount_table[,SOI])[[1]]
write.table(pace_pacethr_deseq, paste0(outfilepathmaster, "validation/", "pace_pacethr/", "pace_pace_v_thr_wcontrol_deseq.csv"), 
            sep = ",", col.names = NA, row.names = TRUE)

# Run the abbreviated function
pace_pace_v_thr_validation_out <- validation_gsea_function(indeseqtable = pace_pacethr_deseq, 
                                                           validationlabel = "pace_pace_v_thr", 
                                                           validationoutfilepath = paste0(outfilepathmaster, "validation/", "pace_pacethr/"))

# Write out the intermediates:
write_out_intermediate_files(rawcounttable = pace_pacethr_rawcount_table, metatable = pace_pacethr_metatable, 
                             compcols = compcols, controlcols = controlcols, 
                             outfilepath = paste0(outfilepathmaster, "validation/", "pace_pacethr/", "intermediate_files/"), 
                             labelparam = "pace_pace_v_thr")

# --------------------------------- MI vs Control from HARP ---------------------------------

## 2 - lets try MI vs Control - have to build in the custom geneset function for this to work...
dir.create(paste0(outfilepathmaster, "validation/", "harp_allmicontrol/"), showWarnings = FALSE, recursive = TRUE)
dir.create(paste0(outfilepathmaster, "validation/", "harp_minocacontrol/"), showWarnings = FALSE, recursive = TRUE)
dir.create(paste0(outfilepathmaster, "validation/", "harp_micadcontrol/"), showWarnings = FALSE, recursive = TRUE)
# Read in
harp_minocacontrol_rawcount_file <- file.path(relative_path, "harp-platelet/output/run2_rmoutliers2/rna_processing/filtrawcounttab.txt")
harp_minocacontrol_rawcount_table <- read.table(harp_minocacontrol_rawcount_file, sep = "\t", header = TRUE, row.names = 1, check.names = FALSE)
harp_minocacontrol_meta_file <- file.path(relative_path, "harp-platelet/output/run2_rmoutliers2/rna_processing/metatable_in.csv")
harp_minocacontrol_metatable <- read.table(harp_minocacontrol_meta_file, sep =",", header = TRUE, row.names = 1, comment.char = "", quote = "")
SOI <- rownames(harp_minocacontrol_metatable)

# HARP needs add on data
library(lubridate)
library(stringr)
harp_minocacontrol_meta_ADDON_file <- file.path(relative_path, "harp-platelet/data/harp-platelet-metadata-ADDON-20220823.csv")
harp_minocacontrol_metatable_ADDON <- read.table(harp_minocacontrol_meta_ADDON_file, sep =",", header = TRUE, row.names = 1, comment.char = "", quote = "", fill = TRUE)

# Create custom Race
race_table <- harp_minocacontrol_metatable_ADDON[SOI,c("American.Indian.or.Alaskan.Native", "Asian", "Black.or.African.American", "White", "Native.Hawaiian.or.Other.Pacific.Islander", "Other")]
race_table[race_table == "No"] <- NA
w <- which(race_table=="Yes",arr.ind=TRUE)
race_table[w] <- names(race_table)[w[,"col"]]
controlcols <- data.frame(race = apply(race_table, 1, function(x) paste0(na.omit(x), collapse = "")))
# Create custom age
harp_minocacontrol_age <- decimal_date(now()) - decimal_date(mdy(str_replace(harp_minocacontrol_metatable_ADDON[SOI,"Date.of.Birth"], '[0-9]+$', '19\\0')))
controlcols[,"age_cat"] <- ifelse(harp_minocacontrol_age >= 74, "74_98",
                           ifelse(harp_minocacontrol_age >= 61, "61_73",
                           ifelse(harp_minocacontrol_age >= 28, "28_60", NA)))


# ## Run all of our DESeq analyses based on the various comps I put in
compcols <- harp_minocacontrol_metatable[SOI,c("comp_angiography__allMI_v_allControl", "comp_angiography__MINOCA_v_allControl", "comp_angiography__MICAD_v_allControl")]
# # Adding in age and sex control - because i think it gives the best results.............
controlcols <- controlcols[SOI,]
harp_minocacontrol_deseq_out <- DESeq_analysis(compcols = compcols[SOI,,drop=FALSE], controlcols = controlcols[SOI,],
                                     rawcounttab = harp_minocacontrol_rawcount_table[,SOI])
harp_allmicontrol_deseq <- harp_minocacontrol_deseq_out[["comp_angiography__allMI_v_allControl"]]
harp_minocacontrol_deseq <- harp_minocacontrol_deseq_out[["comp_angiography__MINOCA_v_allControl"]]
harp_micadcontrol_deseq <- harp_minocacontrol_deseq_out[["comp_angiography__MICAD_v_allControl"]]

write.table(harp_allmicontrol_deseq, paste0(outfilepathmaster, "validation/", "harp_allmicontrol/", "harp_allmicontrol_wcontrol_deseq.csv"), 
            sep = ",", col.names = NA, row.names = TRUE)
write.table(harp_minocacontrol_deseq, paste0(outfilepathmaster, "validation/", "harp_minocacontrol/", "harp_minocacontrol_wcontrol_deseq.csv"), 
            sep = ",", col.names = NA, row.names = TRUE)
write.table(harp_micadcontrol_deseq, paste0(outfilepathmaster, "validation/", "harp_micadcontrol/", "harp_micadcontrol_wcontrol_deseq.csv"), 
            sep = ",", col.names = NA, row.names = TRUE)

# ALL MI
harp_allmicontrol_validation_out <- validation_gsea_function(indeseqtable = harp_allmicontrol_deseq, 
                                                             validationlabel = "harp_allmi_v_control",
                                                             validationoutfilepath = paste0(outfilepathmaster, "validation/", "harp_allmicontrol/"))

# MINOCA
harp_minocacontrol_validation_out <- validation_gsea_function(indeseqtable = harp_minocacontrol_deseq, 
                                                              validationlabel = "harp_minoca_v_control",
                                                              validationoutfilepath = paste0(outfilepathmaster, "validation/", "harp_minocacontrol/"))
# MICAD
harp_micadcontrol_validation_out <- validation_gsea_function(indeseqtable = harp_micadcontrol_deseq, 
                                                             validationlabel = "harp_micad_v_control",
                                                             validationoutfilepath = paste0(outfilepathmaster, "validation/", "harp_micadcontrol/"))


# Write out the intermediates:
write_out_intermediate_files(rawcounttable = harp_minocacontrol_rawcount_table, metatable = harp_minocacontrol_metatable, 
                             compcols = compcols, controlcols = controlcols, 
                             outfilepath = paste0(outfilepathmaster, "validation/", "harp_allmicontrol/", "intermediate_files/"), 
                             labelparam = "harp_allmicontrol")


# --------------------------------- NYU COVID v control ---------------------------------

# 3 COVID v Control
dir.create(paste0(outfilepathmaster, "validation/", "covid_covidcontrol/"), showWarnings = FALSE, recursive = TRUE)

covid_covidcontrol_rawcount_file <- file.path(relative_path, "covid-platelet/output/run3/rna_processing/filtrawcounttab.txt")
covid_covidcontrol_rawcount_table <- read.table(covid_covidcontrol_rawcount_file, sep = "\t", header = TRUE, row.names = 1, check.names = FALSE)
covid_covidcontrol_meta_file <- file.path(relative_path, "covid-platelet/output/run3/rna_processing/metatable_in.csv")
covid_covidcontrol_metatable <- read.table(covid_covidcontrol_meta_file, sep =",", header = TRUE, row.names = 1, comment.char = "", quote = "")
SOI <- rownames(covid_covidcontrol_metatable)

## Run all of our DESeq analyses based on the various comps I put in
compcols <- covid_covidcontrol_metatable[SOI,c("comp_covid__covid_v_control"), drop = FALSE]
controlcols <- covid_covidcontrol_metatable[SOI,c("Sex", "Age_cat", "Race")]
covid_covidcontrol_deseqout <- DESeq_analysis(compcols = compcols[SOI,,drop=FALSE], controlcols = controlcols[SOI,],
                                               rawcounttab = covid_covidcontrol_rawcount_table[,SOI])
covid_covidcontrol_deseq <- covid_covidcontrol_deseqout[[1]]

write.table(covid_covidcontrol_deseq, paste0(outfilepathmaster, "validation/", "covid_covidcontrol/", "harp_minocacontrol_wcontrol_deseq.csv"),
            sep = ",", col.names = NA, row.names = TRUE)
covid_covidcontrol_validation_out <- validation_gsea_function(indeseqtable = covid_covidcontrol_deseq,
                                                              validationlabel = "covid_covidcontrol",
                                                              validationoutfilepath = paste0(outfilepathmaster, "validation/", "covid_covidcontrol/"))
# Write out the intermediates:
write_out_intermediate_files(rawcounttable = covid_covidcontrol_rawcount_table, metatable = covid_covidcontrol_metatable, 
                             compcols = compcols, controlcols = controlcols, 
                             outfilepath = paste0(outfilepathmaster, "validation/", "covid_covidcontrol/", "intermediate_files/"), 
                             labelparam = "covid_covidcontrol")

# --------------------------------- SLE validation ---------------------------------

##### SLE Control and Prot
dir.create(paste0(outfilepathmaster, "validation/", "sle_slecontrol/"), showWarnings = FALSE, recursive = TRUE)
dir.create(paste0(outfilepathmaster, "validation/", "sle_sleprot/"), showWarnings = FALSE, recursive = TRUE)

sle_slecontrol_rawcount_file <- file.path(relative_path, "platelet-sle/output/run10_restart_rm1_20210517/rna_processing/filtrawcounttab.txt")
sle_slecontrol_rawcount_table <- read.table(sle_slecontrol_rawcount_file, sep = "\t", header = TRUE, row.names = 1, check.names = FALSE)
sle_slecontrol_meta_file <- file.path(relative_path, "platelet-sle/output/run10_restart_rm1_20210517/rna_processing/metatable_in.csv")
sle_slecontrol_metatable <- read.table(sle_slecontrol_meta_file, sep =",", header = TRUE, row.names = 1, comment.char = "", quote = "")
SOI <- rownames(sle_slecontrol_metatable)

# Set Control Data
controlcols <- sle_slecontrol_metatable[SOI,"Race",drop=FALSE]
controlcols[,"age_cat"] <- ifelse(sle_slecontrol_metatable[,"Age.VISIT."] >= 47, "47_71",
                           ifelse(sle_slecontrol_metatable[,"Age.VISIT."] >= 36, "36_46",
                           ifelse(sle_slecontrol_metatable[,"Age.VISIT."] >= 19, "19_35", NA)))
compcols <- sle_slecontrol_metatable[SOI,c("comp_sle__sle_v_control", "comp_proteinuria__prot_v_noprot")]
controlcols <- controlcols[SOI,]
sle_slecontrol_deseq_out <- DESeq_analysis(compcols = compcols[SOI,,drop=FALSE], controlcols = controlcols[SOI,],
                                               rawcounttab = sle_slecontrol_rawcount_table[,SOI])
sle_slecontrol_deseq <- sle_slecontrol_deseq_out[[1]]
sle_sleprot_deseq <- sle_slecontrol_deseq_out[[2]]

write.table(sle_slecontrol_deseq, paste0(outfilepathmaster, "validation/", "sle_slecontrol/", "sle_slecontrol_wcontrol_deseq.csv"), 
            sep = ",", col.names = NA, row.names = TRUE)
write.table(sle_sleprot_deseq, paste0(outfilepathmaster, "validation/", "sle_sleprot/", "sle_sleprot_wcontrol_deseq.csv"), 
            sep = ",", col.names = NA, row.names = TRUE)

# SLE Control
sle_slecontrol_validation_out <- validation_gsea_function(indeseqtable = sle_slecontrol_deseq, 
                                                              validationlabel = "sle_slecontrol",
                                                              validationoutfilepath = paste0(outfilepathmaster, "validation/", "sle_slecontrol/"))
# SLE Prot
sle_sleprot_validation_out <- validation_gsea_function(indeseqtable = sle_sleprot_deseq, 
                                                             validationlabel = "sle_sleprot",
                                                             validationoutfilepath = paste0(outfilepathmaster, "validation/", "sle_sleprot/"))

# Write out the intermediates:
write_out_intermediate_files(rawcounttable = sle_slecontrol_rawcount_table, metatable = sle_slecontrol_metatable, 
                             compcols = compcols, controlcols = controlcols, 
                             outfilepath = paste0(outfilepathmaster, "validation/", "sle_slecontrol/", "intermediate_files/"), 
                             labelparam = "sle_slecontrol")

# --------------------------------- Utah COVID Validation ---------------------------------

## EXTERNAL VALIDATION - Utah COVID
dir.create(paste0(outfilepathmaster, "validation/", "utahcovid_covidcontrol/"), showWarnings = FALSE, recursive = TRUE)
utahcovid_covidcontrol_deseq_file <- file.path(relative_path, "covid-platelet/docs/Supplemental_Table_4_ICUandnonICU_v_control.txt")
utahcovid_covidcontrol_deseq <- read.table(utahcovid_covidcontrol_deseq_file, sep = "\t", header = TRUE, row.names = 2)[,c(c(5:10))]
write.table(utahcovid_covidcontrol_deseq, paste0(outfilepathmaster, "validation/", "utahcovid_covidcontrol/", "utahcovid_covidcontrol_deseq.csv"), 
            sep = ",", col.names = NA, row.names = TRUE)
utahcovid_covidcontrol_validation_out <- validation_gsea_function(indeseqtable = utahcovid_covidcontrol_deseq, 
                                                              validationlabel = "utahcovid_covidcontrol",
                                                              validationoutfilepath = paste0(outfilepathmaster, "validation/", "utahcovid_covidcontrol/"))


# --------------------------------- Bonfiovanni MP vs RP validation ---------------------------------

## Another validation - looking at the Bongiovanni 2019 RP dataset, just to see where this lines up?
dir.create(paste0(outfilepathmaster, "validation/", "bongiovanni_rpmp/"), showWarnings = FALSE, recursive = TRUE)
rpdeseq_table_file <- "data/Bongiovanni_2019_RP_DGE_list.csv"
rpdeseq_table_temp <- read.table(rpdeseq_table_file, sep = ",", header = TRUE, row.names = 1, quote = "", strip.white = TRUE)
rpdeseq_table_temp[,"Pvalue"] <- gsub("x10", "e", rpdeseq_table_temp[,"Pvalue"])
rpdeseq_table <- apply(rpdeseq_table_temp, 2, function(x) as.numeric(as.character(x)))
dimnames(rpdeseq_table) <- list(rownames(rpdeseq_table_temp), c("log2FoldChange", "pvalue"))
## Have to turn this into a dummy table to match the format:
bongiovanni_rpmp_deseq <- cbind.data.frame(rpdeseq_table, data.frame("baseMean" = NA, "lfcSE" = NA,"stat" = NA, "padj" = NA))[,c("baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj")]

write.table(bongiovanni_rpmp_deseq, paste0(outfilepathmaster, "validation/", "bongiovanni_rpmp/", "bongiovanni_rpmp_deseq.csv"), 
            sep = ",", col.names = NA, row.names = TRUE)
bongiovanni_rpmp_validation_out <- validation_gsea_function(indeseqtable = bongiovanni_rpmp_deseq, 
                                                                  validationlabel = "bongiovanni_rpmp",
                                                                  validationoutfilepath = paste0(outfilepathmaster, "validation/", "bongiovanni_rpmp/"))



# --------------------------------- PACE event validation ---------------------------------

## Additional validation
# Take our PACE cohort - and do our various event analyses, and then again, run these gsea tests with our custom UP and DOWN, do any enrich in the differences between events?
dir.create(paste0(validationoutfilepath, "pace_event_deseq/"), recursive = TRUE, showWarnings = FALSE)
## Run all of our DESeq analyses based on the various comps I put in
pace_platelet_rawcounttab_file <- "output/run7_rmoutliers2_agesexcontrol_withpaired_20220422/rna_processing/filtrawcounttab.txt"
pace_platelet_rawcounttab <- read.table(pace_platelet_rawcounttab_file, sep = "\t", header = TRUE, row.names = 1, check.names = FALSE)
pace_platelet_metatable_file <- "output/run7_rmoutliers2_agesexcontrol_withpaired_20220422/rna_processing/metatable_in.csv"
pace_platelet_metatable <- read.table(pace_platelet_metatable_file, sep = ",", header = TRUE, row.names = 1, check.names = FALSE, quote = "")
SOI <- rownames(pace_pacethr_metatable)

# grab events for a year (rough) and see for each event how does it do?
EOI <- unique(unlist(strsplit(colnames(pace_platelet_metatable[grepl("PACE_", rownames(pace_platelet_metatable)), 
                                                               grepl("censor_|time_to_", colnames(pace_platelet_metatable))]), 
                              split = "time_to_|censor_")))
EOI <- EOI[!EOI %in% c("")]
## Lets do comparisons of every event within 1 year vs NO events at all for AT least one year
noeventnum <- 365
eventnum <- 365
SOI <- rownames(pace_platelet_metatable)[grepl("PACE_", rownames(pace_platelet_metatable)) & !grepl("-3", rownames(pace_platelet_metatable))]
compoutlist <- list()
for (EOInum in seq_len(length(EOI))) {
  EOIsel <- EOI[EOInum]
  compcolsel <- pace_platelet_metatable[SOI, c(paste0("censor_", EOIsel), paste0("time_to_", EOIsel))]
  compout <- compcolsel
  compout[,paste0("comp_pace", EOIsel, "__", EOIsel, "_v_no", EOIsel)] <- 
    ifelse(compout[,grepl("censor_", colnames(compout))] == 1 & compout[,grepl("time_to_", colnames(compout))] <= eventnum, 1, 
    ifelse(compout[,grepl("censor_", colnames(compout))] == 0 & compout[,grepl("time_to_", colnames(compout))] > noeventnum, 0, NA))
  compoutlist[[EOInum]] <- compout[,3,drop=FALSE]
}
compcols <- do.call(cbind, compoutlist)
write.table(t(apply(compcols, 2, function(x) table(x, useNA="ifany"))), 
            paste0(validationoutfilepath, "pace_event_deseq/", "event_tabulation_table.csv"), sep = ",", col.names = NA, row.names = TRUE)
controlcols <- pace_pacethr_metatable[SOI,c("sex1","age_cat", "race")]

DEseqreslist_pace_events = DESeq_analysis(compcols = compcols[SOI,], controlcols = controlcols, rawcounttab = pace_platelet_rawcounttab[,SOI])

# Write out the intermediates:
write_out_intermediate_files(rawcounttable = pace_platelet_rawcounttab[,SOI], metatable = pace_platelet_metatable[SOI,], 
                             compcols = compcols, controlcols = controlcols, 
                             outfilepath = paste0(outfilepathmaster, "validation/", "pace_event_deseq/", "intermediate_files/"), 
                             labelparam = "pace_event_deseq")

## Write out DEseq analysis
outfilepathDEseq <- paste0(validationoutfilepath, "pace_event_deseq/")
GSEAoutlist_UP <- GSEAoutlist_DN <- GSEAfulloutlist <- list()
for (DEseqrestab in seq_len(length(DEseqreslist_pace_events))) {
  DEseqsel <- DEseqreslist_pace_events[[DEseqrestab]]
  DEseq_label <- names(DEseqreslist_pace_events)[DEseqrestab]
  deseqanalysisoutfolder = paste0(outfilepathDEseq, DEseq_label, "/")
  dir.create(deseqanalysisoutfolder, recursive = TRUE, showWarnings = FALSE)
  deseqanalysisoutfile = paste0(deseqanalysisoutfolder, "deseq_results_", DEseq_label, ".csv")
  # write.table(DEseqsel, deseqanalysisoutfile, quote = FALSE, sep = ",", row.names = TRUE, col.names=NA)
  
  write.table(DEseqsel, deseqanalysisoutfile, quote = FALSE, sep = ",", col.names = NA, row.names = TRUE)
  DEseqsel_validation_out <- validation_gsea_function(indeseqtable = DEseqsel, 
                                                      validationlabel = DEseq_label,
                                                      validationoutfilepath = deseqanalysisoutfolder)
  
  # Save out our GSEA results so we can compare later
  custom_hypergseatest <- DEseqsel_validation_out$hypergsea_out
  custom_fullgseatest <- DEseqsel_validation_out$fullgsea_out
  
  ## Need to write out the GSEA results too to quickly look over:
  GSEAoutlist_UP[[DEseqrestab]] <- custom_hypergseatest$enricherUPout
  GSEAoutlist_DN[[DEseqrestab]] <- custom_hypergseatest$enricherDOWNout
  GSEAfulloutlist[[DEseqrestab]] <- custom_fullgseatest
  names(GSEAoutlist_UP)[DEseqrestab] <- names(GSEAoutlist_DN)[DEseqrestab] <- names(GSEAfulloutlist)[DEseqrestab] <- 
    names(DEseqreslist_pace_events)[DEseqrestab]
  
}

## Now run through the GSEAUP and GSEADOWN and look for my GSOI
# GSEA_UP_summarytable <- do.call(rbind, lapply(GSEAoutlist_UP, function(x) x[c("custom_mgc_hyper_up", "custom_mgc_hyper_down"),])) ## down is never really up (opp), so we dont need to include
GSEA_UP_summarytable <- do.call(rbind, lapply(GSEAoutlist_UP, function(x) x[c("custom_mgc_hyper_up"),])) 
write.table(GSEA_UP_summarytable[order(GSEA_UP_summarytable[,"pvalue"]),c("ID", "GeneRatio", "BgRatio", "pvalue", "p.adjust")], paste0(outfilepathDEseq, "GSEA_UP_mgccustom_summarytable.csv"), sep = ",", col.names = NA, row.names = TRUE)
GSEA_DN_summarytable <- do.call(rbind, lapply(GSEAoutlist_DN, function(x) x[c("custom_mgc_hyper_down"),])) ## up is never really down (opp), so we dont need to include
write.table(GSEA_DN_summarytable[order(GSEA_DN_summarytable[,"pvalue"]),c("ID", "GeneRatio", "BgRatio", "pvalue", "p.adjust")], paste0(outfilepathDEseq, "GSEA_DN_mgccustom_summarytable.csv"), sep = ",", col.names = NA, row.names = TRUE)
# So both the up is enriched up and the down is enriched down
GSEAfull_summarytable <- do.call(rbind, lapply(GSEAfulloutlist, function(x) x[c("custom_mgc_hyper_up", "custom_mgc_hyper_down"),]))
write.table(GSEAfull_summarytable[order(GSEAfull_summarytable[,"pvalue"]),c("ID", "pvalue", "p.adjust", "NES", "core_enrichment")], 
            paste0(outfilepathDEseq, "GSEAfull_mgccustom_summarytable.csv"), sep = ",", col.names = NA, row.names = TRUE)




# --------------------------------- Final validation summary figures ---------------------------------

# Grab all results and summarize as a dot plot:
# "/Users/tosh/Desktop/Ruggles_Lab/projects/platelet-pace/output/hyper_geneset_creation/run12_hyper60_hypo40_AGRCONTROL/validation/pace_pacethr/pace_pace_v_thr_GO_w_custom_fullgsea.csv"

validation_list <- list(run1 = c(label = "pace_pacethr", file = "pace_pace_v_thr_GO_w_custom_fullgsea.csv"),
                        run2 = c(label = "harp_minocacontrol", file = "harp_minoca_v_control_GO_w_custom_fullgsea.csv"),
                        run3 = c(label = "harp_micadcontrol", file = "harp_micad_v_control_GO_w_custom_fullgsea.csv"),
                        run4 = c(label = "harp_allmicontrol", file = "harp_allmi_v_control_GO_w_custom_fullgsea.csv"),
                        run5 = c(label = "covid_covidcontrol", file = "covid_covidcontrol_GO_w_custom_fullgsea.csv"),
                        run6 = c(label = "utahcovid_covidcontrol", file = "utahcovid_covidcontrol_GO_w_custom_fullgsea.csv"),
                        run7 = c(label = "sle_slecontrol", file = "sle_slecontrol_GO_w_custom_fullgsea.csv"),
                        run8 = c(label = "sle_sleprot", file = "sle_sleprot_GO_w_custom_fullgsea.csv"),
                        run9 = c(label = "bongiovanni_rpmp", file = "bongiovanni_rpmp_GO_w_custom_fullgsea.csv"),
                        run10 = c(label = "pace_event_deseq", file = "GSEAfull_mgccustom_summarytable.csv") ## This is an exception - this is the full folder of event analysis
                        )
grabtable_list <- list()
for (validation_sel in validation_list) {
    label_sel <- validation_sel[["label"]]
    file_sel <- validation_sel[["file"]]
    table_sel <- read.table(paste0(outfilepathmaster, "validation/", label_sel, "/", file_sel), sep = ",", header = TRUE, row.names = 1)
    
    grabtable <- table_sel[grepl("mgc", table_sel[,1]), c("NES", "pvalue", "p.adjust")]
    grabtable_list[[label_sel]] <- cbind(validation_run = label_sel, geneset = rownames(grabtable), grabtable)
}
all_mgc_validation_table <- do.call(rbind, grabtable_list)

selected_validation_table <- all_mgc_validation_table[
    c("pace_pacethr.custom_mgc_hyper_up",
    "pace_pacethr.custom_mgc_hyper_down",
    "pace_event_deseq.comp_paceMACE__MACE_v_noMACE.custom_mgc_hyper_up",
    "pace_event_deseq.comp_paceMACE__MACE_v_noMACE.custom_mgc_hyper_down",
    "pace_event_deseq.comp_paceMALE__MALE_v_noMALE.custom_mgc_hyper_up",
    "pace_event_deseq.comp_paceMALE__MALE_v_noMALE.custom_mgc_hyper_down",
    "covid_covidcontrol.custom_mgc_hyper_up",
    "covid_covidcontrol.custom_mgc_hyper_down",
    "sle_slecontrol.custom_mgc_hyper_up",
    "sle_slecontrol.custom_mgc_hyper_down",
    "bongiovanni_rpmp.custom_mgc_hyper_up",
    "bongiovanni_rpmp.custom_mgc_hyper_down"),]


dotplottab <- selected_validation_table[,c("validation_run", "geneset", "NES", "pvalue")]
dotplottab[,"dataset"] <- factor(ifelse(grepl("pace_pacethr", dotplottab[,"validation_run"]), "Atherosclerosis",
                          ifelse(grepl("comp_paceMACE__MACE_v_noMACE", dotplottab[,"geneset"]), "MACE",
                          ifelse(grepl("comp_paceMALE__MALE_v_noMALE", dotplottab[,"geneset"]), "MALE",
                          ifelse(grepl("covid_covidcontrol", dotplottab[,"validation_run"]), "COVID",
                          ifelse(grepl("sle_slecontrol", dotplottab[,"validation_run"]), "Systemic Lupus Erythematosus",
                          ifelse(grepl("bongiovanni_rpmp", dotplottab[,"validation_run"]), "Reticulated Platelets",
                                 NA)))))), levels = rev(c("Atherosclerosis", "MACE", "MALE", "COVID", "Systemic Lupus Erythematosus", "Reticulated Platelets")))
dotplottab[,"geneset"] <- factor(ifelse(grepl("custom_mgc_hyper_up", dotplottab[,"geneset"]), "PRESS UP",
                                 ifelse(grepl("custom_mgc_hyper_down", dotplottab[,"geneset"]), "PRESS DOWN",
                                                                    NA)), levels = c("PRESS UP", "PRESS DOWN"))
dotplottab[,"absNES"] <- abs(dotplottab[,"NES"])
dotplottab[,"signed_pstat"] <- -log10(dotplottab[,"pvalue"]) * sign(dotplottab[,"NES"])
dotplottab[,"pstat"] <- -log10(dotplottab[,"pvalue"]) 

# "#9ac2e6" # light blue
# "#426484" # dark blue
# "#eca3a3" # light red
# "#a20000" # dark red

# read in the old table
library(glue)
dotplottab <- read.csv(glue("{outfilepathmaster}validation/PRESS_set_validation_dotplot_table.csv"))
dotplottab <- dotplottab %>%
  filter(!dataset %in% c('MALE', 'MACE')) %>%
  mutate(dataset = factor(dataset, levels = c('Reticulated Platelets','COVID', 'Systemic Lupus Erythematosus', 'Atherosclerosis')))

   

pout <- ggplot(dotplottab, aes(x=dotplottab[,"geneset"], y=dotplottab[,"dataset"], size=dotplottab[,"pstat"], color=dotplottab[,"NES"]))
# pout <- ggplot(modulesel, aes(x=modulesel[,1], y=modulesel[,3], size=modulesel[,2], color=modulesel[,2]))
pout <- pout + geom_point(alpha = 1, shape = 16)
pout <- pout + scale_color_gradientn(name = "NES", colors = c("#426484", "#9ac2e6", "#eca3a3", "#a20000"),
                                     values = scales::rescale(c(-3.6, -0.0001, 0.0001, 3.6), to=c(0,1)),
                                     limits = c(-3.6,3.6))
pout <- pout + scale_size_continuous(name = "-log(pvalue)", 
                                     range = c(1, 20), 
                                     limits = c(1, 10),
                                     breaks = c(-log10(0.01), -log10(0.001), -log10(0.0001)), 
                                     labels = c(2, 3, 4)
                                    )
pout <- pout + labs(title = "PRESS Set GSEA Results", x = "PRESS Set", y = "Dataset", color = "p. value", size = "NES")
pout <- pout + theme_bw(base_size = 24) + theme(legend.position = "right", axis.text.x = element_text(angle = 0))
pout

write.table(dotplottab, paste0(outfilepathmaster, "validation/", "PRESS_set_validation_dotplot_table.csv"), sep = ",", col.names = NA, row.names = TRUE)
pdf(width = 15, height = 10, paste0(outfilepathmaster, "validation/", "PRESS_set_validation_dotplot.pdf"))
print(pout)
junk <- dev.off()




# --------------------------------- END ---------------------------------





# --------------------------------- Junk Code ---------------------------------

# pout <- ggplot(dotplottab, aes(x=dotplottab[,"geneset"], y=dotplottab[,"dataset"], size=dotplottab[,"absNES"], color=dotplottab[,"signed_pstat"]))
# # pout <- ggplot(modulesel, aes(x=modulesel[,1], y=modulesel[,3], size=modulesel[,2], color=modulesel[,2]))
# pout <- pout + geom_point(alpha = 0.8, shape = 16)
# pout <- pout + scale_color_gradientn(name = "p. value", colors = c("#426484", "#9ac2e6", "#eca3a3", "#a20000"),
#                                      values = scales::rescale(c(-10, -0.0001, 0.0001, 10), to=c(0,1)),
#                                      limits = c(-10,10))
# pout <- pout + scale_size_continuous(name = "NES", range = c(5, 20), limits = c(1,3))
# pout <- pout + labs(title = "PRESS Set GSEA Results", x = "PRESS Set", y = "Dataset", color = "p. value", size = "NES")
# pout <- pout + theme_bw() theme(legend.position = "bottom", axis.text.x = element_text(angle = 90))
# pout



# ## Model to predict epi score from the genes I have
# library(caret)
# library(xgboost)
# fullmodelingdata <- na.omit(merge(t(normcounttable[GOI,]), epi_metatable, by = "row.names"))
# rownames(fullmodelingdata) <- fullmodelingdata[,"FINAL_PACE_ID"]
# fullmodelingdata <- fullmodelingdata[,!colnames(fullmodelingdata) %in% c("FINAL_PACE_ID", "Row.names")]
# 
# # Are we predicting the epi value or the cohort?
# classvariable <- "epi_04um_300s_n"
# outcomelabels <- "hypercohort"
# featurelabels <- GOI
# modeling_indata <- fullmodelingdata[,c(featurelabels, outcomelabels)]
# set.seed(13345)
# 
# # Set test and train
# test_index <- createDataPartition(modeling_indata[,outcomelabels], p = 0.6, list = FALSE)
# testSet <- modeling_indata[-test_index,]
# trainSet <- modeling_indata[test_index,]
# 
# # Split train into features and outcomes
# x_train = as.matrix(trainSet[,!colnames(trainSet) %in% outcomelabels])
# y_train = as.factor(trainSet[,outcomelabels])
# 
# #### Generic control parametrs
# ctrl <- trainControl(method="repeatedcv", 
#                      number=10, 
#                      repeats=5,
#                      savePredictions=TRUE, 
#                      classProbs=TRUE,
#                      summaryFunction = multiClassSummary)
# 
# xgbgrid <- expand.grid(nrounds = 10,
#                        max_depth = 5,
#                        eta = 0.05,
#                        gamma = 0.01,
#                        colsample_bytree = 0.75,
#                        min_child_weight = 0,
#                        subsample = 0.5)
# 
# 
# set.seed(113355123)
# x_train 
# xgb_model = train(x_train, 
#                   y_train,  
#                   trControl = ctrl,
#                   method = "xgbTree",
#                   tuneGrid = xgbgrid)
# xgb_model
# test_predict <- predict(xgb_model, testSet)
# confusionMatrix(test_predict, as.factor(testSet[,outcomelabels]))
# 
# ## Get variable importance
# varimpout <- varImp(xgb_model)$importance
# ggplot(varimpout) + theme_minimal()
# ## Looks like from the quick viz of importance that we want pretty much anything above 1? Cause thats about the top 10% (the rest dropout)
# modelimpGOI <- rownames(varimpout[varimpout[,1] > 1,,drop=FALSE])
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# ## Trying singscore
# library(singscore)
# library(GSEABase)
# 
# ## first rank genes with the rankgene function
# rankData <- rankGenes(normcounttable)
# scoredf <- simpleScore(rankData, upSet = GOIup, downSet = GOIdown)
# wilcox.test(scoredf[grepl("PACE_", rownames(scoredf)),1], scoredf[grepl("THR", rownames(scoredf)),1])
# 
# 
# singscore_stable_bloodgenes <- getStableGenes(n_stable = 10, type = "blood")
# 
# ## get row variances of genes
# scaled_counttable <- data.frame(t(apply(normcounttable, 1, scale)))
# gene_scaleranges <- data.frame(t(apply(scaled_counttable, 1, range)))
# gene_scaleranges[,"range"] <- gene_scaleranges[,2] - gene_scaleranges[,1]
# gene_scaleranges <- gene_scaleranges[order(gene_scaleranges[,"range"]),]
# 
# genestablity_table <- merge(gene_scaleranges, data.frame(rowMeans(normcounttable)), by = "row.names")
# 
# genestablity_table <- Reduce(function(dtf1, dtf2)
#   merge(dtf1, dtf2, by = "gene", all.x = TRUE, sort = FALSE), 
#   list(
#     cbind(gene = rownames(gene_scaleranges), gene_scaleranges), 
#     cbind(gene = rownames(data.frame(rowMeans(normcounttable))), data.frame(rowMeans(normcounttable))), 
#     cbind(gene = rownames(normcounttable), apply(normcounttable, 1, function(x) shapiro.test(x)$p.value)),
#     cbind(gene = rownames(normcounttable), apply(normcounttable, 1, function(x) ks.test(x, pnorm)$p.value))
#   ))
# dimnames(genestablity_table) <- list(genestablity_table[,"gene"], c("gene", "minsd", "maxsd", "rangesd", "genebasemean", "shapiro_pval", "ks_pval"))
# 
# 
# 
# genestablity_table <- genestablity_table[order(genestablity_table[,"rangesd"]),]
# genestablity_table[singscore_stable_bloodgenes,]
# 
# 
# tt1 <- apply(normcounttable, 1, function(x) shapiro.test(x)$p.value)
# tt1 <- apply(normcounttable, 1, function(x) ks.test(x, pnorm)$p.value)
