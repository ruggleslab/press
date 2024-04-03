###########################################################################
#
#                            harp_predictions
#
###########################################################################
# Author: Matthew Muller
# Date: 2023-07-27
# Script Name: harp_predictions
# Output directory:
experiment <- "harp_predictions"
run <- 2
outdir <- file.path('output', paste0(experiment, '__run_', run))
dir.create(outdir, showWarnings = F)

# Notes about the experiment run:
notes <- "Predictions of harp using the PRESS scoring." #nolint

# save notes to file
write(notes, file.path(outdir, "notes.txt"))

#======================== LIBRARIES ========================#
packages <- c("lintr", "httpgd", "languageserver", "devtools", "sys", "dplyr", "tidyverse", "ggplot2", "ggpubr", "SummarizedExperiment", "glue", "readxl")
pkgs <- lapply(packages, function(x) suppressMessages(require(x, character.only=T,quietly=T))) # nolint
print(packages[!sapply(pkgs, isTRUE)])

# LOAD FUNCTIONS
# space reserved for sourcing in functions
source_url('https://raw.githubusercontent.com/mattmuller0/Rtools/main/general_functions.R')
source_url('https://raw.githubusercontent.com/mattmuller0/Rtools/main/plotting_functions.R')
source_url('https://raw.githubusercontent.com/mattmuller0/Rtools/main/converting_functions.R')

#======================== CODE ========================#
# load in the data
harp_rawcounts <- read.csv('data/harp_plt/harp_allmicontrol__rawcounttable.csv', header = TRUE, row.names = 1)
harp_metadata <- read.csv('data/harp_plt/harp_allmicontrol__metatable.csv', header = TRUE, row.names = 1)
harp_metadata_addon <- read.csv('data/harp_plt/harp-platelet-metadata-ADDON-20220823.csv', header = TRUE, row.names = 1)
harp_metadata_addon2 <- read_xlsx('data/harp_plt/Metadata HARP for WB RNAseq 3_14_2022x.xlsx')

# subset the metadata addon to the harp metadata
harp_metadata_addon <- harp_metadata_addon[rownames(harp_metadata),]

# add the addon to the metadata
harp_metadata <- harp_metadata_addon

# gsub the harp metadata names
rownames(harp_metadata) <- gsub("-", "\\.", rownames(harp_metadata))
harp_metadata_addon2$Subject.ID <- gsub("-", "\\.", harp_metadata_addon2$Subject.ID)

# subset to the press genes
press <- read.csv('data/press451_genes.csv', header = FALSE) %>% pull(.)

# add missing rows of press genes
harp_rawcounts %>% add_missing_rows(press) -> harp_rawcounts

# make the se and save it
harp_se <- make_se(harp_rawcounts, harp_metadata)
harp_dds <- DESeqDataSet(harp_se, design = ~ 1) %>% DESeq()
# save_se(harp_dds[press,], file.path(outdir, "harp_se.RData"), normalize = 'mor', log = TRUE)


# Now from here we are going to move to python for predictions.
preds <- read.csv('output/harp_predictions__run_2/harp_hyper_v_hypo_predictions.csv', header = TRUE, row.names = 1)
harp_metadata$press <- scale(preds$harp_preds)

colnames(harp_metadata)

# plot the four groups of harp
p1 <- ggplot(harp_metadata, aes(x = Angiography.Report, y = press, fill = Angiography.Report)) +
  geom_boxplot() +
  labs(x = 'Predicted PRESS', y = 'HARP') +
  theme_matt(18) +
  theme(legend.position = 'none') +
  stat_compare_means()
p1

# list from tessa of who to include
harp_selection <- c(
# MI-CAD   
'HARP-01-0246-1',
'HARP-01-0098-1',
'HARP-01-0052-1',
'HARP-01-0039-1',
'HARP-01-0033-1',
'HARP-01-0035-1',
'HARP-01-109-1',
'HARP-01-0204-1',
'HARP-01-0034-1',
'HARP-01-0135-1',
'HARP-01-0177-1',
'HARP-01-0245-1',
'HARP-01-0236-1',
'HARP-01-0070-1',
'HARP-01-0043-1',
'HARP-01-0055-1',
'HARP-01-0057-1',
'HARP-01-0174-1',
'HARP-01-0241-1',
'HARP-01-0043-1',
# Obstructive
'HARP-01-5061-1',
'HARP-01-5062-1',
'HARP-01-5057-1',
'HARP-01-5012-2',
'HARP-01-5013-1',
'HARP-01-5024-1',
'HARP-01-5027-1',
'HARP-01-5033-1',
'HARP-01-5046-1'
) %>% 
gsub("-", "\\.", .)

# get a list of the IDs we are missing from the metadata
harp_selection[!harp_selection %in% rownames(harp_metadata)]

# change HARP.01.0135.1 to MI-CAD
harp_metadata[rownames(harp_metadata) == 'HARP.01.0135.1', 'Angiography.Report'] <- 'MI-CAD'

# merge the four groups into just MI and Control
harp_metadata <- harp_metadata %>%
    # filter(rownames(.) %in% harp_selection) %>%
    mutate(
        MI_v_Ctrl = ifelse(Angiography.Report == 'MI-CAD', 'MI-CAD',
        ifelse(Angiography.Report != 'MINOCA', 'Control', NA)),
        
        MI_v_Obstr = ifelse(Angiography.Report == 'MI-CAD', 'MI-CAD',
        ifelse(Angiography.Report == 'Obstructive', 'Obstructive', NA)),
        
        MINOCA_v_Ctrl = case_when(
            Angiography.Report == 'MINOCA' ~ 'MINOCA',
            Angiography.Report == 'Obstructive' ~ 'Control',
            Angiography.Report == 'Non-obstructive' ~ 'Control',
            TRUE ~ NA_character_
        
        )
    )
# harp_metadata$MI_v_Ctrl <- factor(harp_metadata$MI_v_Ctrl, levels = c('Control', 'MI-CAD'))
# harp_metadata$MI_v_Obstr <- factor(harp_metadata$MI_v_Obstr, levels = c('Obstructive', 'MI-CAD'))

# add the bmi
harp_metadata$bmi <- harp_metadata_addon2 %>%
  as.data.frame() %>%
  column_to_rownames('Subject.ID') %>%
  .[rownames(harp_metadata), 'BMI']

# plotting
p2 <- ggplot(filter(harp_metadata, !is.na(MI_v_Ctrl)), aes(x = MI_v_Ctrl, y = press, fill = MI_v_Ctrl)) +
  geom_boxplot() +
  labs(x = 'PRESS Score', y = 'HARP') +
  theme_matt(18) +
  theme(legend.position = 'none') +
  stat_compare_means(method = 't.test') +
  labs(y = NULL, x = NULL)
ggsave(file.path(outdir, "harp_predictions_mi_v_control.pdf"), p2)

p3 <- ggplot(filter(harp_metadata, !is.na(MI_v_Obstr)), aes(x = MI_v_Obstr, y = press, fill = MI_v_Obstr)) +
  geom_boxplot() +
  theme_matt(18) +
  theme(legend.position = 'none') +
  stat_compare_means(method = 't.test') +
  labs(y = 'PRESS Score')
ggsave(file.path(outdir, "harp_predictions_mi_v_obstr.pdf"), p3)

p4 <- ggplot(filter(harp_metadata, !is.na(MINOCA_v_Ctrl)), aes(x = MINOCA_v_Ctrl, y = press, fill = MINOCA_v_Ctrl)) +
  geom_boxplot() +
  theme_matt(18) +
  theme(legend.position = 'none') +
  stat_compare_means(method = 't.test') +
  labs(y = NULL, x = NULL)
ggsave(file.path(outdir, "harp_predictions_minoca_v_control.pdf"), p4)

# plot the two plots and save
plots <- plot_grid(p2, p3, ncol = 2)

# save the plot
ggsave(file.path(outdir, "harp_predictions.pdf"), p3)

# make a table of this all
with(harp_metadata[harp_selection,], table(Angiography.Report))

# get an adjusted odds ratio of the two groups
library(oddsratio) # nolint
library(mfx) # nolint

# prep our data
data_df <- harp_metadata %>% 
    filter(!is.na(MI_v_Obstr)) %>%
    mutate(
      race = ifelse(
        White == 'Yes', 'White',
          ifelse(
            Black.or.African.American == 'Yes', 'Black',
            'Other'
          )
      )
    )

# get the year of birth from the date of birth
data_df$yob <- strsplit(data_df$Date.of.Birth, split = '/') %>% 
  sapply(function(x) x[3]) %>% 
  as.numeric()
data_df$age <- 121 - data_df$yob

# refactor the MI_v_Obstr
data_df$MI_v_Obstr <- factor(data_df$MI_v_Obstr, levels = c('Obstructive', 'MI-CAD'))
data_df$Diabetes[data_df$Diabetes == ''] <- 'missing'
data_df$Diabetes <- factor(data_df$Diabetes, levels = c('No', 'missing', 'Yes'))
data_df$Prior.Stroke.TIA <- factor(data_df$Prior.Stroke.TIA, levels = c('Yes', 'No'))
data_df$missing <- ifelse(is.na(data_df$bmi), 'missing', 'not missing')
data_df$bmi <- as.numeric(data_df$bmi)
data_df$bmi[is.na(data_df$bmi) & data_df$MI_v_Obstr == 'Obstructive'] <- mean(data_df$bmi[data_df$MI_v_Obstr == 'Obstructive'], na.rm = T)
data_df$bmi[is.na(data_df$bmi) & data_df$MI_v_Obstr == 'MI-CAD'] <- mean(data_df$bmi[data_df$MI_v_Obstr == 'MI-CAD'], na.rm = T)
colnames(data_df)

data_df <- data_df %>% drop_na(MI_v_Obstr)

# Ethnicity + Prior.Stroke.TIA + bmi + 
# try an mfx model OR
model <- glm(
  MI_v_Obstr ~ press +race + Ethnicity + age + bmi + missing, 
  data = data_df,
  family = binomial(link = 'logit')
  )
summary(model)

# get the odds ratio
or_adj <- cbind("Odds ratio" = coef(model), confint.default(model, level = 0.95))['press',]
or_adj$pvalue <- summary(model)$coefficients['press', 'Pr(>|z|)']
or_adj <- as.data.frame(or_adj)
or_adj$beta <- beta_adj <- summary(model)$coefficients['press', 'Estimate']
or_adj$adjustments <- paste0('Adjusted for: ', paste0('race + age + Ethnicity + Hypertension + Diabetes + bmi + Prior.Stroke.TIA'))
colnames(or_adj) <- c('Beta', '2.5%', '97.5%', 'p-value', 'beta', 'adjustments')

# write the or to file
write.csv(or_adj, file.path(outdir, "beta.csv"), row.names = F)

# the statement of choice
library(glue)
glue('After multivariable adjustment for age, race/ethnicity, BMI, diabetes, and stroke, PRESS was significantly associated with acute MI (β={or_adj$Beta}, 95% CI {or_adj$`2.5%`} to {or_adj$`97.5%`}, p={or_adj$`p-value`}).') # nolint


#======================== END ========================#
save.image(file.path(outdir, "image.RData"))

sink(file.path(outdir, "session.log"), append = TRUE, split = TRUE)
sessionInfo()
sink()

# confidence interval for beta val
confint.default(model, level = 0.95)