Wdir = "Your directory"
setwd(Wdir)

library(phyloseq)
# library(speedyseq)
library(ANCOMBC)
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggvenn)
library(ggh4x)
library(broman)
library(ggtext)
library(tidyverse)
library(patchwork)

###################### 16S ######################################################
rm(list = ls())

a16S = readRDS("./RDSdata/NA_free_THREEmonths_bsthit_pruned_16S_org.rds")

sample_data(a16S) = data.frame(a16S@sam_data)

class(a16S@sam_data)

b16.vv.6 = subset_samples(a16S, Cultivar == "VV" & Month == "June")
b16.vv.8 = subset_samples(a16S, Cultivar == "VV" & Month == "August")
b16.vv.10 = subset_samples(a16S, Cultivar == "VV" & Month == "October")

b16.jb.6 = subset_samples(a16S, Cultivar == "JB" & Month == "June")
b16.jb.8 = subset_samples(a16S, Cultivar == "JB" & Month == "August")
b16.jb.10 = subset_samples(a16S, Cultivar == "JB" & Month == "October")

pdj = "fdr"
# dir.create("./ANCOMBC/B16S", recursive = T)
# dir.create("./ANCOMBC/ITS", recursive = T)
packageVersion("ANCOMBC")
ancom.16S.vv.6 = ANCOMBC::ancombc2(b16.vv.6,
                                    tax_level = "Genus",
                                    fix_formula  = "Chemical",
                                    p_adj_method = pdj,
                                    pseudo = 0,
                                    pseudo_sens = T,
                                    # struc_zero = T,
                                    # neg_lb = T,
                                    iter_control = list(tol = 1e-5, max_iter = 500, verbose = FALSE),
                                    em_control = list(tol = 1e-5, max_iter = 500),
                                    n_cl = 4,
                                    verbose = T)

saveRDS(ancom.16S.vv.6, "./ANCOMBC/B16S/bac_VV_JUNE.rds")

ancom.16S.vv.8 = ANCOMBC::ancombc2(b16.vv.8,
                                   tax_level = "Genus",
                                   fix_formula  = "Chemical",
                                   p_adj_method = pdj,
                                   pseudo = 0,
                                   pseudo_sens = T,
                                   # struc_zero = T,
                                   # neg_lb = T,
                                   iter_control = list(tol = 1e-5, max_iter = 500, verbose = FALSE),
                                   em_control = list(tol = 1e-5, max_iter = 500),
                                   n_cl = 4,
                                   verbose = T)

saveRDS(ancom.16S.vv.8, "./ANCOMBC/B16S/bac_VV_AUGUST.rds")


ancom.16S.vv.10 = ANCOMBC::ancombc2(b16.vv.10,
                                   tax_level = "Genus",
                                   fix_formula  = "Chemical",
                                   p_adj_method = pdj,
                                   pseudo = 0,
                                   pseudo_sens = T,
                                   iter_control = list(tol = 1e-5, max_iter = 500, verbose = FALSE),
                                   em_control = list(tol = 1e-5, max_iter = 500),
                                   n_cl = 4,
                                   verbose = T)

saveRDS(ancom.16S.vv.10, "./ANCOMBC/B16S/bac_VV_OCT.rds")

############ justin 
ancom.16S.jb.6 = ANCOMBC::ancombc2(b16.jb.6,
                                   tax_level = "Genus",
                                   fix_formula  = "Chemical",
                                   p_adj_method = pdj,
                                   pseudo = 0,
                                   pseudo_sens = T,
                                   iter_control = list(tol = 1e-5, max_iter = 500, verbose = FALSE),
                                   em_control = list(tol = 1e-5, max_iter = 500),
                                   n_cl = 4,
                                   verbose = T)

saveRDS(ancom.16S.jb.6, "./ANCOMBC/B16S/bac_JB_JUNE.rds")

ancom.16S.jb.8 = ANCOMBC::ancombc2(b16.jb.8,
                                   tax_level = "Genus",
                                   fix_formula  = "Chemical",
                                   p_adj_method = pdj,
                                   pseudo = 0,
                                   pseudo_sens = T,
                                   iter_control = list(tol = 1e-5, max_iter = 500, verbose = FALSE),
                                   em_control = list(tol = 1e-5, max_iter = 500),
                                   n_cl = 4,
                                   verbose = T)

saveRDS(ancom.16S.jb.8, "./ANCOMBC/B16S/bac_JB_AUGUST.rds")


ancom.16S.jb.10 = ANCOMBC::ancombc2(b16.jb.10,
                                    tax_level = "Genus",
                                    fix_formula  = "Chemical",
                                    p_adj_method = pdj,
                                    pseudo = 0,
                                    pseudo_sens = T,
                                    iter_control = list(tol = 1e-5, max_iter = 500, verbose = FALSE),
                                    em_control = list(tol = 1e-5, max_iter = 500),
                                    n_cl = 4,
                                    verbose = T)

saveRDS(ancom.16S.jb.10, "./ANCOMBC/B16S/bac_JB_OCT.rds")

#################################################################################
rm(list = ls())
ITS = readRDS("./RDSdata/NA_free_THREEmonths_bsthit_pruned_ITS_org.rds")

ITS.vv.6  = subset_samples(ITS, Cultivar == "VV" & Month == "June")
ITS.vv.8  = subset_samples(ITS, Cultivar == "VV" & Month == "August")
ITS.vv.10 = subset_samples(ITS, Cultivar == "VV" & Month == "October")

ITS.jb.6  = subset_samples(ITS, Cultivar == "JB" & Month == "June")
ITS.jb.8  = subset_samples(ITS, Cultivar == "JB" & Month == "August")
ITS.jb.10 = subset_samples(ITS, Cultivar == "JB" & Month == "October")

pdj = "fdr"
# dir.create("./ANCOMBC/B16S", recursive = T)
# dir.create("./ANCOMBC/ITS", recursive = T)

ancom.ITS.vv.6 = ANCOMBC::ancombc2(ITS.vv.6,
                                   tax_level = "Genus",
                                   fix_formula  = "Chemical",
                                   p_adj_method = pdj,
                                   pseudo = 0,
                                   pseudo_sens = T,
                                   iter_control = list(tol = 1e-5, max_iter = 500, verbose = FALSE),
                                   em_control = list(tol = 1e-5, max_iter = 500),
                                   n_cl = 4,
                                   verbose = T)

saveRDS(ancom.ITS.vv.6, "./ANCOMBC/ITS/ITS_VV_JUNE.rds")

#-
ancom.ITS.vv.8 = ANCOMBC::ancombc2(ITS.vv.8,
                                   tax_level = "Genus",
                                   fix_formula  = "Chemical",
                                   p_adj_method = pdj,
                                   pseudo = 0,
                                   pseudo_sens = T,
                                   iter_control = list(tol = 1e-5, max_iter = 500, verbose = FALSE),
                                   em_control = list(tol = 1e-5, max_iter = 500),
                                   n_cl = 4,
                                   verbose = T)

saveRDS(ancom.ITS.vv.8, "./ANCOMBC/ITS/ITS_VV_AUGUST.rds")


ancom.ITS.vv.10 = ANCOMBC::ancombc2(ITS.vv.10,
                                    tax_level = "Genus",
                                    fix_formula  = "Chemical",
                                    p_adj_method = pdj,
                                    pseudo = 0,
                                    pseudo_sens = T,
                                    iter_control = list(tol = 1e-5, max_iter = 500, verbose = FALSE),
                                    em_control = list(tol = 1e-5, max_iter = 500),
                                    n_cl = 4,
                                    verbose = T)

saveRDS(ancom.ITS.vv.10, "./ANCOMBC/ITS/ITS_VV_OCT.rds")

############ justin 
ancom.ITS.jb.6 = ANCOMBC::ancombc2(ITS.jb.6,
                                   tax_level = "Genus",
                                   fix_formula  = "Chemical",
                                   p_adj_method = pdj,
                                   pseudo = 0,
                                   pseudo_sens = T,
                                   iter_control = list(tol = 1e-5, max_iter = 500, verbose = FALSE),
                                   em_control = list(tol = 1e-5, max_iter = 500),
                                   n_cl = 4,
                                   verbose = T)

saveRDS(ancom.ITS.jb.6, "./ANCOMBC/ITS/ITS_JB_JUNE.rds")

ancom.ITS.jb.8 = ANCOMBC::ancombc2(ITS.jb.8,
                                   tax_level = "Genus",
                                   fix_formula  = "Chemical",
                                   p_adj_method = pdj,
                                   pseudo = 0,
                                   pseudo_sens = T,
                                   iter_control = list(tol = 1e-5, max_iter = 500, verbose = FALSE),
                                   em_control = list(tol = 1e-5, max_iter = 500),
                                   n_cl = 4,
                                   verbose = T)

saveRDS(ancom.ITS.jb.8, "./ANCOMBC/ITS/ITS_jb_AUGUST.rds")


ancom.ITS.jb.10 = ANCOMBC::ancombc2(ITS.jb.10,
                                    tax_level = "Genus",
                                    fix_formula  = "Chemical",
                                    p_adj_method = pdj,
                                    pseudo = 0,
                                    pseudo_sens = T,
                                    iter_control = list(tol = 1e-5, max_iter = 500, verbose = FALSE),
                                    em_control = list(tol = 1e-5, max_iter = 500),
                                    n_cl = 4,
                                    verbose = T)

saveRDS(ancom.ITS.jb.10, "./ANCOMBC/ITS/ITS_JB_OCT.rds")


########################################################################

# plot ANCOM results

########################################################################
rm(list = ls())

source("./tools/visual_tools.R")

# load saved data
ancom.16S.vv.6 = readRDS("./ANCOMBC/B16S/bac_VV_JUNE.rds")
ancom.16S.vv.8 = readRDS("./ANCOMBC/B16S/bac_VV_AUGUST.rds")
ancom.16S.vv.10 = readRDS("./ANCOMBC/B16S/bac_VV_OCT.rds")

ancom.16S.jb.6 = readRDS("./ANCOMBC/B16S/bac_JB_JUNE.rds")
ancom.16S.jb.8 = readRDS("./ANCOMBC/B16S/bac_JB_AUGUST.rds")
ancom.16S.jb.10 = readRDS("./ANCOMBC/B16S/bac_JB_OCT.rds")

ancom.16S.vv.6.df =  reshape_ancom(ancom.16S.vv.6, sel_var = "Chemical")
ancom.16S.vv.8.df =  reshape_ancom(ancom.16S.vv.8, sel_var = "Chemical")
ancom.16S.vv.10.df = reshape_ancom(ancom.16S.vv.10, sel_var = "Chemical")

ancom.16S.vv.6.df$cultivar = "Vardar Valley"
ancom.16S.vv.8.df$cultivar = "Vardar Valley"
ancom.16S.vv.10.df$cultivar = "Vardar Valley"

ancom.16S.vv.6.df$month = "June"
ancom.16S.vv.8.df$month = "August"
ancom.16S.vv.10.df$month = "October"


ancom.16S.jb.6.df =  reshape_ancom(ancom.16S.jb.6, sel_var = "Chemical")
ancom.16S.jb.8.df =  reshape_ancom(ancom.16S.jb.8, sel_var = "Chemical")
ancom.16S.jb.10.df = reshape_ancom(ancom.16S.jb.10, sel_var = "Chemical")


ancom.16S.jb.6.df$cultivar = "Justin Brouwers"
ancom.16S.jb.8.df$cultivar = "Justin Brouwers"
ancom.16S.jb.10.df$cultivar = "Justin Brouwers"

ancom.16S.jb.6.df$month = "June"
ancom.16S.jb.8.df$month = "August"
ancom.16S.jb.10.df$month = "October"



ancom.16S.vv.6.df %>% filter(qval < 0.05) # 14
unique(ancom.16S.vv.6.df %>% filter(qval < 0.05 & abs(logFC) >= 0.5) %>% pull(taxon)) #12

ancom.16S.vv.8.df %>% filter(qval < 0.05) # 0
ancom.16S.vv.10.df %>% filter(qval < 0.05) # 0

ancom.16S.jb.6.df %>% filter(qval < 0.05) # 0
ancom.16S.jb.8.df %>% filter(qval < 0.05) # 2
ancom.16S.jb.10.df %>% filter(qval < 0.05) # 46
unique(ancom.16S.jb.10.df %>% filter(qval < 0.05 & abs(logFC) >= 0.5) %>% pull(taxon)) #12

sig_16S.vv.6  = ancom.16S.vv.6.df  %>% filter(qval < 0.05) # 
sig_16S.vv.8  = ancom.16S.vv.8.df  %>% filter(qval < 0.05) #
sig_16S.vv.10 = ancom.16S.vv.10.df %>% filter(qval < 0.05) # 
sig_16S.jb.6  = ancom.16S.jb.6.df  %>% filter(qval < 0.05) # 
sig_16S.jb.8  = ancom.16S.jb.8.df  %>% filter(qval < 0.05) # 
sig_16S.jb.10 = ancom.16S.jb.10.df %>% filter(qval < 0.05) #


sig_16S.vv.6.TF = filter(sig_16S.vv.6, target_var == "TF") %>% pull(taxon)
unique(sig_16S.vv.6.TF) #11

sig_16S.vv.6.VG = filter(sig_16S.vv.6, target_var == "VG") %>% pull(taxon)
unique(sig_16S.vv.6.VG) #0

sig_16S.vv.6.WP = filter(sig_16S.vv.6, target_var == "WP") %>% pull(taxon)
unique(sig_16S.vv.6.WP) #3

# svg("./ANCOMBC/venn_16S_vv_6.svg",width = 5, height = 5, pointsize = 20)
(venn16S.vv.6 = ggvenn(list(TransFilm    = sig_16S.vv.6.TF,
                           `Vapor Gard` = sig_16S.vv.6.VG,
                           `Wilt-Pruf`   = sig_16S.vv.6.WP), show_percentage = F, text_size = 6, set_name_size = 4) + 
    theme(plot.title = element_text(hjust = 0.5)) + labs(title="Vardar Valley (June)")) 
# dev.off()

sig_16S.jb.10.TF = filter(sig_16S.jb.10, target_var == "TF") %>% pull(taxon)
sig_16S.jb.10.VG = filter(sig_16S.jb.10, target_var == "VG") %>% pull(taxon)
sig_16S.jb.10.WP = filter(sig_16S.jb.10, target_var == "WP") %>% pull(taxon)

(venn16S.jb.10 = ggvenn(list(TransFilm    = sig_16S.jb.10.TF,
                            `Vapor Gard` = sig_16S.jb.10.VG,
                            `Wilt-Pruf`   = sig_16S.jb.10.WP), show_percentage = F, text_size = 6,
                        set_name_size = 4) +
    theme(plot.title = element_text(hjust = 0.5)) + labs(title="Justin Brouwers (October)")) 

tiff("./ANCOMBC/venn_16S.tif", width = 4000, height = 3600, res=400)
venn16S.vv.6 + venn16S.jb.10 + plot_annotation(tag_levels = "a") & theme(plot.tag = element_text(face="bold"))
dev.off()

df.16S.acom = rbind(ancom.16S.vv.6.df,
                    ancom.16S.vv.8.df,
                    ancom.16S.vv.10.df,
                    ancom.16S.jb.6.df,
                    ancom.16S.jb.8.df,
                    ancom.16S.jb.10.df)

str(df.16S.acom)

# venn diagram


topBac = rownames(readRDS("./Comp/topTaxa_16S_heatmap.rds"))
topBac


df.16S.acom$cultivar = factor(df.16S.acom$cultivar, levels = c("Vardar Valley", "Justin Brouwers"))
df.16S.acom$month = factor(df.16S.acom$month, levels = c("June", "August", "October"))

tiff("./ANCOMBC/plot_sig_16S_LogFC0_5_withColorLabel.tif", width = 4000, height = 3000, res=400)
(sig.16S = df.16S.acom %>%
  filter(qval < 0.05 & abs(logFC) >= 0.5 ) %>%
  mutate(textColor = paste("<span style = 'color: ",
                             ifelse(taxon %in% topBac, "red", "black"),
                             ";'>",
                             taxon,
                             "</span>", sep=""), 
           colorLabel = fct_reorder(textColor, as.character(taxon))) %>%
  ggplot(aes(x = logFC, y = reorder(colorLabel, -logFC), fill = target_var)) +
  geom_col(position = position_dodge(0.9), alpha=0.7) +
  geom_errorbar(aes(xmin = logFC - se, xmax = logFC + se), position = position_dodge(0.9), width=0.2) +
  geom_vline(aes(xintercept = 0), color="grey50", linetype="dashed") +
    facet_nested(  ~ cultivar + month) +
  scale_fill_brewer(palette = "Dark2", labels=c("TransFilm", "Vapor Gard", "Wilt-Pruf")) +
  theme_classic() +
  theme(panel.grid.major.y = element_line(color="grey80", linetype = "dotted"),
        axis.text.x = element_text(size=10, color="black"),
        axis.text.y = element_markdown(size=10),
        legend.position = "bottom") +
  labs(x="Log fold change", y="Bacterial genus", fill = "Anti-desiccants"))
dev.off()

#saveRDS(sig.16S, "./ANCOMBC/B16S/sig_16S_logFC.rds")

(sig_count_16S = df.16S.acom %>%
  filter(qval < 0.05) %>%
  group_by(target_var, cultivar, month) %>%
  summarise(Enriched = sum(logFC > 0), Suppressed = sum(logFC < 0)) %>%
  tidyr::gather(key="Enrichment", value = "Count", Enriched:Suppressed) %>%
  filter(Count > 0) %>%
  mutate(AD = ifelse(target_var == "TF", "TransFilm", ifelse(target_var == "VG", "Vapor Gard", "Wilt-Pruf"))) %>%
  ggplot(aes(y = Count, x = AD, fill = Enrichment)) +
  geom_col(position = position_stack(), width = 0.5, alpha=0.5) +
  geom_text(aes(label = Count), position = position_stack(vjust = 0.5)) +
  facet_nested(~ cultivar + month) +
  theme_classic() +
  scale_y_continuous(expand = c(0,0)) +
  labs(x = "Anti-desiccant", y="Number of impacted bacterial genera") +
  scale_fill_manual(values=c("cornflowerblue", "brown1")) +
  theme(axis.title.x = element_text(color="black", size=12, margin = margin(t=10)),
        axis.title.y = element_text(color="black", size=12, margin = margin(r=10)),
        axis.text.x = element_text(color="black", size=10, angle = 45, vjust = 1, hjust = 1),
        axis.text.y = element_text(color="black", size=10),
        legend.position = "right"))



# venn diagram and count plot
tiff("./ANCOMBC/plot_sig16S_all.tif", width = 6000, height = 3600, res=400)
((sig_count_16S) / (venn16S.vv.6 + venn16S.jb.10) | (sig.16S))  + plot_annotation(tag_levels = "a") & theme(plot.tag = element_text(face = "bold"))
dev.off()

######### iTS ###############################################
rm(list = ls())
source("./tools/visual_tools.R")


ancom.ITS.vv.6  = readRDS("./ANCOMBC/ITS/ITS_VV_JUNE.rds")
ancom.ITS.vv.8  = readRDS("./ANCOMBC/ITS/ITS_VV_AUGUST.rds")
ancom.ITS.vv.10 = readRDS("./ANCOMBC/ITS/ITS_VV_OCT.rds")

ancom.ITS.jb.6  = readRDS("./ANCOMBC/ITS/ITS_JB_JUNE.rds")
ancom.ITS.jb.8  = readRDS("./ANCOMBC/ITS/ITS_JB_AUGUST.rds")
ancom.ITS.jb.10 = readRDS("./ANCOMBC/ITS/ITS_JB_OCT.rds")


ancom.ITS.vv.6.df =  reshape_ancom(ancom.ITS.vv.6, sel_var = "Chemical")
ancom.ITS.vv.8.df =  reshape_ancom(ancom.ITS.vv.8, sel_var = "Chemical")
ancom.ITS.vv.10.df = reshape_ancom(ancom.ITS.vv.10, sel_var = "Chemical")

ancom.ITS.vv.6.df$cultivar = "Vardar Valley"
ancom.ITS.vv.8.df$cultivar = "Vardar Valley"
ancom.ITS.vv.10.df$cultivar = "Vardar Valley"

ancom.ITS.vv.6.df$month = "June"
ancom.ITS.vv.8.df$month = "August"
ancom.ITS.vv.10.df$month = "October"


########### justin brouwers ###############################################
ancom.ITS.jb.6.df =  reshape_ancom(ancom.ITS.jb.6, sel_var = "Chemical")
ancom.ITS.jb.8.df =  reshape_ancom(ancom.ITS.jb.8, sel_var = "Chemical")
ancom.ITS.jb.10.df = reshape_ancom(ancom.ITS.jb.10, sel_var = "Chemical")

ancom.ITS.jb.6.df$cultivar = "Justin Brouwers"
ancom.ITS.jb.8.df$cultivar = "Justin Brouwers"
ancom.ITS.jb.10.df$cultivar = "Justin Brouwers"

ancom.ITS.jb.6.df$month = "June"
ancom.ITS.jb.8.df$month = "August"
ancom.ITS.jb.10.df$month = "October"

ancom.ITS.vv.6.df %>% filter(qval < 0.05) # 27
ancom.ITS.vv.8.df %>% filter(qval < 0.05) # 5
ancom.ITS.vv.10.df %>% filter(qval < 0.05) # 27

ancom.ITS.jb.6.df %>% filter(qval < 0.05) # 17
ancom.ITS.jb.8.df %>% filter(qval < 0.05) # 12
ancom.ITS.jb.10.df %>% filter(qval < 0.05) # 8

df.ITS.acom = rbind(ancom.ITS.vv.6.df,
                    ancom.ITS.vv.8.df,
                    ancom.ITS.vv.10.df,
                    ancom.ITS.jb.6.df,
                    ancom.ITS.jb.8.df,
                    ancom.ITS.jb.10.df)

str(df.ITS.acom)

df.ITS.acom$cultivar = factor(df.ITS.acom$cultivar, levels = c("Vardar Valley", "Justin Brouwers"))
df.ITS.acom$month = factor(df.ITS.acom$month, levels = c("June", "August", "October"))

# color the top based on the heatmap

topFUNGI = rownames(readRDS("./Comp/topTaxa_ITS_heatmap.rds"))
topFUNGI

# colorTop_fungi = df.ITS.acom %>%
#   filter(qval < 0.05 & abs(logFC) >= 0.5) %>%
#    %>% pull(colors)
# 
# colorTop_fungi
str(df.ITS.acom)

tiff("./ANCOMBC/plot_sig_ITS_logFC0_5_colorLabel.tiff", width = 4000, height = 3500, res=400)
(sig.ITS = df.ITS.acom %>%
    filter(qval < 0.05 & abs(logFC) >= 0.5) %>%
    mutate(textColor = paste("<span style = 'color: ",
                             ifelse(taxon %in% topFUNGI, "red", "black"),
                             ";'>",
                             taxon,
                             "</span>", sep=""), 
           colorLabel = fct_reorder(textColor, as.character(taxon))) %>%
    ggplot(aes(x = logFC, y = reorder(colorLabel, -logFC), fill = target_var)) +
    geom_col(position = position_dodge(0.9), alpha=0.7) +
    geom_errorbar(aes(xmin = logFC - se, xmax = logFC + se), position = position_dodge(0.9), width=0.2) +
    geom_vline(aes(xintercept = 0), color="grey50", linetype="dashed") +
    facet_nested(  ~ cultivar + month) +
    scale_fill_brewer(palette = "Dark2", labels=c("TransFilm", "Vapor Gard", "Wilt-Pruf")) +
    theme_classic() +
    theme(panel.grid.major.y = element_line(color="grey70", linetype = "dotted"),
          axis.text.x = element_text(size=10, color="black"),
          axis.text.y = element_markdown(size=10),
          legend.position = "bottom",
          # strip.background = element_blank(),
          #ggh4x.facet.nestline = element_line(color="black")
    ) +
    labs(x="Log fold change", y="Fungal genus", fill = "Anti-desiccants"))
dev.off()

############################### venn diagram #######################################
sig_ITS.vv.6 = ancom.ITS.vv.6.df %>% filter(qval < 0.05) # 27
sig_ITS.vv.8 = ancom.ITS.vv.8.df %>% filter(qval < 0.05) # 5
sig_ITS.vv.10 = ancom.ITS.vv.10.df %>% filter(qval < 0.05) # 27

sig_ITS.jb.6 = ancom.ITS.jb.6.df %>% filter(qval < 0.05) # 17
sig_ITS.jb.8 = ancom.ITS.jb.8.df %>% filter(qval < 0.05) # 12
sig_ITS.jb.10 = ancom.ITS.jb.10.df %>% filter(qval < 0.05) # 8


################################################################################
sig_ITS.vv.6.TF = filter(sig_ITS.vv.6, target_var == "TF") %>% pull(taxon)
sig_ITS.vv.6.VG = filter(sig_ITS.vv.6, target_var == "VG") %>% pull(taxon)
sig_ITS.vv.6.WP = filter(sig_ITS.vv.6, target_var == "WP") %>% pull(taxon)

svg("./ANCOMBC/venn_ITS_vv_6.svg",width = 5, height = 5, pointsize = 20)
(venn.vv.6 = ggvenn(list(TransFilm = sig_ITS.vv.6.TF,
                         `Vapor Gard` = sig_ITS.vv.6.VG,
                         `Wilt-Pruf` = sig_ITS.vv.6.WP), show_percentage = F, text_size = 10))
dev.off()

sig_ITS.vv.8.TF = filter(sig_ITS.vv.8, target_var == "TF") %>% pull(taxon)
sig_ITS.vv.8.VG = filter(sig_ITS.vv.8, target_var == "VG") %>% pull(taxon)
sig_ITS.vv.8.WP = filter(sig_ITS.vv.8, target_var == "WP") %>% pull(taxon)

svg("./ANCOMBC/venn_ITS_vv_8.svg",width = 5, height = 5, pointsize = 20)
(venn.vv.8 = ggvenn(list(TransFilm =     sig_ITS.vv.8.TF,
                         `Vapor Gard` = sig_ITS.vv.8.VG,
                         `Wilt-Pruf` =   sig_ITS.vv.8.WP), show_percentage = F, text_size = 10))
dev.off()

sig_ITS.vv.10.TF = filter(sig_ITS.vv.10, target_var == "TF") %>% pull(taxon)
sig_ITS.vv.10.VG = filter(sig_ITS.vv.10, target_var == "VG") %>% pull(taxon)
sig_ITS.vv.10.WP = filter(sig_ITS.vv.10, target_var == "WP") %>% pull(taxon)


svg("./ANCOMBC/venn_ITS_vv_10.svg",width = 5, height = 5, pointsize = 20)
(venn.vv.10 = ggvenn(list(TransFilm =     sig_ITS.vv.10.TF,
                         `Vapor Gard` =  sig_ITS.vv.10.VG,
                         `Wilt-Pruf` =    sig_ITS.vv.10.WP), show_percentage = F, text_size = 10))
dev.off()
############################# Justin Brouwers #######################################

sig_ITS.jb.6.TF = filter(sig_ITS.jb.6, target_var == "TF") %>% pull(taxon)
sig_ITS.jb.6.VG = filter(sig_ITS.jb.6, target_var == "VG") %>% pull(taxon)
sig_ITS.jb.6.WP = filter(sig_ITS.jb.6, target_var == "WP") %>% pull(taxon)


svg("./ANCOMBC/venn_ITS_jb_6.svg",width = 5, height = 5, pointsize = 20)
(venn.jb.6 = ggvenn(list(TransFilm =     sig_ITS.jb.6.TF,
                         `Vapor Gard` = sig_ITS.jb.6.VG,
                         `Wilt-Pruf` =   sig_ITS.jb.6.WP), show_percentage = F, fill_color = c("#fc6c85","#ceff1d","#5d76cb"), text_size = 10))
dev.off()

sig_ITS.jb.8.TF = filter(sig_ITS.jb.8, target_var == "TF") %>% pull(taxon)
sig_ITS.jb.8.VG = filter(sig_ITS.jb.8, target_var == "VG") %>% pull(taxon)
sig_ITS.jb.8.WP = filter(sig_ITS.jb.8, target_var == "WP") %>% pull(taxon)

svg("./ANCOMBC/venn_ITS_jb_8.svg",width = 5, height = 5, pointsize = 20)
(venn.jb.8 = ggvenn(list(TransFilm =     sig_ITS.jb.8.TF,
                         `Vapor Gard` = sig_ITS.jb.8.VG,
                         `Wilt-Pruf` =   sig_ITS.jb.8.WP), show_percentage = F, fill_color = c("#fc6c85","#ceff1d","#5d76cb"), text_size = 10))
dev.off()

sig_ITS.jb.10.TF = filter(sig_ITS.jb.10, target_var == "TF") %>% pull(taxon)
sig_ITS.jb.10.VG = filter(sig_ITS.jb.10, target_var == "VG") %>% pull(taxon)
sig_ITS.jb.10.WP = filter(sig_ITS.jb.10, target_var == "WP") %>% pull(taxon)

svg("./ANCOMBC/venn_ITS_jb_10.svg",width = 5, height = 5, pointsize = 20)
(venn.jb.10 = ggvenn(list(TransFilm =      sig_ITS.jb.10.TF,
                          `Vapor Gard` =  sig_ITS.jb.10.VG,
                          `Wilt-Pruf` =    sig_ITS.jb.10.WP), show_percentage = F, fill_color = c("#fc6c85","#ceff1d","#5d76cb"), text_size = 10))
dev.off()


###########################
# how many were impacted
##########################

df.ITS.acom %>%
  filter(qval < 0.05 & abs(logFC)>=0.5) %>%
  group_by(target_var, cultivar) %>%
  summarise(n())
  
# target_var cultivar        `n()`
# <chr>      <fct>           <int>
#   1 TF         Vardar Valley      11
# 2 TF         Justin Brouwers     4
# 3 VG         Vardar Valley      20
# 4 VG         Justin Brouwers     9
# 5 WP         Vardar Valley      24
# 6 WP         Justin Brouwers    16

svg("./ANCOMBC/sig_ITS_counts.svg",width = 10, height = 5, pointsize = 10)
df.ITS.acom %>%
  filter(qval < 0.05) %>%
  group_by(target_var, cultivar, month) %>%
  summarise(Enriched = sum(logFC > 0), Suppressed = sum(logFC < 0)) %>%
  tidyr::gather(key="Enrichment", value = "Count", Enriched:Suppressed) %>%
  filter(Count > 0) %>%
  mutate(AD = ifelse(target_var == "TF", "TransFilm", ifelse(target_var == "VG", "Vapor Gard", "Wilt-Pruf"))) %>%
  ggplot(aes(y = Count, x = AD, fill = Enrichment)) +
  geom_col(position = position_stack(), width = 0.5, alpha=0.5) +
  geom_text(aes(label = Count), position = position_stack(vjust = 0.5)) +
  facet_nested(~ cultivar + month) +
  theme_classic() +
  scale_y_continuous(expand = c(0,0)) +
  labs(x = "Anti-desiccant", y="Number of impacted fungal genera") +
  scale_fill_manual(values=c("cornflowerblue", "brown1")) +
  theme(axis.title.x = element_text(color="black", size=12, margin = margin(t=10)),
        axis.title.y = element_text(color="black", size=12, margin = margin(r=10)),
        axis.text.x = element_text(color="black", size=10, angle = 45, vjust = 1, hjust = 1),
        axis.text.y = element_text(color="black", size=10))
dev.off()  

# tiff("./ANCOMBC/plot_sig_ITS.tiff", width = 4000, height = 3000, res=400)
# (sig.ITS = df.ITS.acom %>%
#     filter(qval < 0.05 & abs(logFC) >= 0.5) %>%
#     ggplot(aes(x = logFC, y = reorder(taxon, -logFC), fill = target_var)) +
#     geom_col(position = position_dodge(0.9), alpha=0.7) +
#     geom_errorbar(aes(xmin = logFC - se, xmax = logFC + se), position = position_dodge(0.9), width=0.2) +
#     geom_vline(aes(xintercept = 0), color="grey70", linetype="dashed") +
#     facet_nested(  ~ cultivar + month) +
#     scale_fill_brewer(palette = "Dark2", labels=c("TransFilm", "Vapor Guard", "Wilt-Pruf")) +
#     theme_classic() +
#     theme(panel.grid.major.y = element_line(color="grey70", linetype = "dotted"),
#           axis.text.x = element_text(size=10, color="black"),
#           axis.text.y = element_text(size=10, color="black"),
#           legend.position = "bottom",
#           # strip.background = element_blank(),
#           #ggh4x.facet.nestline = element_line(color="black")
#           ) +
#     labs(x="Log fold change", y="Fungal genus", fill = "Anti-desiccants"))
# dev.off()

# saveRDS(sig.ITS, "./ANCOMBC/ITS/sig_ITS_logFC.rds")


####################################################
# plot
#####################################################
# 
# sig.16S.plot = readRDS("./ANCOMBC/B16S/sig_16S_logFC.rds")
# sig.ITS.plot = readRDS("./ANCOMBC/ITS/sig_ITS_logFC.rds")
# 
# library(patchwork)
# 
# tiff("./ANCOMBC/plot_sig_micro.tif", width = 4000, height = 5000, res=400)
# sig.16S.plot / sig.ITS.plot + plot_annotation(tag_levels = "a")
# dev.off()
