Wdir = "Your directory"
setwd(Wdir)

library(dplyr)
library(tidyverse)
library(phyloseq)
library(ampvis2)
library(patchwork)
library(ComplexHeatmap)



#==================================================================================
# summary
#==================================================================================
rm(list = ls())

a16S = readRDS("./RDSdata/NA_free_THREEmonths_bsthit_pruned_16S_org.rds")
ITS = readRDS("./RDSdata/NA_free_THREEmonths_bsthit_pruned_ITS_org.rds")

meta16 = sample_data(a16S)
metaITS = sample_data(ITS)

meta16$Cultivar2 = ifelse(meta16$Cultivar == "VV", "Vardar Valley", "Justin Brouwers")
meta16$Cultivar2 = factor(meta16$Cultivar2, levels = c("Vardar Valley", "Justin Brouwers"))
meta16$Cultivar2
sample_data(a16S) = meta16

metaITS$Cultivar2 = ifelse(metaITS$Cultivar == "VV", "Vardar Valley", "Justin Brouwers")
metaITS$Cultivar2 = factor(meta16$Cultivar2, levels = c("Vardar Valley", "Justin Brouwers"))
metaITS$Cultivar2
sample_data(ITS) = metaITS

# phylum level

phy.16S = a16S %>%
  speedyseq::tax_glom("Phylum") %>%
  speedyseq::psmelt()

gen.16S = a16S %>%
  speedyseq::tax_glom("Genus") %>%
  speedyseq::psmelt()

length(unique(phy.16S$Phylum))# 24

length(unique(gen.16S$Phylum)) # 24

length(unique(gen.16S$Class)) # 49
length(unique(gen.16S$Order)) # 105
length(unique(gen.16S$Family)) # 182
length(unique(gen.16S$Genus)) # 327

phy.16S$Phylum

(Phy.16S.pie = phy.16S %>%
  dplyr::group_by(Phylum, Cultivar) %>%
  dplyr::summarise(sum_abund = sum(Abundance)) %>%
  dplyr::group_by(Cultivar) %>%
  dplyr::mutate(total_cul = sum(sum_abund)) %>%
  dplyr::mutate(rel_abund = sum_abund / total_cul, new_label = ifelse(rel_abund< 0.01, "< 1%", Phylum)) %>%
  dplyr::group_by(Cultivar, new_label) %>%
  dplyr::summarise(rel_abund2 = sum(rel_abund), display = paste0(round(100*rel_abund2, 1), "%")) %>%
  dplyr::arrange(desc(rel_abund2)))



Phy.16S.pie.VV = Phy.16S.pie %>% filter(Cultivar == "VV") %>% arrange(desc(rel_abund2))
Phy.16S.pie.JB = Phy.16S.pie %>% filter(Cultivar == "JB") %>% arrange(desc(rel_abund2))

phy.16S %>%
  dplyr::group_by(Phylum, Cultivar) %>%
  dplyr::summarise(sum_abund = sum(Abundance)) %>%
  dplyr::group_by(Cultivar) %>%
  dplyr::mutate(total_cul = sum(sum_abund)) %>%
  dplyr::mutate(rel_abund = sum_abund / total_cul) %>%
  filter(Cultivar == "JB")

Combined16S_phy = unique(as.character(c(Phy.16S.pie.VV$new_label, Phy.16S.pie.JB$new_label)))
Combined16S_phy


pcol = ggsci::pal_npg()(n=length(Combined16S_phy))
pcol

names(pcol) = Combined16S_phy

pcol


Phy.16S.pie.VV$new_label = factor(Phy.16S.pie.VV$new_label, levels = Phy.16S.pie.VV$new_label)
Phy.16S.pie.JB$new_label = factor(Phy.16S.pie.JB$new_label, levels = Phy.16S.pie.VV$new_label)

(pie.16S.plot.VV = ggpubr::ggpie(data = Phy.16S.pie.VV, x = "rel_abund2", label = "display", fill = "new_label", color="white", palette = pcol) + labs(fill = "Bacterial phylum", title = "Vardar Valley") + theme(
  axis.text.x =  element_text(size = 14),
  #axis.text.y = element_text(size=12),
  legend.text = element_text(size=14),
  legend.title = element_text(size=14),
  plot.title = element_text(hjust = 0.5)))

saveRDS(pie.16S.plot.VV, "./Comp/pie_16S_VV.rds")


(pie.16S.plot.JB = ggpubr::ggpie(data = Phy.16S.pie.JB, x = "rel_abund2", label = "display", fill = "new_label", color="white") + labs(fill = "Bacterial phylum", title = "Justin Brouwers") + scale_fill_manual(values = pcol, drop=F) + theme(
  axis.text.x =  element_text(size = 14),
  #axis.text.y = element_text(size=12),
  legend.text = element_text(size=14),
  legend.title = element_text(size=14),
  plot.title = element_text(hjust = 0.5)))

saveRDS(pie.16S.plot.JB, "./Comp/pie_16S_JB.rds")

ggpubr::ggarrange(pie.16S.plot.VV,pie.16S.plot.JB, common.legend = T, legend = "bottom")


#=============================================================================
# ITS
#==============================================================================
# phylum level

phy.ITS = ITS %>%
  speedyseq::tax_glom("Phylum") %>%
  speedyseq::psmelt()

gen.ITS = ITS %>%
  speedyseq::tax_glom("Genus") %>%
  speedyseq::psmelt()

length(unique(phy.ITS$Phylum))# 2

length(unique(gen.ITS$Phylum)) # 2

length(unique(gen.ITS$Class)) # 27
length(unique(gen.ITS$Order)) # 89
length(unique(gen.ITS$Family)) # 198
length(unique(gen.ITS$Genus)) # 328

unique(phy.ITS$Phylum)

(Phy.ITS.pie = phy.ITS %>%
    dplyr::group_by(Phylum, Cultivar) %>%
    dplyr::summarise(sum_abund = sum(Abundance)) %>%
    dplyr::group_by(Cultivar) %>%
    dplyr::mutate(total_cul = sum(sum_abund)) %>%
    dplyr::mutate(rel_abund = sum_abund / total_cul, new_label = ifelse(rel_abund< 0.01, "< 1%", Phylum)) %>%
    dplyr::group_by(Cultivar, new_label) %>%
    dplyr::summarise(rel_abund2 = sum(rel_abund), display = paste0(round(100*rel_abund2, 1), "%")) %>%
    dplyr::arrange(desc(rel_abund2)))



Phy.ITS.pie.VV = Phy.ITS.pie %>% filter(Cultivar == "VV") %>% arrange(desc(rel_abund2))
Phy.ITS.pie.JB = Phy.ITS.pie %>% filter(Cultivar == "JB") %>% arrange(desc(rel_abund2))


CombinedITS_phy = unique(as.character(c(Phy.ITS.pie.VV$new_label, Phy.ITS.pie.JB$new_label)))
CombinedITS_phy


pcol.its = ggsci::pal_gsea()(n=length(CombinedITS_phy))
pcol.its

names(pcol.its) = CombinedITS_phy

pcol.its


Phy.ITS.pie.VV$new_label = factor(Phy.ITS.pie.VV$new_label, levels = Phy.ITS.pie.VV$new_label)
Phy.ITS.pie.JB$new_label = factor(Phy.ITS.pie.JB$new_label, levels = Phy.ITS.pie.JB$new_label)

(pie.ITS.plot.VV = ggpubr::ggpie(data = Phy.ITS.pie.VV, x = "rel_abund2", label = "display", fill = "new_label", color="white") + labs(fill = "Fungal phylum", title = "Vardar Valley") + theme(
  axis.text.x =  element_text(size = 14),
  #axis.text.y = element_text(size=12),
  legend.text = element_text(size=14),
  legend.title = element_text(size=14),
  plot.title = element_text(hjust = 0.5)))

saveRDS(pie.ITS.plot.VV, "./Comp/pie_ITS_VV.rds")

(pie.ITS.plot.JB = ggpubr::ggpie(data = Phy.ITS.pie.JB, x = "rel_abund2", label = "display", fill = "new_label", color="white") + labs(fill = "Fungal phylum", title = "Justin Brouwers") + theme(
  axis.text.x =  element_text(size = 14),
  #axis.text.y = element_text(size=12),
  legend.text = element_text(size=14),
  legend.title = element_text(size=14),
  plot.title = element_text(hjust = 0.5)))

saveRDS(pie.ITS.plot.JB, "./Comp/pie_ITS_JB.rds")

ggpubr::ggarrange(pie.ITS.plot.VV,pie.ITS.plot.JB, common.legend = T, legend = "bottom")




# combine with heatmap
# svg("./Comp/phylum.svg", width = 10, height = 10)
# (pie.16S.plot.VV + 
#     pie.16S.plot.JB + plot_layout(guides = "collect")) /
#   (pie.ITS.plot.VV + 
#      pie.ITS.plot.JB +  plot_layout(guides = "collect")) + plot_layout(ncol = 1) + plot_annotation(tag_levels = "a") & theme(legend.position = "right",plot.tag = element_text(face="bold"))
# dev.off()


#####################

# between treated and nontreated
# core taxa
# share between June, August, October 

####################


########### shared

gen.16S.jun.vv = gen.16S %>% filter(Month == "June" & Cultivar == "VV" & Abundance != 0)
gen.16S.aug.vv = gen.16S %>% filter(Month == "August" & Cultivar == "VV"& Abundance != 0)
gen.16S.oct.vv = gen.16S %>% filter(Month == "October" & Cultivar == "VV"& Abundance != 0)

gen.16S.jun.jb = gen.16S %>% filter(Month == "June" & Cultivar == "JB" & Abundance != 0)
gen.16S.aug.jb = gen.16S %>% filter(Month == "August" & Cultivar == "JB"& Abundance != 0)
gen.16S.oct.jb = gen.16S %>% filter(Month == "October" & Cultivar == "JB"& Abundance != 0)

gen.16S.vv = gen.16S %>% filter(Cultivar == "VV" & Abundance != 0)
gen.16S.jb = gen.16S %>% filter(Cultivar == "JB" & Abundance != 0)


(gen16S.cultivar.venn = ggvenn::ggvenn(list(`Vardar Valley` = gen.16S.vv$Genus,
                    `Justin Brouwers` = gen.16S.jb$Genus),
               show_percentage = F, text_size = 5))

(gen16S.vv.season.venn = ggvenn::ggvenn(list(June = gen.16S.jun.vv$Genus,
                    August = gen.16S.aug.vv$Genus,
                    October = gen.16S.oct.vv$Genus), show_percentage = F, text_size = 6) +
  labs(title="Vardar Valley")+ theme(plot.title = element_text(hjust = 0.5, size = 17)))

# 3 months shared 250
# June and August shared 10
# June and October shared 23
# August and October shared 19
# Unique to June 13, August 1, October 9


(gen16S.jb.season.venn = ggvenn::ggvenn(list(June = gen.16S.jun.jb$Genus,
                    August = gen.16S.aug.jb$Genus,
                    October = gen.16S.oct.jb$Genus), show_percentage = F, text_size = 6) +
  labs(title="Justin Brouwers")+ theme(plot.title = element_text(hjust = 0.5, size = 17)))


# 3 months shared 259
# June and August shared 19
# June and October shared 16
# August and October shared 21
# Unique to June 5, August 3, October 3

gen.ITS.jun.vv = gen.ITS %>% filter(Month == "June" & Cultivar == "VV" & Abundance != 0)
gen.ITS.aug.vv = gen.ITS %>% filter(Month == "August" & Cultivar == "VV"& Abundance != 0)
gen.ITS.oct.vv = gen.ITS %>% filter(Month == "October" & Cultivar == "VV"& Abundance != 0)
gen.ITS.jun.jb = gen.ITS %>% filter(Month == "June" & Cultivar == "JB" & Abundance != 0)
gen.ITS.aug.jb = gen.ITS %>% filter(Month == "August" & Cultivar == "JB"& Abundance != 0)
gen.ITS.oct.jb = gen.ITS %>% filter(Month == "October" & Cultivar == "JB"& Abundance != 0)

gen.ITS.vv = gen.ITS %>% filter(Cultivar == "VV" & Abundance != 0)
gen.ITS.jb = gen.ITS %>% filter(Cultivar == "JB" & Abundance != 0)


(genITS.cultivar.venn = ggvenn::ggvenn(list(`Vardar Valley`   = gen.ITS.vv$Genus,
                                            `Justin Brouwers` = gen.ITS.jb$Genus),
                                       show_percentage = F, text_size = 6, fill_color = c("salmon", "cyan4")))


(genITs.vv.season.venn = ggvenn::ggvenn(list(June =    gen.ITS.jun.vv$Genus,
                    August =  gen.ITS.aug.vv$Genus,
                    October = gen.ITS.oct.vv$Genus), show_percentage = F, text_size = 6, fill_color = c("salmon", "cyan4", "goldenrod1"))+ 
                      theme(plot.title = element_text(hjust = 0.5, size=17))+ 
                      labs(title="Vardar Valley"))

# 3 months shared 243
# June and August shared 4
# June and October shared 20
# August and October shared 36
# Unique to June 7, August 3, October 13


(genITS.jb.season.venn = ggvenn::ggvenn(list(June =    gen.ITS.jun.jb$Genus,
                    August =  gen.ITS.aug.jb$Genus,
                    October = gen.ITS.oct.jb$Genus), show_percentage = F, text_size = 6, fill_color = c("salmon", "cyan4", "goldenrod1"))+
    labs(title="Justin Brouwers") + 
    theme(plot.title = element_text(hjust = 0.5, size=17)))

# 3 months shared 222
# June and August shared 14
# June and October shared 22
# August and October shared 26
# Unique to June 15, August 7, October 10

tiff("./Comp/venn_season.tif", width = 3600, height = 4600, res=400)
gen16S.vv.season.venn + gen16S.jb.season.venn + genITs.vv.season.venn + genITS.jb.season.venn + plot_annotation(tag_levels = "a") & theme(plot.tag = element_text(face="bold"))
dev.off()


#############################################################################
# bar plot
###############################################################################
rm(list = ls())

a16S = readRDS("./RDSdata/NA_free_THREEmonths_bsthit_pruned_16S_org.rds")
ITS = readRDS("./RDSdata/NA_free_THREEmonths_bsthit_pruned_ITS_org.rds")

meta16 = sample_data(a16S)
metaITS = sample_data(ITS)

meta16$Cultivar2 = ifelse(meta16$Cultivar == "VV", "Vardar Valley", "Justin Brouwers")
meta16$Cultivar2 = factor(meta16$Cultivar2, levels = c("Vardar Valley", "Justin Brouwers"))
meta16$Cultivar2
sample_data(a16S) = meta16

metaITS$Cultivar2 = ifelse(metaITS$Cultivar == "VV", "Vardar Valley", "Justin Brouwers")
metaITS$Cultivar2 = factor(meta16$Cultivar2, levels = c("Vardar Valley", "Justin Brouwers"))
metaITS$Cultivar2
sample_data(ITS) = metaITS

source("./tools/visual_tools.R")

# compared between treated and nontreated
gen.16S = a16S %>%
  speedyseq::tax_glom("Genus") %>%
  speedyseq::psmelt()

gen.16S$treat = ifelse(gen.16S$Chemical == "NT", "Nontreated", "Treated")
gen.16S$treat = factor(gen.16S$treat, levels = c("Nontreated", "Treated"))

# 
gen.ITS = ITS %>%
  speedyseq::tax_glom("Genus") %>%
  speedyseq::psmelt()

gen.ITS$treat = ifelse(gen.ITS$Chemical == "NT", "Nontreated", "Treated")
gen.ITS$treat = factor(gen.ITS$treat, levels = c("Nontreated", "Treated"))


(cul16S.bar = plotBars(gen.16S, taxlevel =  "Genus", fontangle = 0, fontHjust = 0.5, abs_abund = F, groupby = c("Cultivar2"), choosecolorset = ggsci::pal_startrek()(7), plotx = "Cultivar2", otherColor = "grey80", plotMargin_r = 40, nshow = 10, facet_group = "", plottitle = "", xlab = "", labfill = "Genus"))


(culITS.bar = plotBars(gen.ITS, taxlevel =  "Genus", fontangle = 0, fontHjust = 0.5, abs_abund = F, groupby = c("Cultivar2"), choosecolorset = ggsci::pal_jco()(10), plotx = "Cultivar2", otherColor = "grey80", plotMargin_r = 40, nshow = 10, facet_group = "", plottitle = "", xlab = "", labfill = "Genus"))


(gen16S.treated.bar = plotBars(gen.16S, taxlevel =  "Genus", fontangle = 0, fontHjust = 0.5, abs_abund = F, groupby = c("Cultivar2", "treat"), plotx = "treat", otherColor = "grey80", plotMargin_r = 40, nshow = 10, facet_group = "~ Cultivar2", plottitle = "", xlab = "", labfill = "Genus"))

saveRDS(gen16S.treated.bar, "./Comp/gen16S_treated.rds")

# with time
(gen16S.treated.T.bar = plotBars(gen.16S, taxlevel =  "Genus", fontangle = 0, fontHjust = 0.5, abs_abund = F, groupby = c("Cultivar2", "Month", "treat"), plotx = "treat", otherColor = "grey80", plotMargin_r = 40, nshow = 10, facet_group = "~ Cultivar2 + Month", plottitle = "", xlab = "", labfill = "Genus"))

################ ITS ############################################################

(genITS.treated.bar = plotBars(gen.ITS, taxlevel =  "Genus", fontangle = 0, fontHjust = 0.5, abs_abund = F, choosecolorset = ggsci::pal_d3()(10), groupby = c("Cultivar2", "treat"), plotx = "treat", otherColor = "grey80", plotMargin_r = 40, nshow = 10, facet_group = "~ Cultivar2", plottitle = "", xlab = "", labfill = "Fungal genus"))

saveRDS(genITS.treated.bar, "./Comp/genITS_treated.rds")

(genITS.treated.T.bar = plotBars(gen.ITS, taxlevel =  "Genus", fontangle = 45, fontHjust = 1, fontVjust = 1, abs_abund = F, choosecolorset = ggsci::pal_d3()(10), groupby = c("Cultivar2", "Month", "treat"), plotx = "treat", otherColor = "grey80", plotMargin_r = 40, nshow = 10, facet_group = "~ Cultivar2 + Month", plottitle = "", xlab = "", labfill = "Fungal genus"))

saveRDS(genITS.treated.T.bar, "./Comp/genITS_treated_wTime.rds")

# (chem16S.bar = plotBars(gen.16S, taxlevel =  "Genus", fontangle = 0, fontHjust = 0.5, abs_abund = F, groupby = c("Chemical", "Month", "Cultivar2"), choosecolorset = ggsci::pal_startrek()(7), plotx = "Chemical", otherColor = "grey80", plotMargin_r = 40, nshow = 10, facet_group = "~ Cultivar2 + Month", plottitle = "", xlab = "", labfill = "Genus"))


gen16S.treated.bar + genITS.treated.bar
#####################################################################################
# heat map for details
########################################################################################


library(ampvis2)
devtools::source_gist("8d0ca4206a66be7ff6d76fc4ab8e66c6")

# a16S = readRDS("./RDSdata/NA_free_THREEmonths_bsthit_pruned_16S_org.rds")
# ITS = readRDS("./RDSdata/NA_free_THREEmonths_bsthit_pruned_ITS_org.rds")
# 
# 
# 
# meta16 = sample_data(a16S)
# metaITS = sample_data(ITS)
# 
# meta16$Cultivar2 = ifelse(meta16$Cultivar == "VV", "Vardar Valley", "Justin Brouwers")
# meta16$Cultivar2 = factor(meta16$Cultivar2, levels = c("Vardar Valley", "Justin Brouwers"))
# meta16$Cultivar2
# sample_data(a16S) = meta16
# 
# metaITS$Cultivar2 = ifelse(metaITS$Cultivar == "VV", "Vardar Valley", "Justin Brouwers")
# metaITS$Cultivar2 = factor(meta16$Cultivar2, levels = c("Vardar Valley", "Justin Brouwers"))
# metaITS$Cultivar2
# sample_data(ITS) = metaITS

a16S@sam_data$treat = ifelse(a16S@sam_data$Chemical == "NT", "Nontreated", "Treated")
a16S@sam_data$treat = factor(a16S@sam_data$treat, levels = c("Nontreated", "Treated"))


ITS@sam_data$treat = ifelse(ITS@sam_data$Chemical == "NT", "Nontreated", "Treated")
ITS@sam_data$treat = factor(ITS@sam_data$treat, levels = c("Nontreated", "Treated"))

a16S.viz = phyloseq_to_ampvis2(a16S)
ITS.viz = phyloseq_to_ampvis2(ITS)



str(meta16)
?amp_heatmap

# (hmp.16S.treat = amp_heatmap(a16S.viz,
#                             group_by = "treat",
#                             facet_by = c("Cultivar2", "Month"),
#                             tax_aggregate = "Genus",
#                             tax_add = "Phylum",
#                             normalise = T,
#                             # plot_colorscale = "sqrt",
#                             min_abundance = 0.1,
#                             tax_show = 15,
#                             plot_values = T,
#                             plot_values_size = 4) + 
#   theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
#         legend.position = "none",
#         legend.text = element_text(size=12)))



# (hmp.16S.treat = amp_heatmap(a16S.viz,
#                        group_by = "treat",
#                        facet_by = c("Cultivar2", "Month"),
#                        tax_aggregate = "Genus",
#                        tax_add = "Phylum",
#                        normalise = T,
#                        # plot_colorscale = "sqrt",
#                        min_abundance = 0.1,
#                        tax_show = 15,
#                        plot_values = F,
#                        plot_values_size = 3) + 
#     theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
#           legend.position = "none",
#           legend.text = element_text(size=12)))
# 
(hmp.16S.treat.val = amp_heatmap(a16S.viz,
                       group_by = "treat",
                       facet_by = c("Cultivar2", "Month"),
                       tax_aggregate = "Genus",
                       tax_add = "Phylum",
                       normalise = T,
                       # plot_colorscale = "sqrt",
                       min_abundance = 0.1,
                       tax_show = 15,
                       plot_values = T,
                       plot_values_size = 4) +
    theme(
          #axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          legend.position = "right",
          legend.text = element_text(size=12)))

(hmp.16.val = amp_heatmap(a16S.viz,
                          group_by = "Chemical",
                          facet_by = c("Cultivar2", "Month"),
                          tax_aggregate = "Genus",
                          tax_add = "Phylum",
                          normalise = T,
                          # plot_colorscale = "sqrt",
                          min_abundance = 0.1,
                          tax_show = 15,
                          plot_values = T,
                          plot_values_size = 4) + 
    theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
          legend.position = "none",
          legend.text = element_text(size=12)))

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
hmp.16.val.textmap = amp_heatmap(a16S.viz,
                                 group_by = "treat",
                                 facet_by = c("Cultivar2", "Month"),
                                 tax_aggregate = "Genus",
                                 # tax_add = "Phylum",
                                 normalise = T,
                                 # plot_colorscale = "sqrt",
                                 min_abundance = 0.1,
                                 tax_show = 15,
                                 plot_values = T,
                                 plot_values_size = 4, textmap = T)


# colnames(hmp.16.val.textmap) = c("CT_JB_8", "CT_JB_6", "CT_JB_10",
#                                  "CT_VV_8", "CT_VV_6", "CT_VV_10", 
#                                  "Trt_JB_8", "Trt_JB_6", "Trt_JB_10", 
#                                  "Trt_VV_8", "Trt_VV_6", "Trt_VV_10")

# hmp.16.val.textmap$taxon = rownames(hmp.16.val.textmap)
# 
# rowAnno = c(rep("Vardar Valley", 3),
#           rep("Justin Brouwers", 3))
# 
# Month = c(rep(c("June", "August", "October"),2))
# Month
# annoDF = data.frame(Cultivars = rowAnno, Month=Month)
# # 
# annoDF
# # 
# colum_ha = HeatmapAnnotation(df = annoDF, col = list(Cultivars=c("Vardar Valley"="orange", "Justin Brouwers"="blue4"), Month = c("June"="#31de5f", "August"="#14662a", "October"="#0a3315")))
# 
# 
# hmp.16.val.textmap %>%
#   mutate(JB_6 = Trt_JB_6 - CT_JB_6,
#          JB_8 = Trt_JB_8 - CT_JB_8,
#          JB_10 = Trt_JB_10 - CT_JB_10,
#          VV_6 = Trt_VV_6 - CT_VV_6,
#          VV_8 = Trt_VV_8 - CT_VV_8,
#          VV_10 = Trt_VV_10 - CT_VV_10) %>%
#   dplyr::select(VV_6, VV_8, VV_10,JB_6, JB_8, JB_10) %>%
#   ComplexHeatmap::Heatmap(cluster_rows = F, 
#                           cluster_columns = F, 
#                           rect_gp = gpar(col="white", lwd=2),
#                           top_annotation = colum_ha,
#                           column_split = factor(c("VV", "VV", "VV", "JB", "JB", "JB"),levels = c("VV", "JB")),
#                           column_names_rot = 45,
#                           column_title = NULL)

# saveRDS(hmp.16.val.textmap, "./Comp/topTaxa_16S_heatmap.rds")
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

(hmp.16S.cultivar = amp_heatmap(a16S.viz,
                          group_by = "Cultivar2",
                          facet_by = "Cultivar2",
                          tax_aggregate = "Genus",
                          tax_add = "Phylum",
                          normalise = T,
                          # plot_colorscale = "sqrt",
                          min_abundance = 0.1,
                          tax_show = 15,
                          plot_values = T,
                          plot_values_size = 4) + 
    theme(axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          legend.position = "right",
          legend.text = element_text(size=12)))

# (hmp.16.val = amp_heatmap(a16S.viz,
#                        group_by = "Chemical",
#                        facet_by = c("Cultivar2", "Month"),
#                        tax_aggregate = "Genus",
#                        tax_add = "Phylum",
#                        normalise = T,
#                        # plot_colorscale = "sqrt",
#                        min_abundance = 0.1,
#                        tax_show = 15,
#                        plot_values = T,
#                        plot_values_size = 3) + 
#     theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
#           legend.position = "none",
#           legend.text = element_text(size=12)))

?amp_heatmap
(hmp.ITS = amp_heatmap(ITS.viz,
                       group_by = "Chemical",
                       facet_by = c("Cultivar2", "Month"),
                       tax_aggregate = "Genus",
                       tax_add = "Phylum",
                       # normalise = F,
                       # plot_colorscale = "sqrt",
                       min_abundance = 0.1,
                       tax_show = 15,
                       plot_values = F,
                       plot_values_size = 3) + 
    theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
            legend.text = element_text(size=12),
          legend.title = element_text(size=10)))

(hmp.ITS.cultivar = amp_heatmap(ITS.viz,
                       group_by = "Cultivar2",
                       facet_by = c("Cultivar2"),
                       tax_aggregate = "Genus",
                       tax_add = "Phylum",
                       # normalise = F,
                       # plot_colorscale = "sqrt",
                       min_abundance = 0.1,
                       tax_show = 15,
                       plot_values = T,
                       plot_values_size = 4) + 
    theme(axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          legend.text = element_text(size=12),
          legend.title = element_text(size=10)))



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
hmp.ITS.textmap = amp_heatmap(ITS.viz,
                      group_by = "Chemical",
                      facet_by = c("Cultivar2", "Month"),
                      tax_aggregate = "Genus",
                      # tax_add = "Phylum",
                      # normalise = F,
                      # plot_colorscale = "sqrt",
                      min_abundance = 0.1,
                      tax_show = 15,
                      plot_values = F,
                      plot_values_size = 3,
                      textmap = T)

saveRDS(hmp.ITS.textmap, "./Comp/topTaxa_ITS_heatmap.rds")
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

(hmp.ITS.treat.val = amp_heatmap(ITS.viz,
                       group_by = "treat",
                       facet_by = c("Cultivar2", "Month"),
                       tax_aggregate = "Genus",
                       tax_add = "Phylum",
                       # normalise = F,
                       # plot_colorscale = "sqrt",
                       min_abundance = 0.1,
                       tax_show = 15,
                       plot_values = T,
                       plot_values_size = 4) + 
    theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
          legend.text = element_text(size=12),
          legend.title = element_text(size=10)))


########## pie chart and heatmap on cultivar effect #########################################


svg("./Comp/pie_and_bar_update.svg", height = 10, width = 15.5, pointsize = 10)
(((pie.16S.plot.VV + pie.16S.plot.JB) + plot_layout(guides = "collect") & theme(legend.position = "bottom")) / ((pie.ITS.plot.VV + pie.ITS.plot.JB) + plot_layout(guides = "collect") & theme(legend.position = "bottom")) | ((hmp.16S.cultivar / hmp.ITS.cultivar) +plot_layout(guides = "collect"))) + plot_layout(widths = c(2.5,1)) + plot_annotation(tag_levels = "a") & theme(plot.tag = element_text(size=20, face = "bold"))
dev.off()

tiff("./Comp/comp_heatmap_treat_val.tiff", width = 4200, height = 3500, res=400)
hmp.16S.treat.val + hmp.ITS.treat.val + plot_layout(ncol = 1, guides = "collect") + plot_annotation(tag_levels = "a") & theme(plot.tag = element_text(size=12, face = "bold"))
dev.off()

tiff("./Comp/comp_heatmap.tiff", width = 4200, height = 3200, res=400)
hmp.16S + hmp.ITS + plot_layout(ncol = 1, guides = "collect") + plot_annotation(tag_levels = "a") &
  theme(plot.tag = element_text(size=12, face = "bold"))
dev.off()



#tiff("./Comp/comp_heatmap.tiff", width = 4200, height = 3200, res=400)
# gen16S.treated.bar + genITS.treated.bar + hmp.16S + hmp.ITS + plot_layout(ncol = 1, guides = "collect") + plot_annotation(tag_levels = "a") &
#   theme(plot.tag = element_text(size=12, face = "bold"))
#dev.off()