Wdir = "Your directory"
setwd(Wdir)

library(microbiome)
library(phyloseq)
library(dplyr)
library(tidyr)
library(vegan)
library(ampvis2)
library(vegan)
library(microbiome)
library(RColorBrewer)
library(microbiomeutilities)
library(ggpubr)
library(rstatix)
library(metagMisc)
library(pwr)
library(lmerTest)


rm(list=ls())


a16S = readRDS("./RDSdata/bsthit_pruned_16S_org.rds")
taxa_names(a16S) = gsub(":.*", "", rownames(a16S@otu_table))
a16S@tax_table = a16S@tax_table[, -ncol(a16S@tax_table)]
colnames(a16S@tax_table) = c("Kingdom", "Phylum",   "Class",    "Order",    "Family",   "Genus",    "Species")
a16S = subset_taxa(a16S, Kingdom == "Bacteria")
a16S

tax.matrix = as.matrix(a16S@tax_table)
unique(tax.matrix[,"Kingdom"]) # Bacteria

sum(tax.matrix == "c__Proteobacteria") #5 change the labeling
tax.matrix[tax.matrix == "c__Proteobacteria"] = "p__Proteobacteria"
tax_table(a16S) = tax.matrix
meta16S = data.frame(sample_data(a16S))
head(meta16S)

meta16S$Location = ifelse(meta16S$Cultivar == "VV", "Loc1", ifelse(meta16S$Cultivar == "JB", "Loc2", "NA"))
meta16S$Chemical = factor(meta16S$Chemical, levels = c("NT", "TransFilm", "Vapor_Gard", "Wilt_Pruf"))
meta16S$Month = factor(meta16S$Month, levels = c("June", "August", "October"))
levels(meta16S$Chemical) = c("NT", "TF", "VG", "WP")
meta16S$Chemical
sample_data(a16S) = meta16S
saveRDS(a16S, "./RDSdata/THREEmonths_bsthit_pruned_16S_org.rds")

########################################
# more cleaning with ITS
########################################

ITS = readRDS("./RDSdata/bsthit_pruned_ITS_org.rds")

taxa_names(ITS) = gsub(":.*", "", rownames(ITS@otu_table))
ITS@tax_table = ITS@tax_table[, -ncol(ITS@tax_table)]
colnames(ITS@tax_table) = c("Kingdom", "Phylum",   "Class",    "Order",    "Family",   "Genus",    "Species")

tax.df.its = as.data.frame(tax_table(ITS))
unique(tax.df.its$Kingdom) # fungi

metaITS = data.frame(sample_data(ITS))
head(metaITS)

# location 1 = VV, location 2 = JB
metaITS$Location = ifelse(metaITS$Cultivar == "VV", "Loc1", ifelse(metaITS$Cultivar == "JB", "Loc2", "NA"))
metaITS$Chemical = factor(metaITS$Chemical, levels = c("NT", "TransFilm", "Vapor_Gard", "Wilt_Pruf"))
levels(metaITS$Chemical) = c("NT", "TF", "VG", "WP")
metaITS$Month
metaITS$Month = factor(metaITS$Month, levels = c("June", "August", "October"))
sample_data(ITS) = metaITS
saveRDS(ITS, "./RDSdata/THREEmonths_bsthit_pruned_ITS_org.rds")

#------------------------------------------------------------------------
#-----------------------------------------------------------------------

rm(list = ls())
a16S = readRDS("./RDSdata/THREEmonths_bsthit_pruned_16S_org.rds")
ITS = readRDS("./RDSdata/THREEmonths_bsthit_pruned_ITS_org.rds")

min(sample_sums(a16S)) #63
min(sample_sums(ITS)) # 7392
# remove NA
a16S = subset_samples(a16S, !is.na(Chemical))
ITS = subset_samples(ITS, !is.na(Chemical))

min(sample_sums(a16S)) #63
min(sample_sums(ITS)) # 7392

meta16S = data.frame(sample_data(a16S))
metaITS = data.frame(sample_data(ITS))

#saveRDS(a16S, "./RDSdata/NA_free_THREEmonths_bsthit_pruned_16S_org.rds")
#saveRDS(ITS, "./RDSdata/NA_free_THREEmonths_bsthit_pruned_ITS_org.rds")

#----------------------------------------------------------------------------
devtools::source_gist("8d0ca4206a66be7ff6d76fc4ab8e66c6")

amp16S = phyloseq_to_ampvis2(a16S)
ampITS = phyloseq_to_ampvis2(ITS)

otu16S.mtx = otu_table(a16S)
class(otu16S.mtx) = "matrix"

rarecurve(t(otu16S.mtx), label = F)


(rarecur.16S.funlab = ampvis2::amp_rarecurve(amp16S, stepsize = 100, color_by = "Chemical") +
    labs(color = "Antidesiccants", x="", y="Number of bacterial OTUs") +
    # geom_vline(aes(xintercept=)) +
    scale_x_continuous(limits = c(0, 2500)) +
    theme(text = element_text(size=12, color="black"),
          axis.text.x = element_text(size=12, color="black"),
          axis.text.y = element_text(size=12, color="black"),
          axis.title.x = element_text(size=12, color="black", margin = margin(t=10)),
          axis.title.y = element_text(size=12, color="black", margin = margin(r=10))))

(rarecur.ITS.funlab = ampvis2::amp_rarecurve(ampITS, stepsize = 100, color_by = "Chemical") +
    labs(color = "Antidesiccants", y="Number of fungal OTUs") +
    # geom_vline(aes(xintercept=10000)) +
    scale_x_continuous(limits = c(0, 180000)) +
    theme(text = element_text(size=12, color="black"),
          axis.text.x = element_text(size=12, color="black"),
          axis.text.y = element_text(size=12, color="black"),
          axis.title.x = element_text(size=12, color="black", margin = margin(t=10)),
          axis.title.y = element_text(size=12, color="black", margin = margin(r=10))))


library(patchwork)

tiff("./alpha_diversity/rarefaction.tif", width = 3000, height = 3000, res = 400)
rarecur.16S.funlab / rarecur.ITS.funlab + plot_layout(guides = "collect") + plot_annotation(tag_levels = "a") + theme(plot.tag = element_text(face='bold', size=14))
dev.off()

# good's coverage

a16S_coverage = phyloseq_coverage(a16S, correct_singletons = T)
ITS_coverage = phyloseq_coverage(ITS, correct_singletons = T)

a16S_coverage

length(unique(a16S_coverage$SampleID))
length(unique(ITS_coverage$SampleID))
mean(a16S_coverage$SampleCoverage) # 0.9000804
mean(ITS_coverage$SampleCoverage) # 0.9989687

#---------------------------------------------------------------------------------------
# get alpha diversity
rm(list = ls())

a16S = readRDS("./RDSdata/NA_free_THREEmonths_bsthit_pruned_16S_org.rds")
ITS = readRDS("./RDSdata/NA_free_THREEmonths_bsthit_pruned_ITS_org.rds")

min(sample_sums(a16S)) # 63
median(sample_sums(a16S)) # 626

min(sample_sums(ITS)) #7392
median(sample_sums(ITS)) #53812

meta16S = data.frame(sample_data(a16S))
metaITS = data.frame(sample_data(ITS))


sum(sample_sums(a16S) >= 200) # 50-96, 100-92, 200-83
length(meta16S$Original_label) #96

sum(sample_sums(ITS) >= 10000) # 30,000-78samples #10,000-93 samples


str(meta16S)

# using Scaling with ranked subsampling (SRS)
# library(SRS)
# ?SRScurve
# SRScurve(data.frame(otu_table(a16S)), metric = "richness", step=100, col = c("red", "blue"), rarefy.comparison = T)
# SRS.shiny.app(data.frame(otu_table(a16S)))

# rarefy samples
rare16S = rarefy_even_depth(a16S, sample.size = 50, rngseed = 123, replace = F)
# rare16S = rarefy_even_depth(a16S, sample.size = min(sample_sums(a16S)), rngseed = 123, replace = F)
dv16S = estimate_richness(rare16S)

rareITS = rarefy_even_depth(ITS, sample.size = 10000, rngseed = 123, replace = F)
# rareITS = rarefy_even_depth(ITS, sample.size = min(sample_sums(ITS)), rngseed = 123, replace = F)
dvITS = estimate_richness(rareITS)


head(meta16S)
head(metaITS)

af16S = merge(dv16S, meta16S, by.x = "row.names", by.y = "SampleID")
afITS = merge(dvITS, metaITS, by.x = "row.names", by.y = "row.names")


af16S.1 = af16S %>%
  gather("alpha_index", "index_measure", c("Observed", "Shannon", "InvSimpson", "Chao1", "ACE", "Simpson", "Fisher"))

head(af16S.1)

afITS.1 = afITS %>%
  gather("alpha_index", "index_measure", c("Observed", "Shannon", "InvSimpson","Chao1", "ACE", "Simpson", "Fisher"))

#saveRDS(af16S.1, "./RDSdata/af16S_jun_aug_oct_50.rds")
#saveRDS(afITS.1, "./RDSdata/afITS_jun_aug_oct_10000.rds")


#----------------------------------------------------------
# plots
#---------------------------------------------------------

rm(list=ls())


source("./tools/Basic_stats.R")
af16S.1 = readRDS("./RDSdata/af16S_jun_aug_oct_50.rds")
afITS.1 = readRDS("./RDSdata/afITS_jun_aug_oct_10000.rds")

af16S.2 = af16S.1 %>%
  mutate(locs = ifelse(Location == "Loc1", "Location 1", "Location 2")) %>%
  mutate(Chemical_lab = ifelse(Chemical == "NT", "Nontreated", ifelse(Chemical == "TF", "TransFilm", ifelse(Chemical == "VG", "Vapor Gard", ifelse(Chemical == "WP", "Wilt-Pruf", "NA")))))


afITS.2 = afITS.1 %>%
  mutate(locs = ifelse(Location == "Loc1", "Location 1", "Location 2")) %>%
  mutate(Chemical_lab = ifelse(Chemical == "NT", "Nontreated", ifelse(Chemical == "TF", "TransFilm", ifelse(Chemical == "VG", "Vapor Gard", ifelse(Chemical == "WP", "Wilt-Pruf", "NA")))))


# the observed OTU richness
data.obs.16s = af16S.2 %>%
  filter(alpha_index == "Observed")

data.obs.16s$treated = ifelse(data.obs.16s$Chemical != "NT", "Treated", "Nontreated")


data.obs.ITS = afITS.2 %>%
  filter(alpha_index == "Observed")

data.obs.ITS$treated = ifelse(data.obs.ITS$Chemical != "NT", "Treated", "Nontreated")


# with disease assessment
ggplot(data.obs.16s, aes(x = Disease_assessment, y=index_measure)) +
  geom_point() +
  geom_smooth(aes(color=Month), method = "lm", formula = "y ~ x")


# basic statistics
sum16S = summarySE(data.obs.16s, measurevar = "index_measure", groupvars = c("Month", "Cultivar", "Chemical_lab"))
sum16S

sum16S.2 = summarySE(data.obs.16s, measurevar = "index_measure", groupvars = c("Month", "Cultivar", "treated"))
sum16S.2

(plot16S = sum16S %>%
  ggplot(aes(x = Month, y=index_measure, color=Chemical_lab, group = Chemical_lab)) +
  geom_point(size=3) +
  geom_line(linetype = "dashed", size=1) +
  geom_errorbar(aes(ymax = index_measure + se, ymin = index_measure - se), width = 0.1, size=0.2) +
  facet_wrap(~Cultivar, labeller = as_labeller(c(`VV`="Vardar Valley", `JB`="Justin Brouwers")),nrow = 2, strip.position = "right") +
  labs(x = "Month", y = "Observed OTU richness", color = "Anti-desiccant", title = "Bacterial community") +
  theme_biome_utils() +
  theme(axis.text.x = element_text(size=12, color="black", angle = 0),
        axis.text.y = element_text(size=12, color="black"),
        axis.title.x = element_text(size=12, margin = margin(t=10), color="black"),
        axis.title.y = element_text(size=12, margin = margin(r=10), color="black")))

# treated and nontreated
(plot16S.2 = sum16S.2 %>%
    ggplot(aes(x = Month, y=index_measure, fill=treated, group = treated)) +
    # geom_point(size=3) +
    geom_col(position = position_dodge(width = 0.9)) +
    geom_errorbar(aes(ymax = index_measure + se, ymin = index_measure - se), width = 0.1, size=0.2, position = position_dodge(width = 0.9)) +
    facet_wrap(~Cultivar, labeller = as_labeller(c(`VV`="Vardar Valley", `JB`="Justin Brouwers")),nrow = 2, strip.position = "right") +
    labs(x = "Month", y = "Observed OTU richness", fill = "Anti-desiccant", title = "Bacterial community") +
    theme_biome_utils() +
    theme(axis.text.x = element_text(size=12, color="black", angle = 0),
          axis.text.y = element_text(size=12, color="black"),
          axis.title.x = element_text(size=12, margin = margin(t=10), color="black"),
          axis.title.y = element_text(size=12, margin = margin(r=10), color="black")))


#saveRDS(plot16S, "./RDSdata/plot16S_obs.rds")
#saveRDS(plot16S.2, "./RDSdata/plot16S_obs_treated_bar.rds")


sumITS = summarySE(data.obs.ITS, measurevar = "index_measure", groupvars = c("Month", "Cultivar", "Chemical_lab"))
sumITS

sumITS.2 = summarySE(data.obs.ITS, measurevar = "index_measure", groupvars = c("Month", "Cultivar", "treated"))
sumITS.2


(plotITS = sumITS %>%
    ggplot(aes(x = Month, y=index_measure, color=Chemical_lab, group = Chemical_lab)) +
    geom_point(size=3) +
    geom_line(linetype = "dashed", size=1) +
    geom_errorbar(aes(ymax = index_measure + se, ymin = index_measure - se), width = 0.1, size=0.2) +
    facet_wrap(~Cultivar, labeller = as_labeller(c(`VV`="Vardar Valley", `JB`="Justin Brouwers")),nrow = 2, strip.position = "right") +
    
    labs(x = "Month", y = "Observed OTU richness", color = "Anti-desiccant",title = "Fungal community") +
    theme_biome_utils() +
    theme(axis.text.x = element_text(size=12, color="black", angle = 0),
          axis.text.y = element_text(size=12, color="black"),
          axis.title.x = element_text(size=12, margin = margin(t=10), color="black"),
          axis.title.y = element_text(size=12, margin = margin(r=10), color="black")))

(plotITS.2 = sumITS.2 %>%
    ggplot(aes(x = Month, y=index_measure, fill=treated, group = treated)) +
    # geom_point(size=3) +
    # geom_line(linetype = "dashed", size=1) +
    geom_col(position = position_dodge(width = 0.9))+
    geom_errorbar(aes(ymax = index_measure + se, ymin = index_measure - se), width = 0.1, size=0.2, position = position_dodge(width = 0.9)) +
    facet_wrap(~Cultivar, labeller = as_labeller(c(`VV`="Vardar Valley", `JB`="Justin Brouwers")),nrow = 2, strip.position = "right") +
    
    labs(x = "Month", y = "Observed OTU richness", fill = "Anti-desiccant",title = "Fungal community") +
    theme_biome_utils() +
    theme(axis.text.x = element_text(size=12, color="black", angle = 0),
          axis.text.y = element_text(size=12, color="black"),
          axis.title.x = element_text(size=12, margin = margin(t=10), color="black"),
          axis.title.y = element_text(size=12, margin = margin(r=10), color="black")))

#saveRDS(plotITS, "./RDSdata/plotITS_obs.rds")
#saveRDS(plotITS.2, "./RDSdata/plotITS_obs_treated_bar.rds")


# inverSimpson
data.invS.16s = af16S.2 %>%
  filter(alpha_index == "InvSimpson")

data.invS.ITS = afITS.2 %>%
  filter(alpha_index == "InvSimpson")

data.invS.16s$treated = ifelse(data.invS.16s$Chemical != "NT", "Treated", "Nontreated")
data.invS.ITS$treated = ifelse(data.invS.ITS$Chemical != "NT", "Treated", "Nontreated")

# basic statistics
sum16S.invS = summarySE(data.invS.16s, measurevar = "index_measure", groupvars = c("Month", "Cultivar", "Chemical_lab"))
sum16S.invS

(plot16S.invS = sum16S.invS %>%
    ggplot(aes(x = Month, y=index_measure, color=Chemical_lab, group = Chemical_lab)) +
    geom_point(size=3) +
    geom_line(linetype = "dashed", size=1) +
    geom_errorbar(aes(ymax = index_measure + se, ymin = index_measure - se), width = 0.1, size=0.2) +
    facet_wrap(~Cultivar, labeller = as_labeller(c(`VV`="Vardar Valley", `JB`="Justin Brouwers")),nrow = 2, strip.position = "right") +
    labs(x = "Month", y = "Inverse Simpson's index", color = "Anti-desiccant",title = "Bacterial community") +
    theme_biome_utils() +
    theme(axis.text.x = element_text(size=12, color="black", angle = 0),
          axis.text.y = element_text(size=12, color="black"),
          axis.title.x = element_text(size=12, margin = margin(t=10), color="black"),
          axis.title.y = element_text(size=12, margin = margin(r=10), color="black")))

sum16S.invS.2 = summarySE(data.invS.16s, measurevar = "index_measure", groupvars = c("Month", "Cultivar", "treated"))
sum16S.invS.2

(plot16S.invS.2 = sum16S.invS.2 %>%
    ggplot(aes(x = Month, y=index_measure, fill = treated, group = treated)) +
    # geom_point(size=3) +
    # geom_line(linetype = "dashed", size=1) +
    geom_col(position = position_dodge(width = 0.9)) +
    geom_errorbar(aes(ymax = index_measure + se, ymin = index_measure - se), width = 0.1, size=0.2,position = position_dodge(width = 0.9)) +
    facet_wrap(~Cultivar, labeller = as_labeller(c(`VV`="Vardar Valley", `JB`="Justin Brouwers")),nrow = 2, strip.position = "right") +
    labs(x = "Month", y = "Inverse Simpson's index", fill = "Anti-desiccant",title = "Bacterial community") +
    theme_biome_utils() +
    theme(axis.text.x = element_text(size=12, color="black", angle = 0),
          axis.text.y = element_text(size=12, color="black"),
          axis.title.x = element_text(size=12, margin = margin(t=10), color="black"),
          axis.title.y = element_text(size=12, margin = margin(r=10), color="black")))

#saveRDS(plot16S.invS, "./RDSdata/plot16S_invS.rds")
#saveRDS(plot16S.invS.2, "./RDSdata/plot16S_invS_treated_bar.rds")


sumITS.invS = summarySE(data.invS.ITS, measurevar = "index_measure", groupvars = c("Month", "Cultivar", "Chemical_lab"))
sumITS.invS



(plotITS.invS = sumITS.invS %>%
    ggplot(aes(x = Month, y=index_measure, color=Chemical_lab, group = Chemical_lab)) +
    geom_point(size=3) +
    geom_line(linetype = "dashed", size=1) +
    geom_errorbar(aes(ymax = index_measure + se, ymin = index_measure - se), width = 0.1, size=0.2) +
    facet_wrap(~Cultivar, labeller = as_labeller(c(`VV`="Vardar Valley", `JB`="Justin Brouwers")),nrow = 2, strip.position = "right") +
    
    labs(x = "Month", y = "Inverse Simpson's index", color = "Anti-desiccant",title = "Fungal community") +
    theme_biome_utils() +
    theme(axis.text.x = element_text(size=12, color="black", angle = 0),
          axis.text.y = element_text(size=12, color="black"),
          axis.title.x = element_text(size=12, margin = margin(t=10), color="black"),
          axis.title.y = element_text(size=12, margin = margin(r=10), color="black")))


sumITS.invS.2 = summarySE(data.invS.ITS, measurevar = "index_measure", groupvars = c("Month", "Cultivar", "treated"))
sumITS.invS.2

(plotITS.invS.2 = sumITS.invS.2 %>%
    ggplot(aes(x = Month, y=index_measure, fill=treated, group = treated)) +
    # geom_point(size=3) +
    # geom_line(linetype = "dashed", size=1) +
    geom_col(position = position_dodge(width = 0.9)) +
    geom_errorbar(aes(ymax = index_measure + se, ymin = index_measure - se), width = 0.1, size=0.2, position = position_dodge(width = 0.9)) +
    facet_wrap(~Cultivar, labeller = as_labeller(c(`VV`="Vardar Valley", `JB`="Justin Brouwers")),nrow = 2, strip.position = "right") +
    
    labs(x = "Month", y = "Inverse Simpson's index", fill = "Anti-desiccant",title = "Fungal community") +
    theme_biome_utils() +
    theme(axis.text.x = element_text(size=12, color="black", angle = 0),
          axis.text.y = element_text(size=12, color="black"),
          axis.title.x = element_text(size=12, margin = margin(t=10), color="black"),
          axis.title.y = element_text(size=12, margin = margin(r=10), color="black")))

#saveRDS(plotITS.invS, "./RDSdata/plotITS_invS.rds")
#saveRDS(plotITS.invS.2, "./RDSdata/plotITS_invS_treated_bar.rds")

sumITS.invS.3 = summarySE(data.invS.ITS, measurevar = "index_measure", groupvars = c("Month", "Cultivar"))






library(patchwork)
tiff("./Alpha_diversity/updated_alpha.tif", width = 4000, height = 4000, res=400)
(diversity_plot = plot16S + plot16S.invS + plotITS  + plotITS.invS + plot_layout(guides = "collect") + plot_annotation(tag_levels = "a"))
dev.off()


tiff("./Alpha_diversity/updated_alpha_treated_bar.tif", width = 4000, height = 4000, res=400)
(diversity_plot = plot16S.2 + plot16S.invS.2 + plotITS.2  + plotITS.invS.2 + plot_layout(guides = "collect") + plot_annotation(tag_levels = "a"))
dev.off()


# make barplots

##### Richness

(obs.16s.bar = data.obs.16s %>%
    mutate(cultivar_lab = ifelse(Cultivar == "VV", "Vardar Valley", "Justin Brouwers")) %>%
    ggbarplot(., x = "Month", y="index_measure", facet.by = "cultivar_lab", add = c("mean_se"), error.plot = "errorbar", ylab =  "Observed OTU richness", xlab = "", title = "Bacterial community"))

stat.test.obs.16s.bar =  data.obs.16s %>%
  mutate(cultivar_lab = ifelse(Cultivar == "VV", "Vardar Valley", "Justin Brouwers")) %>%
  group_by(cultivar_lab) %>%
  t_test(index_measure ~ Month, p.adjust.method = "BH") %>% 
  add_xy_position(x = "Month")

obs.16s.bar = obs.16s.bar + 
  stat_pvalue_manual(stat.test.obs.16s.bar, label = "p.adj.signif", tip.length = 0.01, step.increase = 0.1)+
  theme(axis.title.x = element_text(size=12, colour = "black", margin = margin(t=10)),
        axis.text.x = element_text(size=12, colour = "black"),
        axis.title.y = element_text(size=12, colour = "black", margin = margin(r=10)),
        axis.text.y = element_text(size=12, colour = "black"))


(obs.its.bar = data.obs.ITS %>%
    mutate(cultivar_lab = ifelse(Cultivar == "VV", "Vardar Valley", "Justin Brouwers")) %>%
    ggbarplot(., x = "Month", y="index_measure", facet.by = "cultivar_lab", add = c("mean_se"), error.plot = "errorbar", ylab =  "Observed OTU richness", title = "Fungal community"))

stat.test.obs.its.bar =  data.obs.ITS %>%
  mutate(cultivar_lab = ifelse(Cultivar == "VV", "Vardar Valley", "Justin Brouwers")) %>%
  group_by(cultivar_lab) %>%
  t_test(index_measure ~ Month, p.adjust.method = "BH") %>% 
  add_xy_position(x = "Month")

obs.its.bar = obs.its.bar + 
  stat_pvalue_manual(stat.test.obs.its.bar, label = "p.adj.signif", tip.length = 0.01, step.increase = 0.01)+
  theme(axis.title.x = element_text(size=12, colour = "black", margin = margin(t=10)),
        axis.text.x = element_text(size=12, colour = "black"),
        axis.title.y = element_text(size=12, colour = "black", margin = margin(r=10)),
        axis.text.y = element_text(size=12, colour = "black"))

########## inverse Simpsons

(invS.16s.bar = data.invS.16s %>%
    mutate(cultivar_lab = ifelse(Cultivar == "VV", "Vardar Valley", "Justin Brouwers")) %>%
    ggbarplot(., x = "Month", y="index_measure", facet.by = "cultivar_lab", add = c("mean_se"), error.plot = "errorbar", ylab =  "Inverse Simpson's index", xlab = "", title = "Bacterial community"))

stat.test.invS.16s.bar =  data.invS.16s %>%
  mutate(cultivar_lab = ifelse(Cultivar == "VV", "Vardar Valley", "Justin Brouwers")) %>%
  group_by(cultivar_lab) %>%
  t_test(index_measure ~ Month, p.adjust.method = "BH") %>% 
  add_xy_position(x = "Month")

invS.16s.bar = invS.16s.bar + 
  stat_pvalue_manual(stat.test.invS.16s.bar, label = "p.adj.signif", tip.length = 0.01, step.increase = 0.2)+
  theme(axis.title.x = element_text(size=12, colour = "black", margin = margin(t=10)),
        axis.text.x = element_text(size=12, colour = "black"),
        axis.title.y = element_text(size=12, colour = "black", margin = margin(r=10)),
        axis.text.y = element_text(size=12, colour = "black"))


(invS.its.bar = data.invS.ITS %>%
  mutate(cultivar_lab = ifelse(Cultivar == "VV", "Vardar Valley", "Justin Brouwers")) %>%
 ggbarplot(., x = "Month", y="index_measure", facet.by = "cultivar_lab", add = c("mean_se"), error.plot = "errorbar", ylab =  "Inverse Simpson's index",  title = "Fungal community"))

stat.test.invS.its.bar =  data.invS.ITS %>%
  mutate(cultivar_lab = ifelse(Cultivar == "VV", "Vardar Valley", "Justin Brouwers")) %>%
  group_by(cultivar_lab) %>%
  t_test(index_measure ~ Month, p.adjust.method = "BH") %>% 
  add_xy_position(x = "Month")

invS.its.bar = invS.its.bar + 
  stat_pvalue_manual(stat.test.invS.its.bar, label = "p.adj.signif", tip.length = 0.01, step.increase = 0.01) +
  theme(axis.title.x = element_text(size=12, colour = "black", margin = margin(t=10)),
        axis.text.x = element_text(size=12, colour = "black"),
        axis.title.y = element_text(size=12, colour = "black", margin = margin(r=10)),
        axis.text.y = element_text(size=12, colour = "black"))

tiff("./Alpha_diversity/byMonth.tif", width = 5000, height = 3500, res = 400)
obs.16s.bar + invS.16s.bar +  obs.its.bar + invS.its.bar + plot_annotation(tag_levels = "a") 
dev.off()

####### statistical analysis #######################################################

str(data.obs.16s)

data.obs.16s$Cultivar = as.factor(data.obs.16s$Cultivar)

aov16S.obs = aov(index_measure ~ Chemical * Month * Cultivar, data = data.obs.16s)
aovITS.obs = aov(index_measure ~ Chemical * Month * Cultivar, data = data.obs.ITS)

summary(aov16S.obs)
plot(aov16S.obs)

summary(aovITS.obs)
plot(aovITS.obs)

?compare_means
compare_means(index_measure ~ Chemical, data = data.obs.16s, group.by = c("Cultivar", "Month"), p.adjust.method = "BH", method =  "t.test") %>%
  filter(p.adj < 0.05)


compare_means(index_measure ~ Chemical, data = data.obs.ITS, group.by = c("Cultivar", "Month"), p.adjust.method = "BH", method =  "t.test") %>%
  filter(p.adj < 0.05)

head(data.obs.16s)

?t.test


data.invS.ITS %>%
  filter(Cultivar == "VV") %>%
  t.test(index_measure ~ treated, data = ., paired=F, var.equal = F) 

aovITS.obs = aov(index_measure ~ Chemical * Month * Cultivar, data = data.obs.ITS)





summary(aov(index_measure ~ Chemical * Month * Cultivar, data = data.invS.16s))

summary(aov(index_measure ~ Chemical * Month * Cultivar, data = data.invS.ITS))

compare_means(index_measure ~ Chemical, data = data.invS.16s, group.by = c("Cultivar", "Month"), p.adjust.method = "BH", method =  "t.test") %>%
  filter(p.adj < 0.05)

compare_means(index_measure ~ Chemical, data = data.invS.ITS, group.by = c("Cultivar", "Month"), p.adjust.method = "BH", method =  "wilcox.test")


# month obs 16S
data.obs.16s %>%
  # filter(treated == "Treated") %>%
  compare_means(index_measure ~ Month, data = ., group.by = c("Cultivar"), p.adjust.method = "BH", method =  "t.test") %>%
  filter(p.adj < 0.05)

# Cultivar .y.           group1 group2         p  p.adj p.format p.signif method
# <fct>    <chr>         <chr>  <chr>      <dbl>  <dbl> <chr>    <chr>    <chr> 
#   1 JB       index_measure June   October 0.000815 0.0024 0.00081  ***      T-test
# 2 JB       index_measure August October 0.000207 0.0012 0.00021  ***      T-test
# 3 VV       index_measure August October 0.00882  0.018  0.00882  **       T-test

data.invS.16s %>%
  # filter(treated == "Treated") %>%
  compare_means(index_measure ~ Month, data = ., group.by = c("Cultivar"), p.adjust.method = "BH", method =  "t.test") %>%
  filter(p.adj < 0.05)

# Cultivar .y.           group1 group2          p   p.adj p.format p.signif method
# <chr>    <chr>         <chr>  <chr>       <dbl>   <dbl> <chr>    <chr>    <chr> 
#   1 JB       index_measure June   August  0.000397  0.00079 0.0004   ***      T-test
# 2 JB       index_measure August October 0.000899  0.0013  0.0009   ***      T-test
# 3 VV       index_measure June   August  0.0000730 0.00044 0.000073 ****     T-test
# 4 VV       index_measure August October 0.000204  0.00061 0.0002   ***      T-test

data.obs.ITS %>%
  # filter(treated == "Treated") %>%
  compare_means(index_measure ~ Month, data = ., group.by = c("Cultivar"), p.adjust.method = "BH", method =  "t.test") %>%
  filter(p.adj < 0.05)

# Cultivar .y.           group1 group2              p        p.adj p.format    p.signif method
# <chr>    <chr>         <chr>  <chr>           <dbl>        <dbl> <chr>       <chr>    <chr> 
#   1 VV       index_measure June   August  0.00000138    0.0000041    0.000001383 ****     T-test
# 2 VV       index_measure June   October 0.00000000101 0.0000000061 0.000000001 ****     T-test
# 3 VV       index_measure August October 0.000600      0.0012       0.0006      ***      T-test


data.invS.ITS %>%
  # filter(treated == "Treated") %>%
  compare_means(index_measure ~ Month, data = ., group.by = c("Cultivar"), p.adjust.method = "BH", method =  "t.test") %>%
  filter(p.adj < 0.05)

# Cultivar .y.           group1 group2          p   p.adj p.format p.signif method
# <chr>    <chr>         <chr>  <chr>       <dbl>   <dbl> <chr>    <chr>    <chr> 
#   1 JB       index_measure June   August  0.00209   0.0063  0.0021   **       T-test
# 2 JB       index_measure June   October 0.0000323 0.00019 0.000032 ****     T-test
