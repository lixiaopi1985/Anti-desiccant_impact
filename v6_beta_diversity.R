Wdir = "Your directory"
setwd(Wdir)


library(phyloseq)
library(tidyverse)
library(dplyr)
# library(microeco)
# library(file2meco)
# devtools::install_github("david-barnett/microViz")
# library(microViz)
library(ggplot2)
library(vegan)



rm(list = ls())

a16S = readRDS("./RDSdata/NA_free_THREEmonths_bsthit_pruned_16S_org.rds")

a16S@sam_data$treat = ifelse(a16S@sam_data$Chemical == "NT", "Nontreated", "Treated")
a16S@sam_data$treat = factor(a16S@sam_data$treat, levels = c("Nontreated", "Treated"))



source("./VT_tools/plot_betadiv.R")

# overall by treated and nontreated
a16S.jb = subset_samples(a16S, Cultivar == "JB")
a16S.vv = subset_samples(a16S, Cultivar == "VV")

plot_betadiv(a16S.jb, group_col = "treat", norm_method = "hellinger", plot_title = "Justin Brouwers", Pcoa_legend = F, nudge = 0.3)

plot_betadiv(a16S.vv, group_col = "treat", norm_method = "hellinger", plot_title = "Vardar Valley", Pcoa_legend = F, nudge = 0.3)




#meta16S = data.frame(sample_data(a16S))
#levels(meta16S$Chemical) = c("Nontreated", "TransFilm", "Vapor Gard", "Wilt-Pruf")

#sample_data(a16S) = meta16S

# subset data
a16S.june.jb = subset_samples(a16S, Month == "June" & Cultivar == "JB")
a16S.aug.jb = subset_samples(a16S, Month == "August" & Cultivar == "JB")
a16S.oct.jb = subset_samples(a16S, Month == "October" & Cultivar == "JB")

a16S.june.vv = subset_samples(a16S, Month == "June" & Cultivar == "VV")
a16S.aug.vv = subset_samples(a16S, Month == "August" & Cultivar == "VV")
a16S.oct.vv = subset_samples(a16S, Month == "October" & Cultivar == "VV")




pal.16S = c("brown3", "blue2", "cyan4", "darkgoldenrod3")
newleg = c("NT=Nontreated", "TF=TransFilm", "VG=Vapor Gard", "WP=Wilt-Pruf")





(plot16S.jb.6 = plot_betadiv(a16S.june.jb, group_col = "Chemical", norm_method = "hellinger", plot_title = "", Pcoa_legend = F,palette = pal.16S, nudge = 0.3))
(plot16S.jb.8 = plot_betadiv(a16S.aug.jb, group_col = "Chemical", norm_method = "hellinger", plot_title = "", Pcoa_legend = T, palette = pal.16S, nudge = 0.3, legend_pos = "bottom", legend_title = "Anti-desiccant", ifsetLegendLabel = T, newLegendLabel = newleg, legCol = 2, legRow = 2))
(plot16S.jb.10 = plot_betadiv(a16S.oct.jb, group_col = "Chemical", norm_method = "hellinger", plot_title = "", Pcoa_legend = F, palette = pal.16S, nudge = 0.3))

#----------------------------------------------------------------------
(plot16S.vv.6 = plot_betadiv(a16S.june.vv, group_col = "Chemical", norm_method = "hellinger", plot_title = "", Pcoa_legend = F,palette = pal.16S))
(plot16S.vv.8 = plot_betadiv(a16S.aug.vv, group_col = "Chemical", norm_method = "hellinger", plot_title = "", Pcoa_legend = F,palette = pal.16S, nudge = 0.3))
(plot16S.vv.10 = plot_betadiv(a16S.oct.vv, group_col = "Chemical", norm_method = "hellinger", plot_title = "",  Pcoa_legend = F,palette = pal.16S, nudge = 0.25))

# tiff('./Beta/B16S_beta.tif', width = 6500, height = 4200, res=400)
# plot16S.jb.6 +
#   plot16S.jb.8 +
#   plot16S.jb.10 +
#   plot16S.vv.6 +
#   plot16S.vv.8 +
#   plot16S.vv.10
# dev.off()

svg("./Beta/B16S_beta.svg", width = 18, height = 12, pointsize = 30)
plot16S.vv.6 +
  plot16S.vv.8 +
  plot16S.vv.10 +
  plot16S.jb.6 +
  plot16S.jb.8 +
  plot16S.jb.10 
dev.off()

##################################################################################################

rm(list = ls())
source("./tools/plot_betadiv.R")

ITS = readRDS("./RDSdata/NA_free_THREEmonths_bsthit_pruned_ITS_org.rds")


# ITS = readRDS("./RDSdata/NA_free_THREEmonths_bsthit_pruned_ITS_org.rds")
ITS@sam_data$treat = ifelse(ITS@sam_data$Chemical == "NT", "Nontreated", "Treated")
ITS@sam_data$treat = factor(ITS@sam_data$treat, levels = c("Nontreated", "Treated"))

ITS.jb = subset_samples(ITS, Cultivar == "JB")
ITS.vv = subset_samples(ITS, Cultivar == "VV")

plot_betadiv(ITS.jb, group_col = "treat", norm_method = "hellinger", plot_title = "Justin Brouwers", Pcoa_legend = F, nudge = 0.3, CI_Ellipse = T)

plot_betadiv(ITS.vv, group_col = "treat", norm_method = "hellinger", plot_title = "Vardar Valley", Pcoa_legend = F, nudge = 0.3, CI_Ellipse = T)




#metaITS = data.frame(sample_data(ITS))
#levels(metaITS$Chemical) = c("Nontreated", "TransFilm", "Vapor Gard", "Wilt-Pruf")

#sample_data(ITS) = metaITS



ITS.june.jb = subset_samples(ITS, Month == "June" & Cultivar == "JB")
ITS.aug.jb =  subset_samples(ITS, Month == "August" & Cultivar == "JB")
ITS.oct.jb =  subset_samples(ITS, Month == "October" & Cultivar == "JB")

ITS.june.vv = subset_samples(ITS, Month == "June" & Cultivar == "VV")
ITS.aug.vv =  subset_samples(ITS, Month == "August" & Cultivar == "VV")
ITS.oct.vv =  subset_samples(ITS, Month == "October" & Cultivar == "VV")

#str(data.frame(ITS.june.jb@sam_data))

#source("./plot_betadiv.R")

newleg = c("NT=Nontreated", "TF=TransFilm", "VG=Vapor Gard", "WP=Wilt-Pruf")
pal.fun = c("chartreuse4", "blue3", "darkgoldenrod2", "brown3")
# pal.fun.8 = c("chartreuse4", "lightskyblue", "maroon1", "pink1")
# pal.fun.10 = c("chartreuse4", "royalblue4", "red3", "tan3")


# (plotITS.jb.6 = plot_betadiv(ITS.june.jb, group_col = "treat", norm_method = "hellinger", plot_title = "Justin Brouwers - June", Pcoa_legend = F, nudge = 0.3))
# (plotITS.jb.8 = plot_betadiv(ITS.aug.jb, group_col = "treat", norm_method = "hellinger", plot_title = "Justin Brouwers - June", Pcoa_legend = F, nudge = 0.3))
# (plotITS.jb.10 = plot_betadiv(ITS.oct.jb, group_col = "treat", norm_method = "hellinger", plot_title = "Justin Brouwers - June", Pcoa_legend = F, nudge = 0.3))

(plotITS.jb.6 = plot_betadiv(ITS.june.jb, group_col = "Chemical", norm_method = "hellinger", plot_title = "", Pcoa_legend = F, palette = pal.fun, nudge = 0.3))
(plotITS.jb.8 =plot_betadiv(ITS.aug.jb, group_col = "Chemical", norm_method = "hellinger", plot_title = "", Pcoa_legend = T, palette = pal.fun, nudge = 0.3, legend_pos = "bottom", legend_title = "Anti-desiccant", ifsetLegendLabel = T, newLegendLabel = newleg, legRow  = 2, legCol = 2))
(plotITS.jb.10 =plot_betadiv(ITS.oct.jb, group_col = "Chemical", norm_method = "hellinger", plot_title = "", Pcoa_legend = F, palette = pal.fun, nudge = 0.4))

#########################################################################################

# by treatment
# (plotITS.vv.6 = plot_betadiv(ITS.june.vv, group_col = "treat", norm_method = "hellinger", plot_title = "Vardar Valley - June", Pcoa_legend = F, nudge = 0.3))
# (plotITS.vv.8 = plot_betadiv(ITS.aug.vv, group_col = "treat", norm_method = "hellinger", plot_title = "Vardar Valley - June", Pcoa_legend = F, nudge = 0.3))
# (plotITS.vv.10 = plot_betadiv(ITS.oct.vv, group_col = "treat", norm_method = "hellinger", plot_title = "Vardar Valley - June", Pcoa_legend = F, nudge = 0.3))


(plotITS.vv.6 =plot_betadiv(ITS.june.vv, group_col = "Chemical", norm_method = "hellinger", plot_title = "", Pcoa_legend = F, palette = pal.fun, nudge = 0.3))
(plotITS.vv.8 =plot_betadiv(ITS.aug.vv, group_col = "Chemical", norm_method = "hellinger", plot_title = "", Pcoa_legend = F, palette = pal.fun, nudge = 0.3))
(plotITS.vv.10 =plot_betadiv(ITS.oct.vv, group_col = "Chemical", norm_method = "hellinger", plot_title = "", Pcoa_legend = F, palette = pal.fun, nudge = 0.3))

tiff("./Beta/ITS_beta.tif", width = 6600, height = 4300, res=400)
plotITS.vv.6 +
  plotITS.vv.8 +
  plotITS.vv.10 +
  plotITS.jb.6 +
  plotITS.jb.8 +
  plotITS.jb.10
dev.off()

svg("./Beta/ITS_beta.svg", width = 18, height = 12, pointsize = 30)
plotITS.vv.6 +
  plotITS.vv.8 +
  plotITS.vv.10 +
plotITS.jb.6 +
  plotITS.jb.8 +
  plotITS.jb.10
dev.off()
