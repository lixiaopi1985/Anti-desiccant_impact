plotHeat_matrix = function(psmelt, taxlevel="Genus", groupby = c("period", "cul_type"), ntop=20, nshow=10, printlevel=0){
  # make heatmap of this (top 10?)
  library(tidyr)
  
  top20 = psmelt %>%
    group_by(!!!syms(groupby)) %>%
    mutate(total_abund = sum(Abundance)) %>%
    mutate(new_taxa =  !!sym(taxlevel)) %>%
    group_by(!!!syms(c(groupby, "new_taxa"))) %>%
    mutate(rel_abund = 100*Abundance/total_abund) %>%# change cul_type to groupby
    summarise(accum_abund = sum(rel_abund)) %>%
    slice_max(order_by = accum_abund, n=ntop)
  #----------------------------------------------------------------------


  # overall top abundant genera
  domdf = psmelt %>%
    mutate(total_abund = sum(Abundance), rel_abund = 100*Abundance/total_abund) %>%
    mutate(new_taxa = !!sym(taxlevel)) %>%
    group_by(new_taxa) %>%
    mutate(rel_abund = 100*Abundance/total_abund) %>%
    summarise(accum_abund = sum(rel_abund)) %>%
    arrange(desc(accum_abund))
  
  doms = domdf$new_taxa[1:nshow] #ordered
  
  # some top genus might be not in other group, so find them as well
  top_accum = psmelt %>%
    group_by(!!!syms(groupby)) %>%
    mutate(total_abund = sum(Abundance), rel_abund = 100*Abundance/total_abund) %>%
    mutate(new_taxa = !!sym(taxlevel)) %>%
    group_by(!!!syms(c(groupby, "new_taxa"))) %>%
    mutate(rel_abund = 100*Abundance/total_abund) %>%
    summarise(accum_abund = sum(rel_abund)) %>%
    filter(new_taxa %in% unique(top20$new_taxa)) %>%
    arrange(desc(accum_abund))

  top.meta.2 = top_accum %>%
    unite(col=Samples, !!!syms(groupby)) %>%
    spread(Samples, accum_abund) %>%
    as.data.frame()
  
  rownames(top.meta.2) = top.meta.2$new_taxa
  
  # matrix
  top.matrix = as.matrix(top.meta.2[, -1])
  
  if(printlevel==0){
    return(list(m = top.matrix, topbygroup=top_accum, dom = doms, overal=domdf))
  } else if(printlevel == 1){
    
    prinL.df = psmelt %>%
      dplyr::select(Phylum, !!sym(taxlevel)) %>%
      mutate(taxa_cbind = paste(Phylum, !!sym(taxlevel), sep = "-")) %>%
      filter(!!sym(taxlevel) %in% doms) %>%
      dplyr::distinct()
    
    # order it
    new_doms = prinL.df$taxa_cbind[match(doms, prinL.df[[taxlevel]])]
    
    # keep the order right
    
    return(list(m = top.matrix, topbygroup = top_accum, dom = new_doms, overal=domdf))
  }
  
  
}

#---------------------------------------------------------------------------------------------------

# plotHeat_plot = function(tops, COLOR=NA, new_break = NA,fontsize=12, out=T, outfile="heatmap.png", clusterRows=F, clusterCols=F, Legend=T, anno_legend=T, custom_anno=T, annoInput=NA, annoColor=NA, cellw = 40, cellh = 40, ...){
#   
#   library(stringr)
#   library(RColorBrewer)
#   library(pheatmap)
#   
#   top.matrix = tops$m
#   doms = tops$dom
#   
#   org_colname = colnames(top.matrix)
#   
#   color_set = brewer.pal(9, "BuGn")
#   colorTime = colorRampPalette(color_set)
#   
#   if(custom_anno){
#     if(!is.na(annoInput)){
#       
#       top.matrix.ord2 = top.matrix[doms,]
#       anno_colors = annoColor
#       
#       if(out){
#         g = pheatmap::pheatmap(top.matrix.ord2, 
#                                annotation_col = annoInput,
#                                annotation_colors = anno_colors, 
#                                cellwidth = cellw, 
#                                cellheight = cellh, 
#                                color = colorRampPalette(rev(brewer.pal(n=7, name = "RdYlBu")))(length(new_break)), 
#                                cluster_cols = clusterCols, 
#                                cluster_rows = clusterRows, 
#                                breaks = new_break, 
#                                fontsize = fontsize, 
#                                filename = outfile, 
#                                legend = Legend,
#                                annotation_legend = anno_legend, 
#                                # column_split = annoCol$Fumigation_history, 
#                                ...)
#         
#         
#       } else {
#         g = pheatmap::pheatmap(top.matrix.ord2, 
#                                annotation_col = annoCol,
#                                cellwidth = cellw, 
#                                cellheight = cellh,
#                                # color = colorRampPalette(c("navy", "white", "firebrick3"))(length(new_break)),
#                                color = colorRampPalette(rev(brewer.pal(n=7, name = "RdYlBu")))(length(new_break)), 
#                                annotation_colors = anno_colors, 
#                                cluster_cols = clusterCols, 
#                                cluster_rows = clusterRows, 
#                                breaks = new_break, 
#                                fontsize = fontsize, 
#                                # filename = outfile, 
#                                legend = Legend,
#                                annotation_legend = anno_legend, 
#                                # column_split = annoCol$Fumigation_history, 
#                                ...)
#       }
#       
#     } else {
#       stop("Please provide annotation")
#     }
#     
#   } else {
#     
#     # no annotation
#     top.matrix.ord2 = top.matrix[dom, ]
#     
#     if(out){
#       g = pheatmap::pheatmap(top.matrix.ord2, 
#                              # annotation_col = annoCol,
#                              cellwidth = cellw, 
#                              cellheight = cellh, 
#                              color = colorRampPalette(rev(brewer.pal(n=7, name = "RdYlBu")))(length(new_break)), 
#                              # annotation_colors = anno_colors, 
#                              cluster_cols = clusterCols, 
#                              cluster_rows = clusterRows, 
#                              breaks = new_break, 
#                              fontsize = fontsize, 
#                              filename = outfile, 
#                              legend = Legend,
#                              annotation_legend = anno_legend,
#                              # column_split = annoCol$Fumigation_history, 
#                              ...)
#     } else {
#       g = pheatmap::pheatmap(top.matrix.ord2, 
#                              # annotation_col = annoCol,
#                              cellwidth = cellw, 
#                              cellheight = cellh, 
#                              color = colorRampPalette(rev(brewer.pal(n=7, name = "RdYlBu")))(length(new_break)), 
#                              # annotation_colors = anno_colors, 
#                              cluster_cols = clusterCols, 
#                              cluster_rows = clusterRows, 
#                              breaks = new_break, 
#                              fontsize = fontsize, 
#                              # filename = outfile, 
#                              legend = Legend,
#                              annotation_legend = anno_legend,
#                              # column_split = annoCol$Fumigation_history, 
#                              ...)
#     }
#    
#   }
#   
#   return(g)
# }


# add significant * to the row labels
addSig = function(origlabel, cond1, cond2){
  
  # italicized the organisms
  new_labs = c()
  for(i in origlabel){
    if((i %in% cond1) & (i %in% cond2)){
      sigLab = paste0(i, " */*")
      new_labs = c(new_labs, sigLab)
    } else if((i %in% cond1) & !(i %in% cond2)){
      sigLab = paste0(i, " */")
      new_labs = c(new_labs, sigLab)
    } else if(!(i %in% cond1) & (i %in% cond2)){
      sigLab = paste0(i, " /*")
      new_labs = c(new_labs, sigLab)
    } else {
      sigLab = i
      new_labs = c(new_labs, sigLab)
    }
  }
  
  return(new_labs)
}

scientific_name_formatter <- function(raw_name, plotdevice="ggplot") { 
  library(stringr)
  # containing p__, c__, o__, f__, g__, s__
  # name contain numbers and other 
  
  if(plotdevice=="ggplot"){
    formatted_names = sapply(raw_name, function(x){
      
      print(x)
      nsplit = str_split(x, "__")
      
      if(length(nsplit[[1]]) == 2){
        # if the second name contains numbers, -, length > 1, all capitalized do not italicize
        if(grepl("[0-9]+|-|uncultured", nsplit[[1]][2])){
          spname = paste0(nsplit[[1]][1], "__", nsplit[[1]][2])
        } else {
          spname = paste0(nsplit[[1]][1], "__*", nsplit[[1]][2], "*")
        }
        
      } else if(length(nsplit[[1]] == 1)) {
        
        if(grepl("[0-9]+|-|uncultured", nsplit[[1]][1])){
          spname = nsplit[[1]][1]
        } else {
          spname = paste0("*", nsplit[[1]][1], "*")
        }
      }
      return(spname)
    })
    
    
  } else {
    # other format
    
    formatted_names = sapply(raw_name, function(x){
      
      nsplit = str_split(x, "__")
      print(x)
      
      if(length(nsplit[[1]]) == 2){

        if(grepl("[0-9]+|-|uncultured", nsplit[[1]][2])){
          
          # if the second name contains numbers, -, length > 1, all capitalized do not italicize
          spname = toString(paste0(nsplit[[1]][1], "__", nsplit[[1]][2]))
          
          
          
        } else {
          
          # italicize
          # to extract "*/*" from f__XXXXX */*
          
          if(grepl("\\*/$|\\*/\\*$|/\\*$", nsplit[[1]][2])){
            
            # separate out "/*" from the string first
            nsplit_space = str_split(nsplit[[1]][2], " ")
            
            if(length(nsplit_space[[1]]) == 2){
              # contains "/*" etc
              spname = bquote(.(nsplit[[1]][1]) * "__" * italic(.(nsplit_space[[1]][1])) ~ .(nsplit_space[[1]][2]))
              
              
            } else if(length(nsplit_space[[1]]) == 1) {
              # does not contain "/*
              spname = bquote(.(nsplit[[1]][1]) * "__" * italic(.(nsplit_space[[1]][1])))
              
            }
          } else {
            # if f__XXXX does not contain "*/", "*/*", "/*"
            spname = bquote(.(nsplit[[1]][1]) * "__" * italic(.(nsplit[[1]][2])))
            
          }


        }
        
      } else if(length(nsplit[[1]]) == 1) {
        print("Length 1")
        if(grepl("[0-9]+|-|uncultured", nsplit[[1]][1])){
          spname = nsplit[[1]][1]
        } else {
          
          # situation like XXXX */
          if(grepl("\\*/$|\\*/\\*$|/\\*$", nsplit[[1]][1])){
            
            # separate out "/*" from the string first
            nsplit_space = str_split(nsplit[[1]][1], " ")
            
            if(length(nsplit_space[[1]]) == 2){
              # contains "/*" etc
              spname = bquote(italic(.(nsplit_space[[1]][1])) ~ .(nsplit_space[[1]][2]))
              
            } else if(length(nsplit_space[[1]]) == 1) {
              # does not contain "/*
              spname = bquote(italic(.(nsplit[[1]][1])))
            }
          } else {
            # if f__XXXX does not contain "*/", "*/*", "/*"
            spname = bquote(italic(.(nsplit[[1]][1])))
          }
        }
      }
      return(spname)
    })
  }

  
  return(unname(formatted_names))
}


plotBars = function(glomdf, taxlevel="Genus", groupby = c("Season", "time_label", "funlab"), ntop=20, nshow=10, customColors = NULL, choosecolorset = NULL, otherColor = "grey45", barTransparent = 1, plotx = "funlab", facet_group = "~ Season + time_label", axis_y_unit="K", axis_y_scale=1e-3, xlab="Fungicides", ylab="Relative abundance (%)", labfill="", plottitle="", xlabmargin_t=10, ylabmargin_r=10, fontsize=12, printlevels=0, fontangle=0, fontHjust=1, fontVjust=0.5, DFonly=F, plotMargin_t = 2, plotMargin_r = 2, plotMargin_b = 2, plotMargin_l = 2, abs_abund = FALSE, titlejust = 0, italLabel = F){
  
  library(ggplot2)
  library(ggh4x)
  library(dplyr)
  library(rlang)
  
  top = plotHeat_matrix(glomdf, taxlevel = taxlevel, groupby = groupby, ntop = ntop, nshow = nshow, printlevel = printlevels)

  doms = top$dom[!is.na(top$dom)]
  
  cat("Top taxa selected:\n")
  print(paste(doms, sep = "\n"))

  if(italLabel){
    if(printlevels==0){
      
      Togroup = c(groupby, taxlevel)
      
      # sp. not italic, p__, o__, c__, f__, g__, c__Subgroup 6, c__KD4-96 not italic
      topn = glomdf %>%
        group_by(!!!syms(Togroup)) %>%
        summarise(accum_abund = sum(Abundance)) %>%
        mutate(new_taxa = ifelse(!!sym(taxlevel) %in% doms, !!sym(taxlevel), "Others"))%>%
        mutate(new_taxa2 = ifelse(new_taxa == "Others", "Others", scientific_name_formatter(new_taxa))) %>% # format the name
        mutate(new_taxa_itali = ifelse(new_taxa2 == "Others", 
                                       "Others",
                                       ifelse(grepl(" sp.*", new_taxa2),
                                              paste0(paste0(gsub(" sp.*$", "", new_taxa2), "*"), " spp."),
                                              new_taxa2)))
      
      # mutate(new_taxa_itali = ifelse(new_taxa2 == "Others", "Others",
      #                                ifelse(grepl(" sp.", new_taxa2),  
      #                                       paste(paste0("*", gsub(" sp.$", "", new_taxa2), "*"), " spp."), 
      #                                       new_taxa2)))
      
      
      # if new_taxa == other, no change, than look for sp in the species names
      
      
      #topn$new_taxa = factor(topn$new_taxa, levels = c(doms, "Others"))
      
      # changed it to italic
      doms_italic = scientific_name_formatter(doms)
      
      if("Others" %in% unique(topn$new_taxa_itali)){
        dom_level = c(ifelse(grepl(" sp.*", doms_italic),
                             paste(paste0(gsub(" sp.*$", "", doms_italic), "*"), " spp."),
                             doms_italic), "Others")
      } else {
        dom_level = ifelse(grepl(" sp.*", doms_italic),
                           paste(paste0(gsub(" sp.*$", "", doms_italic), "*"), " spp."),
                           doms_italic)
      }
      
      
      
      # dom_level = c(ifelse(grepl(" sp.", doms), 
      #                    paste(paste0("*", gsub(" sp.$", "", doms), "*"), " spp."), 
      #                    paste0("*", doms, "*")), "Others")
      
      
      topn$new_taxa_itali = factor(topn$new_taxa_itali, levels = dom_level)
      
      
    } else if(printlevels == 1){
      # add Phylum to it
      Togroup = c(groupby, "taxa_cbind", taxlevel)
      
      topn = glomdf %>%
        mutate(taxa_cbind = paste(Phylum, !!sym(taxlevel), sep = "-")) %>%
        group_by(!!!syms(Togroup)) %>%
        summarise(accum_abund = sum(Abundance)) %>%
        mutate(new_taxa = ifelse(taxa_cbind %in% doms, taxa_cbind, "Others")) %>%
        mutate(new_taxa2 = ifelse(new_taxa == "Others", "Others", scientific_name_formatter(new_taxa))) %>%
        mutate(new_taxa_itali = ifelse(new_taxa2 == "Others", "Others",
                                       ifelse(grepl(" sp.*", new_taxa2),
                                              paste0(paste0(gsub(" sp.*$", "", new_taxa2), "*"), " spp."),
                                              new_taxa2)))
      
      # mutate(new_taxa_itali = ifelse(new_taxa == "Others", "Others",
      #                                ifelse(grepl(" sp.", new_taxa),  
      #                                       paste(paste0("*", gsub(" sp.$", "", new_taxa), "*"), " spp."), 
      #                                       paste0("*", new_taxa, "*")))) 
      #mutate(new_taxa_itali = ifelse(new_taxa != "Others", paste0("*", new_taxa, "*"), new_taxa))
      #!!!!!!!!!!!!!!!!!!! new_taxa_itali has not tested
      
      # order new_taxa
      topn$new_taxa = factor(topn$new_taxa, levels = c(doms, "Others"))
      #topn$new_taxa_itali = factor(topn$new_taxa_itali, levels = c(paste0("*", doms, "*"), "Others"))
      doms_italic = scientific_name_formatter(doms)
      
      if("Others" %in% unique(topn$new_taxa_itali)){
        dom_level = c(ifelse(grepl(" sp.*", doms_italic),
                             paste(paste0(gsub(" sp.*$", "", doms_italic), "*"), " spp."),
                             doms_italic), "Others")
      } else {
        dom_level = ifelse(grepl(" sp.*", doms_italic),
                           paste(paste0(gsub(" sp.*$", "", doms_italic), "*"), " spp."),
                           doms_italic)
      }
      
      
      # dom_level = c(ifelse(grepl(" sp.", doms), 
      #                      paste(paste0("*", gsub(" sp.$", "", doms), "*"), " spp."), 
      #                      paste0("*", doms, "*")), "Others")
      
      topn$new_taxa_itali = factor(topn$new_taxa_itali, levels = dom_level)
      
    }
  } else {
    
    Togroup = c(groupby, taxlevel)
    
    # sp. not italic, p__, o__, c__, f__, g__, c__Subgroup 6, c__KD4-96 not italic
    topn = glomdf %>%
      group_by(!!!syms(Togroup)) %>%
      summarise(accum_abund = sum(Abundance)) %>%
      mutate(new_taxa_itali = ifelse(!!sym(taxlevel) %in% doms, !!sym(taxlevel), "Others"))
      # mutate(new_taxa2 = ifelse(new_taxa == "Others", "Others", scientific_name_formatter(new_taxa))) %>% # format the name
      # mutate(new_taxa_itali = ifelse(new_taxa2 == "Others", 
      #                                "Others",
      #                                ifelse(grepl(" sp.*", new_taxa2),
      #                                       paste0(paste0(gsub(" sp.*$", "", new_taxa2), "*"), " spp."),
      #                                       new_taxa2)))
      # 
    dom_level = c(doms, "Others")
    topn$new_taxa_itali = factor(topn$new_taxa_itali, levels = dom_level)
  }
  
  

  

  
  if(is.null(customColors)){
    if(is.null(choosecolorset)){
      colors = ggsci::pal_npg()(10)
    } else {
      colors = choosecolorset
    }
    
    colorPal = colorRampPalette(colors)(length(unique(topn$new_taxa_itali)))
    names(colorPal) = levels(topn$new_taxa_itali)
    colorPal["Others"] = otherColor
    
  } else {
    colorPal = customColors
  }


  
  if(is.null(labfill)){
    labfill = NULL 
  } else if(labfill==""){
    labfill = taxlevel
  }
  
  
  if(DFonly){
    return(topn)
  }
  
  if(facet_group == ""){
    
    if(abs_abund){
      warnings("Plotting Absolute abundance")
      g = topn %>%
        ggplot(aes(x = !!sym(plotx), y=accum_abund, fill=new_taxa_itali))+
        geom_col(position = position_stack(), alpha=barTransparent) +
        scale_fill_manual(values = colorPal) +
        # scale_y_continuous(expand = c(0,0)) +
        scale_y_continuous(expand = c(0,0), labels = scales::unit_format(unit = axis_y_unit, scale = axis_y_scale)) +
        theme_classic() +
        labs(x=xlab, y=ylab, fill=labfill, title=plottitle) +
        theme(
          plot.title = element_text(hjust = titlejust),
          legend.text = ggtext::element_markdown(),
          text = element_text(size=fontsize, color="black"),
          #plot.margin = margin(t=plotMargin_t, r=plotMargin_r, b=plotMargin_b, l=plotMargin_l),
          legend.margin = margin(t=plotMargin_t, r=plotMargin_r, b=plotMargin_b, l=plotMargin_l),
          axis.text.x = element_text(size=fontsize, color="black", angle = fontangle, hjust = fontHjust, vjust = fontVjust),
          axis.text.y = element_text(size=fontsize, color="black"),
          axis.title.x = element_text(size=fontsize, margin = margin(t=xlabmargin_t), color="black"),
          axis.title.y = element_text(size=fontsize, margin = margin(r=ylabmargin_r), color="black"))
    } else {
      g = topn %>%
        ggplot(aes(x = !!sym(plotx), y=accum_abund, fill=new_taxa_itali))+
        geom_col(position = position_fill(), alpha=barTransparent) +
        scale_fill_manual(values = colorPal) +
        # scale_y_continuous(expand = c(0,0)) +
        scale_y_continuous(expand = c(0,0), labels = scales::percent_format(accuracy =1, suffix = "")) +
        theme_classic() +
        labs(x=xlab, y=ylab, fill=labfill, title=plottitle) +
        theme(
          plot.title = element_text(hjust = titlejust),
          legend.text = ggtext::element_markdown(),
          text = element_text(size=fontsize, color="black"),
          #plot.margin = margin(t=plotMargin_t, r=plotMargin_r, b=plotMargin_b, l=plotMargin_l),
          legend.margin = margin(t=plotMargin_t, r=plotMargin_r, b=plotMargin_b, l=plotMargin_l),
          axis.text.x = element_text(size=fontsize, color="black", angle = fontangle, hjust = fontHjust, vjust = fontVjust),
          axis.text.y = element_text(size=fontsize, color="black"),
          axis.title.x = element_text(size=fontsize, margin = margin(t=xlabmargin_t), color="black"),
          axis.title.y = element_text(size=fontsize, margin = margin(r=ylabmargin_r), color="black"))
      
    }

    

  } else {
    if(abs_abund){
      warnings("Plotting Absolute abundance")
      g = topn %>%
        ggplot(aes(x = !!sym(plotx), y=accum_abund, fill=new_taxa_itali))+
        geom_col(position = position_stack(), alpha=barTransparent) +
        facet_nested(formula(facet_group), scales = "free") +
        scale_fill_manual(values = colorPal) +
        # scale_y_continuous(expand = c(0,0)) +
        scale_y_continuous(expand = c(0,0), labels = scales::unit_format(unit = axis_y_unit, scale = axis_y_scale)) +
        theme_classic() +
        labs(x=xlab, y=ylab, fill=labfill, title = plottitle) +
        theme(
          plot.title = element_text(hjust = titlejust),
          legend.text = ggtext::element_markdown(),
          text = element_text(size=fontsize, color="black"),
          #plot.margin = margin(t=plotMargin_t, r=plotMargin_r, b=plotMargin_b, l=plotMargin_l),
          legend.margin = margin(t=plotMargin_t, r=plotMargin_r, b=plotMargin_b, l=plotMargin_l),
          axis.text.x = element_text(size=fontsize, color="black",angle = fontangle, hjust = fontHjust, vjust = fontVjust),
          axis.text.y = element_text(size=fontsize, color="black"),
          axis.title.x = element_text(size=fontsize, margin = margin(t=xlabmargin_t), color="black"),
          axis.title.y = element_text(size=fontsize, margin = margin(r=ylabmargin_r), color="black"))
    } else {
      g = topn %>%
        ggplot(aes(x = !!sym(plotx), y=accum_abund, fill=new_taxa_itali))+
        geom_col(position = position_fill(), alpha=barTransparent) +
        facet_nested(formula(facet_group), scales = "free") +
        scale_fill_manual(values = colorPal) +
        # scale_y_continuous(expand = c(0,0)) +
        scale_y_continuous(expand = c(0,0), labels = scales::percent_format(accuracy =1, suffix = "")) +
        theme_classic() +
        labs(x=xlab, y=ylab, fill=labfill, title = plottitle) +
        theme(
          plot.title = element_text(hjust = titlejust),
          legend.text = ggtext::element_markdown(),
          text = element_text(size=fontsize, color="black"),
          #plot.margin = margin(t=plotMargin_t, r=plotMargin_r, b=plotMargin_b, l=plotMargin_l),
          legend.margin = margin(t=plotMargin_t, r=plotMargin_r, b=plotMargin_b, l=plotMargin_l),
          axis.text.x = element_text(size=fontsize, color="black",angle = fontangle, hjust = fontHjust, vjust = fontVjust),
          axis.text.y = element_text(size=fontsize, color="black"),
          axis.title.x = element_text(size=fontsize, margin = margin(t=xlabmargin_t), color="black"),
          axis.title.y = element_text(size=fontsize, margin = margin(r=ylabmargin_r), color="black"))
    }

  }
  return(g)
}


plotORDS = function(phylo, tranform_method = "hellinger", ordmethod="PCoA", distmethod="bray", fs=12, plottype="sample", colorby="funlab", colorLab = "Fungicides", fontcolor="black", alpha=1, ncolor = 4, colorPal = ggsci::pal_npg(), revCol = F, custom_colors = NULL, CI_Ellipse=F, DATA_Ellipse=F, ell_type="t", ell_linetype=2, returnDist = F, gtitle="", plotDist = T, box_alpha=0.5,pairwise_method = "HSD", neg_correct = "none", rn= NULL, jitter_size = 2, point_size=5,letter_size = 5, plotPC_xlab = "", plotPC1_ylab = "", plotPC2_ylab = "", fontAngle=0, fontVjust=0, fontHjust=0, legendPos = "left"){

  library(patchwork)
  library(agricolae)
  
  # ?transform
  phylo.norm = microbiome::transform(phylo, transform = tranform_method)
  
  phylo.ord = ordinate(phylo.norm, method = ordmethod, distmethod)
  
  meta.phy = data.frame(sample_data(phylo.norm))
  
  col_values = colorPal(ncolor)
  
  if(!is.null(custom_colors)){
    col_values = custom_colors
  }
  
  if(revCol){
    col_values = rev(col_values)
  }
  

  
  
  
  g =  plot_ordination(phylo.norm, phylo.ord, type=plottype, color=colorby) +
    geom_vline(aes(xintercept=0), color = "grey80") +
    geom_hline(aes(yintercept=0), color = "grey80") +
    geom_point(size=5, alpha=alpha) +
    scale_color_manual(values = col_values) +
    theme_bw() +
    labs(color=colorLab, title=gtitle) +
    theme(text = element_text(size=fs, color= fontcolor),
          axis.title.x = element_text(size=fs, color = fontcolor, margin=margin(t=10)),
          axis.title.y = element_text(size=fs, color=fontcolor, margin=margin(r=10)),
          axis.text.x = element_text(size=fs, color=fontcolor),
          axis.text.y = element_text(size=fs, color=fontcolor),
          legend.position = legendPos)
  
  if(CI_Ellipse){
    g = g+ggpubr::stat_conf_ellipse(size=1, linetype = ell_linetype)
    
  } else if(DATA_Ellipse){
    g = g+stat_ellipse(type=ell_type, linetype=ell_linetype)
  }
  
  
  
  if(plotDist){
    # create distance object
    phy.dist = vegan::vegdist(t(otu_table(phylo.norm)), method = distmethod)
    phy.dist
    
    # get pcoa ordination
    pcoa = ape::pcoa(phy.dist, correction = neg_correct, rn=rn)
    pcoa
    
    PC1 = pcoa$vectors[, 1] 
    PC2 = pcoa$vectors[, 2]
    PC3 = pcoa$vectors[, 3]
    

    plotdata <- data.frame(rownames(pcoa$vectors),PC1,PC2,PC3, meta.phy[[colorby]])
    colnames(plotdata) = c("sample", "PC1", "PC2", "PC3", "Group")
    
    name_group = unique(plotdata$Group)
    
    pc1 <-round(pcoa$values$Relative_eig[1]*100,digits = 2)
    pc2 <-round(pcoa$values$Relative_eig[2]*100,digits = 2)
    pc3 <-round(pcoa$values$Relative_eig[3]*100,digits = 2)
    

    if(pairwise_method %in% c('LSD','HSD','duncan','scheffe','REGW','SNK')){
      
      yf <- plotdata
      yd1 <- yf %>% dplyr::group_by(Group) %>% dplyr::summarise(Max = max(PC1))
      yd2 <- yf %>% dplyr::group_by(Group) %>% dplyr::summarise(Max = max(PC2))
      yd3 <- yf %>% dplyr::group_by(Group) %>% dplyr::summarise(Max = max(PC3))
      yd1$Max <- yd1$Max + max(yd1$Max)*0.2
      yd2$Max <- yd2$Max + max(yd2$Max)*0.2
      yd3$Max <- yd3$Max + max(yd2$Max)*0.2
      
      fit1 = aov(PC1 ~ Group, data = plotdata)
      fit2 = aov(PC2 ~ Group, data = plotdata)
      fit3 = aov(PC3 ~ Group, data = plotdata)
      
      fit1_output = group_test(fit = fit1,method = pairwise_method, mapping = meta.phy[[colorby]])
      fit2_output = group_test(fit = fit2,method = pairwise_method, mapping = meta.phy[[colorby]])
      fit3_output = group_test(fit = fit3,method = pairwise_method, mapping = meta.phy[[colorby]])
      
      a=data.frame(groups=fit1_output$model$groups$groups,Gid=rownames(fit1_output$model$groups))
      a$Gid=factor(a$Gid,levels = yd1$Group)
      a=as.data.frame(a[order(a$Gid),])
      
      b=data.frame(groups=fit2_output$model$groups$groups,Gid=rownames(fit2_output$model$groups))
      b$Gid=factor(b$Gid,levels = yd1$Group)
      b=as.data.frame(b[order(b$Gid),])
      
      c=data.frame(groups=fit3_output$model$groups$groups,Gid=rownames(fit3_output$model$groups))
      c$Gid=factor(c$Gid,levels = yd1$Group)
      c=as.data.frame(c[order(c$Gid),])
      
      # test
      
      test <- data.frame(PC1 =a$groups,PC2 = b$groups,PC3 = c$groups,
                         yd1 = yd1$Max,yd2 = yd2$Max,yd3 = yd3$Max,Group = yd1$Group)
      rownames(test)=yd1$Group
      
      test$Group <- factor(test$Group,levels = name_group)
      
      
      plot_pc1 = ggplot(plotdata,aes(Group,PC1)) +
        geom_boxplot(aes(fill = Group),outlier.colour = NA, alpha = box_alpha) +
        scale_fill_manual(values=col_values)+
        scale_color_manual(values = col_values) +
        geom_jitter(aes(color=Group), width = 0.1, size=jitter_size) +
        geom_text(data = test,aes(x = Group,y = yd1,label = PC1),
                  size = letter_size,color = "black") +
        theme_bw()+
        labs(y = plotPC1_ylab, subtitle = "Axis 1", x = plotPC_xlab) + 
        theme(
          axis.ticks.length = unit(0.1,"lines"),
          axis.ticks = element_line(color='black'),
          axis.line = element_line(colour = "black"),
          axis.title.x=element_blank(),
          axis.title.y=element_text(size=12, color="black", margin = margin(r=10)),
          axis.text.y=element_text(colour='black',size=12),
          axis.text.x=element_text(colour='black',size=12, angle = fontAngle, vjust = fontVjust, hjust = fontHjust),
          legend.position = "none")
      
      plot_pc2 = ggplot(plotdata,aes(Group,PC2)) +
        geom_boxplot(aes(fill = Group),outlier.colour = NA, alpha = box_alpha) +
        geom_jitter(aes(color=Group), width = 0.1, size=jitter_size) +
        scale_fill_manual(values=col_values) +
        scale_color_manual(values = col_values) +
        geom_text(data = test,aes(x = Group,y = yd2,label = PC2),
                  size = letter_size,color = "black") +
        labs(y = plotPC2_ylab, subtitle = "Axis 2", x = plotPC_xlab) + 
        theme_bw()+
        theme(
          axis.ticks.length = unit(0.1,"lines"),
          axis.ticks = element_line(color='black'),
          axis.line = element_line(colour = "black"),
          axis.title.x=element_blank(),
          # axis.title.y=element_blank(),
          axis.text.x=element_text(colour='black',size=12, angle = fontAngle, vjust = fontVjust, hjust = fontHjust),
          axis.text.y = element_text(colour='black',size=12),
          legend.position = "none")
    }
    
    g = g / (plot_pc1 + plot_pc2) + plot_layout(guides = "collect", heights = c(2,1))
    
    if(returnDist){
      return(list(dist = phy.dist, bplot=g))
    } else {
      return(g)
    }
  } else {
    
    # no adding PC distance
    
    if(returnDist){
      return(list(dist = phy.dist, bplot=g))
    } else {
      return(g)
    }
    }
  

}

group_test=function(fit,method,mapping){
  
  if (method=='HSD') {
    group_num=as.numeric(table(mapping))
    group_even_check=all(group_num==group_num[1])
    if (group_even_check==T) {
      res<- HSD.test(fit,'Group',unbalanced=F)
      pvalue<- HSD.test(fit,'Group',group = F,unbalanced=F)$comparison
    }else{
      res<- HSD.test(fit,'Group',unbalanced=T)
      pvalue<- HSD.test(fit,'Group',group = F,unbalanced=T)$comparison
      #message('Detected unequal replication, HSD test activated unbalanced mode.')
    }
  }else if(method=='LSD'){
    res<- LSD.test(fit,'Group')
    pvalue<- LSD.test(fit,'Group',group = F)$comparison
  }else if(method=='duncan'){
    res<- duncan.test(fit,'Group')
    pvalue<- duncan.test(fit,'Group',group = F)$comparison
  }else if(method=='scheffe'){
    res<- scheffe.test(fit,'Group')
    pvalue<- scheffe.test(fit,'Group',group = F)$comparison
  }else if(method=='REGW'){
    res<- REGW.test(fit,'Group')
    pvalue<- REGW.test(fit,'Group',group = F)$comparison
  }else if(method=='SNK'){
    res<- SNK.test(fit,'Group')
    pvalue<- SNK.test(fit,'Group',group = F)$comparison
  }else{warning("no method matched in multiple comparisons!")}
  deposit=list()
  deposit$model=res
  deposit$comparison=pvalue
  deposit$letter=data.frame(groups=res$groups$groups,Gid=rownames(res$groups))
  return(deposit)
}


DEanalysis = function(phylo, model, fittype="local", Contrast=NULL, alpha=0.05, abs_log2cutoff=1, Colns = c("log2FoldChange", "padj", "Phylum", "Class", "Genus", "group")){
  # create a significant df
  
  gm_mean = function(x, na.rm=TRUE){
    exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
  }
  
  diagdds = phyloseq_to_deseq2(phylo, model)
  
  geoMeans = apply(counts(diagdds), 1, gm_mean)
  diagdds = estimateSizeFactors(diagdds, geoMeans = geoMeans)
  diagdds = DESeq(diagdds, fitType=fittype)
  print(resultsNames(diagdds))
  res = results(diagdds, contrast = Contrast)
  summary(res)
  res = res[order(res$padj, na.last = NA),]
  
  sigtab = res[(res$padj < alpha),]
  
  if(nrow(sigtab)>0){
    sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(phylo)[rownames(sigtab), ], "matrix"))
    sigtab = sigtab[(abs(sigtab$log2FoldChange)>=abs_log2cutoff), ]
    sigtab$group = paste(Contrast[2], Contrast[3], sep="_vs_")
    sigtab = sigtab[, Colns]
    return(sigtab)
  } else {
    return(data.frame())
  }

}


plotDE = function(df, topn = 10, xlab="log2 Fold Change", ylab = "Taxa (Phylum-Genus)", labfill = "Abundance"){
  
  # df = rbind(...)
  # 
  # df_spread = df %>%
  #   dplyr::mutate(taxa=paste(Phylum, Genus, sep="-")) %>%
  #   dplyr::select(taxa, log2FoldChange, group) %>%
  #   tidyr::spread(key = group, value = log2FoldChange)
  # 
  # df_shared = df_spread[complete.cases(df_spread),]
  # 
  # df_gather = df_shared %>%
  #   tidyr::gather(key=group, value = log2FoldChange, colnames(df_shared)[2:length(colnames(df_shared))])

  g = df %>%
    dplyr::mutate(taxa=paste(Phylum, Genus, sep="-")) %>%
    dplyr::select(taxa, log2FoldChange, group) %>%
    mutate(RAchanges = ifelse(log2FoldChange < 0, "Decreased", "Increased")) %>%
    mutate(RAchanges = factor(RAchanges, levels = c("Increased", "Decreased"))) %>%
    dplyr::select(log2FoldChange, taxa, group, RAchanges) %>%
    group_by(RAchanges) %>%
    slice_max(order_by = abs(log2FoldChange), n=topn) %>%
    ggplot(aes(x=log2FoldChange, y = reorder(taxa, log2FoldChange), fill=RAchanges)) +
    geom_bar(stat = "identity") +
    labs(x=xlab, y=ylab, fill=labfill) +
    theme_classic() +
    scale_fill_manual(values = c("#B03A2E", "#1B4F72")) +
    theme(text=element_text(size=12),
              axis.title.x = element_text(size=12, margin=margin(t=10), color="black"),
              axis.title.y = element_text(size=12, margin=margin(r=10), color="black"),
              axis.text.x = element_text(size=12, color="black"),
              axis.text.y = element_text(size=12, color="black"))
  
  return(g)
}


plotANCOM = function(microM_ancom_out, searchtop=NA, ntop=10, out=T, prefix="ANCOM_", savepath="./", xlab="W", ylab="Taxa", grp1="Group1", grp2="Group2", labfills="Enrich groups", setFactorLevels=F, levels=""){
  
  library(ggpubr)
  library(ggplot2)
  
  
  #output barplot and pieplot
  mrtb = marker_table(microM_ancom_out)
  mrtb.df = data.frame(mrtb)
  
  if(setFactorLevels){
    mrtb.df$enrich_group = factor(mrtb.df$enrich_group, levels=levels)
  }
  
  titlegroup = mrtb.df$enrich_group[order(mrtb.df$enrich_group, decreasing = T)]
  
  if(is.na(unique(titlegroup)[1]) ||  is.na(unique(titlegroup)[2])){
    warning("One of group does not exist")
    if(is.na(unique(titlegroup)[1])){
      print("Group1 is NA")
      grp1 = grp1
      comparison = paste0(unique(titlegroup)[1], "_vs_", grp1)
    } else if(is.na(unique(titlegroup)[2])) {
      print("Group2 is NA")
      grp1 = grp2
      comparison = paste0(grp2, "_vs_", unique(titlegroup)[1])
    } else {
      comparison = paste0(grp2, "_vs_", grp1)
    }
    
    
  } else {
    comparison = paste0(unique(titlegroup)[2], "_vs_", unique(titlegroup)[1])
  }
  
  
  if(!is.na(searchtop)){
    if(length(searchtop) < 0){
      warning("Please input search top vector")
    }
    
    mrtb.df$predom = ifelse(mrtb.df$feature %in% searchtop, paste0("*", mrtb.df$feature), mrtb.df$feature)
    mrtb.df$Ifpredom = ifelse(mrtb.df$feature %in% searchtop, "Y", "N")
    # mrtb.df$fontface = ifelse(mrtb.df$Ifpredom == "Y", "bold", "plain")
  }

  if(out){
    filename = paste0(prefix, "_", comparison, ".csv") 
    write.csv(mrtb.df, paste(savepath, filename, sep="/"), row.names = F, quote = F)
  }

  if("predom" %in% colnames(mrtb.df)){
    bar = mrtb.df %>%
      group_by(enrich_group) %>%
      slice_max(order_by = abs(ef_W), n = ntop) %>%
      ggplot(aes(x=ef_W, y=reorder(predom, ef_W), fill=enrich_group)) +
      geom_bar(stat="identity") +
      scale_fill_manual(values = c("#2874A6", "#EC7063")) +
      theme_classic2() +
      labs(title = comparison, x=xlab, y=ylab, fill=labfills) +
      theme(text=element_text(size=12),
            plot.title = element_text(hjust = 0.5),
            axis.title.x = element_text(size=12, margin=margin(t=10), color="black"),
            axis.title.y = element_text(size=12, margin=margin(r=10), color="black"),
            axis.text.x = element_text(size=12, color="black"),
            axis.text.y = element_text(size=12, color="black"))
  } else {
    bar = mrtb.df %>%
      group_by(enrich_group) %>%
      slice_max(order_by = abs(ef_W), n = ntop) %>%
      ggplot(aes(x=ef_W, y=reorder(feature, ef_W), fill=enrich_group)) +
      geom_bar(stat="identity") +
      scale_fill_manual(values = c("#2874A6", "#EC7063")) +
      theme_classic2() +
      labs(title = comparison, x=xlab, y=ylab, fill=labfills) +
      theme(text=element_text(size=12),
            plot.title = element_text(hjust = 0.5),
            axis.title.x = element_text(size=12, margin=margin(t=10), color="black"),
            axis.title.y = element_text(size=12, margin=margin(r=10), color="black"),
            axis.text.x = element_text(size=12, color="black"),
            axis.text.y = element_text(size=12, color="black"))
  }

  return(list(barp = bar, sig_features = mrtb.df$feature, df=mrtb.df))
}



plotPIE = function(phylo_df, mrtb.df, pie_cutoff=0.01, provideColor = NA, defaultColor="Paired", grp1="Group1", grp2="Group2"){
  
  
  library(RColorBrewer)
  
  getMeta = phylo_df %>%
    dplyr::select(Phylum, Genus) %>%
    dplyr::distinct() %>%
    filter(Genus %in% mrtb.df$feature) %>%
    dplyr::select(Phylum) %>%
    dplyr::group_by(Phylum) %>%
    dplyr::count(sort=T) %>%
    ungroup() %>%
    mutate(perc = n/sum(n)) %>%
    mutate(Phylum2 = ifelse(perc<pie_cutoff, paste0("<", 100*pie_cutoff, "%"), Phylum)) %>%
    group_by(Phylum2) %>%
    summarise(n2=sum(n)) %>%
    mutate(total = sum(n2), perc = n2/total, labels=scales::percent(perc, accuracy = 1)) 
  
  
  if(is.na(provideColor)){

    phylumColors = unique(getMeta$Phylum2)
    phylumColors2 = phylumColors[phylumColors != paste0("<", 100*pie_cutoff, "%")]
    phylumColors2 = c(phylumColors2, paste0("<", 100*pie_cutoff, "%"))
    cols = colorRampPalette(brewer.pal(9, defaultColor))(length(phylumColors2))
    names(cols) = phylumColors2

  } else {
    
    # select used colors
    cols = provideColor[names(provideColor) %in% unique(getMeta$Phylum2)]
  }

  # deal with plot title
  titlegroup = mrtb.df$enrich_group[order(mrtb.df$enrich_group, decreasing = T)]
  
  if(is.na(unique(titlegroup)[1]) ||  is.na(unique(titlegroup)[2])){
    warning("One of group does not exist")
    if(is.na(unique(titlegroup)[1])){
      print("Group1 is NA")
      grp1 = grp1
      comparison = paste0(unique(titlegroup)[1], "_vs_", grp1)
    } else if(is.na(unique(titlegroup)[2])) {
      print("Group2 is NA")
      grp1 = grp2
      comparison = paste0(grp2, "_vs_", unique(titlegroup)[1])
    } else {
      comparison = paste0(grp2, "_vs_", grp1)
    }
    
    
  } else {
    comparison = paste0(unique(titlegroup)[2], "_vs_", unique(titlegroup)[1])
  }
  
  
  pie = ggpie(getMeta, x="n2", label = "labels", fill="Phylum2", palette = cols, color = "white", legend="right") + labs(fill="Phylum", title = comparison) + theme(
    plot.title = element_text(hjust = 0.5)
  )
    

  return(list(df=getMeta, pie=pie))
}



plotDE2 = function(..., topn = 10, xlab="log2 Fold Change", ylab = "Taxa (Phylum-Genus)", labfill = "Abundance"){
  
  df = rbind(...)

  df_spread = df %>%
    dplyr::mutate(taxa=paste(Phylum, Genus, sep="-")) %>%
    dplyr::select(taxa, log2FoldChange, group) %>%
    tidyr::spread(key = group, value = log2FoldChange)

  df_shared = df_spread[complete.cases(df_spread),]

  df_gather = df_shared %>%
    tidyr::gather(key=group, value = log2FoldChange, colnames(df_shared)[2:length(colnames(df_shared))])
  
  g = df_gather %>%
    dplyr::select(taxa, log2FoldChange, group) %>%
    mutate(RAchanges = ifelse(log2FoldChange < 0, "Decreased", "Increased")) %>%
    mutate(RAchanges = factor(RAchanges, levels = c("Increased", "Decreased"))) %>%
    dplyr::select(log2FoldChange, taxa, group, RAchanges) %>%
    group_by(RAchanges) %>%
    slice_max(order_by = abs(log2FoldChange), n=topn) %>%
    ggplot(aes(x=log2FoldChange, y = reorder(taxa, log2FoldChange), fill=RAchanges)) +
    geom_bar(stat = "identity") +
    facet_grid(~ group) +
    labs(x=xlab, y=ylab, fill=labfill) +
    theme_classic() +
    scale_fill_manual(values = c("#B03A2E", "#1B4F72")) +
    theme(text=element_text(size=12),
          panel.grid.major.y = element_line(),
          axis.title.x = element_text(size=12, margin=margin(t=10), color="black"),
          axis.title.y = element_text(size=12, margin=margin(r=10), color="black"),
          axis.text.x = element_text(size=12, color="black"),
          axis.text.y = element_text(size=12, color="black"))
  
  return(g)
}

plotCoef_mcm = function(result_diffTest, phylo, fdr = 0.05, taxLevel="Genus", nlevel = 3,orderList = NULL){
  
  # specific for mcm
  all_fdr_values = result_diffTest$p_fdr[!is.na(result_diffTest$p_fdr)]
  sigTaxa = corncob::otu_to_taxonomy(result_diffTest$significant_taxa, data=phylo, level=taxLevel)
  
  try_species = names(all_fdr_values)
  taxaNames = corncob::otu_to_taxonomy(try_species, data = phylo, level = taxLevel) # OTU name - genus 
  
  all_models = result_diffTest[["all_models"]][!is.na(result_diffTest[["all_models"]])]
  names(all_models) = taxaNames


  if( ! is.null(orderList) ){
    
    if((length(orderList) == length(all_models)) | (length(orderList) == length(taxaNames))){
      all_models = all_models[as.character(orderList)]
      taxaNames = taxaNames[order(match(taxaNames, orderList))]
      all_fdr_values = all_fdr_values[names(taxaNames)]
    } else {
      stop("orderList length is different to the model or taxaNames")
    }
    
  }
  
  
  var_levels = c()
  
  coef = c()
  maxse = c()
  minse = c()
  
  vari = c()
  vari_maxse = c()
  vari_minse = c()
  taxaNames_rep = c()
  fdr_rep = c()
  sigTaxa_w_symbol = c()
  
  
  for(i in 1:length(all_models)){
    
    print(i)
    print(paste("Order list:", orderList[i]))
    print(paste("Model name:", names(all_models[i])))
    
    if(taxaNames[i] %in% sigTaxa){
      sigTaxa_w_symbol = c(sigTaxa_w_symbol, rep("*", nlevel))
    } else {
      sigTaxa_w_symbol = c(sigTaxa_w_symbol, rep("ns", nlevel))
    }
    
    taxaNames_rep = c(taxaNames_rep, rep(taxaNames[i], nlevel))
    fdr_rep = c(fdr_rep, rep(all_fdr_values[i], nlevel))
    
    # how many coefficients
    n.mu = length(grep("mu", rownames(all_models[[i]]$coefficients)))
    n.phi = length(grep("phi", rownames(all_models[[i]]$coefficients)))

    
    if(n.mu > 1 & n.phi > 1){
      
      print("Mu and phi both tested")
      
      var_levels = c(var_levels, gsub("mu.", "", names(all_models[[i]]$coefficients[2:n.mu,1])))
      
      # coefficients
      coef = c(coef, all_models[[i]]$coefficients[2:n.mu,1]) 
      maxse = c(maxse, all_models[[i]]$coefficients[2:n.mu,1]+qnorm(0.975)*all_models[[i]]$coefficients[2:n.mu,2])
      minse = c(minse, all_models[[i]]$coefficients[2:n.mu,1]-qnorm(0.975)*all_models[[i]]$coefficients[2:n.mu,2])
      # coefficients pvalu
      
      # variability
      vari = c(vari, all_models[[i]]$coefficients[6:8,1])
      vari_maxse = c(vari_maxse, all_models[[i]]$coefficients[6:8,1]+qnorm(0.975)*all_models[[i]]$coefficients[6:8,2])
      vari_minse = c(vari_minse, all_models[[i]]$coefficients[6:8,1]-qnorm(0.975)*all_models[[i]]$coefficients[6:8,2])
      
      df = data.frame(coef = coef, facet_var = var_levels, variability = vari, xmin = minse, xmax = maxse, var_xmin = vari_minse, var_xmax = vari_maxse, taxa = taxaNames_rep, siglevels=sigTaxa_w_symbol, fdr = fdr_rep)

      
    } else if (n.mu > 1 & n.phi == 1){
      
      
      print("Only Mu tested")
      var_levels = c(var_levels, gsub("mu.", "", names(all_models[[i]]$coefficients[2:n.mu,1])))
      
      # coefficients
      coef = c(coef, all_models[[i]]$coefficients[2:n.mu,1]) 
      maxse = c(maxse, all_models[[i]]$coefficients[2:n.mu,1]+qnorm(0.975)*all_models[[i]]$coefficients[2:n.mu,2])
      minse = c(minse, all_models[[i]]$coefficients[2:n.mu,1]-qnorm(0.975)*all_models[[i]]$coefficients[2:n.mu,2])
      # coefficients pvalu
      
      # variability
      # vari = c(vari, all_models[[i]]$coefficients[6:8,1])
      # vari_maxse = c(vari_maxse, all_models[[i]]$coefficients[6:8,1]+qnorm(0.975)*all_models[[i]]$coefficients[6:8,2])
      # vari_minse = c(vari_minse, all_models[[i]]$coefficients[6:8,1]-qnorm(0.975)*all_models[[i]]$coefficients[6:8,2])
      
      df = data.frame(coef = coef, facet_var = var_levels, xmin = minse, xmax = maxse, taxa = taxaNames_rep, siglevels=sigTaxa_w_symbol, fdr = fdr_rep)
      

    } else if(n.mu == 1) {

      stop("No test")
    }
    }
    
  return(list(DF=df, mod=all_models))
}



plot_BCAsp = function(sp_df, speciesname, formulas, maxIT=25, SN="ES", padjust="tukey", rds_save = F, rds_path = "", posi_dodge=0.9, alpha=0.9, err_dodge_width=0.2, addColor="none", ncolor=4, fontcolor="black", fontsize=12, labfill="Fungicides", xlab="Sampling time (day)", ylab="Square root of abundance", setRange = NULL, sameYmax=F, nudgeY=0){
  

  cat("\n-------------------------------------\n")
  cat(SN, "\n")
  cat(speciesname)
  cat("\n-------------------------------------\n")
  # this function is hard coded
  select_sp = sp_df %>% filter(Species == speciesname & Season == SN)
  
  # 0 distribution
  print(scales::percent(sum(select_sp$Abundance == 0)/length(select_sp$Abundance)), digits = 4)
  
  # negative binomial
  model = glm.nb(formula(formulas), data=select_sp, control = glm.control(maxit = maxIT))
  
  # comparision within each time
  em = emmeans(model, pairwise ~ funlab | time_label, p.adjust= padjust)
  # cld.em = cld(em, Letters = letters, alpha=0.05, sort=T)
  cld.em = cld(em, Letters = letters)
  cld_df = as.data.frame(cld.em)

  # cld_df$group = gsub("3", "c", gsub("2", "b", gsub("1", "a", stringr::str_trim(cld_df$.group))))
  # 
  cld_df$joint = paste(cld_df$time_label, cld_df$funlab, sep = "_")
  select_sp$joint = paste(select_sp$time_label, select_sp$funlab, sep = "_")
  
  (select_sp.agg = select_sp %>%
      group_by(joint, time_label, funlab) %>%
      mutate(sqrtAbund = sqrt(Abundance+1)) %>%
      summarise(avg = mean(sqrtAbund), se = plotrix::std.error(sqrtAbund)))
  
  plotsp.df.ltrs = merge(select_sp.agg, cld_df, by="joint", no.dups = T)
  
  if(is.null(setRange)){
    rangeMax = sqrt(max(sp_df %>% filter(Species == speciesname) %>% pull(Abundance)))
  } else {
    rangeMax =setRange
  }
  
  
  if(addColor=="none"){
    colorpal =ggsci::pal_jco()(ncolor)
  } else {
    colorpal = addColor
  }
  
  if(sameYmax){
    plotsp_plot = plotsp.df.ltrs %>%
      ggplot(aes(x = time_label.x, y=avg, fill=funlab.x)) +
      geom_bar(stat="identity", position = position_dodge(posi_dodge), alpha=alpha, color="black") +
      geom_text(aes(label=.group, y=rangeMax), position = position_dodge(posi_dodge)) +
      geom_errorbar(aes(ymin=avg-se, ymax = avg+se), width=err_dodge_width, position = position_dodge(posi_dodge)) +
      # scale_y_continuous(limits = c(0, rangeMax)) +
      scale_fill_manual(values = colorpal) +
      theme_classic() +
      labs(fill = labfill, x = xlab, y= ylab, title=speciesname) +
      theme(text=element_text(size=12),
            axis.title.x = element_text(size=fontsize, margin=margin(t=10), color=fontcolor),
            axis.title.y = element_text(size=fontsize, margin=margin(r=10), color=fontcolor),
            axis.text.x =  element_text(size=fontsize, color=fontcolor),
            axis.text.y =  element_text(size=fontsize, color=fontcolor))
  } else {
    plotsp_plot = plotsp.df.ltrs %>%
      ggplot(aes(x = time_label.x, y=avg, fill=funlab.x)) +
      geom_bar(stat="identity", position = position_dodge(posi_dodge), alpha=alpha, color="black") +
      geom_text(aes(label=.group, y=avg+se+nudgeY), position = position_dodge(posi_dodge)) +
      geom_errorbar(aes(ymin=avg-se, ymax = avg+se), width=err_dodge_width, position = position_dodge(posi_dodge)) +
      # scale_y_continuous(limits = c(0, rangeMax)) +
      scale_fill_manual(values = colorpal) +
      theme_classic() +
      labs(fill = labfill, x = xlab, y= ylab, title=speciesname) +
      theme(text=element_text(size=12),
            axis.title.x = element_text(size=fontsize, margin=margin(t=10), color=fontcolor),
            axis.title.y = element_text(size=fontsize, margin=margin(r=10), color=fontcolor),
            axis.text.x =  element_text(size=fontsize, color=fontcolor),
            axis.text.y =  element_text(size=fontsize, color=fontcolor))
  }

  
  if(rds_save){
    if(rds_path == ""){
      rds_path = paste0("./BCA/BCA_plots/", speciesname, "_", SN, ".rds")
    }
    saveRDS(plotsp_plot, rds_path)
  }
  
  return(list(plot=plotsp_plot, summary_model=summary(model), statistics_anova=car::Anova(model)))

}


plotTukey = function(obj, diversity_data, plotx="funlab", ploty="value", facet.by="~ full_season", mod = "pairwise ~ funlab | full_season", padj = "tukey", labfill="", xlab="Fungicides", ylab="", title="", fontcolor="black", fontsize=12, labelsize=6){
  
  mod0 = as.formula(mod)
  tukeyout = emmeans::emmeans(obj, mod0, p.adjust.methods = padj)
  
  tukeydf = as.data.frame(multcomp::cld(tukeyout, Letters = letters))
  
  ym = min(diversity_data$value)
  
  g = ggplot(diversity_data, aes(x= !!sym(plotx), y = !!sym(ploty))) +
    geom_boxplot() +
    facet_grid(facet.by) +
    geom_text(data= tukeydf, aes(y=ym, label=.group), size=labelsize) +
    theme_bw() +
    labs(fill = labfill, x = xlab, y= ylab, title=title) +
    theme(text=element_text(size=12),
          axis.title.x = element_text(size=fontsize, margin=margin(t=10), color=fontcolor),
          axis.title.y = element_text(size=fontsize, margin=margin(r=10), color=fontcolor),
          axis.text.x =  element_text(size=fontsize, color=fontcolor),
          axis.text.y =  element_text(size=fontsize, color=fontcolor))

  
  return(g)
  
}


plot_lefse = function(run_lefse_out, phy=NULL, topn=10, colorSet = NULL,relevel=F, new_levels=NULL, padj_cutoff=0.05, xlab="LDA score(log10)", ylab="Genus", labfill = "Enriched group", plottitle= "", fontcolor="black", fontsize=12, bartrans= 0.9, label_print=0){
  library(dplyr)
  library(ggplot2)
  library(RColorBrewer)
  library(ggtext)
  

  
  df = data.frame(marker_table(run_lefse_out))
  
  df$feature = scientific_name_formatter(df$feature)
  
  df$enrich_group = factor(df$enrich_group)

  if(relevel){
    if(is.null(new_levels)){
      stop("Levels need to be specific")
    } else {
      levels(df$enrich_group) = new_levels
    }
  }
  
  if(is.null(colorSet)){
    colorInput = colorRampPalette(brewer.pal(8, "Dark2"))(length(levels(df$enrich_group)))
  } else {
    colorInput = colorSet
  }
  
  df_order = df%>%
    arrange(.data$enrich_group, .data$ef_lda)
  
  feat = df_order$feature
  df_order$feature = factor(feat, levels = feat)
  
  df_top = df_order %>%
    group_by(enrich_group) %>%
    filter(padj < padj_cutoff) %>%
    slice_max(order_by = ef_lda, n = topn)
  
  if(label_print == 1){
    if(is.null(phy)){
      stop("Need phyloseq object")
    }
    tax_m = as.data.frame(tax_table(phy))
    taxa  = tax_m %>% 
      filter(Genus %in% as.character(df_top$feature)) %>%
      mutate(Phylum = scientific_name_formatter(Phylum), Genus = scientific_name_formatter(Genus)) %>%
      mutate(new_feature = paste(Phylum, Genus, sep = "-")) %>%
      select(new_feature, Genus) %>%
      distinct()
    
    rownames(taxa) = taxa$Genus
    
    df_top_merge = merge(df_top, taxa, by.x = "feature", by.y = "Genus")
    
    df_top_merge = df_top_merge %>%
      arrange(.data$enrich_group, .data$ef_lda)
    
    feat = df_top_merge$new_feature
    df_top_merge$new_feature = factor(feat, levels = feat)
    
    print(df_top_merge)
    
    g = df_top_merge %>%
      group_by(enrich_group) %>%
      mutate(enrich_group_num = as.numeric(enrich_group)) %>%
      ggplot(aes(x=ef_lda, y=reorder(new_feature, -enrich_group_num), fill=enrich_group)) +
      geom_col(alpha=bartrans) +
      scale_x_continuous(expand = c(0,0)) +
      scale_y_discrete(expand = c(0.01, 0.01))+ 
      scale_fill_manual(values = colorInput) + 
      labs(fill=labfill, x=xlab, y=ylab, title = plottitle) +
      theme(
        text=element_text(size=12),
        axis.title.x = element_text(size=fontsize, margin=margin(t=10), color=fontcolor),
        axis.title.y = element_text(size=fontsize, margin=margin(r=10), color=fontcolor),
        axis.text.x =  element_text(size=fontsize, color=fontcolor),
        axis.text.y =  element_markdown(size=fontsize, color=fontcolor),
        axis.ticks.y = element_blank(),
        panel.background = element_blank(),
        panel.grid.major.x = element_line(colour = "black", linetype = 8)) +
      guides(guide_legend())
    

     
  } else {

    
    g = df_top %>%
      mutate(enrich_group_num = as.numeric(enrich_group)) %>%
      ggplot(aes(x=ef_lda, y=reorder(feature, -enrich_group_num), fill=enrich_group)) +
      geom_col(alpha=bartrans) +
      scale_x_continuous(expand = c(0,0)) +
      scale_y_discrete(expand = c(0.01, 0.01))+ 
      scale_fill_manual(values = colorInput) + 
      labs(fill=labfill, x=xlab, y=ylab, title = plottitle) +
      theme(
        text=element_text(size=12),
        axis.title.x = element_text(size=fontsize, margin=margin(t=10), color=fontcolor),
        axis.title.y = element_text(size=fontsize, margin=margin(r=10), color=fontcolor),
        axis.text.x =  element_text(size=fontsize, color=fontcolor),
        axis.text.y =  element_markdown(size=fontsize, color=fontcolor),
        axis.ticks.y = element_blank(),
        panel.background = element_blank(),
        panel.grid.major.x = element_line(colour = "black", linetype = 8)) +
      guides(guide_legend())
  }


  


    
  return(list(graph=g, df=df, colorused = colorInput))
  
}

extractCornCob = function(moddf, outputfile=T, outname=""){
  
  library(dplyr)
  
  startdf = data.frame(matrix(ncol = 7, nrow = 0))
  colnames(startdf) = c("taxa","mu.time_label1", "mu.time_label7", "mu.time_label14", "mu.funlabDL", "mu.funlabBM", "mu.funlabC2")
  
  spnames = names(moddf$mod)
  
  print(spnames)
  
  for(i in spnames){
    print(i)
    
    spdf = as.data.frame(moddf$mod[[i]]$coefficients)[, c("Estimate", "Pr(>|t|)")]
    
    colnames(spdf) = c("Est", "Pval")
    spdf2 = spdf %>%
      mutate(OUTtest = ifelse(Est > 0, "+", ifelse(Est < 0, "-", "0"))) %>%
      dplyr::select(OUTtest, Pval)
    
    spdf2$Pval = round(spdf2$Pval,4)
    tsp_df = t(spdf2)[, c("mu.time_label1", "mu.time_label7", "mu.time_label14", "mu.funlabDL", "mu.funlabBM", "mu.funlabC2")]
    
    tsp_df2 = as.data.frame(tsp_df)
    tsp_df2$taxa = i
    
    print(head(tsp_df2))
    tsp_df2 = tsp_df2[,c("taxa","mu.time_label1", "mu.time_label7", "mu.time_label14", "mu.funlabDL", "mu.funlabBM", "mu.funlabC2")]
    
    startdf = rbind(startdf, tsp_df2)
  }
  
  if(outputfile){
    write.csv(startdf, outname, quote = F, row.names = F)
  }
  
  return(startdf)
}

plotLinda = function(lindaplotobj, ylab = "Taxa", xlab = "Log2 Fold Changes", titles=NULL, pointsize=3, pointalpha=0.8, slice_max_n = 50, dotcolor="blue4", dotlabel="Corrected"){
  library(ggplot2)
  library(ggtext)
  
  figs = list()
  
  for(i in 1:length(lindaplotobj$plot.lfc)){
    
    if( (!is.null(titles)) & (length(titles) > 1)){
        Title = titles[i]
      } else {
        Title = paste0("fig",i)
      }
      
    print(Title)
    
    theData = lindaplotobj$plot.lfc[[i]]$data
    
    fig_n = theData %>%
      mutate(OTU = gsub("(^OTU[0-9]+?:).*", "\\1", Taxa)) %>%
      mutate(only_taxa = gsub("(^OTU[0-9]+?:)(.*)", "\\2", Taxa)) %>%
      mutate(italic_taxa = ifelse(grepl(" sp.", only_taxa),  
                              paste(paste0("*", gsub(" sp.$", "", only_taxa), "*"), " sp."), paste0("*", only_taxa, "*"))) %>%
      mutate(new_taxa = paste0(OTU, italic_taxa)) %>%
      dplyr::filter(bias == "Debiased") %>%
      mutate(sign=ifelse(Log2FoldChange > 0, "pos", "neg")) %>%
      group_by(sign) %>%
      slice_max(order_by = abs(Log2FoldChange), n=slice_max_n) %>%
      ggplot(aes(x = Log2FoldChange, y=reorder(new_taxa, Log2FoldChange))) +
      geom_errorbar(aes(xmin = Log2FoldChange - 1.96 * lfcSE,
                        xmax = Log2FoldChange + 1.96 * lfcSE), width = .2) +
      geom_point(aes(color = bias), size = pointsize, alpha=pointalpha) +
      geom_vline(xintercept = 0, color = 'gray', linetype = 'dashed') +
      geom_vline(xintercept = 2, color = 'gray', linetype = 'dashed') +
      geom_vline(xintercept = -2, color = 'gray', linetype = 'dashed') +
      scale_color_manual(values = dotcolor, labels=dotlabel) +
      labs(title = Title, y=ylab, x=xlab, color="Bias correction") +
      theme_classic() +
      theme(
            legend.position = "none",
            text = element_text(size=12),
            panel.grid.major.y = element_line(),
            axis.title.x = element_text(margin = margin(t=10), color="black", size=12),
            axis.title.y = element_text(margin = margin(r=10), color="black", size=12),
            axis.text.x = element_text(color="black", size=12),
            axis.text.y = element_markdown(color="black", size=12))
     
     figs[[Title]] = fig_n
  }
  
  return(figs)
}

summLinda = function(df, pv=0.05, log2=1, slicemax=5){
  
  df$taxa = gsub("OTU[0-9]+?:", "", rownames(df))
  
  df2 = df %>%
    filter(padj < pv & abs(log2FoldChange) >= log2) %>%
    mutate(signs = ifelse(log2FoldChange > 0, "+", "-")) %>%
    group_by(signs,taxa) %>%
    dplyr::summarise(howmany = n()) %>%
    arrange(desc(howmany), .by_group = T) %>%
    slice_max(order_by = howmany, n=slicemax)
  
  return(df2)
}


ci <- function(dataframe, summary_var, group_var = NULL){

  summary_df <- dataframe %>%
    group_by(!!! syms(group_var)) %>%
    summarise(ci = list(mean_cl_normal(!! syms(summary_var)) %>% 
                          rename(mean = y, lower_ci = ymin, upper_ci = ymax))) %>%
    tidyr::unnest()
}


dfcorn_reformat = function(plotCoef_mcm_data, write_out = F, pathname1="corncob_fun_table.tsv", pathname2="corncob_time_table.tsv"){
  # hardcoded
  # reformat
  tax_labs = c()
  DL_estimates = c()
  BM_estimates = c()
  C2_estimates = c()
  DL_pval = c()
  BM_pval = c()
  C2_pval = c()
  
  D1_estimates = c()
  D7_estimates = c()
  D14_estimates = c()
  D1_pval = c()
  D7_pval = c()
  D14_pval = c()
  
  for(i in names(plotCoef_mcm_data$mod)){
    
    print(paste0("Extract for ", i))
    
    esti_table = plotCoef_mcm_data$mod[[i]]

    if(is.null(esti_table$coefficients)){
      print("Contained null 'coefficients")
      all_estimate = summary(esti_table)$coef[,1]
      all_pvalue = summary(esti_table)$coef[,4]
    } else {
      all_estimate = esti_table$coefficients[,1]
      all_pvalue = esti_table$coefficients[,4]
    }
    
    
    e_DL=all_estimate[2]
    e_BM=all_estimate[3]
    e_C2=all_estimate[4]
    
    e_D1=all_estimate[5]
    e_D7=all_estimate[6]
    e_D14=all_estimate[7]
    
    
    # pvalue
    
    p_DL=all_pvalue[2]
    p_BM=all_pvalue[3]
    p_C2=all_pvalue[4]
    
    p_D1=all_pvalue[5]
    p_D7=all_pvalue[6]
    p_D14=all_pvalue[7]
    
    
    tax_labs = c(tax_labs, i)
    DL_estimates = c(DL_estimates, e_DL)
    DL_pval = c(DL_pval, p_DL)
    
    BM_estimates = c(BM_estimates, e_BM)
    BM_pval = c(BM_pval, p_BM)
    
    C2_estimates = c(C2_estimates, e_C2)
    C2_pval = c(C2_pval, p_C2)
    
    # sampling time
    D1_estimates = c(D1_estimates, e_D1)
    D1_pval      = c(D1_pval, p_D1)
    
    D7_estimates = c(D7_estimates, e_D7)
    D7_pval      = c(D7_pval, p_D7)
    
    D14_estimates = c(D14_estimates, e_D14)
    D14_pval      = c(D14_pval, p_D14)
    
  }
  
  fun_df = data.frame(taxa = tax_labs,
                      BM_est = BM_estimates,
                      DL_est = DL_estimates,
                      C2_est = C2_estimates,
                      BM_pvalues = BM_pval,
                      DL_pvalues = DL_pval,
                      C2_pvalues = C2_pval)
  
  fundf2 = fun_df %>%
    mutate(
           BM_sign = ifelse(BM_pvalues < 0.001, "***", ifelse( BM_pvalues >= 0.001 & BM_pvalues < 0.01, "**", ifelse(BM_pvalues >= 0.01 & BM_pvalues < 0.05, "*", ".") )),
           DL_sign = ifelse(DL_pvalues < 0.001, "***", ifelse( DL_pvalues >= 0.001 & DL_pvalues < 0.01, "**", ifelse(DL_pvalues >= 0.01 & DL_pvalues < 0.05, "*", "."))),
           C2_sign = ifelse(C2_pvalues < 0.001, "***", ifelse( C2_pvalues >= 0.001 & C2_pvalues < 0.01, "**", ifelse(C2_pvalues >= 0.01 & C2_pvalues < 0.05, "*", ".") )))
  

    
  
  
  Time_df = data.frame(taxa = tax_labs,
                                  D1_est = D1_estimates,
                                  D7_est = D7_estimates,
                                  D14_est = D14_estimates,
                                  D1_pvalues = D1_pval,
                                  D7_pvalues = D7_pval,
                                  D14_pvalues = D14_pval)
  
  Timedf2 = Time_df %>%
    mutate(D1_sign = ifelse(D1_pvalues < 0.001, "***", ifelse( D1_pvalues >= 0.001 & D1_pvalues < 0.01, "**", ifelse(D1_pvalues >= 0.01 & D1_pvalues < 0.05, "*", "."))),
           D7_sign = ifelse(D7_pvalues < 0.001, "***", ifelse( D7_pvalues >= 0.001 & D7_pvalues < 0.01, "**", ifelse(D7_pvalues >= 0.01 & D7_pvalues < 0.05, "*", ".") )),
           D14_sign = ifelse(D14_pvalues < 0.001, "***", ifelse( D14_pvalues >= 0.001 & D14_pvalues < 0.01, "**", ifelse(D14_pvalues >= 0.01 & D14_pvalues < 0.05, "*", ".") )))

  
  if(write_out){

      
      write.table(fundf2, file = pathname1, sep = "\t", quote = F, row.names = F)
      write.table(Timedf2, file = pathname2, sep = "\t", quote = F, row.names = F)

  }
  
  return(list(outFUN = fundf2, outTime = Timedf2))
  
}

reshape_ancom = function(ancombcout, sel_var="funlab2"){
  
  df = ancombcout$res
  coln = colnames(df)
  
  sel_col = c("taxon", coln[grepl(sel_var, coln)])
  col_logfc = sel_col[grepl("lfc_", sel_col)]
  col_df_se = sel_col[grepl("se_", sel_col)]
  col_W_value = sel_col[grepl("W_", sel_col)]
  col_p_value = sel_col[grepl("p_", sel_col)]
  col_q_value = sel_col[grepl("q_", sel_col)]
  col_diff_abund = sel_col[grepl("diff_", sel_col)]
 
  new_col = c("taxon",
              col_logfc,
              col_df_se,
              col_W_value,
              col_p_value,
              col_q_value,
              col_diff_abund)
  

  library(stringi)
  
  sel_df = df[, new_col]
  
  reshapedf = sel_df %>%
    pivot_longer(
      -taxon,
      names_to = c("type_var", "target_var"),
      names_pattern = "(lfc_|se_|W_|p_|q_|diff_)(.*)",
      values_to = "measures"
    ) %>%
    mutate(target_var = gsub(sel_var, "", target_var), 
           type_var = stri_replace_all_regex(type_var,
                                             pattern = c("lfc_", "se_", "W_", "p_", "q_", "diff_"),
                                             replacement = c("logFC", "se", "W", "pval", "qval", "diff"),
                                             vectorize_all = F)) %>%
    spread(type_var, measures)
  
  return(reshapedf)
}


runALDEX2 = function(phylo, the_model, glomLevel="Genus", mcSamples = 128, denom = "all", verbose=T, ci=F, useMC.effect = F){
  
  library(dplyr)
  library(ALDEx2)

  phylo_glom = phylo %>% speedyseq::tax_glom(glomLevel) 
  tax_table(phylo_glom) = tax_table(phylo_glom)[, 1:7]
  colnames(tax_table(phylo_glom)) = c("Kingdom","Phylum","Class","Order","Family","Genus","Species")
  taxa_names(phylo_glom) = as.character(tax_table(phylo_glom)[, glomLevel])
  
  # metadf = data.frame(sample_data(phylo_glom))
  # run aldex2
  print("Run aldex2.clr")
  clr_out = aldex.clr(otu_table(phylo_glom), the_model, mc.samples = mcSamples, denom = denom, verbose=verbose)
  # run aldex test on glm
  print("Run aldex2 test using GLM")
  clr.glm.test = aldex.glm(clr_out, the_model)
  # run aldex effect
  print("Run aldex2 effect size on the GLM model")
  clr.glm.effect = aldex.glm.effect(clr_out, CI=ci, useMC = useMC.effect)
  
  # combine
  clr.glm.comb = data.frame(clr.glm.test, clr.glm.effect)
  return(list(aldexTest_out = clr.glm.test, aldexEffect_out = clr.glm.effect, combined = clr.glm.comb))
}

reshape_aldex2 = function(aldex_out, whichlist = "combined", sel_var="funlab2"){
  
  library(dplyr)
  
  cat("ALDEX2 output explanation: \n\n")
  
  cat("∗ we.ep - Expected P value of Welch’s t test",
  "∗ we.eBH - Expected Benjamini-Hochberg corrected P value of Welch’s t test",
  "∗ wi.ep - Expected P value of Wilcoxon rank test",
  "∗ wi.eBH - Expected Benjamini-Hochberg corrected P value of Wilcoxon test",
  "∘ kw.ep - Expected P value of Kruskal-Wallace test",
  "∘ kw.eBH - Expected Benjamini-Hochberg corrected P value of Kruskal-Wallace test",
  "∘ glm.ep - Expected P value of glm test",
  "∘ glm.eBH - Expected Benjamini-Hochberg corrected P value of glm test",
  "♢ rab.all - median clr value for all samples in the feature",
  "♢ rab.win.NS - median clr value for the NS group of samples",
  "♢ rab.win.S - median clr value for the S group of samples",
  "♢ rab.X1_BNS.q50 - median expression value of features in sample X1_BNS if include.item.summary=TRUE",
  "♢ dif.btw - median difference in clr values between S and NS groups",
  "♢ dif.win - median of the largest difference in clr values within S and NS groups",
  "♢ effect - median effect size: diff.btw / max(diff.win) for all instances",
  "♢ overlap - proportion of effect size that overlaps 0 (i.e. no effect)", sep = "\n")
  
  if(whichlist == "combined"){
    df = aldex_out$combined %>%
      tibble::rownames_to_column("taxon")
  } 
  
  
  
  coln = colnames(df)
  sel_col = c("taxon", coln[grepl(sel_var, coln)])
  
  col_est = sel_col[grepl("\\.Est", sel_col)]
  col_p_value = sel_col[grepl("\\.pval$", sel_col)]
  col_q_value = sel_col[grepl("\\.pval\\.holm$", sel_col)]
  col_effect = sel_col[grepl("\\.effect$", sel_col)]
  col_effect_95CI_L = sel_col[grepl("\\.effect\\.low$", sel_col)]
  col_effect_95CI_U = sel_col[grepl("\\.effect\\.high$", sel_col)]
  
  new_col = c("taxon", 
              col_est,
              col_p_value,
              col_q_value,
              col_effect,
              col_effect_95CI_L,
              col_effect_95CI_U)
  
  library(stringi)

  sel_df = df[, new_col]
  
  reshapedf = sel_df %>%
    tidyr::pivot_longer(
      -taxon,
      names_to = c("target_var", "type_var"),
      names_pattern = "(.*)(\\.Est$|\\.pval$|\\.pval\\.holm$|\\.effect$|\\.effect\\.low$|\\.effect\\.high$)",
      values_to = "measures"
    ) %>%
    mutate(target_var = gsub(sel_var, "", target_var),
           type_var = stri_replace_all_regex(type_var,
                                             pattern = c("\\.Est",
                                                         "\\.pval$",
                                                         "\\.pval\\.holm",
                                                         "\\.effect",
                                                         "\\.effect\\.low",
                                                         "\\.effect\\.high"),
                                             replacement = c("Estimate",
                                                             "pval",
                                                             "qval",
                                                             "effect_size",
                                                             "effect_CI_lw",
                                                             "effect_CI_up"),
                                             vectorize_all = F)) %>%
    tidyr::spread(type_var, measures) %>%
    mutate(target_var = gsub("\\.", " ", target_var))
  
  return(reshapedf)
}

PIE_df_generic = function(df, tax_level = "Phylum", rel=T, rel_cut = 0.01, abund_cut = 1000){
  
  if(rel){
    pie_out = df %>%
      group_by(!!sym(tax_level)) %>%
      summarise(tax_abund = sum(Abundance)) %>%
      mutate(total_abund = sum(tax_abund), rel_abund = tax_abund / total_abund, new_label = ifelse(rel_abund < rel_cut, paste0("<", 100*rel_cut, "%"), Phylum)) %>%
      group_by(new_label) %>%
      summarise(rel_abund2 = sum(rel_abund), display = paste0(round(100*rel_abund2, 1), "%")) %>%
      arrange(desc(rel_abund2))
  } else {
    pie_out = df %>%
      group_by(!!sym(tax_level)) %>%
      summarise(tax_abund = sum(Abundance)) %>%
      mutate(new_label = ifelse(tax_abund < abund_cut, "Other", Phylum)) %>%
      group_by(new_label) %>%
      summarise(abund = sum(tax_abund), display = abund) %>%
      arrange(desc(abund))
  }
  
  return(pie_out)
}

plot_PIE_generic = function(pieDF, tax_level = "Phylum", pieBorderCol = "white", ColPal = NULL, xlab = "", ylab ="", filllab =NULL, plot_title="", xFontSize=12, legFontSize = 12, legTitleSize=12,legendRow = 1,legendPos = "bottom"){
  
  library(ggpubr)
  
  # make plot
  
  
  pieDF$new_label = factor(pieDF$new_label, levels=pieDF$new_label)
  
  if(!is.null(ColPal)){
    pie_plot = ggpubr::ggpie(data = pieDF, 
                             x = "rel_abund2", 
                             label = "display", 
                             fill = "new_label",
                             color=pieBorderCol, palette = ColPal) 
  } else {
    pie_plot = ggpubr::ggpie(data = pieDF, 
                             x = "rel_abund2", 
                             label = "display", 
                             fill = "new_label",
                             color=pieBorderCol) 
  }
  
  if(is.null(filllab)){
    filllab = tax_level
  }
  
  pie_plot = pie_plot + 
    labs(fill = filllab, title = plot_title) + 
    theme(
    axis.text.x =  element_text(size = xFontSize, margin = margin(t=5)), 
    legend.text = element_text(size= legFontSize),
    legend.title = element_text(size= legTitleSize),
    legend.position = legendPos,
    plot.title = element_text(hjust = 0.5)) +
    guides(fill=guide_legend(nrow = legendRow, byrow = T))
  
  
  return(pie_plot)
}


list2df = function(...){
  tlst = tibble::lst(...)
  tlst = lapply(tlst, sort)
  tlst_df = data.frame(lapply(tlst, `length<-`, max(lengths(tlst))))
  return(list(lst=tlst, df=tlst_df))
}

###########################################################################################################

# from the microbiome package


plot_core_mod <- function(x, prevalences=seq(.1, 1, 0.1), detections=20,
                      plot.type="lineplot", colours=NULL, # gray(seq(0, 1, length=5)),
                      min.prevalence=NULL, taxa.order=NULL, horizontal=FALSE) {
  
  if (length(detections) == 1) {
    detections <- 10^seq(log10(0.001), log10(max(abundances(x),
                                                 na.rm=TRUE)), length=detections)
  }
  
  if (!is_compositional(x)) {
    warning("The plot_core function is typically used with compositional 
                data. The data is not compositional. Make sure that you
                intend to operate on non-compositional data.")
  }
  
  if (plot.type == "lineplot") {
    
    # Calculate the core matrix (prevalences x abundance thresholds)
    coremat <- core_matrix(x, prevalences, detections)
    res <- core_lineplot(coremat)
    
  } else if (plot.type == "heatmap") {
    
    # Here we use taxon x abundance thresholds table
    #  indicating prevalences
    res <- core_heatmap(
      abundances(x),
      dets=detections,
      cols=colours, 
      min.prev=min.prevalence,
      taxa.order=taxa.order)
  }
  
  p <- res$plot
  
  if (horizontal) {
    p <- p + coord_flip() + theme(axis.text.x=element_text(angle=90))
  }
  
  p
  
}


core_matrix <- function(x, prevalences=seq(0.1, 1, 1), detections=NULL) {
    
    # Pick abundances
    data <- microbiome::abundances(x)
    
    # Convert prevalences from percentages to sample counts
    p.seq <- 0.01 * prevalences * ncol(data)
    
    ## Intensity vector
    if (is.null(detections)) {
        detections <- seq(min(data), max(data), length=10)
    }
    i.seq <- detections
    
    coreMat <- matrix(NA, nrow=length(i.seq), ncol=length(p.seq),
        dimnames=list(i.seq, p.seq))
    
    n <- length(i.seq) * length(p.seq)
    cnt <- 0
    for (i in i.seq) {
        for (p in p.seq) {
            # Number of OTUs above a given prevalence
            coreMat[as.character(i), as.character(p)] <-
            sum(rowSums(data > i) >= p)                
        }
    }
    
    # # Convert Prevalences to percentages
    colnames(coreMat) <- as.numeric(colnames(coreMat))/ncol(data)
    rownames(coreMat) <- as.character(as.numeric(rownames(coreMat)))
    
    coreMat
    
}

core_heatmap <- function(x, dets, cols, min.prev, taxa.order)
{

    data <- x
    DetectionThreshold <- Taxa <- Prevalence <- NULL
    
    # Prevalences with varying dets
    
    prev <- lapply(dets, function(th) {
        microbiome::prevalence(data, detection=th)
    })
    prev <- do.call("cbind", prev)
    colnames(prev) <- as.character(dets)

    # Exclude rows and cols that never exceed the given prevalence
    if (!is.null(min.prev)) {
        rinds <- rowMeans(prev > min.prev) > 0
	      cinds <- colMeans(prev > min.prev) > 0
        prev <- prev[rinds, cinds, drop=FALSE]	
    }

    df <- as.data.frame(prev)
    if (nrow(df) == 0) {stop("Too few taxa fulfil the criteria on detection and prevalence. Apply less conservative limits.")}
    
    df$ID <- rownames(prev)

    df <- melt(df, "ID")
    names(df) <- c("Taxa", "DetectionThreshold", "Prevalence")
    
    df$DetectionThreshold <- as.numeric(as.character(df$DetectionThreshold))
    df$Prevalence <- as.numeric(as.character(df$Prevalence))
    df$DetectionThreshold <- factor(df$DetectionThreshold)

    if (is.null(taxa.order)) {
        o <- names(sort(rowSums(prev)))
    } else {
        o <- taxa.order
    }
    
    df$Taxa <- factor(df$Taxa, levels=o)

    outPlot <- ggplot2::ggplot(df, aes(x=DetectionThreshold, y=Taxa, fill=Prevalence)) +
            geom_tile()
      
    
    if (is_compositional(x)) {

        outPlot <- outPlot + scale_x_discrete(labels = function(x)round(100*as.numeric(x),3))

        if (!is.null(cols)) {
          outPlot <- outPlot + scale_fill_gradientn(
                    "Prevalence (%)",
                    # breaks=seq(from=0, to=1, by=0.1),
                    labels=scales::percent_format(suffix = ""),
                    colours=cols,
                    limits=c(0, 1))
        }

    } else {

        if (!is.null(cols)) {

          outPlot <- outPlot + scale_fill_gradientn(
                "Prevalence",
                breaks=seq(from=0, to=1, by=0.1),
                colours=cols,
                limits=c(0, 1))
        }

    }
    outPlot <- outPlot + labs(x = "Detection Threshold")
        
    return(list(plot=outPlot, data=df))
    
}


core_lineplot <- function(x, 
    xlabel="Abundance", ylabel="Core size (N)") {

    Abundance <- Prevalence <- Count <- NULL

    df <- as.data.frame(x)
    df$ID <- rownames(x)
    df <- melt(df, "ID")    
    names(df) <- c("Abundance", "Prevalence", "Count")
    
    df$Abundance <- as.numeric(as.character(df$Abundance))
    df$Prevalence <- as.numeric(as.character(df$Prevalence))
    df$Count <- as.numeric(as.character(df$Count))
    
    p <- ggplot(df, aes(x=Abundance, y=Count,
        color=Prevalence, group=Prevalence))
    
    p <- p + geom_line()
    p <- p + geom_point()
    p <- p + scale_x_log10()
    p <- p + xlab(xlabel)
    p <- p + ylab(ylabel)
    
    list(plot=p, data=x)
}



####################################################################################