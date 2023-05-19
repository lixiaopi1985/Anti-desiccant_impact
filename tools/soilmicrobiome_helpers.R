deseq_phylo = function(phylo, model = "~ State + Season", contrast = c("YearSeason", "2018_LF", "2018_ES"), alpha=0.01, ifout = T, nfold=2, output="./output/differential/sig2017_4fold_order.csv"){
  
  library(DESeq2)
  
  gm_mean = function(x, na.rm=TRUE){
    exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
  }
  
  evalmodel = formula(model)
  
  diagdds = phyloseq_to_deseq2(phylo, evalmodel)
  geoMeans = apply(counts(diagdds), 1, gm_mean)
  diagdds = estimateSizeFactors(diagdds, geoMeans = geoMeans)
  diagdds = DESeq2::DESeq(diagdds, test="Wald", fitType = "parametric")
  # resultsNames(diagdds)
  
  res0 = results(diagdds, contrast = contrast)
  res1 = lfcShrink(diagdds, contrast = contrast,
                   res = res0, type = "normal")
  
  
  
  sigTab = res1[which(res1$padj < alpha),]
  sigTab = cbind(as(sigTab, "data.frame"), as(tax_table(phylo)[rownames(sigTab), ], "matrix"))
  
  sigTab_sig = sigTab%>%
    filter(padj < alpha & abs(log2FoldChange) >= nfold)
  
  if(ifout){
      sigTab_sig %>%
      write.csv(file=output, quote = F)
  }
  
  return(sigTab_sig)
  
}




plotComp = function(df, taxlevel = "Genus", colorPalette = "none", plotx = "Season", ploty="rel_abund", plotfill = "genus_lab", facet = "~ State", fill_legend = "Genus", plot_title="", xLabs = "", ylabs = "Relative abundance", textsize=12, axis_x_vjust = -10){
  library(lazyeval)
  
  
  if(length(facet) > 1){
    facet_express = interp(facet)
    comp = df %>%
      ggplot(aes_string(x=plotx, y=ploty, fill=plotfill)) +
      geom_bar(stat = "identity", position = position_fill()) +
      facet_grid(facet_express)+
      scale_fill_manual(values= colorPalette)+
      scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
      theme_minimal() +
      labs(fill=fill_legend, x = xLabs, y = ylabs, title = plot_title) +
      theme(legend.position = "right",
            plot.title = element_text(hjust = 0.5),
            text = element_text(size=textsize),
            axis.text.x = element_text(color="black", margin = margin(t=axis_x_vjust)),
            axis.title.x = element_text(margin = margin(t=5)),
            panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank())
  } else {
    comp = df %>%
      ggplot(aes_string(x=plotx, y=ploty, fill=plotfill)) +
      geom_bar(stat = "identity", position = position_fill()) +
      # facet_grid(facet_express)+
      scale_fill_manual(values= colorPalette)+
      scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
      theme_minimal() +
      labs(fill=fill_legend, x = xLabs, y = ylabs, title = plot_title) +
      theme(legend.position = "right",
            plot.title = element_text(hjust = 0.5),
            panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank(),
            text = element_text(size=textsize),
            axis.text.x = element_text(color="black", margin = margin(t=axis_x_vjust)),
            axis.title.x = element_text(margin = margin(t=5))) +
      guides(fill=guide_legend(ncol=1))
  }
  
  print(comp)
  
  return(comp)
}

colorPats = function(..., col = "Paired", ncol=12, special=""){
  colPatte = colorRampPalette(brewer.pal(ncol, col))
  
  dflist = list(...)

  if(length(dflist)>1){
    uniqueGen = unique(c(dflist[[1]]$A$genus_lab, dflist[[2]]$A$genus_lab))
  } else {
    uniqueGen = unique(dflist[[1]]$A$genus_lab)
  }
  
  print(paste("Total colors", length(uniqueGen)))
  gen_col = colPatte(length(uniqueGen))
  print(gen_col)
  names(gen_col) = uniqueGen
  
  gen_col[special] = "grey80"
  pos = grep(special, names(gen_col))
  pos
  gen_col.2 = gen_col[-pos]
  orders = c(special, names(gen_col.2))
  orders
  gen_col = gen_col[orders]
  print(paste("After reordering", length(gen_col)))
  
  return(gen_col)
}

kw_test = function(orig_df, bac_vector, taxlevel="Order", taxfilter="Bacteria", filterYear = 2017, groups = c("Season", "State"), out=F, outdir="", test1="Season", test2="State", testOn="RA", speciesCol="Species2", padjust_method="fdr"){
  
  
  # top_bac_figs = list()
  top_bac_data = list()
  test_result = list()

  for(n in 1:length(bac_vector) ){
    
    i = bac_vector[n]
    cat(n, " ", i, "\n")
    
    xdata = orig_df %>%
      dplyr::filter(Year == filterYear & Domain == taxfilter) %>%
      group_by(!!!syms(groups)) %>%
      dplyr::mutate(total_abund = sum(Abundance)) %>%
      group_by(!!!syms(c(groups, speciesCol))) %>%
      dplyr::mutate(rel_abund = 100*(Abundance / total_abund)) %>%
      dplyr::select(Season, State, !!sym(taxlevel), Abundance, rel_abund) %>%
      dplyr::filter(!!sym(taxlevel) == i)

    
    top_bac_data[[i]] = xdata

    if(testOn == "RA"){
      testSeason = kruskal.test(xdata$rel_abund ~ xdata[[test1]])
      testState = kruskal.test(xdata$rel_abund ~ xdata[[test2]])
    } else if(testOn == "TA"){
      testSeason = kruskal.test(xdata$Abundance ~ xdata[[test1]])
      testState = kruskal.test(xdata$Abundance ~ xdata[[test2]])
    }

    
    test_result[[i]]["Season_chisq"] = round(testSeason$statistic,4)
    test_result[[i]]["Season_pvalue"] = round(testSeason$p.value,4)
    test_result[[i]]["State_chisq"] = round(testState$statistic, 4)
    test_result[[i]]["State_pvalue"] = round(testState$p.value,4)
    
    
    if(out){
      write.csv(rbind(testSeason, testState), paste0(outdir, "/kw_test", filterYear, "_", i, ".csv"), quote = F)
    }
  }
  
  # adjust p value
  names(test_result) = gsub(" ", "_", names(test_result))
  df0 = t(data.frame(test_result))
  df = data.frame(df0)
  
  df$Season.padj = p.adjust(df$Season_pvalue, method = padjust_method)
  df$State.padj = p.adjust(df$State_pvalue, method=padjust_method)
  
  df = df[, c("Season_chisq","Season_pvalue","Season.padj", "State_chisq", "State_pvalue","State.padj")]
  

  return(list(data = top_bac_data, testR = test_result, df=df))
  
}

plot_kw = function(kw_out, bac_vector, taxlevel = "Order", facet_ = "~ Order", plotx="State", ploty="rel_abund", fillvar="Season", filllab="Season", xlab="State", ylab="Relative abundance (%)", annoX = 2, gap=1, geomtype="box", Ylimits = NULL, fontsize=12){
  
  # process the data input
  # rbind
  library(lazyeval)
  library(rstatix)
  
  figure_list = list()
  
  xdata = Reduce(function(x, y)rbind(x, y), kw_out$data)
  xdata[[taxlevel]] = gsub(" ", "_", xdata[[taxlevel]])
  xtest = kw_out$df
  xtest[[taxlevel]] = gsub("_", " ", rownames(xtest))

  # xtest = Reduce(function(x, y)rbind(x, y), kw_out$testR)
  
  # row.names(xtest) = names(kw_out$testR)
  # print(xtest)
  # 
  # colnames(xtest) = c("stats_Season" ,"pSeason", "stats_State","pState")
  
  merged = merge(xdata, xtest, by.x = taxlevel, by.y = "row.names")
  merged[[taxlevel]] = gsub("_", " ", merged[[taxlevel]])
  
  
  
  # set scale 
  if(is.null(Ylimits)){
    if(ploty == "rel_abund"){
      ymin = min(merged$rel_abund)
      ymax = max(merged$rel_abund)
    } else if(ploty=="Abundance"){
      ymin = min(merged$Abundance)
      ymax = max(merged$Abundance)
    }

  } else {
    ymin = Ylimits[1]
    ymax = Ylimits[2]
  }

  
  for(i in bac_vector){
    

    AnnotateSeason = xtest[which(xtest[[taxlevel]] == i),"Season.padj"]
    AnnotateState = xtest[which(xtest[[taxlevel]] == i), "State.padj"]
    
    selectdf = merged %>%
      filter(!!sym(taxlevel) == i)

    if(geomtype == "box"){
      figure =
        selectdf %>%
        ggplot(., aes_string(x = plotx, y = ploty, fill = fillvar, color = fillvar)) +
        geom_boxplot(alpha = 0.9) +
        scale_y_continuous(limits = c(ymin, ymax)) +
        facet_wrap(interp(facet_)) +
        annotate("text", x = annoX, y = ymax, label = paste0("Season:", "P.adj=", round(AnnotateSeason, 4)), hjust=0, size=5) +
        annotate("text", x = annoX, y = ymax - gap,  label = paste0("State:",  "P.adj=", round(AnnotateState, 4)), hjust=0, size=5) +
        theme_bw() +
        labs(x= xlab, y=ylab, fill=filllab) +
        theme(
              # text = element_text(size=fontsize),
              axis.text.x =  element_text(colour = "black", size = fontsize),
              axis.text.y = element_text(colour = "black", size=fontsize),
              axis.title.x = element_text(margin = margin(t=0.2), size=fontsize),
              axis.title.y = element_text(margin = margin(r=0.2), size=fontsize))
      
      figure_list[[i]] = figure
      
    } else if(geomtype == "bar"){
      figure =
        selectdf %>%
        group_by(State, Season, !!sym(taxlevel)) %>%
        summarise(accum_abund = sum(rel_abund), sd = sd(rel_abund))  %>%
        ggplot(., aes_string(x = plotx, y = "accum_abund", fill = fillvar)) +
        geom_bar(stat="identity", position=position_dodge(0.9), alpha = 0.9) +
        geom_errorbar(aes(ymin=accum_abund - sd, ymax=accum_abund + sd), width=.2, position = position_dodge(0.9)) +
        scale_y_continuous(limits = c(ymin, ymax)) +
        facet_wrap(interp(facet_)) +
        annotate("text", x = annoX, y = ymax,  label = paste0("Season:", "P.adj=", round(AnnotateSeason, 4)), hjust=0) +
        annotate("text", x = annoX, y = ymax - gap,  label = paste0("State:",  "P.adj=", round(AnnotateState, 4)), hjust=0) +
        theme_bw() +
        labs(x= xlab, y=ylab, fill=filllab) +
        theme(
              # text = element_text(size=12),
              axis.text.x =  element_text(colour = "black", size = fontsize),
              axis.text.y = element_text(colour = "black", size = fontsize),
              axis.title.x = element_text(margin = margin(t=0.2), size = fontsize),
              axis.title.y = element_text(margin = margin(r=0.2), size = fontsize))
      
      
      figure_list[[i]] = figure
    }
    
  }
  
  return(list(figlist = figure_list, mergedDF = merged))
}

aggreDF = function(df, yearfilter = 2017, taxfilter = "Bacteria", ntop = 10, groups=c("Season", "State"), perc = 1, taxlevel = "Genus", sigStart = F, kw_out = NA, alpha = 0.05){
  
  thres = paste0("< ", perc, "%")
  
  if(sigStart){
    
    
    if(!is.na(kw_out)){
      
      xtest_df = kw_out$df
      rownames(xtest_df) = gsub("_", " ", rownames(xtest_df))
      # xtest = Reduce(function(x, y)rbind(x, y), kw_out$testR)
      # row.names(xtest) = names(kw_out$testR)
      # xtest_df = as.data.frame(xtest)
      
      xtest.season = xtest_df$Season.padj
      names(xtest.season) = rownames(xtest_df)
      xtest.state = xtest_df$State.padj
      names(xtest.state) = rownames(xtest_df)
      
      
      cleandf = df %>%
        dplyr::filter(Year == yearfilter & Domain == taxfilter) %>%
        dplyr::group_by(!!!syms(groups)) %>%
        dplyr::mutate(total_abund = sum(Abundance)) %>%
        dplyr::group_by(!!!syms(c(groups, taxlevel))) %>%
        dplyr::mutate(rel_abund = 100*Abundance / total_abund) %>%
        dplyr::select(Season, State, !!sym(taxlevel), Abundance, total_abund, rel_abund) %>%
        mutate(genus_lab_1 = ifelse(rel_abund < perc, thres, as.character(!!sym(taxlevel))))
      
      newlab = c()
      
      for(i in cleandf$genus_lab_1){
        
        
        i = gsub("_", " ", i)
        print(i)
        
        if(i %in% rownames(xtest_df)){
          
          if( (xtest.season[i] < alpha & xtest.state[i] < alpha)){
            lstars = paste(i, "*/*")
            print(lstars)
          } else if( (xtest.season[i] < alpha & xtest.state[i] > alpha)){
            lstars = paste(i, "*/")
            print(lstars)
          } else if((xtest.season[i] > alpha & xtest.state[i] < alpha)){
            lstars = paste(i, "/*")
            print(lstars)
          } else if((xtest.season[i] > alpha & xtest.state[i] > alpha)){
            lstars = i
          } 
        } else {
          lstars = i
        }
        
        newlab = c(newlab, lstars)
      }
      
      cleandf$genus_lab = newlab
      
    } else {
      
      errorCondition("kw data needed")
    }
  } else {
    
    cleandf = df %>%
      dplyr::filter(Year == yearfilter & Domain == taxfilter) %>%
      dplyr::group_by(!!!syms(groups)) %>%
      dplyr::mutate(total_abund = sum(Abundance)) %>%
      dplyr::group_by(!!!syms(c(groups, taxlevel))) %>%
      dplyr::mutate(rel_abund = 100* (Abundance / total_abund)) %>%
      dplyr::select(Season, State, !!sym(taxlevel), Abundance, total_abund, rel_abund) %>%
      mutate(genus_lab = ifelse(rel_abund < perc, thres, as.character(!!sym(taxlevel))))
  }
  
  
  return(list("A"=cleandf, "B"=thres))
}


plotTukeyDive = function(inputDF, aovout,  gap=100, plotx="State", ploty="diversity", xlabs="State", ylabs="Diversity", fillvar="State", facet_by = "~ Season", LetterTerm = "State:Season", groupbys = c("State, Season"), sigfun = multcompLetters4, textsize=12){
  
  library(lazyeval)
  library(stringr)
  
  tukeytest = TukeyHSD(aovout)
  
  tukeyletters = sigfun(aovout, tukeytest)
  print(tukeyletters)
  
  A1 = str_to_title(str_split(LetterTerm, ":")[[1]][1])
  B1 = str_to_title(str_split(LetterTerm, ":")[[1]][2])
  
  joincol = paste0(A1, B1)

  df = as.data.frame(tukeyletters[[LetterTerm]]$Letters)
  colnames(df) = "letters"
  df[[joincol]] = rownames(df)
  
  df_sep = df %>%
    tidyr::separate(!!sym(joincol), c('state', 'season'), remove = F)
  
  
  
  inputDF[[joincol]] = paste(inputDF[[A1]], inputDF[[B1]],  sep = ":")
  
  input.merge = merge(inputDF, df_sep, by.x=joincol, by.y = joincol, all.x = T)

 
  
  # find the max value then keep it there
  
  maxYdf = input.merge %>%
    group_by(!!!syms(groupbys)) %>%
    summarise(maxY = max(diversity)) 
  
  maxYdf[[joincol]] = paste(maxYdf[[A1]], maxYdf[[B1]],  sep = ":")
  
  
  dfwMaxY = merge(df_sep, maxYdf, by.x = joincol, by.y = joincol)
  
  
  #( plot )
  
  g = ggplot(input.merge, aes_string(x = plotx, y = ploty, fill=fillvar, color=fillvar)) +
    geom_boxplot(alpha=0.5) +
    geom_point(size=2, alpha=0.5) +
    geom_text(data = dfwMaxY, aes(x=state, y=maxY+gap, label=letters), color="black") +
    facet_grid(interp(facet_by)) +
    theme_biome_utils() +
    scale_fill_brewer(palette = "Dark2") +
    scale_color_brewer(palette = "Dark2") +
    labs(x=xlabs, y=ylabs) +
    theme(legend.position = "none",
          text = element_text(size=textsize),
          axis.text.x = element_text(color="black"),
          axis.text.y = element_text(color="black"),
          axis.title.x = element_text(margin = margin(t=5)),
          axis.title.y = element_text(margin = margin(r=5)))
  
  return(list(fig = g, tuk = tukeytest))
}

getDiff = function(df){
  
  (df_ES = df %>%
     select(State, Season,diversity) %>%
     filter(Season == "ES") %>%
     group_by(State, Season) %>%
     summarise(avg_es = mean(diversity)))
  
  (df_LF = df %>%
      select(State, Season,diversity) %>%
      filter(Season == "LF")%>%
      group_by(State, Season) %>%
      summarise(avg_lf = mean(diversity)))
  
  
  df_2 = cbind(df_ES, df_LF)
  print(df_2 %>%
    mutate(dff = avg_es - avg_lf))
}


plotBray = function(phylosq_object, group, threshold = 0.05, method = "bh", xlabs = "Comparisons", ylabs = "Distance (Bray Curtis)", fontsize=15, textfontsize=10, gap=0){
  
  library(FSA)
  library(rcompanion)
  
  df = phyloseq_group_dissimilarity(phylosq_object, group = group, between_groups = T, justDF = T)
  df$Group = gsub("-", "_", df$Group)
  kw = kruskal.test(Dist ~ Group, data = df)
  print(kw)
  
  
  dun = dunnTest(Dist ~ Group, data=df, method = "bh")
  print(dun)
  letterdf = cldList(P.adj ~ Comparison, data = dun$res, threshold = threshold)

  print(letterdf)

  sumdf = df %>%
      group_by(Group) %>%
      dplyr::summarise(maxY = max(Dist))
  
  sumdf.merge = merge(sumdf, letterdf, by.x = "Group",by.y = "Group")
  
  g = ggplot(df, aes(x = Group, y = Dist, fill = Group)) +
    geom_boxplot() +
    geom_text(data=sumdf.merge, aes(label = Letter, y = maxY+gap), size=textfontsize) +
    theme_classic() +
    labs(x = xlabs, y = ylabs, fill = xlabs) +
    theme(text = element_text(size=fontsize),
          axis.title.x = element_text(margin = margin(t=20)),
          axis.title.y = element_text(margin = margin(r=20)),
          axis.text.x = element_text(angle=40, hjust = 1))
  
  return(g)
}


plotDendro = function(psq,dist_method = "bray", coL = "State", yLab = "Bray curtis distance", xLab="", textsize = 4){
  
  library(ggdendro)
  
  # psq.sub = subset_samples(psq, Year == year & Season == season)
  psq.sub = psq

  Dist =   phyloseq::distance(psq.sub, method = dist_method)

  
  meta = data.frame(sample_data(psq))
  site_label = meta[[coL]]
  names(site_label) = rownames(meta)
  
  hc = hclust(Dist)
  hc.newlab = site_label[hc$labels]
  hc.newlab
  dend = as.dendrogram(hc)
  dend.data = dendro_data(dend, type = "rectangle")

  dend.data$labels$label = as.character(hc.newlab)

  
  g = ggplot(dend.data$segments) +
    geom_segment(aes(x = x, y=y, xend=xend, yend = yend))+
    geom_text(data = dend.data$labels, aes(x, y, label = label, color=label), hjust = 1, angle = 90, size = textsize) +
    theme_minimal() +
    theme(
      text = element_text(size=12),
      axis.text.y = element_text(colour = "black"),
      axis.title.x = element_text(margin = margin(t=5)),
      axis.title.y = element_text(margin = margin(r=5)),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.text.x = element_blank(),
      legend.position = "none") +
    labs(y = yLab, x=xLab)
  
  return(g)
}
