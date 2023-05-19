

plot_betadiv = function(phyObj, norm = T, 
                        norm_method = "identity", 
                        dist_method = "bray", 
                        neg_correct = "none", 
                        rn= NULL, 
                        group_col = NULL, 
                        pairwise_method = "HSD", 
                        palette=NULL,
                        jitter_size = 2, 
                        box_alpha=0.5, 
                        CI_Ellipse=F, 
                        DATA_Ellipse=F, 
                        ell_type="t", 
                        ell_linetype=2,
                        point_size=5,
                        letter_size = 5,
                        legend_title = "Group",
                        Pcoa_legend = T,
                        legend_pos = "left",
                        plot_title = "",
                        yside_angle = 90,
                        ifsetLegendLabel = F,
                        newLegendLabel = NULL,
                        legCol = 1,
                        legRow = 1,
                        nperm = 1000, nudge=0.2, seeds=123, plotmode = 1, ...){
  
  library(agricolae)
  library(patchwork)
  library(ggside)
  
  if(is.null(palette)){
    palette = c("#E64B35FF","#4DBBD5FF","#00A087FF","#3C5488FF","#F39B7FFF","#8491B4FF",
                "#B2182B","#E69F00","#56B4E9","#009E73","#F0E442","#0072B2","#D55E00","#CC79A7","#CC6666")
  }
  
  # normalization
  if(norm){
    phy.norm = microbiome::transform(phyObj, transform = norm_method)
  } else {
    phy.norm = phyObj
  }
  
  if(is.null(group_col)){
    warnings("Automatic selected the first column as the variable of interest")
    group_col = 1
  } else {
    group_col = group_col
  }
  

  
  
  meta.phy = data.frame(sample_data(phy.norm))
  
  if(ifsetLegendLabel){
    if(is.null(newLegendLabel)){
      warning("No new legend label provided, use Default\n")
      setLegendLabel = levels(meta.phy[[group_col]])
    }
    setLegendLabel = newLegendLabel
  } else {
    cat("Used Default legend labels\n")
    setLegendLabel = levels(meta.phy[[group_col]])
  }
  

  
  # if(ifrelevel){
  #   meta.phy[[group_col]] = factor(meta.phy[[group_col]], levels = newlevel)
  # }
  # 
  # create distance object
  phy.dist = vegan::vegdist(t(otu_table(phy.norm)), method = dist_method)
  phy.dist
  
  # get pcoa ordination
  pcoa = ape::pcoa(phy.dist, correction = neg_correct, rn=rn)
  pcoa
  
  PC1 = pcoa$vectors[, 1] 
  PC2 = pcoa$vectors[, 2]
  PC3 = pcoa$vectors[, 3]
  
  plotdata <- data.frame(rownames(pcoa$vectors),PC1,PC2,PC3, meta.phy[[group_col]])
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
    
    fit1_output = group_test(fit = fit1,method = pairwise_method, mapping = meta.phy[[group_col]])
    fit2_output = group_test(fit = fit2,method = pairwise_method, mapping = meta.phy[[group_col]])
    fit3_output = group_test(fit = fit3,method = pairwise_method, mapping = meta.phy[[group_col]])
    
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
    

    p1 = ggplot(plotdata,aes(Group,PC1)) +
      geom_boxplot(aes(fill = Group),outlier.colour = NA, alpha = box_alpha) +
      scale_fill_manual(values=palette)+
      scale_color_manual(values = palette) +
      geom_jitter(aes(color=Group), width = 0.1, size=jitter_size) +
      geom_text(data = test,aes(x = Group,y = yd1,label = PC1),
                size = letter_size,color = "black") +
      coord_flip() +
      theme_bw()+
      theme(
            axis.ticks.length = unit(0.1,"lines"),
            axis.ticks = element_line(color='black'),
            axis.line = element_line(colour = "black"),
            axis.title.x=element_blank(),
            axis.title.y=element_blank(),
            axis.text.y=element_text(colour='black',size=12),
            axis.text.x=element_blank(),
            legend.position = "none")

    p2 = ggplot(plotdata,aes(Group,PC2)) +
      geom_boxplot(aes(fill = Group),outlier.colour = NA, alpha = box_alpha) +
      geom_jitter(aes(color=Group), width = 0.1, size=jitter_size) +
      scale_fill_manual(values=palette) +
      scale_color_manual(values = palette) +
      geom_text(data = test,aes(x = Group,y = yd2,label = PC2),
                size = letter_size,color = "black") +
      theme_bw()+
      theme(
            axis.ticks.length = unit(0.1,"lines"),
            axis.ticks = element_line(color='black'),
            axis.line = element_line(colour = "black"),
            axis.title.x=element_blank(),
            axis.title.y=element_blank(),
            axis.text.x=element_text(colour='black',size=12, angle = 45, vjust = 1, hjust = 1),
            axis.text.y = element_blank(),
            legend.position = "none")
    
  }else{
    warning('The pairwise comparison method not in the list')
  }
  
  # add PERMONVA result to the plot
  set.seed(seeds)
  fml = paste(quote(phy.dist), "~", group_col)
  fml2 = stats::as.formula(fml)
  
  pmanova = adonis2(fml2, data = meta.phy, method = dist_method, permutations = nperm, ...)
  print(pmanova)
  
  
  R2 = round(pmanova$R2[1], 3)
  pvalue = round(pmanova$`Pr(>F)`[1],4)
  
  cat("\n")
  print(paste("Permanova", fml, " nperm:", nperm))
  print(paste("R2:", R2))
  print(paste("pvalue:", pvalue))

  label_manova = c("'PERMANOVA:\n'", 
                   sprintf("~R^2=='%.3f'\n", R2), 
                   sprintf("~italic(P)-value=='%.4f'\n", pvalue))
  
  
  if(Pcoa_legend){
    
    p12<-ggplot(plotdata, aes(PC1, PC2)) +
      geom_point(aes(color=Group), size=point_size, alpha = box_alpha) +
      annotate("text", 
               x=max(ggplot_build(p1)$layout$panel_scales_y[[1]]$range$range)-max(ggplot_build(p1)$layout$panel_scales_y[[1]]$range$range)*nudge, 
               # three y values
               y=c(max(ggplot_build(p2)$layout$panel_scales_y[[1]]$range$range)-max(ggplot_build(p2)$layout$panel_scales_y[[1]]$range$range)*nudge
                   ,
                   max(ggplot_build(p2)$layout$panel_scales_y[[1]]$range$range)-1*max(ggplot_build(p2)$layout$panel_scales_y[[1]]$range$range)*nudge,
                   max(ggplot_build(p2)$layout$panel_scales_y[[1]]$range$range)-1.5*max(ggplot_build(p2)$layout$panel_scales_y[[1]]$range$range)*nudge
                   ),
               label = label_manova, parse=T) +
      scale_color_manual(values=palette, labels = setLegendLabel)+
      xlab(paste("PCoA1 ( ",pc1,"%"," )",sep="")) +
      ylab(paste("PCoA2 ( ",pc2,"%"," )",sep=""))+
      xlim(ggplot_build(p1)$layout$panel_scales_y[[1]]$range$range) +
      ylim(ggplot_build(p2)$layout$panel_scales_y[[1]]$range$range) +
      theme(text=element_text(size=12))+
      geom_vline(aes(xintercept = 0),linetype="dotted")+
      geom_hline(aes(yintercept = 0),linetype="dotted")+
      theme(panel.background = element_rect(fill='white', colour='black'),
            panel.grid=element_blank(),
            axis.title = element_text(color='black',size=12),
            axis.ticks.length = unit(0.4,"lines"), axis.ticks = element_line(color='black'),
            axis.line = element_line(colour = "black"),
            axis.title.x=element_text(colour='black', size=12),
            axis.title.y=element_text(colour='black', size=12),
            axis.text=element_text(colour='black',size=12),
            legend.title=element_text(size = 12),
            legend.text=element_text(size=12),
            legend.key=element_blank(),
            legend.position = c(legend_pos),
            # legend.background = element_rect(colour = "black"),
            legend.key.height=unit(1,"cm")) +
      guides(fill = guide_legend(ncol = legCol, nrow = legRow, override.aes = list(alpha = box_alpha)), 
             colour = guide_legend(ncol = legCol, nrow = legRow, override.aes = list(alpha = box_alpha)))
  } else {
    
    p12<-ggplot(plotdata, aes(PC1, PC2)) +
      geom_point(aes(color=Group), size=point_size, alpha = box_alpha) +
      scale_color_manual(values=palette, labels = setLegendLabel)+
      xlab(paste("PCoA1 ( ",pc1,"%"," )",sep="")) +
      ylab(paste("PCoA2 ( ",pc2,"%"," )",sep=""))+
      xlim(ggplot_build(p1)$layout$panel_scales_y[[1]]$range$range) +
      ylim(ggplot_build(p2)$layout$panel_scales_y[[1]]$range$range) +
      annotate("text", 
               x=max(ggplot_build(p1)$layout$panel_scales_y[[1]]$range$range)-max(ggplot_build(p1)$layout$panel_scales_y[[1]]$range$range)*nudge, 
               # three y values
               y=c(max(ggplot_build(p2)$layout$panel_scales_y[[1]]$range$range)-max(ggplot_build(p2)$layout$panel_scales_y[[1]]$range$range)*nudge
                   ,
                   max(ggplot_build(p2)$layout$panel_scales_y[[1]]$range$range)-1*max(ggplot_build(p2)$layout$panel_scales_y[[1]]$range$range)*nudge,
                   max(ggplot_build(p2)$layout$panel_scales_y[[1]]$range$range)-1.5*max(ggplot_build(p2)$layout$panel_scales_y[[1]]$range$range)*nudge
               ),
               label = label_manova, parse=T) +
      theme(text=element_text(size=12))+
      geom_vline(aes(xintercept = 0),linetype="dotted")+
      geom_hline(aes(yintercept = 0),linetype="dotted")+
      theme(panel.background = element_rect(fill='white', colour='black'),
            panel.grid=element_blank(),
            axis.title = element_text(color='black',size=12),
            axis.ticks.length = unit(0.4,"lines"), axis.ticks = element_line(color='black'),
            axis.line = element_line(colour = "black"),
            axis.title.x=element_text(colour='black', size=12),
            axis.title.y=element_text(colour='black', size=12),
            axis.text=element_text(colour='black',size=12),
            # legend.title=element_text(size = 12),
            # legend.text=element_text(size=12),
            # legend.key=element_blank(),
            legend.position = "none",
            # legend.background = element_rect(colour = "black"),
            legend.key.height=unit(1,"cm"))
      # guides(fill = guide_legend(ncol = 1))
  }
  
  if(CI_Ellipse){
    p12 = p12+ggpubr::stat_conf_ellipse(mapping = aes(color=Group), linewidth=1, linetype = ell_linetype)
    
  } else if(DATA_Ellipse){
    p12 = p12+ggpubr::stat_ellipse(mapping = aes(color=Group), type=ell_type, linetype=ell_linetype)
  }
  if(plotmode==1){
    plot_all = p12 +
      geom_xsideboxplot(data = plotdata, aes(xfill=Group, y=Group), orientation = "y", alpha = box_alpha) +
      geom_xsidetext(data = test,aes(y = Group,x = yd1, label = PC1),
                     size = letter_size,color = "black") +
      geom_xsidepoint(aes(color=Group, y=Group, x = PC1), size=jitter_size) +
      
      geom_ysideboxplot(aes(yfill = Group, x = Group), orientation = "x", alpha = box_alpha) +
      geom_ysidetext(data = test,aes(x = Group, y = yd2, label = PC2),
                     size = letter_size,color = "black") +
      geom_ysidepoint(aes(color=Group, x=Group, y = PC2), size=jitter_size) +
      
      scale_xfill_manual(values  = palette, labels = setLegendLabel) +
      scale_yfill_manual(values  = palette, labels = setLegendLabel) +
      scale_xsidey_discrete() +
      scale_ysidex_discrete(guide = guide_axis(angle=yside_angle)) +
      theme(
        ggside.panel.scale = 0.2
      ) +
      labs(title = plot_title, fill=legend_title, color=legend_title, xfill = legend_title, yfill = legend_title, xcolor = legend_title, ycolor = legend_title)
  } else {
    plot_all = p1 + plot_spacer() + p12 + p2 + plot_layout(heights = c(1,5), widths=c(5,1), ncol=2, nrow = 2) + plot_annotation(title=plot_title)
  }

  return(plot_all)
  
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