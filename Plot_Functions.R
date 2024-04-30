
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Functions for Plot Generation
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


#> Load Libraries ----

library(Seurat)
library(scales)
library(tidyselect)
library(ggpubr)
library(tidyverse)
library(SeuratWrappers)
library(purrr)
library(ggplot2)
library(cowplot)
library(grid)
library(gridExtra)
library(gtable)
library(reshape2)
library(circlize)
library(tidyverse)
library(Biobase)
library(gginnards)
library(dplyr)
library(igraph)
library(sctransform)
library(glmGamPoi)
library(patchwork)



#Global Functions

#loads an RData file, and returns it
loadRData <- function(fileName){
  load(fileName)
  get(ls()[ls() != "fileName"])
}

`%!in%` = Negate(`%in%`)






#> Stacked Violin Plot Functions ----

#' Generate stacked violin plots from a list of genes with appropriate formatting
#'
#' @param obj Seurat object
#' @param features vector of genes for each vln plot
#' @param ident Metadata column to group by (x axis)
#' @param colors vector of colors for each entity of ident
#' @param xlabs x axis labels (leveled)
#' @param negNum Reduces space between stacked plots, ideally = 8, can tinker
#' @param posNum increases space between stacked plots, ideally = 12, can tinker
#' @param split.by metadata column to split the data by
#' @param ylabs can submit a vector of alternate names for y labels, otherwise uses feature vector
#' @param font font size
#' @param title Y axis label position, can be "left" or "right"
#' @param legend include legend, can be NULL, "left" or "right"
#' @param legendLabs labels in legend as they will display on the plot
#'
#' @return A stacked set of violin plots 
#'
#' @importFrom patchwork 

StackedVlnPlot <- function(obj, 
                           features,
                           ident,
                           colors,
                           xlabs,
                           negNum = 8,
                           posNum = 12,
                           split.by = NULL,
                           ylabs = features,
                           font=12,
                           title = "left",
                           legend = NULL,
                           legendLabs = NULL,
                           ...) {
  
  # Require Packages
  requireNamespace("patchwork", quietly = TRUE)
  
  # Force RNA assay
  DefaultAssay(obj) <- "RNA"
  
  # Set margins & alignment
  plot.margin = unit(c(-0.75, 0, -0.75, 0), "cm")
  ylabs <- str_pad(ylabs, 
                   side="left",
                   width=max(nchar(ylabs)), pad=" ")
  
  
  # Create plots, allow for L/R title position
  plot_list<- purrr::map2(features, ylabs, function(x,y) modify_vlnplot(obj = obj,
                                                                        feature = x, 
                                                                        ident = ident,
                                                                        colors = colors,
                                                                        split.by = split.by,
                                                                        ylabs = y,
                                                                        font = font,
                                                                        position = title))
  
  # Add x-axis labels to bottom plot
  if(!is.null(xlabs)){
    plot_list[[length(plot_list)]]<- plot_list[[length(plot_list)]] +
      theme(axis.ticks.x = element_line(),
            axis.text.x = element_text(angle=90,hjust=0.95,vjust=0.2,size=font)) + 
      scale_x_discrete(labels = xlabs)
  } else
  {
    p <- VlnPlot(obj, 
                 features = features[1], 
                 group.by = ident,
                 split.by = split.by)
    
    # Extract x-axis labels as a vector
    x_labels <- levels(p$data$ident)
    
    if(is.null(x_labels)){
      stop("Order your x axis ident before running")
    }
    
    plot_list[[length(plot_list)]]<- plot_list[[length(plot_list)]] +
      theme(axis.ticks.x = element_line(),
            axis.text.x = element_text(angle=90,hjust=0.95,vjust=0.2,size=font)) + 
      scale_x_discrete(labels = x_labels)
  }
  
  if(!is.null(legend)){
    plot_list[[1]]<- plot_list[[1]] +
      theme(legend.position = legend) +
      scale_fill_manual(values = colors,
                        labels = legendLabs)
  }
  
  
  
  #Wrap the plots, with null plot spacers for alignment
  p <- patchwork::wrap_plots(plotlist = plot_list, ncol = 1)
  
  plot_list_cp <- list()
  nums_nulls <- seq(2,((length(plot_list)*2)-1),2)
  nums_plots <- seq(1,(length(plot_list)*2),2)
  
  for(x in 1:(length(plot_list))){
    plot_list_cp[[nums_plots[x]]] <- plot_list[[x]]
    if(x < length(plot_list)){
      plot_list_cp[[nums_nulls[x]]] <- plot_spacer()
    }
  }
  
  p2 <- patchwork::wrap_plots(plot_list_cp, ncol=1)+ plot_layout(heights = c(rep(c(posNum,-negNum),(length(plot_list)-1)),posNum),guides = "collect")
  
  return(p2)
}


#Sub-Function to construct individual violin plots for the stacked vln plot
modify_vlnplot <- function(obj, 
                           feature, 
                           pt.size = 0, 
                           ident,
                           colors,
                           plot.margin = unit(c(-1.25, 0, -0.75, 0), "cm"),
                           split.by = NULL,
                           ylabs,
                           font,
                           position, 
                           ...) {
  # Set default font
  if(is.null(font)){font <- 14}
  
  #Construct left titled plot
  p<- VlnPlot(obj, 
              features = feature, 
              pt.size = pt.size, 
              group.by = ident,
              cols = colors,
              split.by = split.by)  + 
    xlab("") + 
    ylab(ylabs) + 
    ggtitle("") + 
    theme(legend.position = "none", 
          axis.text.x = element_blank(), 
          axis.ticks.x = element_blank(),
          axis.ticks.y = element_blank(), 
          axis.title.y = element_text(size =font, angle = 0, hjust = 1), 
          axis.text.y = element_blank(), 
          plot.margin = unit(c(0,0,0,0), "cm"),
          axis.line.y.left = element_blank())
  
  # Modify in case of right title
  if(position == "right"){
    p<- p + 
      scale_y_continuous(position = 'right') + 
      theme(axis.title.y.right = element_text(angle=0, vjust = 1,size=font),
            axis.line.y.right = element_blank(),
            axis.line.y.left = element_blank())
  }
  
  return(p)
}


  




    
#> Bar Plot for Categorical Proportions ----

#' Proportion Barplot with Multiple Comparisons
#'
#' @param Object Seurat object
#' @param Category_Ident Variable being compared, typically Disease State
#' @param ClusterIdent Metadata column to calculate porportions for, typically clustering
#' @param Cat_List List format, entities of Category_Ident in desired order
#' @param remove_position Vector of ClusterIdent entities to remove, if not interested in plotting
#' @param my_comparisons Comparisons for significance testing
#' @param colors Vector of colors for Category_Ident plotting
#' @param yhi ymax value
#' @param GroupLabels Leveling for Metadata column to calculate porportions for, typically clustering
#' @param stepval Step distance between signficance bars
#' @param xlab_angle x labels angle
#' @param font Font size
#' @param NewNames Names of Category_Ident entities, as will appear on the plot
#' @param ypos Starting y position for significance bars
#' @param hjust hjust value for x labels
#' @param test Test to run, can be "Wilcox" or "t-test"
#' @param bracket_length Length of brackets for significance bars
#' @param labelmeans T/F for displaying the mean value of each bar on the barplot
#' @param jitter_width Width for the jittering of individual proportion points
#' @param clean_stats T/F for whether to show blanks for NS results, or to show "NS" over the comparison bar
#'
#' @return A barplot with multiple comparisons
#'
#' @import tidyverse
#' @import rstatix
#' @import stringr
#' @import ggpubr
#' @import gginnards

Proportion_Barplot_withMultipleComparisons <- function(Object,
                                                       Category_Ident,
                                                       ClusterIdent,
                                                       Cat_List,
                                                       remove_position = NULL,
                                                       my_comparisons,
                                                       colors,
                                                       yhi,
                                                       GroupLabels,
                                                       stepval,
                                                       xlab_angle,
                                                       font,
                                                       NewNames,
                                                       ypos,
                                                       hjust,
                                                       test = "Wilcox",
                                                       bracket_length,
                                                       labelmeans = FALSE,
                                                       jitter_width,
                                                       clean_stats = TRUE) {
  
  # Require packages
  packages <- c("tidyverse", "rstatix", "stringr", "ggpubr", "gginnards")
  
  for (pkg in packages) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      install.packages(pkg)
    } else {
      library(pkg, character.only = TRUE)
    }
  }
  
  #Set other parameters
  Name_List <- unlist(Cat_List)
  
  # Generate proportion tables across identities for each patient, per condition
  MD <- Object@meta.data
  Table_List <- list()
  
  for(x in 1:length(Cat_List)){
    # Filter per listed condition
    Category <- Cat_List[[x]]
    MD_filtered <- MD[MD[,Category_Ident] == Category,]
    
    # Generate proportion table per patient sample across identities
    CellsperIdent <- as.data.frame(as.data.frame.matrix(prop.table(table(MD_filtered[,ClusterIdent], MD_filtered[,"Label"]), margin = 2))[,unique(MD_filtered$Label)])
    
    # Remove Identities if we are not interested in plotting them 
    # [Note: proportions are still calculated out of all cells, not out of cells in the kept identities. Subset first if you wish to plot this calculation.]
    if(!is.null(remove_position)){
      CellsperIdent <- CellsperIdent %>% filter(!row.names(CellsperIdent) %in% remove_position)
    }
    
    # Account for situations where there is only one sample in a requested condition
    if(length(unique(MD_filtered$Label)) == 1){
      colnames(CellsperIdent) <- unique(MD_filtered$Label)
      rownames(CellsperIdent) <- GroupLabels
    }
    
    #Append condition df to list of dfs
    Table_List[[x]] <- CellsperIdent
  }
  
  
  # Generate a df of each condition, proportion, Mean, and SEM
  
  ## Calculation functions
  sem <- function(x){
    sd(x)/sqrt(length(x))
  }
  
  repCount <- length(unique(MD[,ClusterIdent]))
  if(!is.null(remove_position)){
    repCount <- repCount - length(remove_position) 
  }
  
  ## Build base df for plotting & stats
  
  Base_Data <- c()
  
  for(h in 1:repCount){ #Loop through ClusterIdent variable
    for(x in 1:length(Table_List)){ #Sub-loop through conditions (Category_Ident variable)
      
      # Pull current condition (Category_Ident) proportion table
      CurrentTable <- Table_List[[x]]
      
      # Pull current ClusterIdent variable from the table
      Base_df <- data.frame(CurrentTable[h, , drop = FALSE])
      
      # Calculate Mean, SEM, Max, store in df
      columns <- c("Category", "Proportion", "Patient","XAxis_Identity","Mean","SEM", "Max","GroupLabel")
      Base_Data_Component <- data.frame(matrix(nrow = ncol(Base_df), ncol = length(columns))) %>%
        setNames(columns) %>%
        mutate(Patient = colnames(Base_df),
               Proportion = t(Base_df)[,1],
               Category = Name_List[x],
               XAxis_Identity = rownames(Base_df),
               Mean = mean(Proportion),
               SEM = sem(Proportion),
               Max = max(Proportion),
               GroupLabel= GroupLabels[h])
      
      #Append df for every ClusterIdent of each CategoryIdent
      Base_Data <- rbind(Base_Data,Base_Data_Component)
    }
  }
  
  
  # Calculate significance for indicated comparisons, with Benjamini-Hochberg (BH) adjustment for multiple comparisons
  
  if(test == "t-test"){
    stat.test_fdr <- Base_Data %>%
      group_by(XAxis_Identity) %>%
      t_test(Proportion ~ Category,
             comparisons = my_comparisons) %>%
      adjust_pvalue(method = "BH") %>%
      add_significance() %>%
      mutate(p.adj = format(round(p.adj, 2), nsmall = 2))  %>%
      mutate(XAxis_Identity  = factor(XAxis_Identity, levels =GroupLabels)) %>%
      mutate(group2 = factor(group2, levels =Name_List))
  }
  
  if(test == "Wilcox"){
    for(x in 1:length(my_comparisons)){
      stat.test_fdr_new <- Base_Data %>%
        filter(Category %in% my_comparisons[x][[1]]) %>%
        group_by(XAxis_Identity) %>%
        do(w = wilcox.test(Proportion ~ Category, data = ., paired = FALSE, exact = F)) %>%
        dplyr::summarise(XAxis_Identity, Wilcox = w$p.value) %>%
        adjust_pvalue(method = "BH", p.col = "Wilcox") %>%
        add_significance(p.col = "Wilcox.adj") %>%
        mutate(Wilcox.adj = format(round(Wilcox.adj, 2), nsmall = 2))  %>% 
        mutate(group1 = my_comparisons[x][[1]][1]) %>% 
        mutate(group2 = my_comparisons[x][[1]][2])  %>%
        mutate(XAxis_Identity  = factor(XAxis_Identity, levels =GroupLabels)) %>% 
        mutate(group2 = factor(group2, levels = Name_List)) %>%
        mutate(yposcol = ypos)
      
      if(x ==1){
        stat.test_fdr <- stat.test_fdr_new
      }
      if(x > 1){
        stat.test_fdr <- rbind(stat.test_fdr,stat.test_fdr_new)
      }
    }
  }
  
  ## Order each ident
  stat.test_fdr <- stat.test_fdr[order(factor(stat.test_fdr$group2, levels=levels(stat.test_fdr$group2))),]
  stat.test_fdr <- stat.test_fdr[order(factor(stat.test_fdr$XAxis_Identity, levels=levels(stat.test_fdr$XAxis_Identity))),]
  
  ## Determines whether NS comparisons are shown, or left blank
  if(isTRUE(clean_stats)){
    stat.test_fdr <- stat.test_fdr %>%
      filter(!Wilcox.adj.signif == "ns")
  }
  
  # Confirm factoring
  Base_Data$Category <- factor(Base_Data$Category, levels = Name_List)
  Ordering <- GroupLabels
  if(!is.null(remove_position)){
    Ordering <- Ordering[! Ordering %in% remove_position] 
  }
  Base_Data$XAxis_Identity <- factor(Base_Data$XAxis_Identity, levels = Ordering)
  
  
  # Set variables for plotting
  ylab = "Proportion"
  labels <- unique(Base_Data$GroupLabel)
  names(labels) <- unique(Base_Data$XAxis_Identity)
  ylo=0
  
  
  # Generate bar plot
  plot <- Base_Data %>%
    ggplot(aes(x=Category, y=Proportion, fill=Category)) +
    geom_bar(position = 'dodge', stat = 'summary', fun = 'mean') +
    facet_wrap(~XAxis_Identity, scales ="fixed",nrow=1, strip.position = "bottom") +
    geom_errorbar(aes(ymin=Mean-SEM, ymax=Mean+SEM), width=.2) +
    geom_jitter(color="black", size=0.4, width = .001) +
    scale_fill_manual(values=colors, labels = NewNames)+
    theme_bw() +
    theme(axis.text.x = element_blank(),
          axis.text.y = element_text(size=font),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.line.x.bottom = element_line(color="black"),
          axis.line.y.left =  element_line(color="black"),
          panel.border =  element_blank(),
          strip.background = element_blank(),
          strip.placement = "outside",
          strip.text = element_text(size=font),
          legend.text=element_text(size=font),
          legend.title=element_text(size=(font+2)),
          axis.line.x= element_line(),
          axis.title=element_text(size=(font+2)),
          panel.spacing = unit(0, "lines"),
          strip.text.x = element_text(angle=xlab_angle, hjust=hjust)) +
    ggtitle("") +
    xlab("") + 
    ylab(paste(ylab))+
    ylim(ylo, yhi) + 
    scale_y_continuous(expand = c(0, 0),limits=c(ylo, yhi)) + 
    guides(fill=guide_legend(title="Disease State"))    
  
  
  #Add means as values on the plot if indicated
  if(isTRUE(labelmeans)){
    
    plot <- plot + 
      stat_summary(aes(label=round(..y..,2)), fun.y=mean, geom="text", size=6,
                   vjust = -0.5) 
  }
  
  # Confirm plot factoring
  plot$data$Category <- factor(plot$data$Category, levels = Name_List)
  plot$data$XAxis_Identity <- factor(plot$data$XAxis_Identity, levels =  unique(Base_Data$XAxis_Identity))
  
  # Add significance to the plot
  if(test == "t-test"){
    plot1 <- plot + stat_pvalue_manual(stat.test_fdr, label = "p.adj.signif", y.position = ypos, step.increase = stepval, label.size = (font*0.4)) 
  }
  
  if(test == "Wilcox"){  
    plot1 <- plot + stat_pvalue_manual(stat.test_fdr, 
                                       label = "Wilcox.adj.signif", 
                                       y.position = stat.test_fdr$yposcol, 
                                       step.increase = stepval, 
                                       label.size = (font*0.4), 
                                       tip.length = bracket_length) 
  }
  
  #Return plot object
  return(plot1)
}

  





#> Bar Plot for Continuous Variables across Cells ----

#' Expression Barplot with Multiple Comparisons
#'
#' @param Object Seurat object
#' @param Feature Feature being examined
#' @param FeatureLocation Location Location of feature. Can be "gene" or "meta.data"
#' @param Category_Ident Variable being compared, typically Disease State
#' @param ClusterIdent Metadata column to calculate porportions for, typically clustering
#' @param Cat_List Vector format, entities of Category_Ident in desired order
#' @param my_comparisons Comparisons for significance testing
#' @param colors Vector of colors for Category_Ident plotting
#' @param yTitle y axis title as displayed on plot
#' @param LegendTitle title of legend as displayed on plot
#' @param yhi ymax value
#' @param ylo ymin value
#' @param GroupLabels Leveling for Metadata column to calculate porportions for, typically clustering
#' @param PatientLabel Metadata column containing patient/sample identification
#' @param stepval Step distance between signficance bars
#' @param xlab_angle x labels angle
#' @param font Font size
#' @param NewNames Names of Category_Ident entities, as will appear on the plot
#' @param ypos Starting y position for significance bars
#' @param hjust hjust value for x labels
#' @param test Test to run, will be "Wilcox" for current purposes
#' @param bracket_length Length of brackets for significance bars
#' @param jitter_width Width for the jittering of individual proportion points
#' @param clean_stats T/F for whether to show blanks for NS results, or to show "NS" over the comparison bar
#'
#' @return A box plot with multiple comparisons
#'
#' @import tidyverse
#' @import rstatix
#' @import stringr
#' @import ggpubr
#' @import gginnards

Expression_BoxPlot <- function(Object,
                               Feature,
                               FeatureLocation = NULL,
                               ClusterIdent = "Clusters",
                               Category_Ident = "DiseaseState",
                               Cat_List,
                               my_comparisons = list(c("EoE_Biopsy", "Healthy_Control")),
                               colors,
                               yTitle,
                               LegendTitle,
                               yhi,
                               ylo = 0,
                               GroupLabels,
                               PatientLabel = "Patient_Region",
                               stepval = 0,
                               xlab_angle = 90,
                               font = 14,
                               NewNames = NULL,
                               ypos,
                               hjust,
                               test = "Wilcox",
                               bracket_length = 0.01,
                               jitter_width = 0.1,
                               clean_stats = TRUE) {
  
  ## Import libraries
  
  packages <- c("tidyverse", "rstatix", "stringr", "ggpubr", "gginnards")
  
  for (pkg in packages) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      install.packages(pkg)
    } else {
      library(pkg, character.only = TRUE)
    }
  }
  
  # Check if feature is a gene or meta.data parameter (gene signature / pseudotime, etc)
  if(is.null(FeatureLocation)){
    if(Feature %in% rownames(Object)){
      FeatureLocation = "gene"
    } else if(Feature %in% colnames(Object@meta.data)){
      FeatureLocation = "meta.data"
    }
  }
  
  ## Extract average expression from either meta.data or normalized RNA assay
  
  
  df <- CalculateFeatureExpression(Input)
  
  if(FeatureLocation == "meta.data"){
    Exp_df <- Object@meta.data %>% 
      dplyr::select(!!as.name(Category_Ident), !!as.name(PatientLabel), !!as.name(ClusterIdent), !!as.name(Feature)) %>%
      group_by(!!as.name(ClusterIdent), !!as.name(PatientLabel)) %>%
      dplyr::mutate(Mean = mean(!!as.name(Feature))) %>%
      ungroup() %>%
      distinct(!!as.name(ClusterIdent), !!as.name(PatientLabel), .keep_all = T) 
  }
  
  if(FeatureLocation == "gene"){
    Object$CombCol <- paste0(Object@meta.data[[Category_Ident]],
                             "#",
                             Object@meta.data[[PatientLabel]],
                             "#",
                             Object@meta.data[[ClusterIdent]])
    Exp_df <- AverageExpression(Object, 
                                features = Feature,
                                assay = "RNA",
                                slot = "data",
                                group.by = "CombCol")
    
    Exp_df <- data.frame(
      Group = names(as.data.frame(Exp_df$RNA)),
      Mean = as.numeric(Exp_df$RNA[1, ])
    )  %>%
      separate(Group, into = c(Category_Ident, PatientLabel, ClusterIdent), sep = "#") %>%
      mutate(!!as.name(ClusterIdent) := factor(!!as.name(ClusterIdent), levels = GroupLabels),
             !!as.name(Category_Ident) := factor(!!as.name(Category_Ident), levels = levels(Object@meta.data[[Category_Ident]])))  
  }
  
  # Run statistics across comparisons & correct for multiple comparisons
  
  stat.test_df <- Exp_df %>% 
    filter(!!as.name(Category_Ident) %in% my_comparisons[[1]]) %>%
    group_by(!!as.name(ClusterIdent)) %>%
    do({
      w = wilcox.test(Mean ~ !!as.name(Category_Ident), data = ., paired = FALSE, exact = F)
      data.frame(W = w$statistic, p_value = w$p.value)
    }) %>%
    adjust_pvalue(method = "BH", p.col = "p_value") %>%
    add_significance(p.col = "p_value.adj") %>%
    mutate(p_value = format(round(p_value, 2), nsmall = 2))  %>%
    mutate(!!as.name(Category_Ident) :=  my_comparisons[[1]][2],
           group1 = my_comparisons[[1]][1],
           group2 = my_comparisons[[1]][2],
           yposcol = 2,
           p_final = ifelse(p_value < 0.05, as.character(p_value.adj.signif),
                            ifelse(p_value < 0.061, as.numeric(p_value),"ns"))) # inclue value of p_final if approaching significance
  
  if(isTRUE(clean_stats)){
    stat.test_df <- stat.test_df %>% filter(p_final %!in% "ns")
  }
  
  # Generate box plot
  
  if(is.null(Cat_List)){
    Cat_List <- levels(Object@meta.data[[Category_Ident]])
  }
  
  plot <- Exp_df %>%
    ggplot(aes(x=!!as.name(Category_Ident), 
               y=Mean, 
               fill=!!as.name(Category_Ident))) +
    geom_jitter(color="black", 
                size=0.4,
                width = jitter_width) +
    geom_boxplot(outlier.shape=NA) +
    facet_wrap(vars(!!as.name(ClusterIdent)), 
               scales ="free_x",nrow=1, 
               strip.position = "bottom") +
    scale_fill_manual(values = colors, 
                      labels = Cat_List) +
    theme_bw() +
    theme(axis.text.x = element_blank(),
          axis.text.y = element_text(size=font),
          axis.title.y = element_text(size=font),
          panel.grid.major  = element_blank(),
          panel.grid.minor = element_blank(),
          axis.line.x.bottom = element_line(color="black"),
          axis.line.y.left =  element_line(color="black"),
          panel.border =  element_blank(),
          strip.background = element_blank(),
          strip.placement = "outside",
          strip.text = element_text(size=font),
          legend.text=element_text(size=font),
          legend.title=element_text(size=(font)),
          axis.line.x= element_line(),
          panel.spacing = unit(0, "lines"),
          strip.text.x = element_text(angle = xlab_angle, hjust = 1),
          axis.ticks.x = element_blank()) +
    ggtitle("") +
    xlab("") + 
    ylab(paste(yTitle)) +
    scale_y_continuous(expand = c(0, 0),limits=c(ylo, yhi)) + 
    guides(fill=guide_legend(title=paste(LegendTitle))) + 
    stat_pvalue_manual(stat.test_df, 
                       label = "p_final", 
                       y.position = ypos, 
                       step.increase = stepval, 
                       label.size = (4), 
                       tip.length = bracket_length) 
  
  return(plot)
}


