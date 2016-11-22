#!/usr/bin/env Rscript

#LSBench R statistical utility
#Authors:
#    Roberto Belli <rbelli@inf.ethz.ch>
#    Salvatore Di Girolamo <digirols@inf.ethz.ch>
#    Torsten Hoefler <htor@inf.ethz.ch>

################################################################################
# Helper functions 
################################################################################
#install.packages("data.table")
library(data.table)


Create.Plot.Violin <- function(dataset=NA, limits=NULL, aes.var , aes.box, legend.label, x.axes.label, y.axes.label, ...){
  viol <- ggplot(data=dataset,aes.var)+
    geom_violin()+
    geom_boxplot(notch=TRUE, width=.1, aes.var,inherit.aes=FALSE)+
    #adding summary point (mean)
    stat_summary(fun.y = mean, fun.ymin = min, fun.ymax = max, geom="point", colour = "black", size = 4, shape=18)+
    #facet_grid(. ~ ID, scales="free", space="free")+
    #facet_wrap(~ ID, nrow = 1, scales="free") +
    #scale_color_discrete(name="Synchronization method", labels=label)+
    #scale_shape_discrete(name="Synchronization method",labels=label)+
    scale_fill_discrete(name=legend.label)+
    scale_x_discrete(name=x.axes.label)+
    scale_y_continuous(name=y.axes.label, limits=limits)+ #adding limits here discards points outside range
    #coord_cartesian(ylim=limits)+ #adding limits in that way do not discard anything
    theme_bw(na_theme_size) + 
    theme(legend.position="bottom",legend.title=element_blank(),legend.key.height=na_legend_kew_h, legend.key.width=na_legend_kew_w, legend.text = na_legend_text,legend.background = element_rect(fill="transparent"),axis.ticks = element_blank(), axis.text.x = element_blank())
  return(viol)
}

CreateBoxPlot <- function(dataset, limits=NULL, aes.var , aes.box, x.labels=NULL, legend.label, x.axes.label, y.axes.label, ...){
  box <- ggplot(data=dataset, aes.var)+
    geom_boxplot(notch=TRUE, aes.box, inherit.aes=FALSE,outlier.shape = NA)+
    #adding summary point (mean)
    stat_summary(fun.y = mean, fun.ymin = min, fun.ymax = max, geom="point",  size = 4, shape=18)+
    #facet_grid(. ~ ID, scales="free", space="free")+
    #facet_wrap(as.formula(paste("~",facet.var)), nrow = facet.nrow, ncol = facet.ncol, scales="free") +
    #scale_color_discrete(name="Synchronization method", labels=label)+
    #scale_shape_discrete(name="Synchronization method",labels=label)+
    #scale_fill_discrete(name=legend.label, labels=summary.plot.labels)+
    scale_x_discrete(name=x.axes.label, labels=x.labels)+
    scale_y_continuous(name=y.axes.label)+#, limits=limits)+ #adding limits here discards points outside range
    coord_cartesian(ylim=limits)+ #adding limits in that way do not discard anything
    theme_bw(na_theme_size) + 
    theme(legend.position="none",axis.text.x = element_blank())
  return(box)
}

Create.Plot.HistoDensity <- function(dataset, aes.var, aes.vline, intercept.var, x.limits=NULL, x.axes.label, y.axes.label, ...){
  plotDensity <- ggplot(dataset, aes.var) + 
    geom_histogram(aes(y=..density..),      # Histogram with density instead of count on y-axis
                   binwidth=.1, colour="black", fill="white", freq=TRUE) +
    geom_density(alpha=.2, fill="#FF6666")+  # Overlay with transparent density plot
    geom_vline(aes.vline,   # Ignore NA values for mean
               color="red", linetype="dashed", size=1)+
    scale_x_continuous( name=x.axes.label)+# , limits=x.limits)+
    scale_y_continuous( name=y.axes.label, limits=c(0, NA))+
    coord_cartesian(xlim=x.limits)+ #adding limits in that way do not discard anything
    #xlim(x.limits)+ #adding limits in that way do not discard anything
    theme_bw(na_theme_size)
  return(plotDensity)
}

Create.GGPlot.Density <- function(dataset, aes.var, data.vline, aes.vline, intercept.var, x.limits=NULL, x.axes.label, y.axes.label, ...){
  plotDensity <- ggplot(dataset, aes.var) +
    geom_density(alpha=.2, fill="#FF6666")+  # Overlay with transparent density plot
    geom_vline(data=data.vline, aes.vline,   # Ignore NA values for mean
                linetype="dashed", size=1)+
    scale_x_continuous( name=x.axes.label, limits=x.limits)+ #danger, discards poin
    scale_y_continuous( name=y.axes.label )+
    #xlim(x.limits)+ #adding limits in that way do not discard anything
    theme_bw(na_theme_size)
  return(plotDensity)
}
  
#TODO: fix it.
CreateQQPlot <- function (ds, col) # argument: vector of numbers
{
  #vec <- x[[Time]]
  # following four lines from base R's qqline()
  y <- quantile(ds[[col]], c(0.25, 0.75))
  x <- qnorm(c(0.25, 0.75))
  slope <- diff(y)/diff(x)
  int <- y[1L] - slope * x[1L]
  
  #d <- data.frame(resids = vec)
  
  plot <- ggplot(ds, aes_string(sample=col))+#ggplot(d, aes(sample = resids))+ 
    #facet_wrap(ID~PEs, nrow = 2, scales="free") +
    stat_qq() + 
    geom_abline(slope = slope, intercept = int) +
    scale_x_continuous(name="Theoretical Quantiles")+
    scale_y_continuous(name="Sample Quantiles")+
    #scale_fill_discrete(name="Synchronization method", labels=label)+
    theme_bw(na_theme_size)
    #theme(legend.position="bottom",legend.title=element_blank(),legend.key.height=na_legend_kew_h, legend.key.width=na_legend_kew_w, legend.text = na_legend_text,legend.background = element_rect(fill="transparent"))
  return(plot)
  
}

CreateSummaryPlot <- function(summary, aes.var, legend.label, x.axes.label, y.axes.label, custom.breaks, ...){
  #custom.breaks <- trans_breaks("log2", function(x) 2^x)
  #merging aes (x parameter passed as a string)
  #aes1 <- aes_string(x=x.string)
  #aes2 <- aes(y=median,color=ID,fill=ID, shape=ID, ymin=median-ciStudT, ymax=median+ciStudT)
  #aes.merge <- c(aes1,aes2)
  #class(aes.merge) <- "uneval"
  #creating ggplot
  plot <- ggplot(data=summary , aes.var) +
    #setting points geometry
    geom_point(size=na_geopoint_size)+
    geom_ribbon(alpha=na_ribbon_alpha, show_guide  = F, aes.var)+
    geom_line(size=na_geomline_size)+
    scale_color_discrete(name=legend.label, labels=summary.plot.labels)+
    scale_shape_discrete(name=legend.label,labels=summary.plot.labels)+
    scale_fill_discrete(name=legend.label,labels=summary.plot.labels)+
    # log X scale
    #scale_x_continuous(name=x.axes.label ,trans=log2_trans(), breaks = custom.breaks)+
    scale_x_continuous(name=x.axes.label)+
    scale_y_continuous(name=y.axes.label)+
    #limits=c(0,150),
    guides(color=guide_legend(override.aes=list(fill=NA)))+
    #setting theme black and white an legend position + remove the title of the legend
    theme_bw(na_theme_size) + theme(legend.position=c(na_bench_leg_pos_x,na_bench_leg_pos_y),legend.title=element_blank(),legend.key.height=na_legend_kew_h, legend.key.width=na_legend_kew_w, legend.text = na_legend_text,legend.background = element_rect(fill="transparent"))
  return(plot)
}

CreatePlot <- function(dataset=NA, x_string=NA, ...){
  plot <- ggplot(data=dataset, aes_string(x=Size, y="Time", color="ID",fill="ID", shape="ID", group="ID"))+
    #setting points geometry
    stat_summary(fun.y = mean, fun.ymin = min, fun.ymax = max, geom="point", colour = "black", size = 4, shape=18)+
    stat_smooth(level=.99)+
    #geom_point(size=na_geopoint_size)+
    #geom_line(size=na_geomline_size)+
    scale_color_discrete(name="Synchronization method", labels=summary.plot.labels)+
    scale_shape_discrete(name="Synchronization method",labels=summary.plot.labels)+
    scale_fill_discrete(name="Synchronization method",labels=summary.plot.labels)+
    #breaks set the xtics seq=sequence from 4 to 48 by 4
    #scale_x_continuous(name="Number of Processes",trans=log2_trans(), breaks = trans_breaks("log2", function(x) 2^x))+
    scale_y_continuous(name="Completion Time (us)")+
    guides(color=guide_legend(override.aes=list(fill=NA)))+
    #setting theme black and white an legend position + remove the title of the legend
    theme_bw(na_theme_size) + theme(legend.position=c(na_bench_leg_pos_x,na_bench_leg_pos_y),legend.title=element_blank(),legend.key.height=na_legend_kew_h, legend.key.width=na_legend_kew_w, legend.text = na_legend_text,legend.background = element_rect(fill="transparent"))
  return(plot)
}

PrintGGPlotOnPDF <- function(plot, file.name, w=na_plot_width, h=na_plot_height){
  pdf(file.name, height=h, width=w)
  print(plot)
  dev.off()
}

PrintGGPlotGrobOnPDF <- function(grob, file.name){
  pdf(file.name, height=na_plot_height, width=na_plot_width)
  grid.draw(grob)
  dev.off()
}

PrintQRegressionOnPDF <- function(regr, file.name){
  pdf(file.name, height=na_plot_height, width=na_plot_width)
  plot(regr)
  dev.off()
}
#######################################################################################
# Merge Tables
#######################################################################################
twoWayMerge <- function(table1, table2, names){
  table <- rbind(table1,table2)
  table$ID <- rep(names, c(nrow(table1),nrow(table2)))
  return(table)
}

#
# Merges a list of data frames in a single table.
# it adds a column for distinguish rows coming from
# different tables
# in order to specify distinguishing values names(list) or tables.IDs is used
#
nWayMerge <-function(list, tables.IDs=names(list)){
  library(data.table)
  if(length(list)!=length(tables.IDs)){
    stop("table.names has to be NULL (names(list) will be used) or it has to be the same length of list")
  }
  #creating column with names
  for(index in 1:length(list)){
    list[[index]]$ID_str <- tables.IDs[[index]]
  }
  table <- rbindlist(l=list)
  return(table)
}
#######################################################################################
# File readers section
#######################################################################################
#
# Read all the files in a specific directory aggregating rows in the same
# table
#
ReadAllFilesInDir.Aggregate <- function(dir.path=NA, col=NA, pheader=TRUE){
  dir.files <- list.files(path=dir.path)
  table <- NULL
  for(file in dir.files){
    file.path <- paste( dir.path, file , sep="")
    in.table <- read.table( file.path , header=pheader, col.names=col)
    table <- rbind(table, in.table) 
  }
  return(table)
}

#
# Read all the files in a specific directory aggregating rows in the same
# table
# Deletes first N entries from table ( useful for excluding warmup outliers)
#
ReadAllFilesInDir.AggregateDN <- function(dir.path=NA, col=NA, del.num=1){
  dir.files <- list.files(path=dir.path)
  table <- NULL
  for(file in dir.files){
    file.path <- paste( dir.path, file , sep="")
    in.table <- read.table( file.path , header=TRUE, col.names=col)
    in.table <- in.table[-(del.num),]
    rownames(in.table) <- NULL
    table <- rbind(table, in.table) 
  }
  return(table)
}

#
# Read all the files in a specific directory
# returns a list of data.frame 
# names(list) returnsthe file name
#
ReadAllFilesInDir.List <- function(dir.path=NA, col=NA){
  list <- list()
  list.names <- list()
  dir.files <- list.files(path=dir.path)
  table <- NULL
  for(index in 1:length(dir.files)){
    file <- dir.files[[index]]
    file.path <- paste( dir.path, file , sep="")
    in.table <- read.table( file.path , header=TRUE, col.names=col)
    list.names[[index]] <- index #file
    list[[index]] <- in.table 
  }
  names(list) <- list.names
  return(list)
}


#######################################################################################
# Plotting Utilities
#######################################################################################
# Multiple plot function
#
# ggplot objects can be passed in ..., or to plotlist (as a list of ggplot objects)
# - cols:   Number of columns in layout
# - layout: A matrix specifying the layout. If present, 'cols' is ignored.
#
# If the layout is something like matrix(c(1,2,3,3), nrow=2, byrow=TRUE),
# then plot 1 will go in the upper left, 2 will go in the upper right, and
# 3 will go all the way across the bottom.
#
multiplot.ggplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}

multiplot.grob <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    grid.draw(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      pushViewport(viewport(layout.pos.row = matchidx$row, layout.pos.col = matchidx$col) )
      grid.draw(plots[[i]])
      popViewport(1)
    }
  }
}
