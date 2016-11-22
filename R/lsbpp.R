#!/usr/bin/env Rscript

# LSBench R statistical utility
#
# Copyright (c) 2015 ETH-Zurich. All rights reserved.
# Use of this source code is governed by a BSD-style license that can be
# found in the LICENSE file.
#
# Authors:
#    Salvatore Di Girolamo <digirols@inf.ethz.ch>
#    Roberto Belli <rbelli@inf.ethz.ch>
#    Torsten Hoefler <htor@inf.ethz.ch>

library(ggplot2)
library(docopt)

source("utils.R")
source("aes.R")
source("stats.R")

#Load the data from the 'dir' folder applying the filter specified by 'expr'
loaddata <- function(dir, expr){
    data_unfolded <- ReadAllFilesInDir.Aggregate(dir.path=dir)
    if (!is.null(expr)) data_unfolded = subset(data_unfolded, eval(parse(text=expr)))
    return(data_unfolded)
}

#Shapiro-Wilk normality test
lsb_shapiro <- function(data, col){
    sampled <- as.numeric(data[[col]])
   
    #shapiro.test takes min 3 val
    if (length(sampled)<3) {
        print("The dataset contains too few rows (min. 3 are required)")
        return()
    }

    #shapiro.test takes max 5000 vals
    if (length(sampled)>5000) sampled <- sample(sampled, 5000)
    shapiro.test(sampled)
}

#ANOVA test
lsb_anova <- function(data, xcol, ycol){
    res <- aov(eval(parse(text=ycol)) ~ eval(parse(text=xcol)), data=data)
    print(res)
    summary(res)
}

#Kruskal–Wallis one-way analysis of variance
lsb_kruskal <- function(data, xcol, ycol){
    kruskal.test(eval(parse(text=ycol))~eval(parse(text=xcol)), data=data)
}

#Histogram + density plot
lsb_histo_density <- function(data, col, xlbl){
    str <- paste("mean(", col, ")")
    plot <- Create.Plot.HistoDensity(data, aes_string(col), aes_string(xintercept=str),,, xlbl, "Density")
    return(plot)
}

#Density plot
lsb_density <- function(data, col, xlbl){
  str <- paste("mean(", col, ")")
  plot <- Create.GGPlot.Density(data, aes_string(col), data, aes_string(xintercept=str),,, xlbl, "Density")
  return(plot)
}

#QQ-plot
lsb_qqplot <- function(data, col){
  plot <- CreateQQPlot(data, col)
  return(plot)
}

#Box Plot
lsb_boxplot <- function(d, cols, xlbl, ylbl){
  plot <- CreateBoxPlot(d,, cols , cols,, "", xlbl, ylbl)
  return(plot)
}

#Violin Plot
lsb_violinplot <- function(d, cols, xlbl, ylbl){
  plot <- Create.Plot.Violin(d,, cols , cols, "", xlbl, ylbl)
  return(plot)
}

#log-normalization on data[datacol]
lsb_lognormalize <- function(data, datacol){
  val <- data[[datacol]]
  val <- log(val)
  val <- data.frame(matrix(unlist(val),ncol=1) )
  data[]
  names(val) <- datacol
  return(val)
}

#k-normalization on data[datacol]
lsb_knormalize <- function(data, datacol, k=1){
  var <- data[[datacol]]
  max <- k
  x <- seq_along(var)
  d1 <- split(var, ceiling(x/max))
  df <- numeric(length(d1))
  i <- 1
  for(chunk in d1){
    len <- length(chunk) 
    s<-sum(chunk)
    df[i] <- s / len
    i = i +1
  }
  df <- data.frame(matrix(unlist(df),ncol=1) )
  #data[[datacol]] <- unlist(df)
  names(df) <- datacol
  return(df)
}

#Quantile Regression Plot between data1[col1] and data2[col2] 
lsb_qr <- function(data1, col1, data2, col2){
    library(quantreg)   
    data_merged = data.frame(data1[[col1]], data2[[col2]])
    colnames(data_merged) <- c("a", "b")
    pdf(file="out.pdf")
    plot(summary(rq(a~b,tau = 1:49/50,data=data_merged)))
}


#### Main 

#Parse arguments
"Usage:
 lsbpp.R <datadir> [--filter=<expr>] [--lognormalize|--knormalize=<k>] [(--summarize <groupvars> <meauservar> [--ci=<ci>] [--qi=<qi>])] (shapiro) <datacol>
 lsbpp.R <datadir> [--filter=<expr>] [--lognormalize|--knormalize=<k>] [(--summarize <groupvars> <meauservar> [--ci=<ci>] [--qi=<qi>])] (qqplot) <datacol> [--plotout=<path>]
 lsbpp.R <datadir> [--filter=<expr>] [--lognormalize|--knormalize=<k>] [(--summarize <groupvars> <meauservar> [--ci=<ci>] [--qi=<qi>])] (histodensity|density) <datacol>  [--xlabel=<string>] [--plotout=<path>]
 lsbpp.R <datadir> [--filter=<expr>] [(--summarize <groupvars> <meauservar> [--ci=<ci>] [--qi=<qi>])] (anova|kruskal) <xcol> <ycol>
 lsbpp.R <datadir> [--filter=<expr>] [(--summarize <groupvars> <meauservar> [--ci=<ci>] [--qi=<qi>])] (boxplot|violinplot) <xcol> <ycol> [--xlabel=<string>] [--ylabel=<string>] [--plotout=<path>]
 lsbpp.R <datadir> [--filter=<expr>] <data_cmp_dir> [--filter_cmp=<expr>] qr <datacol> <datacol_cmp> [--plotout=<path>] 

Options:
 --xlabel=<string> x-axis label [default: x]
 --ylabel=<string> y-axis label [default: y]
 --plotout=<path> Path where to write the plot [default: out.pdf]
 --ci=<ci> Confidence Interval [default: 0.95]
 --qi=<qi> Quantile Interval [default: 0.95]
" -> doc

opts <- docopt(doc)

#print(opts["--summarize"])
#print(opts[["datadir"]])
#print(opts[["datacol"]])
#print(opts[["--filter"]])
#print(opts[["shapiro"]])
#print(opts[["--lognormalize"]])
#print(opts[["--knormalize"]])
#print(opts[["--summarize"]])


datacol = opts[["datacol"]]
xcol = opts[["xcol"]]
ycol = opts[["ycol"]]
xlbl = opts[["--xlabel"]]
ylbl = opts[["--ylabel"]]
plotout = opts[["--plotout"]]

#print("Reading data...")
#Load the data
data <- loaddata(opts[["datadir"]], opts[["--filter"]])
#print("Data has been read")

#Summarize it if requested

#for (i in 1::opts[["--summarize"]]){

if (opts[["--summarize"]]){
    #print(opts[["meauservar"]])
    #print(opts[["groupvars"]])
    #print(opts[["--ci"]])
    #print(opts[["--qi"]])

    
    #print("Summarize")
    data <- CalculateDataSummary(data=data, measurevar=opts[["meauservar"]], groupvars=unlist(strsplit(opts[["groupvars"]], split=",")), conf.interval=as.numeric(opts[["--ci"]]), quantile.interval=as.numeric(opts[["--qi"]]))
    #print("done!")
}

#Normalize it if requested
if (opts[["--lognormalize"]]) {data <- lsb_lognormalize(data, datacol)}
if (!is.null(opts[["--knormalize"]])) {data <- lsb_knormalize(data, datacol, as.numeric(opts[["--knormalize"]]))}


#Execute the command
plot = NULL
if (opts[["shapiro"]]) {lsb_shapiro(data, datacol)
}else if (opts[["histodensity"]]) {plot <- lsb_histo_density(data, datacol, xlbl)
}else if (opts[["density"]]) {plot <- lsb_density(data, datacol, xlbl)
}else if (opts[["qqplot"]]) { plot <- lsb_qqplot(data, datacol)
}else if (opts[["boxplot"]]) { plot <- lsb_boxplot(data, aes(factor(eval(parse(text=xcol))), eval(parse(text=ycol))), xlbl, ylbl)
}else if (opts[["violinplot"]]) { plot <- lsb_violinplot(data, aes(factor(eval(parse(text=xcol))), eval(parse(text=ycol))), xlbl, ylbl)
}else if (opts[["anova"]]) { lsb_anova(data, xcol, ycol) 
}else if (opts[["kruskal"]]) { lsb_kruskal(data, xcol, ycol) 
}else if (opts[["qr"]]) { 
    data_cmp <- loaddata(opts[["data_cmp_dir"]], opts[["--filter_cmp"]])
    lsb_qr(data, opts[["datacol"]], data_cmp, opts[["datacol_cmp"]]) 
}


#switch(op, 
#    shapiro = lsb_shapiro(data, xcol),
#    anova = lsb_anova(data, xcol, ycol),
#    kruskal = lsb_kruskal(data, xcol, ycol),
#    histodensity = plot<-lsb_histo_density(data, xcol),
#    density = plot<-lsb_density(data, xcol),
#    qqplot = plot <- lsb_qqplot(data, xcol),
#    boxplot = plot <- lsb_boxplot(data, aes(factor(eval(parse(text=xcol))), eval(parse(text=ycol)))),
#    violinplot = plot <- lsb_violinplot(data, aes(factor(eval(parse(text=xcol))), eval(parse(text=ycol))))
#)

#Print the plot (if any)
if (!is.null(plot)){
    f.name <- plotout
    pdf(f.name, height=na_plot_height, width=na_plot_width)
    print(plot)
}
    
unlink(data)


