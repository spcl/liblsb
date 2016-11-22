library(ggplot2)

#install.packages("ggplot2")
#install.packages("data.table")


setwd("/home/salvo/ETH/benchmarking/R/")

source("utils.R")
source("stats.R")
source("aes.R")


data.bh_study <-  ReadAllFilesInDir.Aggregate(dir.path="sample_data/", col=c("type", "rank", "csize", "htsize", "memsize", "iter", "mybodies", "id", "time", "overhead"))
re.bh_study <- CalculateDataSummary(data=data.bh_study, measurevar="time", groupvars=c("type", "htsize", "memsize", "csize"), conf.interval=.95, quantile.interval=.95)
re.bh_study.toplot <- subset(re.bh_study, memsize==1048576)

summary.plot.labels <- c("CLaMPI (Adaptive)", "CLaMPI (Fixed)", "Native")
aes.var <- aes(x=htsize/1000, y=median, ymin=CI.NNorm.high, ymax=CI.NNorm.low, color=factor(type), shape=factor(type))
plot_lat <- ggplot(data=re.bh_study.toplot, aes.var) +
  geom_point(size=na_geopoint_size)+  
  geom_line(size=na_geomline_size)+
  geom_errorbar(width = na_geom_errorbar_width, lwd=na_geom_errorbar_width, color="black")+
  scale_color_discrete(name="", labels=summary.plot.labels)+
  scale_shape_discrete(name="",labels=summary.plot.labels)+
  scale_x_continuous(name="Hash Table Size x 10^3", breaks=c(1, 15, 30))+
  scale_y_continuous(name="Force Comp. (us)")+
  theme_bw(na_theme_size) +
  theme(legend.position="bottom", legend.direction = "horizontal", legend.key = element_rect(fill = "transparent", colour = "transparent"), legend.background = element_rect(fill="transparent", colour="transparent")) + 
  theme(plot.margin=na_plot_margin) + theme(text = element_text(size=na_theme_text_size)) + theme(legend.key.height=na_theme_legend_key_height) +
  theme(axis.title.y=element_text(vjust=na_theme_y_vjust)) +
  theme(axis.title.x=element_text(vjust=na_theme_x_vjust))

print(plot_lat)

PrintGGPlotOnPDF(plot_lat, "test.pdf")
