library(ggplot2)
library(reshape2)
library(scales)
library(plyr)
library(grid)

na_geopoint_size <- 7
#na_geopoint_size <- 5
na_geomline_size <- 2
#na_geomline_size <- 1
na_ribbon_alpha <- 0.2
na_perc_size <- 6
#na_theme_size <- 16
na_theme_size <- 22
na_plot_height <- 7
na_plot_width <- 10
#na_legend_key_size <- unit(1, "char")
na_legend_key_size <- unit(3, "char")
#na_legend_text <- element_text(size = na_theme_size)
na_legend_text_size <- 22
na_legend_text <- element_text(size = na_legend_text_size)
na_legend_kew_h <- unit(3, "char")
na_legend_kew_w <- unit(3, "char")
#na_legend_kew_h <- unit(2, "char")
#na_legend_kew_w <- unit(2, "char")
na_bench_leg_pos_x <- 0.25
na_bench_leg_pos_y <- 0.85
na_legend_position <- c(na_bench_leg_pos_x,na_bench_leg_pos_y)
na_viewport_border_size <- 2


na_mbench_plot_x_breaks <- c(8,64,512,4096,32768,262144) 
  #trans_breaks("log2", function(x) 2^x, n=12)
na_mbench_plot_x_labels <- c("8","64","512","4096","32768","262144") 
  #trans_format("log2", math_format(2^.x))
na_mbench_plot_x_trans <- log2_trans()
na_mbench_plot_x_name <- "Number of Transferred Bytes"


na_mbench_plot_y_breaks <- c(0.1,0.3 ,1, 3,10, 30)
#  trans_breaks("log10", function(x) 10^x)
na_mbench_plot_y_labels <-  c(0.1,0.3 ,1, 3,10, 30)
  #trans_format("log10", math_format(10^.x))
na_mbench_plot_y_trans <- log10_trans()                                       
na_mbench_plot_y_name <- "Latency (us)" ; 
#expression(paste("Latency (", "mu" ,"s)" )) not editable with illustrator


na_mbench_subplot_x_limits <- c(8,128)
na_mbench_subplot_x_breaks <- c(8,16,32,64,128)

na_mbech_viewport_width <- 0.4
na_mbech_viewport_height <- 0.35


na_geom_errorbar_width <- 0.2
na_geom_errorbar_lwd <- 0.7

na_plot_margin <- unit(c(0.1,0.5,0.3,0.5), "cm")
na_theme_text_size <- 20
na_theme_x_vjust <- -0.4
na_theme_y_vjust <- 1.2
na_theme_legend_key_height <- unit(1.5,"line")

