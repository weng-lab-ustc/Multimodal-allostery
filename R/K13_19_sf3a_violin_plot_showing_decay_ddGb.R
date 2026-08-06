library(data.table)
library(krasddpcams)

anno <- fread("C:/Users/36146/OneDrive - USTC/DryLab/base_information_for_K13_K19_project/anno_final_for_8.csv")

rects_sheet <- data.frame(xstart = c(3,38,51,77,109,139),
                          xend = c(9,44,57,84,115,143),
                          col = c("b1","b2","b3","b4","b5","b6"))  

#rects<-data.frame(xstart=c(2,37,49,77,111,141),
#                  xend=c(10,46,58,83,116,143), 
#                  col=c("b1","b2","b3","b4","b5","b6"))


######################   use previously weng data
ddG_file <- "C:/Users/36146/OneDrive/桌面/z7z8/KRAS_nature_fuxian/DATA/weights_Binding_RAF.txt"


ddG <- fread(ddG_file)
ddG[, `:=`(Pos = Pos_ref + 1)]  
data_plot_mutation <- merge(ddG, anno, by = "Pos", all = T)

rects_dt <- as.data.table(rects_sheet)

data_plot_mutation[, `:=`(colors_type = "others")]  
data_plot_mutation[Pos >= rects_dt[col == "b1", xstart] & 
                     Pos <= rects_dt[col == "b1", xend], `:=`(colors_type = "b1")]
data_plot_mutation[Pos >= rects_dt[col == "b2", xstart] & 
                     Pos <= rects_dt[col == "b2", xend], `:=`(colors_type = "b2")]
data_plot_mutation[Pos >= rects_dt[col == "b3", xstart] & 
                     Pos <= rects_dt[col == "b3", xend], `:=`(colors_type = "b3")]
data_plot_mutation[Pos >= rects_dt[col == "b4", xstart] & 
                     Pos <= rects_dt[col == "b4", xend], `:=`(colors_type = "b4")]
data_plot_mutation[Pos >= rects_dt[col == "b5", xstart] & 
                     Pos <= rects_dt[col == "b5", xend], `:=`(colors_type = "b5")]
data_plot_mutation[Pos >= rects_dt[col == "b6", xstart] & 
                     Pos <= rects_dt[col == "b6", xend], `:=`(colors_type = "b6")]

data_plot_mutation[, `:=`(colors_type2 = "others")]
data_plot_mutation[colors_type != "others", `:=`(colors_type2 = "beta sheet")]

data_plot_mutation_beta <- data_plot_mutation[colors_type2 == "beta sheet", ]
data_plot_mutation_beta$colors_type <- factor(data_plot_mutation_beta$colors_type, 
                                              levels = c("b2", "b3", "b1", "b4", "b5", "b6")) 

ggplot2::ggplot(data_plot_mutation_beta, 
                ggplot2::aes(x = colors_type, y = `mean_kcal/mol`)) + 
  ggplot2::geom_violin() + 
  ggplot2::geom_jitter(size = 1 * 0.35, height = 0) + 
  ggplot2::ylab("Binding free energy change \n(kcal/mol)") + 
  ggplot2::xlab("beta sheet") + 
  ggplot2::theme_classic() + 
  ggplot2::theme(text = ggplot2::element_text(size = 10), 
                 axis.text = ggplot2::element_text(size = 10), 
                 axis.text.x = ggplot2::element_text(angle = 90, hjust = 1, vjust = 0.5),  
                 legend.text = ggplot2::element_text(size = 10))


ggplot2::ggsave("C:/Users/36146/OneDrive - USTC/Manuscripts/K13_K19/figures/20260521_start_Updating/result3_figure3/figure3e_ddG_betasheet weng RAF1 energy data 2.pdf", device = cairo_pdf,height = 4,width=4)





######################   use new RAF1 energy data
ddG_file <- "C:/Users/36146/OneDrive - USTC/DryLab/MoCHI_8binders_l2_e6_RA_old_new_merge_at_mochi_20260104_lr_0.025_2048/task_901/weights/weights_Binding_RAF.txt"


ddG <- fread(ddG_file)
ddG[, `:=`(Pos = Pos_ref + 1)]  
data_plot_mutation <- merge(ddG, anno, by = "Pos", all = T)

rects_dt <- as.data.table(rects_sheet)

data_plot_mutation[, `:=`(colors_type = "others")]  
data_plot_mutation[Pos >= rects_dt[col == "b1", xstart] & 
                     Pos <= rects_dt[col == "b1", xend], `:=`(colors_type = "b1")]
data_plot_mutation[Pos >= rects_dt[col == "b2", xstart] & 
                     Pos <= rects_dt[col == "b2", xend], `:=`(colors_type = "b2")]
data_plot_mutation[Pos >= rects_dt[col == "b3", xstart] & 
                     Pos <= rects_dt[col == "b3", xend], `:=`(colors_type = "b3")]
data_plot_mutation[Pos >= rects_dt[col == "b4", xstart] & 
                     Pos <= rects_dt[col == "b4", xend], `:=`(colors_type = "b4")]
data_plot_mutation[Pos >= rects_dt[col == "b5", xstart] & 
                     Pos <= rects_dt[col == "b5", xend], `:=`(colors_type = "b5")]
data_plot_mutation[Pos >= rects_dt[col == "b6", xstart] & 
                     Pos <= rects_dt[col == "b6", xend], `:=`(colors_type = "b6")]

data_plot_mutation[, `:=`(colors_type2 = "others")]
data_plot_mutation[colors_type != "others", `:=`(colors_type2 = "beta sheet")]

data_plot_mutation_beta <- data_plot_mutation[colors_type2 == "beta sheet", ]
data_plot_mutation_beta$colors_type <- factor(data_plot_mutation_beta$colors_type, 
                                              levels = c("b2", "b3", "b1", "b4", "b5", "b6")) 

ggplot2::ggplot(data_plot_mutation_beta, 
                ggplot2::aes(x = colors_type, y = `mean_kcal/mol`)) + 
  ggplot2::geom_violin() + 
  ggplot2::geom_jitter(size = 1 * 0.35, height = 0) + 
  ggplot2::ylab("Binding free energy change(RAF1) \n(kcal/mol)") + 
  ggplot2::xlab("beta sheet") + 
  ggplot2::theme_classic() + 
  ggplot2::theme(text = ggplot2::element_text(size = 10), 
                 axis.text = ggplot2::element_text(size = 10), 
                 axis.text.x = ggplot2::element_text(angle = 90, hjust = 1, vjust = 0.5),  
                 legend.text = ggplot2::element_text(size = 10))


ggplot2::ggsave("C:/Users/36146/OneDrive - USTC/Manuscripts/K13_K19/figures/20260521_start_Updating/result3_figure3/ddG_betasheet decay new RAF1 energy data2.pdf", device = cairo_pdf,height = 4,width=4)

