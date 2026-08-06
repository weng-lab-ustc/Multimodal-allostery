library(data.table)
library(ggVennDiagram)
library(ggplot2)

RAF1_VS_K27 <- fread("C:/Users/36146/OneDrive - USTC/Manuscripts/K13_K19/figures/20260521_start_Updating/result4_figure5_2/RAF1_vs_K27_anticorrelated_mutations.csv")
RAF1_VS_K19 <- fread("C:/Users/36146/OneDrive - USTC/Manuscripts/K13_K19/figures/20260521_start_Updating/result4_figure5_2/RAF1_vs_K19_anticorrelated_mutations.csv")
K27_VS_K19 <- fread("C:/Users/36146/OneDrive - USTC/Manuscripts/K13_K19/figures/20260521_start_Updating/result4_figure5_2/K27_vs_K19_anticorrelated_mutations.csv")

# 提取每个数据框的 mutation 列
mutations_RAF1_VS_K27 <- RAF1_VS_K27$mutation
mutations_RAF1_VS_K19 <- RAF1_VS_K19$mutation
mutations_K27_VS_K19 <- K27_VS_K19$mutation

# 创建列表
mutation_list <- list(
  RAF1_VS_K27 = mutations_RAF1_VS_K27,
  RAF1_VS_K19 = mutations_RAF1_VS_K19,
  K27_VS_K19 = mutations_K27_VS_K19
)

# 绘制 Venn 图并保存为 PDF
venn_plot <- ggVennDiagram(mutation_list, 
                           label_alpha = 0,
                           edge_size = 0,  # 去掉圆圈黑色线 (设置为0)
                           label_size = 5,  # 增大数字标签大小
                           set_size = 5) +  # 增大集合名称大小
  scale_fill_gradient(low = "white", high = "#75C2F6") +
  scale_color_manual(values = rep("transparent", 3)) +  # 确保圆圈边框透明
  theme_void() +  # 使用简洁主题，移除多余背景
  theme(legend.position = "right",
        #plot.title = element_text(hjust = 0.5, size = 16),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 12)) 
#labs(title = "Mutation Venn Diagram")

# 查看绘图窗口尺寸
print(venn_plot)

# 保存为 PDF（设置足够大的尺寸）
ggsave("C:/Users/36146/OneDrive - USTC/Manuscripts/K13_K19/figures/20260521_start_Updating/result4_figure5_2/RAF1 vs K27_K27 vs K19_RAF1 vs K19_anticorrelated muts_overlap_vennplot.pdf", 
       plot = venn_plot,
       width = 5.5,   # 宽度（英寸）
       height = 5,   # 高度（英寸）
       dpi = 300,
       units = "in",
       device = "pdf")

