library(data.table)
library(ggVennDiagram)
library(ggplot2)

RAF1_VS_K27 <- fread("C:/Users/36146/OneDrive - USTC/Manuscripts/K13_K19/figures/20260521_start_Updating/result4_figure5_2/RAF1_vs_K27_anticorrelated_mutations.csv")
RAF1_VS_K19 <- fread("C:/Users/36146/OneDrive - USTC/Manuscripts/K13_K19/figures/20260521_start_Updating/result4_figure5_2/RAF1_vs_K19_anticorrelated_mutations.csv")
K27_VS_K19 <- fread("C:/Users/36146/OneDrive - USTC/Manuscripts/K13_K19/figures/20260521_start_Updating/result4_figure5_2/K27_vs_K19_anticorrelated_mutations.csv")

# Extract the 'mutation' column from each data frame.
mutations_RAF1_VS_K27 <- RAF1_VS_K27$mutation
mutations_RAF1_VS_K19 <- RAF1_VS_K19$mutation
mutations_K27_VS_K19 <- K27_VS_K19$mutation

# Create list
mutation_list <- list(
  RAF1_VS_K27 = mutations_RAF1_VS_K27,
  RAF1_VS_K19 = mutations_RAF1_VS_K19,
  K27_VS_K19 = mutations_K27_VS_K19
)

# Draw a Venn diagram and save it as a PDF.
venn_plot <- ggVennDiagram(mutation_list, 
                           label_alpha = 0,
                           edge_size = 0,  
                           label_size = 5, 
                           set_size = 5) +  
  scale_fill_gradient(low = "white", high = "#75C2F6") +
  scale_color_manual(values = rep("transparent", 3)) + 
  theme_void() +  
  theme(legend.position = "right",
        #plot.title = element_text(hjust = 0.5, size = 16),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 12)) 
#labs(title = "Mutation Venn Diagram")

# 
print(venn_plot)

# 
ggsave("C:/Users/36146/OneDrive - USTC/Manuscripts/K13_K19/figures/20260521_start_Updating/result4_figure5_2/RAF1 vs K27_K27 vs K19_RAF1 vs K19_anticorrelated muts_overlap_vennplot.pdf", 
       plot = venn_plot,
       width = 5.5,   
       height = 5,  
       dpi = 300,
       units = "in",
       device = "pdf")

