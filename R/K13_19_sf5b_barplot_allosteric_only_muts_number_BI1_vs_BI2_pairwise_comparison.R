# =========================================================
# Create a bar chart showing the counts of "allosteric-only" instances for BI1 vs. BI2.
# =========================================================

library(ggplot2)
library(dplyr)
library(tidyr)
library(data.table)

allosteric_number_pairwise_comparison <- fread("C:/Users/36146/OneDrive - USTC/Manuscripts/K13_K19/figures/20260521_start_Updating/result4_figure5/mutation_classification_summary_8binder.csv")

names(allosteric_number_pairwise_comparison)

BI1 <- c("RAF1", "RALGDS", "PI3KCG", "SOS1", "K55", "K27")
BI2 <- c("K13", "K19")



# Define the pairing order to be preserved.
pair_order <- c(
  "RAF1 vs K13", "RAF1 vs K19",
  "RALGDS vs K13", "RALGDS vs K19",
  "PI3KCG vs K13", "PI3KCG vs K19",
  "SOS1 vs K13", "SOS1 vs K19",
  "K55 vs K13", "K55 vs K19",
  "K27 vs K13", "K27 vs K19"
)

# Filter data
filtered_data <- allosteric_number_pairwise_comparison %>%
  filter(Pair %in% pair_order) %>%
  mutate(Pair = factor(Pair, levels = pair_order)) %>%
  select(Pair, Allosteric_only_in_X, Allosteric_only_in_Y)

# Convert the data format to long format to facilitate plotting.
plot_data <- filtered_data %>%
  pivot_longer(cols = c(Allosteric_only_in_X, Allosteric_only_in_Y),
               names_to = "Type",
               values_to = "Count") %>%
  mutate(Type = factor(Type, 
                       levels = c("Allosteric_only_in_X", "Allosteric_only_in_Y"),
                       labels = c("Allosteric only in X", "Allosteric only in Y")))

# Define color
custom_colors <- c("Allosteric only in X" = "#1B38A6", 
                   "Allosteric only in Y" = "#75C2F6")

# Create a bar chart
p <- ggplot(plot_data, aes(x = Pair, y = Count, fill = Type)) +
  geom_bar(stat = "identity", position = position_dodge(0.7), width = 0.6, 
           color = "white", size = 0.3) + 
  scale_fill_manual(values = custom_colors) +
  labs(title = "Allosteric Mutations Exclusive to Each Binder Pair",
       x = "Binder Pair",
       y = "Number of Allosteric Mutations",
       fill = "Type") +
  theme_classic() +  
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
        axis.text.y = element_text(size = 10),
        axis.title = element_text(size = 12),
        axis.line = element_line(color = "black", size = 0.5),  
        axis.ticks = element_line(color = "black", size = 0.5),  
        plot.title = element_text(size = 14, hjust = 0.5),
        legend.position = "top",
        legend.title = element_text(size = 11),
        legend.text = element_text(size = 10),
        panel.background = element_rect(fill = "white"), 
        plot.background = element_rect(fill = "white"))  

# 
print(p)

# 
output_file <- "C:/Users/36146/OneDrive - USTC/Manuscripts/K13_K19/figures/20260521_start_Updating/result4_figure5_2/allosteric_only_barplot in pairwise comparison between BI1 and BI2 2.pdf"
ggsave(output_file, plot = p, width = 5, height = 3, dpi = 300)

cat("\n========== Plot Complete ==========\n")
cat("Plot saved to:", output_file, "\n")
