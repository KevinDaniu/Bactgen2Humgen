####设置工作空间####
setwd("C:/Users/yuxip/Desktop/LZH/Tool/")
rm(list = ls())
gc()

####导入R包####
library(dplyr)
library(plotly)
library(readr)
library(htmlwidgets)
library(ggplot2)
library(forcats)

####可视化####
#----导入数据----#
data <- read.csv("Data_Show/data/high_conserved_summary_annotated.csv")
data <- data[,-3]

#----可视化数据----#
colnames(data) <- c("strain", "n_conserved",
                    "species", "family", "genus")

data <- data |>
  mutate(
    species = ifelse(is.na(species) | species == "", "Unclassified", species),
    genus   = ifelse(is.na(genus)   | genus   == "", "Unclassified", genus),
    family  = ifelse(is.na(family)  | family  == "", "Unclassified", family)
  )

# 1. Species level
species_df <- data |>
  group_by(taxon = species) |>
  summarise(
    n_strain = n(),
    total    = sum(n_conserved),
    mean     = mean(n_conserved),
    .groups  = "drop"
  ) |>
  arrange(desc(total)) |>
  slice_head(n = 30) |>
  mutate(taxon = factor(taxon, levels = taxon[order(-total)]))

p1 <- plot_ly(
  species_df,
  x = ~taxon, y = ~total,
  type = "bar",
  text = ~paste0("Strains: ", n_strain,
                 "<br>Total clusters: ", total,
                 "<br>Mean: ", round(mean, 1)),
  hoverinfo = "text",
  marker = list(color = "#4C78A8")
) |>
  layout(
    title = "Top Species by Conserved Cluster Count",
    xaxis = list(title = "", tickangle = -45, categoryorder = "array",
                 categoryarray = ~taxon),
    yaxis = list(title = "Total conserved clusters"),
    margin = list(b = 160)
  )
p1

# 2. Genus level
genus_df <- data |>
  group_by(taxon = genus) |>
  summarise(
    n_strain = n(),
    total    = sum(n_conserved),
    mean     = mean(n_conserved),
    .groups  = "drop"
  ) |>
  arrange(desc(total)) |>
  slice_head(n = 30) |>
  mutate(taxon = factor(taxon, levels = taxon[order(-total)]))

p2 <- plot_ly(
  genus_df,
  x = ~taxon, y = ~total,
  type = "bar",
  text = ~paste0("Strains: ", n_strain,
                 "<br>Total clusters: ", total,
                 "<br>Mean: ", round(mean, 1)),
  hoverinfo = "text",
  marker = list(color = "#54A24B")
) |>
  layout(
    title = "Top Genera by Conserved Cluster Count",
    xaxis = list(title = "", tickangle = -45, categoryorder = "array",
                 categoryarray = ~taxon),
    yaxis = list(title = "Total conserved clusters"),
    margin = list(b = 160)
  )
p2

# 3. Family level
family_df <- data |>
  group_by(taxon = family) |>
  summarise(
    n_strain = n(),
    total    = sum(n_conserved),
    mean     = mean(n_conserved),
    .groups  = "drop"
  ) |>
  arrange(desc(total)) |>
  slice_head(n = 30) |>
  mutate(taxon = factor(taxon, levels = taxon[order(-total)]))

p3 <- plot_ly(
  family_df,
  x = ~taxon, y = ~total,
  type = "bar",
  text = ~paste0("Strains: ", n_strain,
                 "<br>Total clusters: ", total,
                 "<br>Mean: ", round(mean, 1)),
  hoverinfo = "text",
  marker = list(color = "#E45756")
) |>
  layout(
    title = "Top Families by Conserved Cluster Count",
    xaxis = list(title = "", tickangle = -45, categoryorder = "array",
                 categoryarray = ~taxon),
    yaxis = list(title = "Total conserved clusters"),
    margin = list(b = 160)
  )
p3

saveWidget(p1, "Data_Show/result/species_bar.html", selfcontained = TRUE)
saveWidget(p2, "Data_Show/result/genus_bar.html",   selfcontained = TRUE)
saveWidget(p3, "Data_Show/result/family_bar.html",  selfcontained = TRUE)

#----图片绘制----#
# 统一主题
cns_theme <- theme_classic(base_size = 12) +
  theme(
    axis.line        = element_line(linewidth = 0.4, color = "black"),
    axis.ticks       = element_line(linewidth = 0.4, color = "black"),
    axis.ticks.length = unit(-0.08, "cm"),          # 刻度朝内
    axis.text.x      = element_text(angle = 90, hjust = 1, vjust = 0.5,
                                    size = 12, color = "black"),
    axis.text.y      = element_text(size = 12, color = "black"),
    axis.title       = element_text(size = 12, color = "black"),
    plot.title       = element_text(size = 14, face = "bold", hjust = 0,
                                    margin = margin(b = 4)),
    plot.margin      = margin(5, 8, 5, 5),
    legend.position  = "none",
    panel.background = element_rect(fill = "white", color = NA),
    plot.background  = element_rect(fill = "white", color = NA)
  )

# ---- 1. Species level ----
p1_gg <- ggplot(species_df, aes(x = fct_reorder(taxon, -total), y = total)) +
  geom_col(fill = "#4C78A8", width = 0.7) +
  labs(title = "Top Strains by Conserved Cluster Count",
       x = NULL, y = "Total conserved clusters") +
  cns_theme

# ---- 2. Genus level ----
p2_gg <- ggplot(genus_df, aes(x = fct_reorder(taxon, -total), y = total)) +
  geom_col(fill = "#54A24B", width = 0.7) +
  labs(title = "Top Genera by Conserved Cluster Count",
       x = NULL, y = "Total conserved clusters") +
  cns_theme

# ---- 3. Family level ----
p3_gg <- ggplot(family_df, aes(x = fct_reorder(taxon, -total), y = total)) +
  geom_col(fill = "#E45756", width = 0.7) +
  labs(title = "Top Families by Conserved Cluster Count",
       x = NULL, y = "Total conserved clusters") +
  cns_theme

# ---- 展示 ----
p1_gg
p2_gg
p3_gg

# ---- 保存：PDF 矢量 + TIFF 高分辨率（期刊常用格式）----
ggsave("Data_Show/result/species_bar_cns.pdf",  p1_gg,
       width = 7, height = 4, units = "in")
ggsave("Data_Show/result/genus_bar_cns.pdf",    p2_gg,
       width = 7, height = 4, units = "in")
ggsave("Data_Show/result/family_bar_cns.pdf",   p3_gg,
       width = 7, height = 4, units = "in")

