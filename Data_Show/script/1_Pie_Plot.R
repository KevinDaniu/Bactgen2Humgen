####设置工作空间####
setwd("C:/Users/yuxip/Desktop/LZH/Tool/")
rm(list = ls())
gc()

####导入R包####
library(dplyr)
library(plotly)
library(readr)
library(htmlwidgets)
library(ggsci)
library(ggplot2)

####可视化####
#----导入数据----#
fam <- read_tsv("Data_Show/data/stat_level/family_count.tsv",
                col_names = c("taxon", "count"))
gen <- read_tsv("Data_Show/data/stat_level/genus_count.tsv",
                col_names = c("taxon", "count"))
re_fam <- read_tsv("Data_Show/data/stat_level_rep0.25/family_count.tsv",
                   col_names = c("taxon", "count"))
re_gen <- read_tsv("Data_Show/data/stat_level_rep0.25/genus_count.tsv",
                   col_names = c("taxon", "count"))

#----配色：拼接多个 ggsci 调色板，得到离散不重复的大量颜色----#
big_pal <- unique(c(
  pal_npg("nrc")(10),      # NPG
  pal_lancet("lanonc")(9), # Lancet
  pal_nejm("default")(8),      # NEJM
  pal_jco("default")(10),  # JCO
  pal_d3("category20")(20),
  pal_igv("default")(51),
  pal_uchicago("default")(9),
  pal_aaas("default")(10),
  pal_jama("default")(7),
  pal_bmj("default")(9),
  pal_flatui("default")(10),
  pal_frontiers("default")(10),
  pal_simpsons("springfield")(16),
  pal_startrek("uniform")(7),
  pal_tron("legacy")(7),
  pal_futurama("planetexpress")(12)
))

# 若类别数超过拼接色数，循环复用
pal <- function(n) rep(big_pal, length.out = n)

#----交互式饼图----#
p1 <- plot_ly(fam, labels = ~taxon, values = ~count, type = "pie",
              textinfo = "none", hoverinfo = "label+value+percent",
              marker = list(colors = pal(nrow(fam)),
                            line = list(color = "#FFFFFF", width = 0.5))) |>
  layout(title = "Family distribution")
p2 <- plot_ly(gen, labels = ~taxon, values = ~count, type = "pie",
              textinfo = "none", hoverinfo = "label+value+percent",
              marker = list(colors = pal(nrow(gen)),
                            line = list(color = "#FFFFFF", width = 0.5))) |>
  layout(title = "Genus distribution")

p3 <- plot_ly(re_fam, labels = ~taxon, values = ~count, type = "pie",
              textinfo = "none", hoverinfo = "label+value+percent",
              marker = list(colors = pal(nrow(re_fam)),
                            line = list(color = "#FFFFFF", width = 0.5))) |>
  layout(title = "Family distribution")
p4 <- plot_ly(re_gen, labels = ~taxon, values = ~count, type = "pie",
              textinfo = "none", hoverinfo = "label+value+percent",
              marker = list(colors = pal(nrow(re_gen)),
                            line = list(color = "#FFFFFF", width = 0.5))) |>
  layout(title = "Genus distribution")

saveWidget(p1, "Data_Show/result/family_pie.html", selfcontained = TRUE)
saveWidget(p2, "Data_Show/result/genus_pie.html",   selfcontained = TRUE)
saveWidget(p3, "Data_Show/result/re_family_pie.html",  selfcontained = TRUE)
saveWidget(p4, "Data_Show/result/re_genus_pie.html",  selfcontained = TRUE)

#----饼图----#
####配色：拼接多个 ggsci 调色板，得到大量离散色####
big_pal <- unique(c(
  pal_npg("nrc")(10),      # NPG
  pal_lancet("lanonc")(9), # Lancet
  pal_nejm("default")(8),      # NEJM
  pal_jco("default")(10),  # JCO
  pal_d3("category20")(20),
  pal_igv("default")(51),
  pal_uchicago("default")(9),
  pal_aaas("default")(10),
  pal_jama("default")(7),
  pal_bmj("default")(9),
  pal_flatui("default")(10),
  pal_frontiers("default")(10),
  pal_simpsons("springfield")(16),
  pal_startrek("uniform")(7),
  pal_tron("legacy")(7),
  pal_futurama("planetexpress")(12)
))
pal <- function(n) rep(big_pal, length.out = n)

####CNS 风格饼图：全类别，图例只标 TopN####
cns_pie <- function(df, title, file, legend_top = 15) {
  d <- df %>%
    arrange(desc(count)) %>%
    mutate(
      taxon = factor(taxon, levels = taxon),
      # 只保留前 legend_top 个的图例名，其余设为 NA（不显示图例）
      lab = ifelse(row_number() <= legend_top, as.character(taxon), NA)
    )
  
  p <- ggplot(d, aes(x = 1, y = count, fill = taxon)) +
    geom_col(width = 1, color = "white", linewidth = 0.3) +
    coord_polar(theta = "y") +
    scale_fill_manual(values = pal(nrow(d)),
                      breaks = d$taxon[d$taxon %in% d$lab[!is.na(d$lab)]],
                      labels = d$lab[!is.na(d$lab)],
                      na.value = "grey80") +
    labs(title = title) +
    theme_void(base_family = "sans") +
    theme(
      plot.title = element_text(size = 14, hjust = 0.02),
      legend.position = "right",
      legend.title = element_blank(),
      legend.text = element_text(size = 8),
      legend.key.size = unit(0.35, "cm"),
      plot.margin = margin(10, 10, 10, 10)
    )
  
  ggsave(file, p, width = 12, height = 9, dpi = 300)
}

####生成并保存####
cns_pie(fam,    "Family distribution",           "Data_Show/result/family_pie.pdf")
cns_pie(gen,    "Genus distribution",            "Data_Show/result/genus_pie.pdf")
cns_pie(re_fam, "Family distribution (rep 0.25)","Data_Show/result/re_family_pie.pdf")
cns_pie(re_gen, "Genus distribution (rep 0.25)", "Data_Show/result/re_genus_pie.pdf")
