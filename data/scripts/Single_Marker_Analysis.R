#Packages


library(tidyverse)
library(readr)
library(janitor)


# Load data

df <- read_csv("Single_Marker_Analysis_Results.csv", show_col_types = FALSE) %>%
  clean_names() %>%
  mutate(
    marker  = as.numeric(marker),
    f_ratio = as.numeric(f_ratio),
    pos_mb  = marker / 1e6
  ) %>%
  filter(!is.na(f_ratio), !is.na(pos_mb)) %>%
  arrange(pos_mb)

# Identify max F-ratio SNP(s)

maxF <- max(df$f_ratio, na.rm = TRUE)

df <- df %>%
  mutate(
    is_max = f_ratio == maxF,
    snp_id = factor(marker, levels = marker),
    mb_lab = sprintf("%.2f", pos_mb)
  )


annot_df <- df %>%
  filter(is_max) %>%
  summarise(
    snp_id = snp_id[ceiling(n() / 2)],  # center over red bars
    y = max(f_ratio) * 1.02,
    label = "Yield QTL:\n36.69 Mb"
  )


# Plot

p <- ggplot(df, aes(x = snp_id, y = f_ratio)) +
  geom_col(
    aes(color = is_max),
    fill = "white",
    width = 0.9,
    linewidth = 0.8
  ) +
  scale_color_manual(
    values = c(`FALSE` = "black", `TRUE` = "red"),
    guide = "none"
  ) +
  scale_x_discrete(labels = df$mb_lab) +
  geom_text(
    data = annot_df,
    aes(x = snp_id, y = y, label = label),
    inherit.aes = FALSE,
    color = "black",
    size = 3.5
    ,
    vjust = 0
  ) +
  labs(
    x = "SNP position (Mb)",
    y = "F-ratio"
  ) +
  theme_classic(base_size = 16) +
  theme(
    axis.text.x = element_text(
      angle = 90, vjust = 0.5, hjust = 1, size = 12
    ),
    plot.margin = margin(10, 10, 25, 10)
  ) +
  coord_cartesian(clip = "off")

print(p)


# Save TIFF (300 dpi)
ggsave(
  filename = "F_ratio_barplot_QTL.tiff",
  plot = p,
  width = 9,
  height = 7,
  units = "in",
  dpi = 300,
  compression = "lzw"
)


# Results table

results_table <- df %>%
  select(
    marker,        # SNP position (bp)
    pos_mb,        # SNP position (Mb)
    f_ratio,
    is_max
  ) %>%
  arrange(pos_mb)
``
write_csv(
  results_table,
  "Single_Marker_Analysis_Results_with_Fratio.csv"
)
names(df)
