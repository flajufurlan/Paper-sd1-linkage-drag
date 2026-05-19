
# Packages

library(ggplot2)
library(dplyr)
library(readxl)
library(stringr)
library(rstatix)
library(ggpubr)

# Read data

data <- read_excel("New_Haps.xlsx")
names(data) <- trimws(names(data))

# Convert Yield to t/ha

data <- data %>%
  mutate(Yield_tha = ALL_LA_AR_Yield * 0.00112085)


# Select region

region_col  <- "34.22—38.11"
facet_title <- "Sub-region 2 (34.22–38.11 Mb)"

df <- data %>%
  transmute(
    Yield_tha = Yield_tha,
    Haplotype = .data[[region_col]]
  ) %>%
  filter(!is.na(Haplotype) & Haplotype != "") %>%
  mutate(
    Haplotype = factor(Haplotype),
    Region    = facet_title
  )


# Tukey HSD

tukey <- df %>%
  tukey_hsd(Yield_tha ~ Haplotype) %>%
  arrange(p.adj) %>%
  mutate(
    y.position = max(df$Yield_tha) +
      seq(0.15, by = 0.15, length.out = n()),
    p.adj = sprintf("%.1e", p.adj)
  )


# Plot

p <- ggplot(df, aes(x = Haplotype, y = Yield_tha)) +
  geom_boxplot(
    fill = NA,
    color = "red",
    linewidth = 1,
    outlier.colour = "black"
  ) +
  facet_wrap(~Region) +
  theme_bw(base_size = 10) +
  theme(
    panel.grid.major = element_line(color = "grey85", linewidth = 0.4),
    panel.grid.minor = element_line(color = "grey92", linewidth = 0.25),
    
    strip.background = element_rect(fill = "grey85", color = "grey35"),
    strip.text = element_text(size = 10),
    
    axis.text.x  = element_text(size = 10, angle = 45, hjust = 1),
    axis.text.y  = element_text(size = 10),
    axis.title.x = element_text(size = 10),
    axis.title.y = element_text(size = 10)
  ) +
  labs(
    x = "Haplotype",
    y = expression(Yield~(t~ha^{-1}))
  ) +
  stat_pvalue_manual(
    tukey,
    label = "p.adj",
    tip.length = 0.01,
    size = 6 / .pt, vjust = 0.3
  )


print(p)

# Save 
ggsave(
  filename = "Subregion2_Haplotype_Yield.tiff",
  plot     = p,
  width    = 3.15,   # 80 mm
  height   = 3.4,
  units    = "in",
  device   = "jpeg",
  dpi      = 300,
)
``

