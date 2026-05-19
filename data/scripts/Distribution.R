
#Packages
library(readxl)
library(dplyr)
library(ggplot2)
library(ggrepel)
library(patchwork)


# Read data
df <- read.csv("Blues_Within_Across_Enviroments.csv") %>%
  rename(
    yield_kg_ha = `Mean Yield - Across LA & AR`,
    height_cm   = `Mean Plant Height - Across LA & AR`
  ) %>%
  mutate(
    yield_t_ha = yield_kg_ha / 1000
  )


# Cultivar mapping

cultivar_map <- c(
  "CPRS_MP"  = "Cypress",
  "CL111_MP" = "CL111",
  "CL153_MP" = "CL153",
  "CTHL_MP"  = "Catahoula",
  "CL172_MP" = "CL172",
  "PSDO_MP"  = "Presidio",
  "ROYJ_MP"  = "RoyJ",
  "LKST_MP"  = "Lakast"
)

label_df <- df %>%
  mutate(Cultivar = cultivar_map[Genotype]) %>%
  filter(!is.na(Cultivar))



dens_y <- function(x, v) {
  d <- density(v)
  approx(d$x, d$y, xout = x)$y * length(v)
}

label_height <- label_df %>%
  mutate(y_lab = dens_y(height_cm, df$height_cm) * 1.10)

label_yield <- label_df %>%
  mutate(y_lab = dens_y(yield_t_ha, df$yield_t_ha) * 1.10)


# Panel A — Plant Height

p_height <- ggplot(df, aes(height_cm)) +
  geom_histogram(
    binwidth = 2,
    fill = "white",
    color = "black",
    linewidth = 0.3
  ) +
  geom_density(
    aes(y = after_stat(count)),
    linewidth = 0.35,
    color = "grey60",
    alpha = 0.6
  ) +
  geom_text_repel(
    data = label_height,
    aes(x = height_cm, y = y_lab, label = Cultivar),
    angle = 45,
    size = 2.2,
    box.padding = 0.25,
    point.padding = 0.15,
    force = 1.2,
    min.segment.length = Inf,
    segment.color = NA,
    max.overlaps = Inf
  ) +
  annotate(
    "text",
    x = -Inf, y = Inf, label = "A",
    fontface = "bold", size = 3,
    hjust = -0.15, vjust = 1.1
  ) +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.08))) +
  coord_cartesian(clip = "off") +
  labs(x = "Plant Height (cm)", y = "Count") +
  theme_classic(base_size = 10)

# Panel B — Yield

p_yield <- ggplot(df, aes(yield_t_ha)) +
  geom_histogram(
    binwidth = 0.3,
    fill = "white",
    color = "black",
    linewidth = 0.3
  ) +
  geom_density(
    aes(y = after_stat(count)),
    linewidth = 0.35,
    color = "grey60",
    alpha = 0.6
  ) +
  geom_text_repel(
    data = label_yield,
    aes(x = yield_t_ha, y = y_lab, label = Cultivar),
    angle = 45,
    size = 2.2,
    box.padding = 0.25,
    point.padding = 0.15,
    force = 1.2,
    min.segment.length = Inf,
    segment.color = NA,
    max.overlaps = Inf
  ) +
  annotate(
    "text",
    x = -Inf, y = Inf, label = "B",
    fontface = "bold", size = 3,
    hjust = -0.15, vjust = 1.1
  ) +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.08))) +
  coord_cartesian(clip = "off") +
  labs(
    x = expression("Yield (t ha"^{-1}*")"),
    y = "Count"
  ) +
  theme_classic(base_size = 10)


# Combine and save

final_plot <- p_height + p_yield + plot_layout(ncol = 2)

tiff(
  "Distributions_AllMainCultivars_300dpi.tiff",
  width = 7,
  height = 3.5,
  units = "in",
  res = 300,
  compression = "lzw"
)

final_plot
dev.off()
