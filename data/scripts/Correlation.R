#Packages

library(readxl)
library(ggplot2)
library(dplyr)

# Read data
dat <- read.csv("Blues_Within_Across_Enviroments.csv")

# Rename for convenience and convert yield to t ha-1
dat <- dat %>%
  rename(
    Yield_kg_ha = `Mean Yield - Across LA & AR`,
    PlantHeight_cm = `Mean Plant Height - Across LA & AR`
  ) %>%
  mutate(
    Yield_t_ha = Yield_kg_ha / 1000
  )

# Correlation test
cor_test <- cor.test(dat$Yield_t_ha, dat$PlantHeight_cm, method = "pearson")

r_val <- round(cor_test$estimate, 2)
p_val <- signif(cor_test$p.value, 2)

# Build plot
p <- ggplot(dat, aes(x = Yield_t_ha, y = PlantHeight_cm)) +
  geom_point(color = "black", size = 1.6, alpha = 0.8) +
  geom_smooth(method = "lm", se = FALSE, color = "red", linewidth = 0.8) +
  labs(
    x = expression("Yield (t ha"^-1*")"),
    y = "Plant Height (cm)"
  ) +
  annotate(
    "text",
    x = Inf, y = Inf,
    label = paste0("r = ", r_val, "\n", "P = ", p_val),
    hjust = 1.1, vjust = 1.3,
    size = 3.5
  ) +
  theme_bw(base_size = 10) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )

# Save as TIFF 
ggsave(
  filename = "Yield_vs_PlantHeight_corr.tiff",
  plot = p,
  width = 80,        # mm 
  height = 80,       # mm
  units = "mm",
  dpi = 300,
  compression = "lzw"
)
``
