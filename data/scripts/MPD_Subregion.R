#Packages

library(dplyr)
library(ggplot2)
library(ggpubr)
library(rstatix)
library(stringr)
library(readxl)

# ==== Load data
dat <- read_excel("MPD_Analysis.xlsx", sheet = "Sheet1")

# ==== Convert Yield (lbs/A → t/ha) ====
dat <- dat %>%
  mutate(Yield_tha = Yield * 0.00112085)

# ==== Define columns ====
geno_col   <- "Line"
yield_col  <- "Yield_tha"
region_cols <- c(
  "Sub-region 33.26-33.97 Mb",
  "Sub-region 34.27-37.69 Mb"
)

plot_names <- c(
  "Sub-region (33.26–33.97 Mb)",
  "Sub-region (34.27–37.69 Mb)"
)

# ==== Text sizes (ALL = 10) ====
BASE_SIZE  <- 10
TICK_SIZE  <- 10
TITLE_SIZE <- 10
STRIP_SIZE <- 10

# ==== Output directory ====
outdir <- "MPD_boxplots"
dir.create(outdir, showWarnings = FALSE)

# ==== Plot loop ====
plots <- list()

for (i in seq_along(region_cols)) {
  
  region <- region_cols[i]
  facet_title <- plot_names[i]
  
  df_long <- dat %>%
    select(
      !!sym(geno_col),
      !!sym(yield_col),
      !!sym(region)
    ) %>%
    rename(
      Genotype  = !!sym(geno_col),
      Yield_tha = !!sym(yield_col),
      Haplotype = !!sym(region)
    ) %>%
    filter(!is.na(Haplotype) & Haplotype != "") %>%
    mutate(
      Region = facet_title,
      Haplotype = factor(Haplotype)
    )
  
  # Skip if < 2 haplotypes
  if (nlevels(df_long$Haplotype) < 2) next
  
  # ==== Tukey HSD ====
  tukey <- df_long %>%
    tukey_hsd(Yield_tha ~ Haplotype) %>%
    arrange(p.adj)
  
  y_max   <- max(df_long$Yield_tha, na.rm = TRUE)
  y_range <- diff(range(df_long$Yield_tha, na.rm = TRUE))
  
  tukey <- tukey %>%
    mutate(
      y.position = y_max + seq_len(n()) * 0.1 * y_range,
      p.label = formatC(as.numeric(p.adj),
                        format = "e", digits = 1)
    )
  
  # ==== Plot ====
  p <- ggplot(df_long, aes(x = Haplotype, y = Yield_tha)) +
    geom_boxplot(
      fill = NA,
      color = "red",
      linewidth = 1,
      outlier.colour = "black"
    ) +
    facet_wrap(~Region, scales = "free_x") +
    theme_bw(base_size = BASE_SIZE) +
    theme(
      axis.text.x  = element_text(size = TICK_SIZE, angle = 45, hjust = 1),
      axis.text.y  = element_text(size = TICK_SIZE),
      axis.title.x = element_text(size = TITLE_SIZE),
      axis.title.y = element_text(size = TITLE_SIZE),
      strip.text   = element_text(size = STRIP_SIZE),
      legend.position = "none"
    ) +
    labs(
      x = "Haplotype",
      y = "Yield (t ha⁻¹)"
    ) +
    stat_pvalue_manual(
      tukey,
      label = "p.label",
      tip.length = 0.01,
      size = 6/ .pt
    )
  
  print(p)
  
  plots[[facet_title]] <- p
  
  # ==== Save TIFF (300 dpi) ====
  file_stub <- str_replace_all(facet_title, "[^A-Za-z0-9]+", "_")
  
  ggsave(
    filename = file.path(outdir, paste0(file_stub, ".tiff")),
    plot = p,
    width = 3.15,
    height = 3.5,
    units = "in",
    dpi = 300,
    compression = "lzw"
  )
}
