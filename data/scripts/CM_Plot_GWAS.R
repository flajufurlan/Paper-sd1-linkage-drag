
#Packages
library(CMplot)
library(dplyr)
library(stringr)
library(tiff)
library(grid)


fname <- "GAPIT.Association.GWAS_Results.FarmCPU.Plant.height.LSU_01.0000132.csv"

memo_tag <- fname %>%
  basename() %>%
  str_remove("\\.csv$") %>%
  str_replace_all("[^A-Za-z0-9._-]", "_")

gwas_df <- read.csv(fname, check.names = FALSE)

# Column mapping
col_map <- list(
  SNP = c("SNP", "SNP.", "Marker", "rs"),
  CHR = c("Chromosome", "Chr", "chr"),
  BP  = c("Position", "Pos", "BP"),
  PV  = c("P.value", "p.value", "P")
)

get_col <- function(x, nms) intersect(x, nms)[1]

cm <- gwas_df %>%
  transmute(
    SNP        = .data[[get_col(col_map$SNP, names(gwas_df))]],
    Chromosome = gsub("^chr", "", .data[[get_col(col_map$CHR, names(gwas_df))]], ignore.case = TRUE),
    Position   = as.numeric(.data[[get_col(col_map$BP, names(gwas_df))]]),
    P.value    = as.numeric(.data[[get_col(col_map$PV, names(gwas_df))]])
  ) %>%
  filter(!is.na(P.value), !is.na(Position))

# SIGNIFICANCE THRESHOLDS

bonf <- 0.05 / nrow(cm)
thr_vec <- c(1e-5, bonf)

highlight_snps <- c("182171b124800000060d")

# MANHATTAN PLOT 

CMplot(cm,
       type = "p",
       plot.type = "m",
       LOG10 = TRUE,
       threshold = thr_vec,
       threshold.lwd = c(1.5, 1.5),
       threshold.lty = c(2, 1),
       
       # ---- TEXT SIZE ≈ 10 pt ----
       axis.cex = 2,
       lab.cex  = 2,
       main.cex = 2,
       
       mar = c(5, 6, 4, 2),
       
       highlight = highlight_snps,
       highlight.col = "red",
       highlight.pch = 22,
       highlight.cex = 2.0,
       highlight.text = "Chr1:38.33 Mb",
       highlight.text.cex = 2,
       highlight.text.font = 3,
       
       file = "tiff",
       file.name = paste0(memo_tag, "_Manhattan_10pt"),
       dpi = 600,
       width = 14,
       height = 7
)

# QQ PLOT 

CMplot(cm,
       type = "p",
       plot.type = "q",
       LOG10 = TRUE,
       conf.int = TRUE,
       main = "",
       
       axis.cex = 1.5,
       lab.cex  = 1.5,
       
       file = "tiff",
       file.name = paste0(memo_tag, "_QQ_10pt"),
       dpi = 600,
       width = 6,
       height = 6
)

# READ TIFFS

plantheight_manhattan <- readTIFF("Rect_Manhtn.GAPIT.Association.GWAS_Results.FarmCPU.Plant.height.LSU_01.0000132_Manhattan_10pt.tiff")
plantheight_qq        <- readTIFF("QQplot.GAPIT.Association.GWAS_Results.FarmCPU.Plant.height.LSU_01.0000132_QQ_10pt.tiff")

yield_manhattan <- readTIFF("Rect_Manhtn.GAPIT.Association.GWAS_Results.FarmCPU.Yield.LSU_01.0000138_Manhattan_10pt.tiff")
yield_qq        <- readTIFF("QQplot.GAPIT.Association.GWAS_Results.FarmCPU.Yield.LSU_01.0000138_QQ_10pt.tiff")

# FINAL MULTI‑PANEL FIGURE (TPG SIZE)

tiff(
  "Figure_GWAS_PlantHeight_Yield_FINAL.tiff",
  width = 7.0,
  height = 7.6,
  units = "in",
  res = 600,
  compression = "lzw"
)

grid.newpage()
pushViewport(
  viewport(layout = grid.layout(
    nrow = 2,
    ncol = 2,
    widths  = unit(c(0.65, 0.35), "npc"),
    heights = unit(c(0.5, 0.5), "npc")
  ))
)

# Row A: Plant height
grid.raster(plantheight_manhattan,
            vp = viewport(layout.pos.row = 1, layout.pos.col = 1))
grid.raster(plantheight_qq,
            vp = viewport(layout.pos.row = 1, layout.pos.col = 2))

# Row B: Yield
grid.raster(yield_manhattan,
            vp = viewport(layout.pos.row = 2, layout.pos.col = 1))
grid.raster(yield_qq,
            vp = viewport(layout.pos.row = 2, layout.pos.col = 2))

# ---- PANEL LABELS (14 pt) ----
grid.text("A",
          x = unit(0.01, "npc"),
          y = unit(0.98, "npc"),
          just = c("left", "top"),
          gp = gpar(fontsize = 12, fontface = "bold"))

grid.text("B",
          x = unit(0.01, "npc"),
          y = unit(0.48, "npc"),
          just = c("left", "top"),
          gp = gpar(fontsize = 12, fontface = "bold"))

dev.off()

#182171b124800000060d
