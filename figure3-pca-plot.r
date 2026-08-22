#!/usr/bin/env Rscript
# ============================================================

INPUT_DIR  <- "/Users/easonlu/Desktop/group5r2=0.5-revised/ld/LDpruneR2-05/pca_fullpc_20260624_222942_2/snps_0095_LDpruned"   # working directory
GROUP_FILE <- file.path(INPUT_DIR, "group.txt")                       # group info file (header: ID,Rename,Country,Host,MtDNA,Group1,Group2)
EVEC_FILE  <- file.path(INPUT_DIR, "pca.eigenvec")                    # eigenvector file from PLINK
EVAL_FILE  <- file.path(INPUT_DIR, "pca.eigenval")                    # eigenvalue file from PLINK
OUTPUT_DIR <- INPUT_DIR                                                # output directory (can be changed)

POINT_SIZE      <- 5
POINT_STROKE    <- 0.25
POINT_ALPHA     <- 0.5
ELLIPSE_ALPHA   <- 0.18
LABEL_SIZE      <- 1.0
REPEL_SEGMENT_SIZE <- 0.25
PLOT_WIDTH      <- 7
PLOT_HEIGHT     <- 5.5
PLOT_DPI        <- 1200

SHAPE_MAP <- c("Ac" = 24,        # triangle (filled)
               "Am" = 21,        # circle (filled)
               "other" = 21)     # default circle
SHAPE_LABELS <- c("Ac" = "host Ac", "Am" = "host Am", "other" = "Other")

host_colors <- c(
  "Ac" = "#0066CC", 
  "Am" = "#FF3333"  
)

country_colors <- c(
  "JPN" = "#e0a7fa",
  "USA" = "#ff6f63",
  "CHN" = "#f96e16",
  "ZAF" = "#8c6f67",
  "RUS" = "#f5bfad",
  "VNM" = "#06b6d4",
  "KOR" = "#8b5cf6",
  "CHE" = "#ffff95",
  "THA" = "#a3ddd0"
)

mtdna_colors <- c(
  "C" = "#0066CC",
  "K1" = "#FF3333",
  "V1" = "#33CC33",
  "J1" = "#FFCC00"
)

DRAW_ELLIPSE <- FALSE

SHOW_LEGEND <- FALSE
if (!require("data.table", quietly = TRUE)) install.packages("data.table")
if (!require("ggplot2", quietly = TRUE)) install.packages("ggplot2")
if (!require("scales", quietly = TRUE)) install.packages("scales")
if (!require("ggrepel", quietly = TRUE)) install.packages("ggrepel")
if (!require("dplyr", quietly = TRUE)) install.packages("dplyr")
if (!require("ggExtra", quietly = TRUE)) install.packages("ggExtra")   # for marginal plots
if (!require("svglite", quietly = TRUE)) install.packages("svglite")   # for SVG output

library(data.table)
library(ggplot2)
library(scales)
library(ggrepel)
library(dplyr)
library(ggExtra)
library(svglite)

evec <- fread(EVEC_FILE, header = FALSE)
colnames(evec) <- c("FID", "IID", paste0("PC", 1:(ncol(evec)-2)))

eval <- fread(EVAL_FILE, header = FALSE)
colnames(eval) <- c("Eigenvalue")
eval[, Variance := Eigenvalue / sum(Eigenvalue)]
eval[, PC := 1:.N]

if (!file.exists(GROUP_FILE)) {
  stop("Missing group.txt file.")
}
group <- fread(GROUP_FILE, header = TRUE)
required_cols <- c("ID", "Rename", "Country", "Host", "MtDNA", "Group1", "Group2")
missing <- setdiff(required_cols, colnames(group))
if (length(missing) > 0) {
  stop("group.txt missing required columns: ", paste(missing, collapse=", "))
}

evec <- merge(evec, group, by.x = "IID", by.y = "ID", all.x = TRUE)
evec[, Shape := fifelse(Host == "Ac", "Ac",
                        fifelse(Host == "Am", "Am", "other"))]

theme_pub <- function(base_size = 14) {
  legend_pos <- if (SHOW_LEGEND) "right" else "none"
  theme_bw(base_size = base_size) +
    theme(
      panel.grid.major = element_line(color = "gray", linewidth = 0.5), 
      panel.grid.minor = element_blank(),                              
      axis.line = element_line(color = "black", linewidth = 0.5),
      axis.ticks = element_line(color = "black"),
      axis.title = element_text(face = "bold", size = 13),
      axis.text = element_text(size = 11, color = "black"),
      plot.title = element_text(face = "bold", hjust = 0.5, size = 14),
      legend.title = element_blank(),
      legend.position = legend_pos,
      legend.text = element_text(size = 9),
      legend.key.size = unit(0.4, "cm"),
      panel.border = element_rect(color = "black", linewidth = 1.5, fill = NA) 
    )
}

generate_pca_plots <- function(data, group_col, prefix, color_vector) {
  
  cat("\n Generating plots for", group_col, "(output prefix:", prefix, ")...\n")
  
  p <- ggplot(data, aes(x = PC1, y = PC2, 
                        fill = !!sym(group_col), 
                        shape = Shape))
  
  if (DRAW_ELLIPSE) {
    counts <- data %>% group_by(!!sym(group_col)) %>% summarise(n = n())
    eligible <- counts %>% filter(n >= 3) %>% pull(!!sym(group_col))
    if (length(eligible) > 0) {
      p <- p + stat_ellipse(
        data = data %>% filter(!!sym(group_col) %in% eligible),
        type = "norm", level = 0.95,
        alpha = ELLIPSE_ALPHA, geom = "polygon",
        color = NA, show.legend = FALSE
      )
    } else {
      cat(" No group with >=3 samples for ellipses.\n")
    }
  }
  
  p <- p +
    geom_point(size = POINT_SIZE, stroke = POINT_STROKE, 
               color = "black", alpha = POINT_ALPHA) +
    scale_fill_manual(values = color_vector, na.value = "grey70") +
    scale_shape_manual(
      values = SHAPE_MAP,
      labels = SHAPE_LABELS,
      drop = FALSE
    ) +
    labs(
      x = paste0("PC1 (", percent(eval$Variance[1]), ")"),
      y = paste0("PC2 (", percent(eval$Variance[2]), ")")
    ) +
    theme_pub() +
    guides(
      fill  = guide_legend(override.aes = list(shape = 21, size = 3, 
                                               color = "black", stroke = POINT_STROKE, alpha = 1)),
      shape = guide_legend(override.aes = list(size = 3, color = "black", 
                                               stroke = POINT_STROKE, alpha = 1))
    )
  
  # ---- Version with sample labels (no title) ----
  labeled <- p +
    geom_text_repel(aes(label = Rename),
                    size = LABEL_SIZE,
                    color = "black",
                    max.overlaps = Inf,
                    segment.color = "grey70",
                    segment.size = REPEL_SEGMENT_SIZE)
  
  ggsave(filename = file.path(OUTPUT_DIR, paste0("pca_scatter_", prefix, "_labeled.png")),
         plot = labeled, width = PLOT_WIDTH, height = PLOT_HEIGHT, dpi = PLOT_DPI)
  clean <- p
  ggsave(filename = file.path(OUTPUT_DIR, paste0("pca_scatter_", prefix, "_clean.png")),
         plot = clean, width = PLOT_WIDTH, height = PLOT_HEIGHT, dpi = PLOT_DPI)
  
  cat("  - Adding marginal density plots...\n")
  
  data_marginal <- p$data %>% filter(!is.na(!!sym(group_col)))
  
  p_main <- ggplot(data_marginal, aes(x = PC1, y = PC2, 
                                      fill = !!sym(group_col), 
                                      shape = Shape)) +
    geom_point(size = POINT_SIZE, stroke = POINT_STROKE, 
               color = "black", alpha = POINT_ALPHA) +
    scale_fill_manual(values = color_vector, na.value = "grey70") +
    scale_shape_manual(
      values = SHAPE_MAP,
      labels = SHAPE_LABELS,
      drop = FALSE
    ) +
    labs(
      x = paste0("PC1 (", percent(eval$Variance[1]), ")"),
      y = paste0("PC2 (", percent(eval$Variance[2]), ")")
    ) +
    theme_pub() +
    guides(
      fill  = guide_legend(override.aes = list(shape = 21, size = 3, 
                                               color = "black", stroke = POINT_STROKE, alpha = 1)),
      shape = guide_legend(override.aes = list(size = 3, color = "black", 
                                               stroke = POINT_STROKE, alpha = 1))
    )
  
  p_marginal <- ggMarginal(
    p_main,
    type = "density",
    groupFill = TRUE,
    alpha = 0.5,
    size = 8,          
    margins = "both",
    xparams = list(adjust = 1),
    yparams = list(adjust = 1)
  )
  
  marginal_file <- file.path(OUTPUT_DIR, paste0("pca_scatter_", prefix, "_marginal.png"))
  png(marginal_file, width = PLOT_WIDTH + 2, height = PLOT_HEIGHT + 2, units = "in", res = PLOT_DPI)
  print(p_marginal)
  dev.off()
  
  cat("  - pca_scatter_", prefix, "_labeled.png\n", sep="")
  cat("  - pca_scatter_", prefix, "_clean.png\n", sep="")
  cat("  - pca_scatter_", prefix, "_marginal.png\n", sep="")
  
  # ---- Generate separate SVG legend for this grouping ----
  cat("  - Generating SVG legend...\n")
  
  legend_plot <- ggplot(data_marginal, aes(x = 1, y = 1, 
                                           fill = !!sym(group_col), 
                                           shape = Shape)) +
    geom_point(alpha = 0) +  
    scale_fill_manual(values = color_vector, na.value = "grey70", drop = FALSE) +
    scale_shape_manual(values = SHAPE_MAP, labels = SHAPE_LABELS, drop = FALSE) +
    theme_void() +
    theme(
      legend.position = "right",
      legend.title = element_blank(),
      legend.text = element_text(size = 9),
      legend.key.size = unit(0.4, "cm"),
      legend.box = "vertical",         
      legend.margin = margin(0, 0, 0, 0)
    ) +
    guides(
      fill = guide_legend(
        override.aes = list(shape = 21, size = 3, alpha = 1, 
                            color = "black", stroke = POINT_STROKE),
        order = 1
      ),
      shape = guide_legend(
        override.aes = list(size = 3, alpha = 1, color = "black", stroke = POINT_STROKE),
        order = 2
      )
    )
  
  n_fill <- length(unique(data_marginal[[group_col]]))
  n_shape <- length(unique(data_marginal$Shape))
  legend_height <- (n_fill + n_shape) * 0.35 + 0.5
  legend_width <- 2.2  
  
  legend_file <- file.path(OUTPUT_DIR, paste0("pca_legend_", prefix, ".svg"))
  ggsave(legend_file, plot = legend_plot, 
         width = legend_width, height = legend_height, 
         dpi = PLOT_DPI, device = "svg")
  
  cat("  - pca_legend_", prefix, ".svg\n", sep="")
}

generate_pca_plots(evec, "Host",    "Host",   host_colors)

generate_pca_plots(evec, "MtDNA",   "Group1", mtdna_colors)

generate_pca_plots(evec, "Country", "Group2", country_colors)

p_scree <- ggplot(eval, aes(x = PC, y = Variance)) +
  geom_line(color = "black", linewidth = 0.7) +
  geom_point(size = 0.5, color = "black") +
  scale_x_continuous(breaks = 1:nrow(eval)) +
  labs(
    title = "Scree Plot",
    x = "Principal Component",
    y = "Proportion of Variance Explained"
  ) +
  theme_pub()

ggsave(filename = file.path(OUTPUT_DIR, "pca_scree.png"),
       plot = p_scree, width = 6, height = 4.5, dpi = PLOT_DPI)