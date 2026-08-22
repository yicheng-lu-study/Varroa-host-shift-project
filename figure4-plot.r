library(ggplot2)
library(reshape2)
library(ggsci)
library(patchwork)
library(dplyr)

setwd("/Users/easonlu/Desktop/admixture-random-seed/admixture-seed-1008")

file_prefix <- "QC"
order_file <- "samplelist_v1221.txt"

aspect_ratio_width <- 11
aspect_ratio_height <- 1

custom_colors <- list(
  k2 = c("#0066CC", "#EDE5F9"),
  k3 = c("#FFCC00", "#0066CC", "#EDE5F9"),
  k4 = c("#0066CC", "#FFCC00", "#33CC33", "#EDE5F9"),
  k5 = c("#0066CC", "#33CC33", "#FF3333", "#EDE5F9", "#FFCC00"),
  k6 = c("#FF3333", "#0066CC", "#33CC33", "#9933CC", "#FFCC00", "#EDE5F9"),
  k7 = c("#9933CC", "#2A2A2A", "#EDE5F9", "#FFCC00", "#FF3333", "#0066CC", "#33CC33"),
  k8 = c("#FF6600", "#33CC33", "#00CCCC", "#EDE5F9", "#9933CC", "#FF3333", "#0066CC", "#FFCC00"),
  k9 = c("#00CCCC", "#9933CC", "#FF6600", "#0066CC", "#2A2A2A", "#33CC33", "#EDE5F9", "#FF3333", "#FFCC00")
)


if (file.exists(order_file)) {
  order_df <- read.table(order_file, header = FALSE, sep = "\t", 
                         stringsAsFactors = FALSE, col.names = "sample_order")
  sample_order <- order_df$sample_order
  message(paste(length(sample_order), "samples"))
} else {
  sample_order <- NULL
  warning("file for sample order donnot exist: ", order_file)
}

process_k_value <- function(k) {
  infile <- paste0(file_prefix, ".", k, ".Q.tab")
  
  if (!file.exists(infile)) {
    stop("File donnot exist: ", infile)
  }
  
  tbl <- read.table(infile)
  colnames(tbl) <- c("strain", paste("Population", seq(ncol(tbl)-1), sep = ''))
  
  for (i in 2:ncol(tbl)) {
    tbl[, i] <- as.numeric(as.character(tbl[, i]))
  }
  
  if (!is.null(sample_order)) {

    missing_samples <- setdiff(tbl$strain, sample_order)
    if (length(missing_samples) > 0) {
      warning("The following samples are missing from the sequential file and will be appended to the end: ", 
              paste(missing_samples, collapse = ", "))
      sample_order <- c(sample_order, missing_samples)
    }
    
    tbl$strain <- factor(tbl$strain, levels = sample_order)
    tbl <- tbl[order(tbl$strain), ]
  } else {
    tbl$strain <- factor(tbl$strain, levels = unique(tbl$strain))
  }
  
  main_pop <- apply(tbl[, -1, drop = FALSE], 1, 
                    function(t) colnames(tbl)[which.max(t) + 1])
  tbl$pop <- main_pop
  
  data_long <- melt(
    data = tbl,
    id.vars = c("strain", "pop"),
    variable.name = "Population",
    value.name = "Ancestry"
  )
  
  data_long$Ancestry <- as.numeric(as.character(data_long$Ancestry))
  
  k_colors <- custom_colors[[paste0("k", k)]]
  names(k_colors) <- paste0("Population", 1:k)
  
  return(list(data = data_long, colors = k_colors, k = k))
}

create_plot <- function(processed_data, show_labels = FALSE) {
  data_long <- processed_data$data
  k_colors <- processed_data$colors
  k <- processed_data$k
  
  n_samples <- length(unique(data_long$strain))
  plot_width <- n_samples * 0.1
  plot_height <- plot_width * aspect_ratio_height / aspect_ratio_width
  
  p <- ggplot(data_long, aes(x = strain, y = Ancestry, fill = Population)) +
    geom_bar(stat = "identity", position = "fill", width = 1) +
    scale_fill_manual(values = k_colors) +
    scale_y_continuous(expand = c(0, 0)) +
    labs(title = paste0("K = ", k), x = NULL, y = "Ancestry Proportion") +
    theme_minimal() +
    theme(
      axis.text.x = if(show_labels) element_text(angle = 90, hjust = 1, vjust = 0.5, size = 6) else element_blank(),
      axis.ticks.x = element_blank(),
      legend.position = "none",  
      panel.grid = element_blank(),
      plot.title = element_text(hjust = 0.5, size = 12, face = "bold"),
      plot.margin = margin(5, 5, 5, 5),
      axis.text.y = element_text(size = 6), 
      axis.title.y = element_text(size = 8)  
    )
  
  return(p)
}

k_values <- 2:7
all_plots <- list()

for (k in k_values) {
  processed_data <- process_k_value(k)
  all_plots[[paste0("k", k)]] <- processed_data
}

plots_with_labels <- list()
plots_without_labels <- list()

for (k in k_values) {
  processed_data <- all_plots[[paste0("k", k)]]
  plots_with_labels[[k-1]] <- create_plot(processed_data, show_labels = TRUE)
  plots_without_labels[[k-1]] <- create_plot(processed_data, show_labels = FALSE)
}

combined_plot_with_labels <- wrap_plots(plots_with_labels, ncol = 1)
combined_plot_without_labels <- wrap_plots(plots_without_labels, ncol = 1)

n_samples <- length(unique(all_plots[[1]]$data$strain))
single_plot_width <- n_samples * 0.1
single_plot_height <- single_plot_width * aspect_ratio_height / aspect_ratio_width
combined_height <- single_plot_height * length(k_values)

ggsave(
  filename = "Admixture_K2-5_combined_with_labels-seed-1008.png",
  plot = combined_plot_with_labels,
  width = single_plot_width,
  height = combined_height,
  dpi = 700,
  limitsize = FALSE
)

ggsave(
  filename = "Admixture_K2-5_combined_with_labels-seed-1008.pdf",
  plot = combined_plot_with_labels,
  width = single_plot_width,
  height = combined_height,
  device = "pdf",
  limitsize = FALSE
)

ggsave(
  filename = "Admixture_K2-5_combined_without_labels-seed-1008.png",
  plot = combined_plot_without_labels,
  width = single_plot_width,
  height = combined_height,
  dpi = 700,
  limitsize = FALSE
)

ggsave(
  filename = "Admixture_K2-5_combined_without_labels-seed-1008.pdf",
  plot = combined_plot_without_labels,
  width = single_plot_width,
  height = combined_height,
  device = "pdf",
  limitsize = FALSE
)