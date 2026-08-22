# setwd
setwd("/Users/easonlu/Desktop/1127_version/figure1/1211version/Ac")
output_dir <- "./output"
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# define colors
colors <- c("#0e2c82", "#346695", "#9DBFDC", "#E8EEF6", "#2ca02c", "#5BABA0", "#1C5A41", "#96CA00", "#FFEA00", "#b6b51f", "#da1f18","#FF2400")
names(colors) <- c("K1-1", "K1-3", "K1-5", "K1-7", "C1-1", "C1-2", "C2-1", "C3-1", "J1-3", "J1-7", "V1-2", "V1-4")

# read data
data <- read.table("input.txt", header = FALSE, stringsAsFactors = FALSE)

# pie
for (i in 1:nrow(data)) {
  row_name <- data[i, 1]
  values <- as.numeric(data[i, 2:13])  # col2 - col13
  output_file <- file.path(output_dir, paste0(row_name, ".png"))
  png(output_file, width = 2000, height = 2000, res = 500, bg = "transparent")
  par(mar = c(0, 0, 0, 0), oma = c(0, 0, 0, 0))
  pie(values, 
      col = colors,
      labels = NA,
      border = "white", 
      lwd = 2,
      main = "")
  dev.off()
}

# legend
legend_file <- file.path(output_dir, "legend.png")
png(legend_file, width = 3000, height = 3000, res = 500, bg = "transparent")
par(mar = c(2, 2, 2, 2), oma = c(1, 1, 1, 1))

plot(1, type = "n", axes = FALSE, xlab = "", ylab = "", 
     xlim = c(0, 1), ylim = c(0, 1))
legend("center", 
       legend = names(colors),
       fill = colors,
       border = "white",
       bty = "n",
       cex = 2.5,
       ncol = 2,
       x.intersp = 0.5,
       y.intersp = 1.5,
       text.width = 0.3)
dev.off()

cat("Success!", output_dir, "filepath\n")