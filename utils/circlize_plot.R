### description:
# Plot genome‑wide DMR distributions as circos tracks:
#  - Initialize a circular ideogram using hg38 cytobands
#  - Add a genomic track of DMR points (hypermethylated in orange, hypomethylated in blue)
#  - Draw baseline and threshold lines for visual reference
#  - Output to PDF with specified title

library(circlize)

plot_circlize <- function (DMR, file_path, title_name, cytoband_hg38_path) {
  pdf(file_path, width = 8, height = 8)
  par(mar = c(1, 1, 2, 1), cex.main = 1.5)

  # Configure circos layout: gaps between chromosomes, start at 12 o'clock
  circos.par( gap.after = 0.5, start.degree = 90, gap.after = c(rep(1, 22), 5, 5), cell.padding = c(0, 0, 0, 0),
                track.margin = c(0, 0), clock.wise = TRUE )
  
  # Draw ideogram with axis and labels
  circos.initializeWithIdeogram( cytoband = cytoband_hg38_path, plotType = c("axis", "labels"), 
                                  ideogram.height = 0.05, axis.labels.cex = 0.6  )
  if (nrow(DMR) != 0) {
      circos.genomicTrackPlotRegion(
          DMR,
          ylim = c(-1, 1),
          bg.col = "#F0F0F0",
          track.height = 0.6,
          track.margin = c(0, 0),
          bg.border = NA,
          panel.fun = function(region, value, ...) {
              for (y in c(-0.5, 0.5)) {
                circos.lines(CELL_META$cell.xlim, c(y, y),
                            col = adjustcolor("lightgrey", alpha.f = 0.5),
                            lty = 1, lwd = 0.8)
              }
              circos.lines(CELL_META$cell.xlim, c(0, 0), col = adjustcolor("red", alpha.f = 0.3), lty = 1, lwd = 0.8)
              cols <- ifelse(value[[1]] > 0,
                            adjustcolor("orange", alpha.f = 0.8),
                            adjustcolor("blue", alpha.f = 0.8))
              circos.genomicPoints(region, value, col = cols, pch = 16, cex = 0.5)
          }
      )
  }
  circos.clear()
  title(main = title_name, line = -1)
  dev.off()
  message(title_name, " Done!")
}