### description:
# Compute genome-wide feature/window densities for one or more GRanges-like inputs:
#  - Normalize input to data frame(s) of genomic intervals
#  - Slide non‑overlapping or half‑overlapping windows of specified size across each chromosome
#  - Count or sum features per window (by "number" or other metric)
#  - Pad missing chromosomes/windows with zero values
#  - Return a list of per‑input density data frames with columns: chr, start, end, value

Calculate_Genomic_Density <- 
  function (data, window.size = NULL, overlap = TRUE, count_by = "number") {
    data = circlize:::normalizeToDataFrame(data)
    if (!circlize:::is.dataFrameList(data)) {
        data = list(data)
    }

    # Compute raw genomic densities for each dataset
    df = vector("list", length = length(data))
    for (i in seq_along(data)) {
        df[[i]] = genomicDensity(data[[i]], window.size = window.size, 
            overlap = overlap, count_by = count_by)
    }

    # Determine all chromosome sectors from circlize canvas
    all_chrs <- get.all.sector.index()
    step <- if (overlap) window.size / 2 else window.size

    # For each density table, pad missing chromosomes/windows with zeros
    for (i in seq_along(df)) {
      present <- unique(df[[i]]$chr)
      missing <- setdiff(all_chrs, present)
      for (chr in missing) {
        xr <- get.cell.meta.data("xrange", sector.index = chr)[1]
        starts <- seq(from = 1, to = step * (floor(xr / step) - 1), by = step)
        ends <- pmin(starts + window.size - 1, xr)
        pad_df <- data.frame(chr = chr,
                              start = starts,
                              end = ends,
                              value = 0)
        df[[i]] <- rbind(df[[i]], pad_df)
      }
      
      # Sort by chromosome and start coordinate
      df[[i]] <- df[[i]][order(df[[i]]$chr, df[[i]]$start), ]
      row.names(df[[i]]) <- paste(
        df[[i]]$chr,
        sprintf("%.0f", df[[i]]$start),
        sprintf("%.0f", df[[i]]$end),
        sep = "_"
      )
    }
    df
  }