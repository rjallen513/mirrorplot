# This funtion is for creating a mirror plot

mirrorplot <- function(DF, CHR, START, END, SENTINEL = "sentinel",
                       GENE_START = NA, GENE_END = NA,
                       T1THRESH1 = NA, T1THRESH2 = NA, T2THRESH1 = NA, T2THRESH2 = NA,
                       TITLE = NA,
                       COLOURS = c("blue", "red3", "orange", "yellow", "lemonchiffon2", "gray60"),
                       GENE_COL = "greenyellow",
                       SHAPES = c(21,24)){

  # Check dataframe has the required fields
  if (!("chr" %in% names(DF))) stop(paste("Column chr not found"))
  if (!("pos" %in% names(DF))) stop(paste("Column pos not found"))
  if (!("rsid" %in% names(DF))) stop(paste("Column rsid not found"))
  if (!("trait1_p" %in% names(DF))) stop(paste("Column trait1_p not found"))
  if (!("trait2_p" %in% names(DF))) stop(paste("Column trait2_p not found"))

  # Create r2 and highlight columns if not already included
  if (!("r2" %in% names(DF))) DF$r2 <- 0
  if (!("highlight" %in% names(DF))) DF$highlight <- 0

  # Subset data to just that needed for the plot
  mirrorplot_temp <- subset(DF, chr==CHR)
  mirrorplot_temp <- subset(mirrorplot_temp, pos >= START)
  mirrorplot_temp <- subset(mirrorplot_temp, pos <= END)

  # Calculate -log10 p values for the two traits
  mirrorplot_temp$t1logp <- -log10(mirrorplot_temp$trait1_p)
  mirrorplot_temp$t2logp <- -log10(mirrorplot_temp$trait2_p)

  # Set the colours
  mirrorplot_temp$col <- "white"
  mirrorplot_temp$col <- ifelse(mirrorplot_temp$r2 >= 0.8, COLOURS[2], mirrorplot_temp$col)
  mirrorplot_temp$col <- ifelse(mirrorplot_temp$r2 < 0.8 & mirrorplot_temp$r2 >= 0.6, COLOURS[3], mirrorplot_temp$col)
  mirrorplot_temp$col <- ifelse(mirrorplot_temp$r2 < 0.6 & mirrorplot_temp$r2 >= 0.4, COLOURS[4], mirrorplot_temp$col)
  mirrorplot_temp$col <- ifelse(mirrorplot_temp$r2 < 0.4 & mirrorplot_temp$r2 >= 0.2, COLOURS[5], mirrorplot_temp$col)
  mirrorplot_temp$col <- ifelse(mirrorplot_temp$r2 < 0.2, COLOURS[6], mirrorplot_temp$col)
  mirrorplot_temp$col <- ifelse(mirrorplot_temp$rsid==SENTINEL, COLOURS[1], mirrorplot_temp$col)
  mirrorplot_temp$col <- ifelse(is.na(mirrorplot_temp$col)==T, "white", mirrorplot_temp$col)

  # Set the shape to be plotted
  mirrorplot_temp$shape <- ifelse(mirrorplot_temp$highlight==1, SHAPES[2], SHAPES[1])

  # Calculate the scale for the y axis
  a1 <- max(mirrorplot_temp$t1logp, na.rm=T) / 10
  scale1 <- ceiling(a1)

  a2 <- max(mirrorplot_temp$t2logp, na.rm=T) / 10
  scale2 <- ceiling(a2)

  lab1 <- 10 * scale2
  lab2 <- 5 * scale2
  lab3 <- 5 * scale1
  lab4 <- 10 * scale1

  # Order the dataframe so that it plots the highlighted variants on top and then the most significant variants for trait 1 on top
  mirrorplot_temp <- mirrorplot_temp[order(mirrorplot_temp$t1logp),]
  mirrorplot_temp <- mirrorplot_temp[order(mirrorplot_temp$highlight),]

  # Make the plot
  plot(mirrorplot_temp$pos, mirrorplot_temp$t1logp/scale1, ylim=c(-10, 10), pch=mirrorplot_temp$shape, bg=mirrorplot_temp$col, col="black", cex=2, axes=F, xlim=c(START, END), xlab="", ylab="-log(p value)", cex.lab=2, main = TITLE, cex.main = 2)
  abline(h=-log10(T1THRESH1)/scale1, col="red", lty=3)
  abline(h=-log10(T1THRESH2)/scale1, col="blue", lty=3)
  abline(h=log10(T2THRESH1)/scale2, col="red", lty=3)
  abline(h=log10(T2THRESH2)/scale2, col="blue", lty=3)
  points(mirrorplot_temp$pos, -mirrorplot_temp$t2logp/scale2, pch=mirrorplot_temp$shape, bg=mirrorplot_temp$col, col="black",cex=2)
  abline(h=0, col="black", lwd=5)
  axis(2, at=c(-10, -5, 0, 5, 10), labels=c(lab1, lab2, "0", lab3, lab4), cex.axis=2, tck=-0.015)
  rect(xleft=GENE_START, xright=GENE_END, ybottom=-0.5, ytop=0.5, col="greenyellow", density = 50, border="black")

}


