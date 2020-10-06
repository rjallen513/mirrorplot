#' Generating plots to compare results from two genetic analyses
#' This code will generate a mirror plot to compare the results in a genetic region from two genetic analyses. For example this could be between two genome-wide analyses or between a genome-wide analysis and a gene expression analysis such as eQTL/pQTL data.
#'
#' @param DF A data.frame with columns "chr", "pos", "rsid", "trait1_p", "trait2_p" and optionally, "r2" and "highlight". Each row of the data.frame should be for a single variant. The "chr"" column contains the chromosome the variant is on, "pos" is the chromosomal position of the variant, "rsid" is the rsid, "trait1_p" is the p value for association with the first trait (which will be plotted above the x axis) and "trait2_p" is the p value for the association with the second trait (which will be plotted below the x axis). "r2" is a measure of how much LD the variant is with the variant of interest. This column is used to colour the variants. If r2 isn't in the dataframe then the code sets this to 0 for everyone. "highlight" is an optional column equal to 1 if you want the variant to have a different shape and 0 if not. The code will also make sure the variants with a 1 are plotted above the variants with a 0 so they won't be hidden. If the dataframe doesn't have this column then it sets everyone to 0 (i.e. it doesn't highlight any variants and each variant is plotted with the same shape).
#'
#' @param CHR The chromosome you want to plot
#' @param START The position you want to start the plot from
#' @param END The position you want to plot up to
#'
#' There are also additional variables you can input but don't have to:
#' @param SENTINEL The rsid of a variant you want to highlight in a different colour (for example the variant that the LD was calculated with respect to)
#' @param GENE_START If you want to plot where a gene is (as shown by a green box on the x axis) give the start position of the gene
#' @param GENE_END The end position of the gene
#' @param T1THRESH1 A threshold to plot as a horizontal dotted line in red for the first trait. The default is to not plot any threshold lines.}
#' @param T1THRESH2 A second horizontal line to plot for the first trait in blue.
#' @param T2THRESH1 A threshold line for the second trait (in red).
#' @param T2THRESH2 A threshold for the second trait (in blue).
#' @param TITLE A main title for the plot. The default is to not include a title.
#' @param COLOURS A vector of six colours. The first being the colour to highlight the SENTINEL variant and then 5 colours for the rest of the variants (with the colour for those in highest LD first). The default colours are used in the plot at the top.}
#' @param GENE_COL The colour for the gene box. Default is "greenyellow".
#' @param SHAPES A vector of two shapes. The first shape determines the baseline shape and the second shape for the shape of the variants with a 1 in the highlight column of the data frame. For circles use 21, squares 22, diamonds 23, triangles 24 and upside down triangles 25. The default is for most variants to be circles and the highlighted variants to be triangles (i.e. SHAPES = c(21,24)).}
#'
#' @examples{
#' mirrorplot(mirrorplot_testdata1, CHR=1, START=1, END=50, SENTINEL="rs20")
#'
#' mirrorplot(mirrorplot_testdata2, CHR=2, START=1, END=50, SENTINEL="rs65",
#'   GENE_START=10, GENE_END=20,
#'   T1THRESH1=5e-8, T2THRESH2=5e-8,
#'   TITLE="Example plot",
#'   COLOURS=c("red", "grey20", "grey40", "grey60", "grey80", "grey90"))
#' }


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

  # Check the input chromosome is valid and subset data to just that chromosome
  if (!(CHR %in% DF$chr)) stop(paste("The supplied chromosome number '", CHR, "' was not found in the dataset", sep=""))
  mirrorplot_temp <- subset(DF, chr==CHR)

  # check the input start and end values are valid for the specified chromosome
  if (END <= START) stop(paste("The end value must be higher than the start value"))
  if (START > max(mirrorplot_temp$pos)) stop(paste("The supplied start position '", START, "' is not within the range of positions in the dataset for chromosome number '", CHR, "'", sep=""))
  if (END < min(mirrorplot_temp$pos)) stop(paste("The supplied end position '", END, "' is not within the range of positions in the dataset for chromosome number '", CHR, "'", sep=""))

  # subset the data to just that needed for the plot (between start and end)
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
  par(mfrow=c(1,1), mar=c(6,7,4,2), bty="l",yaxs="i",xaxs="i")
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


