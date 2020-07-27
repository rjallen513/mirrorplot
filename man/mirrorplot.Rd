\name{mirrorplot}
\alias{mirrorplot}
\title{Generating plots to compare results from two genetic analyses}
\usage{
mirrorplot(DF, CHR, START, END, SENTINEL = "sentinel",
                       GENE_START = NA, GENE_END = NA,
                       T1THRESH1 = NA, T1THRESH2 = NA, T2THRESH1 = NA, T2THRESH2 = NA,
                       TITLE = NA,
                       COLOURS = c("blue", "red3", "orange", "yellow", "lemonchiffon2", "gray60"),
                       GENE_COL = "greenyellow",
                       SHAPES = c(21,24))
}
\description{
This code will generate a mirror plot to compare the results in a genetic region from two genetic analyses. For example this could be between two genome-wide analyses or between a genome-wide analysis and a gene expression analysis such as eQTL/pQTL data.
}
\arguments{
\item{DF}{A data.frame with columns "chr", "pos", "rsid", "trait1_p", "trait2_p" and optionally, "r2" and "highlight". Each row of the data.frame should be for a single variant. The "chr"" column contains the chromosome the variant is on, "pos" is the chromosomal position of the variant, "rsid" is the rsid, "trait1_p" is the p value for association with the first trait (which will be plotted above the x axis) and "trait2_p" is the p value for the association with the second trait (which will be plotted below the x axis). "r2" is a measure of how much LD the variant is with the variant of interest. This column is used to colour the variants. If r2 isn't in the dataframe then the code sets this to 0 for everyone. "highlight" is an optional column equal to 1 if you want the variant to have a different shape and 0 if not. The code will also make sure the variants with a 1 are plotted above the variants with a 0 so they won't be hidden. If the dataframe doesn't have this column then it sets everyone to 0 (i.e. it doesn't highlight any variants and each variant is plotted with the same shape).}

\item{CHR}{The chromosome you want to plot}
\item{START}{The position you want to start the plot from}
\item{END}{The position you want to plot up to}

There are also additional variables you can input but don't have to:

\item{SENTINEL}{The rsid of a variant you want to highlight in a different colour (for example the variant that the LD was calculated with respect to)}
\item{GENE_START}{If you want to plot where a gene is (as shown by a green box on the x axis) give the start position of the gene, and}
\item{GENE_END}{The end position of the gene}
\item{T1THRESH1}{A threshold to plot as a horizontal dotted line in red for the first trait. The default is to not plot any threshold lines.}
\item{T1THRESH2}{A second horizontal line to plot for the first trait in blue.}
\item{T2THRESH1}{A threshold line for the second trait (in red).}
\item{T2THRESH2}{A threshold for the second trait (in blue).}
\item{TITLE}{A main title for the plot. The default is to not include a title.}
\item{COLOURS}{A vector of six colours. The first being the colour to highlight the SENTINEL variant and then 5 colours for the rest of the variants (with the colour for those in highest LD first). The default colours are used in the plot at the top.}
\item{GENE_COL}{The colour for the gene box. Default is "greenyellow".}
\item{SHAPES}{A vector of two shapes. The first shape determines the baseline shape and the second shape for the shape of the variants with a 1 in the highlight column of the data frame. For circles use 21, squares 22, diamonds 23, triangles 24 and upside down triangles 25. The default is for most variants to be circles and the highlighted variants to be triangles (i.e. SHAPES = c(21,24)).}
}
\examples{

##### Example 1 #####

# Example where the two traits colocalise. The r2 column shows the linkage disequilibrium with the sentinel variant rs20.

mirrorplot(mirrorplot_testdata1, CHR=1, START=1, END=50, SENTINEL="rs20")


##### Example 2 #####

# Example where the two traits do not colocalise. The r2 column shows the linkage disequilibrium with the sentinel variant rs65. There is no highlight column included in the dataframe.

# Plot the data adding a genome-wide significance threshold at 5x10-8 for both traits (in red for the first trait and in blue for the second trait), a title, changing the colours and plotting a gene position.

mirrorplot(mirrorplot_testdata2, CHR=2, START=1, END=50, SENTINEL="rs65",
  GENE_START=10, GENE_END=20,
  T1THRESH1=5e-8, T2THRESH2=5e-8,
  TITLE="Example plot",
  COLOURS=c("red", "grey20", "grey40", "grey60", "grey80", "grey90"))
}
