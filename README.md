# mirrorplot

This R package is designed to generate plots to compare the region plots from 2 genetic analyses. This could be for example the results from two genome-wide association studies (GWAS) or from a GWAS and an expression analysis such as those from eQTL studies. These plots are good for displaying results alongside a colocalisation analysis.

Below is an example of a mirror plot comparing the results from a genome-wide analysis of idiopathic pulmonary fibrosis (IPF, above the x axis) to results from lung eQTL from the GTEx consortium for the *DEPTOR* gene (below the x axis). Each point is a genetic variant. Points above the x axis show the strength of association with the first trait (in this case risk of IPF) and points below the x axis show the strength of association with the second trait (in this case gene expression). Variants are coloured by their linkage disequilibrium with the GWAS sentinel and the location of the *DEPTOR* gene is shown by the green box on the x axis. From this plot it appears the GWAS and eQTL signal are driven by the same variants as the plot below the x axis mirrors the plot above the x axis. This plot was taken from [Allen et al (2020) Genome-Wide Association Study of Susceptibility to Idiopathic Pulmonary Fibrosis, *AJRCCM*](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7047454/).

![Example mirror plot](https://github.com/rjallen513/images/blob/master/DEPTOR_gtex_Lung.png?raw=true)

The basic use of this function is:

    mirrorplot(DF, CHR = ?, START = ?, END = ?)

To make the plot above we would use the code:

    mirrorplot(df, CHR = 2, START = 54970822, END = 56970822, SENTINEL = "rs148616939", GENE_START = 55401927, GENE_END = 55459699, T1THRESH1 = 5e-8, TITLE = "DEPTOR - GTEx (Lung)")

More details on how to use the function including for how to change the appearance of the plot see the manual mirrorplot.Rd in the "man" folder. Two example datasets are available in the "data" folder.


# How to install

    library(devtools)
    install_git("https://github.com/rjallen513/mirrorplot.git")
    library(mirrorplot)

# Input Data

The R function above uses a dataframe that must include the following columns:

    chr: Chromosome
    pos: Chromosomal position
    rsid: Variant rsid
    trait1_p: A p value for the first trait (that will be plotted above the x axis)
    trait2_p: A p value for the second trait (that will be plotted below the x axis)

The dataframe can also have the following optional columns:

    r2: A measure of how much LD the variant is with the variant of interest. This column is used to colour the variants. If this isn't in the dataframe then the code sets this to 0 for everyone.
    highlight: A column equal to 1 if you want the variant to have a different shape and 0 if not. For example, in the image at the top of the page, the variants in the credible set are shown by triangles. The code will also make sure the variants with a 1 are plotted above the variants with a 0 so they won't be hidden. If the dataframe doesn't have this column then it sets everyone to 0 (i.e. it doesn't highlight any variants and each variant is plotted with the same shape).

# Using the function

To use the function copy and paste it in the terminal/console. You can then generate the plot using the following command:

    mirrorplot(DF, CHR = ?, START = ?, END = ?)

To run it must have the following input variables:

    DF: The name of the input dataframe
    CHR: The chromosome you want to plot
    START: The position you want to start from
    END: The position you want to plot up to

There are also additional variables you can input but don't have to:

    SENTINEL: The rsid of a variant you want to highlight in blue (for example the variant that the LD was calculated with respect to)
    GENE_START: If you want to plot where a gene is (as shown by a green box on the x axis) give the start position of the gene, and
    GENE_END: The end position of the gene
    T1THRESH1: A threshold to plot as a dotted line in red for the first trait (for example in the plot above there is a threshold plotted at genome-wide significance of p = 5x10-8). The default is to not plot any threshold lines.
    T1THRESH2: A second line to plot for the first trait in blue (for example a suggestive significance line)
    T2THRESH1: A threshold for the second trait (in red)
    T2THRESH2: A threshold for the second trait (in blue)
    TITLE: A main title for the plot. The default is to not include a title.
    COLOURS: A vector of six colours. The first being the colour to highlight the SENTINEL variant and then 5 colours for the rest of the variants (with the colour for those in highest LD first). The default colours are used in the plot at the top.
    GENE_COL: The colour for the gene box. Default is "greenyellow".
    SHAPES: A vector of two shapes. The first shape determines the baseline shape and the second shape for the shape of the variants with a 1 in the highlight column of the data frame. For circles use 21, squares 22, diamonds 23, triangles 24 and upside down triangles 25. The default is for most variants to be circles and the highlighted variants to be triangles (i.e. SHAPES = c(21,24)).

# Acknowledgements
Thank you to Nick Shrine for helping me to create this package.
