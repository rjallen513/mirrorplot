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


# Acknowledgements
Thank you to Nick Shrine for helping me to create this package.
