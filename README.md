# seekerBio
This is an R package that have functions to find gene and SNP information in [GWAS catalog](https://www.ebi.ac.uk/gwas/), [ENSEMBL](https://www.ensembl.org/info/docs/api/index.html), [Reactome](https://reactome.org) 

## Installation

```
if (!require(devtools)) {
    install.packages("devtools")
}
devtools::install_github("Erickcufe/seekerBio")
```

## Source
[Reactome](https://reactome.org)

Antonio Fabregat, Steven Jupe, Lisa Matthews, Konstantinos Sidiropoulos, Marc Gillespie, Phani Garapati, Robin Haw, Bijay Jassal, Florian Korninger, Bruce May, Marija Milacic, Corina Duenas Roca, Karen Rothfels, Cristoffer Sevilla, Veronica Shamovsky, Solomon Shorser, Thawfeek Varusai, Guilherme Viteri, Joel Weiser, Guanming Wu, Lincoln Stein, Henning Hermjakob, Peter D’Eustachio, The Reactome Pathway Knowledgebase, Nucleic Acids Research, Volume 46, Issue D1, 4 January 2018, Pages D649–D655, https://doi.org/10.1093/nar/gkx1132

[Ensembl REST](https://rest.ensembl.org)

Andrew Yates, Kathryn Beal, Stephen Keenan, William McLaren, Miguel Pignatelli, Graham R. S. Ritchie, Magali Ruffier, Kieron Taylor, Alessandro Vullo, Paul Flicek, The Ensembl REST API: Ensembl Data for Any Language, Bioinformatics, Volume 31, Issue 1, 1 January 2015, Pages 143–145, https://doi.org/10.1093/bioinformatics/btu613

[GWAS Catalog](https://www.ebi.ac.uk/gwas/)

Buniello A, MacArthur JAL, Cerezo M, Harris LW, Hayhurst J, Malangone C, McMahon A, Morales J, Mountjoy E, Sollis E, Suveges D, Vrousgou O, Whetzel PL, Amode R, Guillen JA, Riat HS, Trevanion SJ, Hall P, Junkins H, Flicek P, Burdett T, Hindorff LA, Cunningham F and Parkinson H.
The NHGRI-EBI GWAS Catalog of published genome-wide association studies, targeted arrays and summary statistics 2019.
Nucleic Acids Research, 2019, Vol. 47 (Database issue): D1005-D1012.


## Examples:
```{r}
library(seekerBio)
```
```{r}
seeker_snp_ld_plot("rs7412", population_study="1000GENOMES:phase_3:FIN", color_select = "green")
```

```{r}

seeker_snp_freq("rs7412", "rs12952")

seeker_gwas("Obesity")
```
```{r}
snps <- seeker_snp_freq("rs7412")
ordered_snps <- seeker_snp_freq_format(snps)

```
