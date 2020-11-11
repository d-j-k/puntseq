![Scripting](https://img.shields.io/badge/Language-Python,_R_&_Bash-yellow.svg) [![Copyright](https://img.shields.io/badge/Copyright-(c)_2020_PuntSeq-green.svg)](License) [![DOI](https://img.shields.io/badge/DOI-10.1101/2020.02.06.936302-blue.svg)](https://doi.org/10.1101/2020.02.06.936302)

![alt text](https://github.com/d-j-k/puntseq/blob/master/figure1.png)

# Freshwater monitoring by nanopore sequencing

*Lara Urban, Andre Holzer, J Jotautas Baronas, Michael Hall, Philipp Braeuninger-Weimer, Michael J Scherm, Daniel J Kunz, Surangi N Perera, Daniel E Martin-Herranz, Edward T Tipper, Susannah J Salter, and Maximilian R Stammnitz*

Clean freshwater lies at the heart of human society. While most traditional water monitoring approaches test for specific chemicals or pathogens, the direct tracing of aquatic DNA poses a more holistic alternative which has hitherto been underappreciated due to challenges in logistics and investment. Here we present a simple, fast, inexpensive and reliable freshwater diagnostics workflow centred around portable nanopore DNA sequencing. Using bacterial mock communities and spatiotemporal microbiata from an example river in Cambridge (UK), our study shows how nanopore sequencing can be readily integrated for the assessment of aquatic bacterial diversity and pollution. We provide a computational benchmark that features more than ten taxonomic classification tools to derive guidelines for bacterial DNA analyses with nanopore data. Through complementary physicochemical measurements, we find that nanopore metagenomics can depict fine temporal gradients along the main hydrological axis of an urban-rural interface, in addition to yielding high-resolution pathogen maps that address concerns of public health.

![alt text](https://github.com/d-j-k/puntseq/blob/master/figure2.png)
*Figure: Overview of the experimental design of this study*

## Description

Here, we provide additional environmental data, the final classifications (using [Minimap2](https://github.com/lh3/minimap2), -k = 15) of all nanopore sequencing reads from three sampling dates (April, June, and August 2018) across nine sampling locations, including rarefied datasets. 
We additionally provide a Snakemake framework that integrates all data pre-processing steps and a Singularity that contains all necessary software. We further provide scripts for the downstream analyses (written in R and python, integrated in a markdown file) and an appropriate conda environment.

Using this platform, the user will be able to replicate all results presented in the corresponding study [preprint](https://www.biorxiv.org/content/10.1101/2020.02.06.936302).

Download the raw data from our [ENA repository](https://www.ebi.ac.uk/ena/data/view/PRJEB34900).

See [here](https://github.com/d-j-k/puntseq/tree/master/analysis) for the detailed description of our raw nanopore data pre-processing steps.

Check out the [project website](https://www.puntseq.co.uk/) and follow us on [Twitter](https://twitter.com/puntseq) for more updates!


## Citation

If you use any part or modified version of the here provided code in your work, you should:

1. state explicitly that you've used the PuntSeq code (or a modified version of it, if this is the case).
2. read and cite the following paper:

- Urban L, Holzer A, Baronas JJ, Hall M, Braeuninger-Weimer P, Scherm MJ, Kunz DJ, Perera SN, Martin-Herranz DE, Tipper ET, Salter SJ and Stammnitz MR (2020), **Freshwater monitoring by nanopore sequencing**, bioRxiv 2020.02.06.936302; doi: [https://doi.org/10.1101/2020.02.06.936302](https://doi.org/10.1101/2020.02.06.936302)
