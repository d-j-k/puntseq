![alt text](https://github.com/d-j-k/puntseq/blob/master/figure1.png)


# Freshwater monitoring with nanopore sequencing

*[Lara Urban](https://github.com/LaraUrban), [Andre Holzer](https://github.com/AndreHolzer), J Jotautas Baronas, [Michael Hall](https://github.com/mbhall88), Philipp Braeuninger-Weimer, [Michael J Scherm](https://github.com/MScherm), [Daniel J Kunz](https://github.com/d-j-k), Surangi N Perera, [Daniel E Martin-Herranz](https://github.com/demh), Edward T Tipper, Susannah J Salter, and [Maximilian R Stammnitz](https://github.com/MaximilianStammnitz)*

Clean freshwater lies at the heart of human society. While most traditional water monitoring approaches test for specific chemicals or pathogens, the direct tracing of aquatic DNA poses a more holistic alternative which has hitherto been underappreciated due to challenges in logistics and investment. Here we present a simple, fast, inexpensive and reliable freshwater diagnostics workflow centred around portable nanopore DNA sequencing. Using bacterial mock communities and spatiotemporal microbiata from an example river in Cambridge (UK), our study shows how nanopore sequencing can be readily integrated for the assessment of aquatic bacterial diversity and pollution. We provide a computational benchmark that features more than ten taxonomic classification tools to derive guidelines for bacterial DNA analyses with nanopore data. Through complementary physicochemical measurements, we find that nanopore metagenomics can depict fine temporal gradients along the main hydrological axis of an urban-rural interface, in addition to yielding high-resolution pathogen maps that address concerns of public health.

Here, we provide additional environmental data, the final classifications (using [Minimap2](https://github.com/lh3/minimap2), -k = 15) of all nanopore sequencing reads from three sampling dates (April, June, and August 2018) across nine sampling locations (see overview figure below). This includes rarefied datasets, scripts for the downstream analyses (written in R and python, integrated in a markdown file) and an appropriate conda environment.

Using this platform, the user will be able to replicate all results presented in our preprint [link will be added soon].

See [here](https://github.com/d-j-k/puntseq/tree/master/analysis) for the detailed description of our raw nanopore data pre-processing steps.

See [here](https://www.puntseq.co.uk/) for more updates about the PuntSeq project!


Overview of the experimental design:

![alt text](https://github.com/d-j-k/puntseq/blob/master/figure2.png)
