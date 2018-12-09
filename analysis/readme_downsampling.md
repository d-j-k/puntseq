# Downsampling

Downsampling of the reads to assess whether we sequenced deeply enough.
The scripts used for the analysis are in `analysis/scripts`.

`get_taxinfo.R` takes all centrifuge outputs and queries the taxIDs in NCBI. 
The results of this query are saved in a table `taxID_info.rds` which can be queried by taxID. This allows to get the species/genus name and taxIDs/names of ranks above the queried taxID, i.e. this would allow you to convert the whole mapped reads to a level of phylum from the current mixed output.

`ran.R` is used to read the taxID table. `assess_downsampling.R` runs the downsampling (relatively slow) and generates the plots.

## Ideas
- Good–Turing frequency estimation
