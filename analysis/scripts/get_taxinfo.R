#######################
## Get Taxonomy Info ##
#######################

library(taxize)
# demh NCBI key
# usethis::edit_r_environ()
# ENTREZ_KEY='12345...'
library(readr)

files_centrifuge = list.files("analysis/data/", patter="tab$", recursive=T, full.names=T)
taxID_unique = c()
taxID_info = c()

# get list of all taxIDs
for (file_centrifuge in files_centrifuge){
  out_centrifuge = read_tsv(file_centrifuge, col_types="cciiiiii")
  taxID_unique = unique(c(taxID_unique, out_centrifuge$taxID))
}


# remove 0 i.e. no taxonomy
taxID_unique= taxID_unique[!(taxID_unique == 0)]

# query taxIDs in NCBI
taxID_unique_chunks = split(taxID_unique, ceiling(seq_along(taxID_unique)/20))

i = 0
for (taxID_chunk in taxID_unique_chunks){
  i = i + 1
  print(i)
  rank_info = NULL
  attempt = 1
  while(is.null(rank_info) && attempt <= 10) {
    attempt <- attempt + 1
    try(rank_info <- classification(taxID_chunk, db = 'ncbi'))
  }
  taxID_info = c(taxID_info, rank_info)
}

# save query results
saveRDS(taxID_info, "analysis/data/taxID_info.rds")
