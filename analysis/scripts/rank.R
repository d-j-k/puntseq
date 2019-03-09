##############
## Taxonomy ##
##############

library(readr)
library(tidyr)
library(dplyr)


taxID_info = readRDS("analysis/data/taxID_info.rds")


get_rank <- function(taxID, taxID_info){
  # return rank for taxID
  if (taxID %in% names(taxID_info)){
    info = taxID_info[[paste0(taxID)]]
    rank = info$rank[info$id == taxID]
  } else {
    rank = NA
  }
  # WORKAROUND bad annotation
  if (identical(rank, character(0))){
    rank = NA
  }
  return(rank)
}

get_taxID_phylum <- function(taxID, taxID_info){
  # return rank for taxID
  if (taxID %in% names(taxID_info)){
    info = taxID_info[[paste0(taxID)]]
    if ("phylum" %in% info$rank){
      taxID_phylum = info$id[info$rank == "phylum"]
    } else {
      taxID_phylum = NA
    }
  } else {
    taxID_phylum = NA
  }
  # WORKAROUND bad annotation
  if (identical(taxID_phylum, character(0))){
    taxID_phylum = NA
  }
  return(taxID_phylum)
}


get_taxID_rank <- function(taxID, rank, taxID_info){
  # return taxID for rank above taxID
  if (taxID %in% names(taxID_info)){
    info = taxID_info[[paste0(taxID)]]
    if (rank %in% info$rank){
      taxID_rank = info$id[info$rank == rank]
    } else {
      taxID_rank = NA
    }
  } else {
    taxID_rank = NA
  }
  # WORKAROUND bad annotation
  if (identical(taxID_rank, character(0))){
    taxID_rank = NA
  }
  return(taxID_rank)
}


get_rank_counts <-function(results_centrifuge, taxID_info){
  results_centrifuge = results_centrifuge %>%
                        mutate(rank = unlist(sapply(results_centrifuge$taxID, function(taxID) get_rank(taxID, taxID_info))))
  results_centrifuge = results_centrifuge %>%
                        mutate(taxID_species = as.integer(unlist(sapply(results_centrifuge$taxID, function(taxID) get_taxID_rank(taxID, "species", taxID_info)))))
  results_centrifuge = results_centrifuge %>%
                        mutate(taxID_genus = as.integer(unlist(sapply(results_centrifuge$taxID, function(taxID) get_taxID_rank(taxID, "genus", taxID_info)))))
  results_centrifuge = results_centrifuge %>%
                        mutate(taxID_family = as.integer(unlist(sapply(results_centrifuge$taxID, function(taxID) get_taxID_rank(taxID, "family", taxID_info)))))
  results_centrifuge = results_centrifuge %>%
                        mutate(taxID_phylum = as.integer(unlist(sapply(results_centrifuge$taxID, function(taxID) get_taxID_rank(taxID, "phylum", taxID_info)))))

  # filter for species and unique reads
  df_unique_species = results_centrifuge %>%
    unite(readID_species, readID, taxID_species, remove=F) %>%
    filter(!(duplicated(readID_species) | duplicated(readID_species, fromLast = TRUE))) %>%
    filter(!duplicated(taxID_species))

  df_unique_genus = results_centrifuge %>%
    unite(readID_genus, readID, taxID_genus, remove=F) %>%
    filter(!(duplicated(readID_genus) | duplicated(readID_genus, fromLast = TRUE))) %>%
    filter(!duplicated(taxID_genus))

  df_unique_family = results_centrifuge %>%
    unite(readID_family, readID, taxID_family, remove=F) %>%
    filter(!(duplicated(readID_family) | duplicated(readID_family, fromLast = TRUE))) %>%
    filter(!duplicated(taxID_family))

    df_unique_phylum = results_centrifuge %>%
    unite(readID_phylum, readID, taxID_phylum, remove=F) %>%
    filter(!(duplicated(readID_phylum) | duplicated(readID_phylum, fromLast = TRUE))) %>%
    filter(!duplicated(taxID_phylum))

  rank_counts = c(n_species = dim(df_unique_species)[1],
                  n_genus = dim(df_unique_genus)[1],
                  n_family = dim(df_unique_family)[1],
                  n_phylum = dim(df_unique_phylum)[1])
  return(rank_counts)
}
