Wdir = "your directory"
setwd(Wdir)


library(tidyr)
library(tidyverse)
library(dplyr)
library(stringr)
library(taxonomizr)
library(ape)
library(phyloseq)

rm(list=ls())

source("./tools/paf2count.R")
# taxonomir database

# this needs you to set it up

db = "F:/Taxonomir/taxonomir_database.sqlite"
#-----------------------------------------------------------------------------------------------
# 16S from centrifuge
#-----------------------------------------------------------------------------------------------
ant16S_otu = read.delim("../../OTU_tables/OTU_antidesc_16S_boxwoodchloroplast.tsv", row.names = 1, sep = "\t")
head(ant16S_otu)
nrow(ant16S_otu)
# modify columns
coln1 = gsub("assignment.nochim.|.tsv", "", colnames(ant16S_otu))
coln1
colnames(ant16S_otu) = coln1

# dir.create("./tsvfiles/")
taxaLin16S = centr2taxtable(ant16S_otu, db = db, writeout = T, outputname = "./tsvfiles/centri_ANTI_16S_taxa.tsv")


head(taxaLin16S)

# convert to matrix
ant16S_taxm = as.matrix(taxaLin16S)

# change row names to feature 1
rownames(ant16S_otu) = paste0("bOTU", seq_along(rownames(ant16S_otu)))
head(ant16S_otu)

# form OTU table
otu16S.ant = otu_table(ant16S_otu, taxa_are_rows = T)
rownames(ant16S_taxm) = rownames(ant16S_otu)


# form taxa table
tax16S.ant = tax_table(ant16S_taxm)

# first phyloseq object
phylo16S.ant = phyloseq(otu16S.ant, tax16S.ant)
phylo16S.ant

# phyloseq-class experiment-level object
# otu_table()   OTU Table:         [ 2549 taxa and 99 samples ]
# tax_table()   Taxonomy Table:    [ 2549 taxa by 7 taxonomic ranks ]

# updated
# otu_table()   OTU Table:         [ 2244 taxa and 99 samples ]
# tax_table()   Taxonomy Table:    [ 2244 taxa by 7 taxonomic ranks ]

# make tree
rntree16S.ant = rtree(ntaxa(phylo16S.ant), rooted = T, tip.label = taxa_names(phylo16S.ant))

# load meta data
meta16_ant = read.csv("../../MetaData/metadata_antid_16S.csv")
head(meta16_ant)
nrow(meta16_ant)

meta16_ant$SampleID = gsub("-", "_", meta16_ant$SampleID)
head(meta16_ant)
rownames(meta16_ant) = meta16_ant$SampleID

meta16_ant_sam = sample_data(meta16_ant)
coln = tools::toTitleCase(colnames(meta16_ant_sam))
coln
colnames(meta16_ant_sam) = coln
meta16_ant_sam$Month = ifelse(meta16_ant_sam$Collect_date == "6/16/2021", "June", ifelse(meta16_ant_sam$Collect_date == "8/26/2021", "August", ifelse(meta16_ant_sam$Collect_date == "10/18/2021", "October", "NA")))

psq16S.ant = merge_phyloseq(phylo16S.ant, meta16_ant_sam, rntree16S.ant)
psq16S.ant

# updated:
# phyloseq-class experiment-level object
# otu_table()   OTU Table:         [ 2244 taxa and 99 samples ]
# sample_data() Sample Data:       [ 99 samples by 13 sample variables ]
# tax_table()   Taxonomy Table:    [ 2244 taxa by 7 taxonomic ranks ]
# phy_tree()    Phylogenetic Tree: [ 2244 tips and 2243 internal nodes ]

# phyloseq-class experiment-level object
# otu_table()   OTU Table:         [ 2549 taxa and 99 samples ]
# sample_data() Sample Data:       [ 99 samples by 11 sample variables ]
# tax_table()   Taxonomy Table:    [ 2549 taxa by 7 taxonomic ranks ]
# phy_tree()    Phylogenetic Tree: [ 2549 tips and 2548 internal nodes ]

unique(tax16S.ant[,"Kingdom"])

# OTU1   "Bacteria" 
# OTU121 "Eukaryota"
# OTU821 NA  

#dir.create("./RDSdata")
saveRDS(psq16S.ant, "./RDSdata/original_phyloseq_16S_ANT.rds")



#############################################################################
# fungi
#############################################################################

rm(list = ls())

ITS_otu = read.delim("../../OTU_tables/antITS_otu.tsv", row.names = 1)
head(ITS_otu)
colnames(ITS_otu)

ITS_taxa = read.delim("../../OTU_tables/antITS_taxa.tsv", row.names = 1)
head(ITS_taxa)

ITS_taxa_clean = ITS_taxa %>%
  separate(Taxon, c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"), sep = ";")
head(ITS_taxa_clean)


ITS_taxa_clean$Kingdom = gsub("k__", "", ITS_taxa_clean$Kingdom)
ITS_taxa_clean$Phylum = gsub("p__", "",  ITS_taxa_clean$Phylum)
ITS_taxa_clean$Class = gsub("c__", "",   ITS_taxa_clean$Class)
ITS_taxa_clean$Order = gsub("o__", "",   ITS_taxa_clean$Order)
ITS_taxa_clean$Family = gsub("f__", "",  ITS_taxa_clean$Family)
ITS_taxa_clean$Genus = gsub("g__", "",   ITS_taxa_clean$Genus)
ITS_taxa_clean$Species = gsub("s__", "", ITS_taxa_clean$Species)


ITS_taxam = as.matrix(ITS_taxa_clean)
ITS_taxam[1:5, 1:7]

# phyloseq
phyloITS_otu = otu_table(ITS_otu, taxa_are_rows = T)

phyloITS_taxa = tax_table(ITS_taxam)

phyloITS = phyloseq(phyloITS_otu, phyloITS_taxa)

rntreeITS = rtree(ntaxa(phyloITS), rooted = T, tip.label = taxa_names(phyloITS))


metaITS = read.csv("../../MetaData/metadata_antid_ITS.csv", row.names = 1)

head(metaITS)
coln = tools::toTitleCase(colnames(metaITS))
coln
colnames(metaITS) = coln
metaITS.sam = sample_data(metaITS)

metaITS.sam$Month = ifelse(metaITS.sam$Collect_date == "6/16/2021", "June", ifelse(metaITS.sam$Collect_date == "8/26/2021", "August", ifelse(metaITS.sam$Collect_date == "10/18/2021", "October", "NA")))
metaITS.sam

ncol(metaITS.sam)

psqITS = merge_phyloseq(phyloITS, 
                        metaITS.sam, 
                        rntreeITS)

psqITS

# phyloseq-class experiment-level object
# otu_table()   OTU Table:         [ 911 taxa and 99 samples ]
# sample_data() Sample Data:       [ 99 samples by 10 sample variables ]
# tax_table()   Taxonomy Table:    [ 911 taxa by 7 taxonomic ranks ]
# phy_tree()    Phylogenetic Tree: [ 911 tips and 910 internal nodes ]

saveRDS(psqITS, "./RDSdata/original_ITS_antides_with_minimap2.rds")


###############################################################################
# data cleaning
###############################################################################

# cleanning of the data
rm(list=ls())

devtools::source_gist("8d0ca4206a66be7ff6d76fc4ab8e66c6")

an16S = readRDS("./RDSdata/original_phyloseq_16S_ANT.rds")
anITS = readRDS("./RDSdata/original_ITS_antides_with_minimap2.rds")

m_tax16 = as(tax_table(an16S), "matrix")
m_tax16[1:4,1:7]

m_taxITS = as(tax_table(anITS), "matrix")
m_taxITS[1:5, 1:7]
m_taxITS[m_taxITS=="unidentified"] = NA
tax_table(anITS) = m_taxITS



meta16 = data.frame(sample_data(an16S))
head(meta16)
metaITS = data.frame(sample_data(anITS))
head(metaITS)


min(sample_sums(an16S)) #491
min(sample_sums(anITS)) #53

library(vegan)

class(otu_table(an16S)) # integer


otu16S.mtx = otu_table(an16S)
class(otu16S.mtx) = "matrix"

otuITS.mtx = otu_table(anITS)
class(otuITS.mtx) = "matrix"

rarecurve(t(otu16S.mtx), step = 100)
abline(v = 1000, col="red")

rarecurve(t(otuITS.mtx), step = 100)
abline(v = 1000, col="red")

# prune sample
an16S_prune_sam = subset_samples(an16S, sample_sums(an16S) >= 1000)
min(sample_sums(an16S_prune_sam)) #1351

anITS_prune_sam = subset_samples(anITS, sample_sums(anITS) >= 1000)
min(sample_sums(anITS_prune_sam)) #7402


# prune taxa ========== keep Cyanobacteria but remove Boxwood chloroplasts
an16S_prune = subset_taxa(an16S_prune_sam, !is.na(Phylum) & Kingdom == "Bacteria")
anITS_prune = subset_taxa(anITS_prune_sam, !is.na(Phylum))

min(sample_sums(an16S_prune)) #71
min(taxa_sums(an16S_prune)) #1

min(sample_sums(anITS_prune)) #7392
min(taxa_sums(anITS_prune)) #1

# further prunning
an16S_prune2 = subset_taxa(an16S_prune, taxa_sums(an16S_prune) >= 10)
anITS_prune2 = subset_taxa(anITS_prune, taxa_sums(anITS_prune) >= 10)


an16S_prune2

# updated 3/2/2023
# phyloseq-class experiment-level object
# otu_table()   OTU Table:         [ 599 taxa and 98 samples ]
# sample_data() Sample Data:       [ 98 samples by 13 sample variables ]
# tax_table()   Taxonomy Table:    [ 599 taxa by 7 taxonomic ranks ]
# phy_tree()    Phylogenetic Tree: [ 599 tips and 598 internal nodes ]

# =========== old =======================================
# otu_table()   OTU Table:         [ 997 taxa and 96 samples ]
# sample_data() Sample Data:       [ 96 samples by 11 sample variables ]
# tax_table()   Taxonomy Table:    [ 997 taxa by 7 taxonomic ranks ]
# phy_tree()    Phylogenetic Tree: [ 997 tips and 996 internal nodes ]

# phyloseq-class experiment-level object
# otu_table()   OTU Table:         [ 997 taxa and 98 samples ]
# sample_data() Sample Data:       [ 98 samples by 13 sample variables ]
# tax_table()   Taxonomy Table:    [ 997 taxa by 7 taxonomic ranks ]
# phy_tree()    Phylogenetic Tree: [ 997 tips and 996 internal nodes ]


anITS_prune2

# updated 3/2/2023
# phyloseq-class experiment-level object
# otu_table()   OTU Table:         [ 472 taxa and 96 samples ]
# sample_data() Sample Data:       [ 96 samples by 12 sample variables ]
# tax_table()   Taxonomy Table:    [ 472 taxa by 7 taxonomic ranks ]
# phy_tree()    Phylogenetic Tree: [ 472 tips and 471 internal nodes ]

# phyloseq-class experiment-level object
# otu_table()   OTU Table:         [ 472 taxa and 96 samples ]
# sample_data() Sample Data:       [ 96 samples by 10 sample variables ]
# tax_table()   Taxonomy Table:    [ 472 taxa by 7 taxonomic ranks ]
# phy_tree()    Phylogenetic Tree: [ 472 tips and 471 internal nodes ]

# phyloseq-class experiment-level object
# otu_table()   OTU Table:         [ 472 taxa and 96 samples ]
# sample_data() Sample Data:       [ 96 samples by 12 sample variables ]
# tax_table()   Taxonomy Table:    [ 472 taxa by 7 taxonomic ranks ]
# phy_tree()    Phylogenetic Tree: [ 472 tips and 471 internal nodes ]

library(microbiomeutilities)
library(microbiome)

# saveRDS(an16S_prune2, "./RDSdata/pruned_16S_org.rds")
# saveRDS(anITS_prune2, "./RDSdata/pruned_ITS_org.rds")


bst16 = microbiomeutilities::format_to_besthit(an16S_prune2)
bstITS = microbiomeutilities::format_to_besthit(anITS_prune2)

# saveRDS(bst16, "./RDSdata/bsthit_pruned_16S_org.rds")
# saveRDS(bstITS, "./RDSdata/bsthit_pruned_ITS_org.rds")

#########################################################
# ends here
#########################################################
