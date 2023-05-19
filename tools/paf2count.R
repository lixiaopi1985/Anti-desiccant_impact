paf2countdf2 = function(paf_path, taxadb="/media/swaggyp1985/HDD6T/VT_Projects_2020/Database/16S/SILVA_138.1/curated_16S_SILVA_138-1_ssu_nr99/curated_taxa_silva_138-1_ssu_nr99.tsv", dlim="\t", taxcol = "Taxon", mergetaxaby = "Feature.ID", filterbyblock=T, blocksize=1000, writeout=T, outdir="paf_count_out"){
  
  
  # for(i in c("readr", "taxonomizr", "doBy", "plyr", "dplyr", "stringr", "data.table")){
  #   if(!require(i)){
  #     install.packages(i)
  #   }
  # }
  # works with know database curated locally
  library(readr)
  # library(taxonomizr)
  # library(doBy)
  library(plyr)
  library(dplyr)
  library(stringr)
  library(data.table)
  
  # load paf file
  paf = readr::read_delim(paf_path, dlim, escape_double = F, col_names = F, trim_ws = T)
  coln = c("q_name", "q_len", "q_start", "q_end", "strand", "target_name", "target_len", "target_start", "target_end", "n_res_match", "align_block", "map_Q", "NM", "ms", "AS", "nn", "tp")
  
  df<-paf[,1:17] # keep useful columns
  colnames(df) = coln
  # AS score - Smith Waterman alignment score
  df$AS = as.numeric(gsub("AS:i:", "", df$AS))
  df2 = cbind.data.frame(df$q_name, df$target_name, df$n_res_match, df$align_block, df$AS, df$map_Q)
  colnames(df2) = c("Query", "Target", "n_res_match", "align_block", "AS", "MapQ")
  df2 = as.data.table(df2)
  # n base match percentage
  # df2$per_match = 100*(df2$n_res_match/df2$align_block)
  
  # filter align block size
  if(filterbyblock){
    df2 = df2[df2$align_block>=blocksize,]
  }
  # aggregate data by query, taxID, per match, mapQ, take max of AS of each group
  # why can't we filter just by max AS grouped by query and taxID
  # print("Aggregate dataframe AS ~ Query + Target + per_match + MapQ, get max AS for each group")
  # aggdf = summaryBy(AS ~ Query + Target + per_match + MapQ , data = df2, FUN=max)
  
  # check if there is a db used by taxonomizr, if not build
  #-----------------------------------------------------------------
  # print("Checking if the database, if not, build database by taxonomizr")
  # print("This will take a while")
  # 
  # if(! db %in% list.files()){
  #   warning(paste0("<<< ", db, " not in the file directory, making a new database with the name >>> "))
  #   prepareDatabase(db)
  # }
  
  # taxonomy table for the accession ID
  # acessID = unique(df$target_name)
  # taxIDs = accessionToTaxa(acessID, "accessionTaxa.sql")
  # head(taxIDs)
  # linage = getTaxonomy(taxIDs, "accessionTaxa.sql")
  # taxaT = data.frame(accID = acessID, taxID = taxIDs, lin = linage, row.names = NULL)
  # taxaT = as.data.table(taxaT)
  taxadf = read.delim(taxadb)
  head(taxadf)
  
  #----------------------------------------------------------------------------------
  print("Merge aggregated df with taxonomic information")
  dfmerge = merge(df2, taxadf, by.x="Target", by.y=mergetaxaby, all.x=T) # change aggdf > df2
  # remove NA taxID, which indicates this entry has been removed in the NCBI database
  print("Remove NA taxID, which which indicates this entry has been removed in the NCBI database")
  dfmerge2 = dfmerge[!is.na(dfmerge[[taxcol]]),]
  
  print(paste0("Removed ", nrow(dfmerge)-nrow(dfmerge2), " reads with NA tax ID"))
  # dfmerge2$taxID = as.character(dfmerge2$taxID)
  
  # sort by AS max score and remove duplicate rows
  print("Group by read, sort by AS score, take the first entry")
  taxID = sym(taxcol)
  dfcount = dfmerge2 %>%
    group_by(Query) %>%
    arrange(desc(abs(AS))) %>% # AS.max to AS
    slice(which.max(AS)) %>%
    dplyr::select(Query, !!taxID) %>%
    dplyr::group_by(!!taxID) %>%
    dplyr::summarise(count = n())
  
  
  # transform to count data
  # keep 2 columns now
  # df_keep = dfsort[, c("taxID", "Query")]
  # df_taxa_lin = dfsort[, c("taxID", "lin.superkingdom", "lin.phylum", "lin.class", "lin.order", "lin.family", "lin.genus", "lin.species")]
  # group by Target
  # print("Generate read counts based on taxon")
  # dfcount = df_keep %>%
  #   group_by(taxID) %>%
  #   dplyr::summarise(count = n())
  # 
  
  # rowN = dfcount$taxID
  # dfout = dfcount[, "count", keep=T]
  # dfout_m = as.data.frame(dfout)
  # row.names(dfout_m) = rowN
  
  SampleName = gsub("-nochim-", "", str_extract(basename(paf_path), "-.+_barcode[0-9]{1,2}"))
  dfcount = as.data.table(dfcount)
  colnames(dfcount) = c(taxcol, SampleName)
  # dfout_m = as.data.frame(dfout_m)
  
  # output file
  if(writeout){
    if(!dir.exists(outdir)){
      dir.create(outdir)
    }
    
    pathout = file.path(outdir, paste0(SampleName, ".tsv"))
    write.table(dfcount, file=pathout, sep = "\t", quote=F, row.names=F, col.names = T)
  }
  
  return(dfcount)
}

split_table = function(df, outfile){
  # change labels for each row
  
  taxDF = df[, "Taxon"]
  rownames(taxDF) = paste0("OTU", seq(1, nrow(taxDF)))
  taxDF$OTU = rownames(taxDF)
  taxDF = taxDF[, c("OTU", "Taxon")]
  
  
  countDF = df[, 2:ncol(df)]
  rownames(countDF) = paste0("OTU", seq(1, nrow(countDF)))
  
  cln = colnames(countDF)
  # print(cln)
  countDF$OTU = rownames(countDF)
  # print(head(countDF))
  newcoln = as.character(c("OTU", cln))
  # print(newcoln)
  countDF = countDF[, ..newcoln]
  
  
  write.table(taxDF, file=paste0(outfile, "_taxa.tsv"), sep = "\t", col.names = T, row.names = F, quote = F)
  write.table(countDF, file=paste0(outfile, "_otu.tsv"), sep="\t", col.names = T, row.names = F, quote = F)
}


centr2taxtable = function(centriout, db="accessionTaxa.sqlite", writeout=T, outputname="centrifuge_nonchime_taxonomy_table.tsv"){
  
  if(!file.exists(db)){
    taxonomizr::prepareDatabase(db)
  }
  
  taxaID = rownames(centriout)
  taxaLin = getTaxonomy(taxaID, db)
  print(head(taxaLin))
  taxaLin_rownames = trimws(rownames(taxaLin))
  print(identical(taxaID, taxaLin_rownames)) # true
  
  rownames(taxaLin) = taxaLin_rownames
  colnames(taxaLin) = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
  
  if(writeout){
    write.table(taxaLin, outputname, sep = "\t", quote = F, row.names = T)
  }
  
  return(taxaLin)
}
