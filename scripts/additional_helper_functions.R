

#if (!require("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install("Homo.sapiens")
# BiocManager::install("GenomicRanges")

library(Homo.sapiens)
library(dplyr)


# Code is from bNMF for T2D pipeline file format_bNMF_results.Rmd, slightly adapted for the purposes of this work
# snps is a vector of variants in the form "chr:pos" or "chr:pos:ref:alt"
# Returns a data frame with variant name, gene name and ENTREZID. 
# If a gene name was not found, defaults to variant name

query_locus_names <- function(snps){
  
  geneRanges <- function(db, column="ENTREZID"){
    g <- genes(db, columns=column)
    col <- mcols(g)[[column]]
    genes <- granges(g)[rep(seq_along(g), elementNROWS(col))]
    mcols(genes)[[column]] <- as.character(unlist(col))
    genes
  }
  
  gns = geneRanges(Homo.sapiens, column="SYMBOL")
  # need chr | start | end, with gene names for row names
  df_gns <- data.frame(gns)
  gtf.gene <- df_gns %>%
    mutate(chr = gsub("chr","",seqnames)) %>%
    dplyr::select(chr, start, end, SYMBOL)
  
  # need chr | start | end, with gene names for row names
  df_entrez <- geneRanges(Homo.sapiens, column="ENTREZID") %>%
    data.frame() %>%
    mutate(chr = gsub("chr","",seqnames)) %>%
    mutate(ChrPos = paste(chr,start,sep=":")) %>% 
    dplyr::select(ChrPos, ENTREZID)
  
  #' Convert from string to range
  #' 
  #' @param pos A vector of strings ex. chr1 2938302 2938329
  #' @param delim Delimiter for string splitting
  #' @param region Boolean of whether region or just one position
  #'
  #' @returns Dataframe of ranges
  #' 
  string2range <- function(pos, delim=' ', region=TRUE) {
    posp <- as.data.frame(do.call(rbind, strsplit(pos, delim)))
    posp[,1] <- posp[,1]
    posp[,2] <- as.numeric(as.character(posp[,2]))
    if(region) {
      posp[,3] <- as.numeric(as.character(posp[,3]))
    } else {
      posp[,3] <- posp[,2]
    }
    return(posp)
  }
  
  #' Convert from ranges to GRanges
  #' 
  #' @param df Dataframe with columns as sequence name, start, and end
  #' 
  #' @returns GRanges version 
  #' 
  range2GRanges <- function(df) {
    require(GenomicRanges)
    require(IRanges)
    gr <- GenomicRanges::GRanges(
      seqnames = df[,1],
      ranges=IRanges(start = df[,2], end = df[,3])
    )
    return(gr)
  }
  
  # convert SNPs to GRanges
  snps.ranges <- string2range(snps, delim=":", region=FALSE)
  snps.granges <- range2GRanges(snps.ranges)
  names(snps.granges) <- snps
  
  # convert genes to GRanges
  gtf.granges <- range2GRanges(gtf.gene)
  names(gtf.granges) <-  gtf.gene$SYMBOL  #gene.names
  
  hits <- GenomicRanges::nearest(snps.granges, gtf.granges)
  # make vector of SNPs to gene
  
  df_gns <- df_gns %>% 
    mutate(chr = gsub("chr","",seqnames)) %>%
    mutate(ChrPos = paste(chr,start,sep=":")) %>%
    dplyr::select(chr, start, end, ChrPos, SYMBOL) %>%
    merge(df_entrez, by="ChrPos")
  
  df_hits <- data.frame(gene=names(gtf.granges)[hits]) %>%
    mutate(variant = names(snps.granges)) %>%
    merge(df_gns[,c("SYMBOL","ENTREZID")], by.x="gene", by.y="SYMBOL") %>%
    filter(!duplicated(variant))
  
  # Make the duplicated genes have unique names.
  #If no gene is found for a variant, then the name of the variant is returned
  if (nrow(df_hits) > 0) {
    dup_genes <- df_hits$gene[!is.na(df_hits$gene)]
    duplicates <- data.frame(duplicated = table(dup_genes) > 1)
    duplicates$gene <- rownames(duplicates)
    
    for (r in 1:nrow(duplicates)) {
      row <- duplicates[r,]
      if (!is.na(row$duplicated) && row$duplicated) {
        idx <- which(df_hits$gene == row$gene)
        df_hits$gene[idx] <- paste0(row$gene, "_", seq_along(idx))
      }
    }
  }
  
  res <- data.frame(variant = snps) %>%
    left_join(df_hits, by="variant") %>%
    mutate(gene = dplyr::coalesce(gene, variant)) # coalese fills in empty genes w/ variant
  
  return(res)
}