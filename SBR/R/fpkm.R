"fpkm" <- function(counts,lengths)  {  # from https://support.bioconductor.org/p/91218/
          x = t(t(counts)/colSums(counts))
          x*1e9/lengths
} 
