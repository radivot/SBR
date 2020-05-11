"cpm" <- function(counts)  {  # from https://support.bioconductor.org/p/91218/
  t(t(counts)*1e6/colSums(counts))
} 
