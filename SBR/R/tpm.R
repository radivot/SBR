"tpm" <- function(counts,lengths)  {  # from https://support.bioconductor.org/p/91218/
        x = counts/lengths   # lengths units of kb or b do not matter
        t(t(x)*1e6/colSums(x))
} 
