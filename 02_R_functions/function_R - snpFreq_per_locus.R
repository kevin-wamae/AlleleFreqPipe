#function for checking snp freq per column/locus
Func_SnpFreq <- function(x) {
  cbind(n = table(x), 
    freq = prop.table(table(x))*100)
}