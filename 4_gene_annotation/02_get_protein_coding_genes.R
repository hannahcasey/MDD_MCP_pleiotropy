library(rtracklayer)

gtf <- import("./resources/gencode.v49.annotation.gtf") 
pc <- subset(gtf, type == "gene" & mcols(gtf)$gene_type == "protein_coding")
gene_id <- sort(unique(mcols(pc)$gene_id))
ens <- sub("\\..*$", "", gene_id) 
writeLines(ens, "./resources/toppfun_background_protein_coding_symbols.txt")
