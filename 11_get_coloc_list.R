library(rtracklayer)

genes <- read.table("~/Desktop/PhD/projects/MCP_MDD_colocalisation//output/FUMA_job661883/genes.txt", header = T)
proteins <- read.table("~/Desktop/PhD/projects/CP_MDD_inflammation_MR/resources/3k_proteins.txt", header = T)



table(genes$symbol %in% proteins$AARSD1)


gtex_genes <- import("~/Desktop/PhD/projects/MCP_MDD_colocalisation/resources/gencode.v26.GRCh38.genes.gtf")
table(genes$symbol %in% gtex_genes$transcript_name)

