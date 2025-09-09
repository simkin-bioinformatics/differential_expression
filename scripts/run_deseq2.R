# run DESeq2 on the salmon quant files
library("tximport")
library("tximportData")
library("DESeq2")

TRANSCRIPTOME <- snakemake@input[['transcriptome']]
samples <- read.csv(snakemake@input[['samples_csv']], header=TRUE)
RUNS <- snakemake@params[['RUNS']]
output_folder <- snakemake@params[['output_folder']]
files <- file.path(output_folder, "quants", RUNS, "quant.sf")
tx2gene <- read.csv(snakemake@input[['tx2gene']])

# import salmon quant data
txi <- tximport(files, type="salmon", tx2gene=tx2gene)
ddsTxi <- DESeqDataSetFromTximport(txi,
                                   colData = samples,
                                   design = ~ condition)

# relevel based on chosen reference
reference <- snakemake@params[['reference']]
comparison <- snakemake@params[['comparison']]
name <- paste0("condition_", comparison, "_vs_", reference)
ddsTxi$condition <- relevel(ddsTxi$condition, ref = reference)
dds <- DESeq(ddsTxi)

# get the unshrunk results for the chosen comparison
res <- results(dds, name=name)

# get the log fold change shrinkage results
resLFC <- lfcShrink(dds, coef=name, type="apeglm")

# order the results
resOrdered <- resLFC[order(res$pvalue),]

# write results to csv
write.csv(resOrdered, snakemake@output[['difex']], row.names = TRUE)


