# run DESeq2 on the salmon quant files
library("tximport")
library("DESeq2")

samples <- read.csv(snakemake@input[['samples_csv']], header=TRUE)
RUNS <- snakemake@params[['RUNS']]
output_folder <- snakemake@params[['output_folder']]
files <- file.path(output_folder, "salmon_quants", RUNS, "quant.sf")
tx2gene <- read.csv(snakemake@input[['tx2gene']], sep="\t")

# import salmon quant data
txi <- tximport(files, type="salmon", tx2gene=tx2gene, dropInfReps = TRUE)
ddsTxi <- DESeqDataSetFromTximport(txi,
                                   colData = samples,
                                   design = ~ condition)

# relevel based on chosen reference
comparison <- snakemake@params[['comparison']]
reference <- snakemake@params[['reference']]

name <- paste0("condition_", comparison)
ddsTxi$condition <- relevel(ddsTxi$condition, ref = reference)
keep <- rowSums(counts(ddsTxi) >= 10) >= min(table(ddsTxi$condition))
dds <- DESeq(ddsTxi[keep, ])

# get the unshrunk results for the chosen comparison
# res <- results(dds, name=name)

# get the log fold change shrinkage results
resLFC <- lfcShrink(dds, coef=name, type="apeglm")

# order the results
resOrdered <- resLFC[order(resLFC$pvalue),]
# arc_summary <- resOrdered['Arc',]
# rownames(arc_summary)['Arc'] <- name

# write results to csv
write.csv(resOrdered, snakemake@output[['difex']], row.names = TRUE)