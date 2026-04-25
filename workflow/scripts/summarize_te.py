import os
import pandas as pd
import copy

files = snakemake.input.difex
comparisons = [os.path.basename(file).replace("_deseq.csv","") for file in files]
empty_comparison_dict = {}
for comparison in comparisons:
	empty_comparison_dict[comparison] = ""

all_significant_TE = {}
comparison_dict = {}
for file in files:
	comparison = os.path.basename(file).replace("_deseq.csv","")
	comparison_dict[comparison] = {}
	deseq = pd.read_csv(file, index_col=0)
	deseq = deseq.dropna()
	deseq['padj'] = deseq['padj'].astype(float)
	deseq = deseq[deseq['padj'] < snakemake.params.p_val]
	for gene, row in deseq.iterrows():
		if gene not in all_significant_TE:
			all_significant_TE[gene] = copy.deepcopy(empty_comparison_dict)
		all_significant_TE[gene][comparison] = row['log2FoldChange']


df = pd.DataFrame.from_dict(all_significant_TE, orient='index')
df.to_csv(snakemake.output.gene_summary, index = True)

# create a filtered dataframe of only the TE/genes of interest
filtered_df = df[df.index.isin(snakemake.params.genes_of_interest)]
filtered_df.to_csv(snakemake.output.filtered_gene_summary, index = True)
