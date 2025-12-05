import pandas as pd
import os

difex = snakemake.input.difex
gene_summary = snakemake.output.gene_summary
gene = snakemake.params.gene

combined_df = pd.DataFrame({})
for csv in difex:
	df = pd.read_csv(csv)
	df.rename(columns = {'Unnamed: 0': 'Comparison'}, inplace = True)
	df = df[df['Comparison'] == gene]
	df['Comparison'] = df['Comparison'].replace(gene, gene+' in ' + os.path.basename(csv).replace('_deseq.csv', ''))
	if combined_df.empty:
		combined_df = df
	combined_df = pd.concat([combined_df, df], ignore_index = True)

pvalue_sorted_df = combined_df.sort_values(by = 'pvalue').drop_duplicates()
pvalue_sorted_df.to_csv(gene_summary, index = False)
