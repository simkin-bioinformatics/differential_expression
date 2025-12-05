import os
import pandas as pd
import copy
deseq_output = "/home/charlie/git/differential_expression/output/PRJNA880828_vs_mus_te_from_dfam_uncombined/deseq_output"

files = snakemake.input.difex
comparisons = [os.path.basename(file).replace("_deseq.csv","") for file in files]
# comparisons = [file.replace("_deseq.csv","") for file in files]
empty_comparison_dict = {}
for comparison in comparisons:
	empty_comparison_dict[comparison] = ""

all_significant_TE = {}
comparison_dict = {}
for file in files:
	comparison = os.path.basename(file).replace("_deseq.csv","")
	comparison_dict[comparison] = {}
	deseq = pd.read_csv(file)
	# print(file)
	deseq = deseq.dropna()
	deseq['padj'] = deseq['padj'].astype(float)
	deseq = deseq[deseq['padj'] < snakemake.params.p_val]
	for index, row in deseq.iterrows():
		if row['Unnamed: 0'] not in all_significant_TE:
			all_significant_TE[row['Unnamed: 0']] = copy.deepcopy(empty_comparison_dict)
		all_significant_TE[row['Unnamed: 0']][comparison] = row['log2FoldChange']
		# all_significant_TE.add(row['Unnamed: 0'])
		# comparison_dict[comparison][row['Unnamed: 0']] = row['log2FoldChange']
		# if row['Unnamed: 0'] == 'RMER17B':
			# print(file)
			# print(f"TE: {row['Unnamed: 0']}, padj: {row['padj']}")
	# print(deseq)


df = pd.DataFrame.from_dict(all_significant_TE, orient='index')
df.to_csv(snakemake.output.gene_summary, index = True)
# print(df)

# create a filtered dataframe of only the TE/genes of interest
filtered_df = df[df.index.isin(snakemake.params.genes_of_interest)]
filtered_df.to_csv(snakemake.output.filtered_gene_summary, index = True)
# print(filtered_df)