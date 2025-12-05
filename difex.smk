configfile: "config.yaml"
import os
import pandas as pd
import ast
transcriptome = config['transcriptome']
transcriptome_name = os.path.splitext(os.path.basename(transcriptome))[0]
output_folder = os.path.join(
	config['output_folder'], 
	f"{config['experiment_id']}_vs_{transcriptome_name}"
)
index = os.path.join(
	config['output_folder'],
	"salmon_indices",
	f"{transcriptome_name}_index"
)
samples_csv = config['samples']
fastqs_dir = os.path.dirname(samples_csv)

df = pd.read_csv(samples_csv)
RUNS = df['Run'].tolist()

# reference = config['reference']
comparisons = config['comparisons']
# print(reference)
# print(comparison)
print(transcriptome)

rule all:
	input: 
		copied_config = os.path.join(output_folder, "config.yaml"),
		gene_summary = os.path.join(output_folder, "significant_change_summaries", "full_summary.csv"),
		filtered_gene_summary = os.path.join(output_folder, "significant_change_summaries", "filtered_summary.csv")
		

rule create_index:
	input: transcriptome = transcriptome
	output: index = directory(index)
	shell:
		"salmon index -t {input.transcriptome} -i {output.index}"

rule quantify_with_salmon:
	input: 
		index = index, 
		run_1 = os.path.join(fastqs_dir, "{current_run}_1.fastq"),
		run_2 = os.path.join(fastqs_dir, "{current_run}_2.fastq")
	output: quant = directory(os.path.join(output_folder, "salmon_quants", "{current_run}"))
	shell:
		'''
			salmon quant -i {input.index} \
			-l A \
			-1 {input.run_1} \
			-2 {input.run_2} \
			-p 8 \
			--validateMappings \
			--gcBias \
			-o "{output.quant}"
		'''

rule run_deseq2:
	input:
		quantified_runs = expand(
			os.path.join(output_folder, "salmon_quants", "{current_run}"), 
			current_run = RUNS
		),
		transcriptome = transcriptome,
		tx2gene = config['tx2gene'],
		samples_csv = samples_csv,
	output:
		difex = os.path.join(output_folder, "deseq_output", "{comparison}"+"_deseq.csv")
	params:
		output_folder = output_folder,
		RUNS = RUNS,
		comparison = lambda wildcards: wildcards.comparison,
		reference = lambda wildcards: wildcards.comparison.split('_')[-1]
	script:
		"scripts/run_deseq2.R"

rule summarize_genes_of_interest:
	input:
		difex = expand(os.path.join(output_folder, "deseq_output", "{comparison}"+"_deseq.csv"), comparison = comparisons),
	output:
		gene_summary = os.path.join(output_folder, "significant_change_summaries", "full_summary.csv"),
		filtered_gene_summary = os.path.join(output_folder, "significant_change_summaries", "filtered_summary.csv")
	params:
		genes_of_interest = config['genes_or_te_of_interest'],
		p_val = config['p-value_cutoff']
	script:
		"scripts/summarize_te.py"

rule copy_config:
	output:
		copied_config = os.path.join(output_folder, "config.yaml")
	shell:
		"cp -f config.yaml {output.copied_config}"