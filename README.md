Inputs

experimental samplesheet downloaded from ncbi with runs and a user added "condition" column
	The condition column serves to concisely indicate the condition of each run
	for example SI04 indicates the treatment (S) the tissue (I) and the day (04)

experimental fastqs from ncbi datasets download with bash and the fasterq-dump tool
runs=("SRR21585405" "SRR21585406")
for run in "${runs[@]}"; do
  fasterq-dump $run
done

transcriptome of species corresponding to experiment downloaded from ncbi
tx2gene.tsv file for the transcriptome from ncbi containing tx_id and gene_id columns
	I wrote a quick script that pulled transcript id and gene id directly from the downloaded transcriptome fasta

transcriptome of transposons for species pulled from dfam (I used the parse_embl pipeline)
tx2gene.tsv file the simply maps each transposon id to itself so that the pipeline will run

Run the pipeline once with the transcriptome of the full species and a second time with the transcriptome of transposons pulled from dfam