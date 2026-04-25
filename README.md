# Inputs
1. experimental metadata sheet downloaded from ncbi with runs and a user added "condition" column. The condition column serves to concisely indicate the condition of each run. For example SI04 indicates the treatment (S) the tissue (I) and the day (04)
2. experimental fastqs from ncbi datasets download with bash and the fasterq-dump tool.  There is a download_runs.smk pipeline that can be run if the metadata sheet is put into an empty directory.  The fastq files from the experiment will be downloaded into that same directory
3. transcriptome of species corresponding to experiment downloaded from ncbi
4. tx2gene.tsv file for the transcriptome from ncbi containing tx_id and gene_id columns.  I wrote a quick script that pulled transcript id and gene id directly from the downloaded transcriptome fasta
5. transcriptome of transposons for species pulled from dfam (I used the parse_embl pipeline)
6. tx2gene.tsv file the simply maps each transposon id to itself so that the pipeline will run

# Running the pipeline
1. Edit the config file to point to the sample sheet created from the metadata downloaded from the SRA website.  Be sure that you have added a condition column to it.  If you haven't downloaded the fastq files from the experiment go ahead and run the download_runs pipeline.
2. Edit the transcriptome and tx2gene in the config to point to the transcript downloaded for the organism from SRA
3. Choose an output folder
4. Edit the comparisons for all of the two-way comparisons between conditions that you wish to make. An example comparison is KI03_vs_SI03
5. Choose a p-value cutoff
6. run the pipeline, don't worry about the genes_or_te of interest for now.  You can go back and add some later
7. rerun the pipeline with the transcriptome and tx2gene file from dfam.  This will who the differential expression of transposons for the experimental dataset
8. After examining the full summary outputs you may want to go back and edit the genes_or_te of interest to focus on certain genes or transposons.
