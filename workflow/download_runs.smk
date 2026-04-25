configfile: "config/config.yaml"
import pandas as pd
from pathlib import Path
output_folder = Path(config['samples']).parent

df = pd.read_csv(config['samples'])
RUNS = df['Run'].tolist()

rule all:
    input: 
        expand(output_folder / "{run}_1.fastq.gz", run=RUNS),
        expand(output_folder / "{run}_2.fastq.gz", run=RUNS)

rule download_runs:
  input: 
    samples = config['samples']
  output: 
    run_1 = output_folder / "{run}_1.fastq",
    run_2 = output_folder / "{run}_2.fastq",
  params:
    output_folder = output_folder
  shell: 
    '''
    fasterq-dump {wildcards.run} -O {params.output_folder}
    '''

rule compress_runs:
    input:
        run_1 = output_folder / "{run}_1.fastq",
        run_2 = output_folder / "{run}_2.fastq",
    output:
        run_1_gz = output_folder / "{run}_1.fastq.gz",
        run_2_gz = output_folder / "{run}_2.fastq.gz",
    shell: 
        '''
        pigz {input.run_1}
        pigz {input.run_2}
        ''' 