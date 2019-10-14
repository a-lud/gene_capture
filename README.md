# Gene Capture Assembly Pipeline

This is a simple pipeline that aligns gene capture data to a transcriptome reference file (Trinity), generating a consensus 'gene sequence' at the other end.

## Requirements

This pipeline requires Python and a few of its libraries. I've included an example of how to install this pipeline in a contained environment below with all required software.

```
$ module load Anaconda3/5.0.1
$ conda create -n geneCapture_env argparse pandas mosdepth
$ conda activate geneCapture_env
```

Once you've created the environment, you'd just need to activate the environment in your `SLURM` submission script, like so:

```
#!/bin/bash
#SBATCH --job-name=gene_capture
#SBATCH -p batch
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 8
#SBATCH --time=8:00:00
#SBATCH --mem=16GB
#SBATCH -o /home/a1645424/fastdir/path/to/geneCapture_slurm/%x_%j.out
#SBATCH -e /home/a1645424/fastdir/path/to/geneCapture_slurm/%x_%j.err
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=firstName.lastName@adelaide.edu.au

## Load modules
module load Anaconda3/5.0.1

PIPELINE=/path/to/gene_capture/pipeline
READS=/path/to/reads
REFERENCES=/path/to/references
OUTDIR=/path/to/output_dir

## Running the pipeline
${PIPELINE}/getConsensus.py -kv /path/to/keyValue.csv -rd ${READS} -fd ${REFRENCES} -o ${OUTDIR} -af /path/to/adapters.fa -j ${SLURM_CPUS_PER_TASK}
```

The above is a slurm submission script and should be saved in a shell script file: e.g. `geneCapture_submission.sh`

**NOTE**: MosDepth is only available on Linux systems! It will not install on Mac.

## Usage

The pipeline has the following arguments:

**Key-value file (-kv)**: This is a two columned CSV file that has the columns `sample` and `reference`. For each sample, match it to a **unique** basename of the reference it belongs to. Do not provide the file extension.
    -kv /path/to/key_value.csv

**Reads file directory (-rd)**: A directory path to where the gene capture reads are located.
  
- -rd /path/to/reads_dir

**Reference file directory (-fd)**: A directory path to where the transcriptomic (or any reference) files are located.
  
- -fd /path/to/references_dir

**Output directory (-o)**: Provide the directory path to where you'd like the output of the pipeline to be stored.
  
- -o /path/to/output_dir

**Adapter file (-af)**: Provide the file path to a fasta file containing the adapter sequences used for your sequencing. If you do not have them, I've included a comprehensive adapter file as part of this repo.

- -af /path/to/adapters.fa

**Threads**: If you have a multi-threaded machine, provide the number of threads you'd like to use in multi-threaded processes (e.g. FastQC, Trimmomatic and BWA)

- -j 8
