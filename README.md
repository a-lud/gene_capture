# Gene Capture Assembly Pipeline

This is a simple pipeline that aligns gene capture data to a transcriptome reference file (Trinity), generating a consensus 'gene sequence' at the other end.

## Requirements

This pipeline requires Python and a few of its libraries. I've included an example of how to install this pipeline in a contained environment below with all required software.

```
$ module load Anaconda3/5.0.1
$ conda create -n geneCapture_env argparse pandas mosdepth trimmomatic bwa fastqc

## Can activate the environment using the command below
$ conda activate geneCapture_env
```

If all the software installs correctly, you should be met with a message demonstrating how to activate the conda environment. See the section `Running the pipeline on Phoenix` below to see how to run the pipeline.

**NOTE**: MosDepth is only available on Linux systems! It will not install on Mac.

## The pipeline: getConsensus.py

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

Below is the full help page

```
(CODEML_env) [13:08:06] alastairludington:gene_capture $ ./getConsensus.py -h
usage: getConsensus.py [-h] -kv KEYVALUEPATH -rd READSDIR [-fd REFDIR] -o
                       OUTPUTDIR -af ADAPTERFILE [-j THREADS]

Pipeline to assemble gene capture data using closely related references

optional arguments:
  -h, --help            show this help message and exit
  -kv KEYVALUEPATH, --keyValuePath KEYVALUEPATH
                        File path to key-value table
  -rd READSDIR, --readsDir READSDIR
                        Directory path to gene capture data
  -fd REFDIR, --refDir REFDIR
                        Directory path to reference fasta files
  -o OUTPUTDIR, --outputDir OUTPUTDIR
                        Directory path for output location
  -af ADAPTERFILE, --adapterFile ADAPTERFILE
                        File path to adapter sequence file for Trimmomatic
  -j THREADS, --threads THREADS
                        Number of threads for different tasks (default: 2)

Developed: Alastair Ludington
Institution: University of Adelaide
Date: 11/10/2019
```

## Running the pipeline on Phoenix

Once you've created your conda environment with all the required software, you need to create a `SLURM` submission script to run the pipeline. I've generated a mock template below of what needs to be in the script.

```
#!/bin/bash
#SBATCH --job-name=gene_capture
#SBATCH -p batch
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 8
#SBATCH --time=8:00:00
#SBATCH --mem=16GB
#SBATCH -o /home/a123456/fastdir/path/to/geneCapture_slurm/%x_%j.out
#SBATCH -e /home/a123456/fastdir/path/to/geneCapture_slurm/%x_%j.err
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

Copy the above contents into a file called `geneCapture_submission.sh`, or something similar.

This would then be executed from Phoenix's head node using the command `$ sbatch geneCapture_submission.sh`, which
would submit the job to the Phoenix execution queue.
